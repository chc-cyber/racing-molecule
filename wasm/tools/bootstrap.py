"""Workspace-local toolchain and pinned upstream sources; no global installs."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
LOCAL = WORKSPACE / ".local" / "wasm"
PINS = json.loads((ROOT / "dependencies.json").read_text())
SOURCE_CHECKSUMS = json.loads((ROOT / "source-checksums.json").read_text())


def download(url: str, name: str) -> Path:
    path = LOCAL / "downloads" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.is_file():
        print(f"Downloading {name}", flush=True)
        request = urllib.request.Request(url, headers={"User-Agent": "molecular-car-wasm-build"})
        with urllib.request.urlopen(request, timeout=120) as source, path.with_suffix(path.suffix + ".part").open("wb") as target:
            while data := source.read(1024 * 1024):
                target.write(data)
        path.with_suffix(path.suffix + ".part").replace(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != SOURCE_CHECKSUMS[name]:
        raise RuntimeError(f"Dependency archive does not match pinned SHA-256: {name}")
    receipt = LOCAL / "downloads" / "receipts.json"
    entries = json.loads(receipt.read_text()) if receipt.is_file() else {}
    previous = entries.get(name)
    if previous and previous != {"url": url, "sha256": digest}:
        raise RuntimeError(f"Downloaded dependency changed: {name}")
    entries[name] = {"url": url, "sha256": digest}
    receipt.write_text(json.dumps(entries, indent=2) + "\n")
    return path


def extract(archive: Path, target: Path):
    if (target / ".extracted").is_file():
        return
    target.mkdir(parents=True, exist_ok=True)
    resolved = target.resolve()
    with tarfile.open(archive) as bundle:
        members = bundle.getmembers()
        prefix = members[0].name.split("/")[0] + "/"
        for member in members:
            if member.name == prefix[:-1]:
                continue
            if not member.name.startswith(prefix):
                raise RuntimeError("Unexpected archive layout")
            member.name = member.name[len(prefix):]
            if not (target / member.name).resolve().is_relative_to(resolved):
                raise RuntimeError("Archive member escapes dependency directory")
            if member.issym() or member.islnk():
                # Upstream documentation has links; Windows extraction must not
                # create links or follow them outside the workspace. None are
                # required by our selected C++ build targets.
                continue
            bundle.extract(member, target, filter="data")
    (target / ".extracted").write_text(hashlib.sha256(archive.read_bytes()).hexdigest() + "\n")


def main():
    LOCAL.mkdir(parents=True, exist_ok=True)
    sources = [
        ("emsdk", f"https://codeload.github.com/emscripten-core/emsdk/tar.gz/{PINS['emsdk_commit']}", "emsdk.tar.gz"),
        ("lammps", f"https://codeload.github.com/lammps/lammps/tar.gz/{PINS['lammps_git_object']}", "lammps.tar.gz"),
        ("rdkit", f"https://codeload.github.com/rdkit/rdkit/tar.gz/{PINS['rdkit_commit']}", "rdkit.tar.gz"),
        ("boost", "https://archives.boost.io/release/1.85.0/source/boost_1_85_0.tar.gz", "boost.tar.gz"),
        ("catch2", f"https://codeload.github.com/catchorg/Catch2/tar.gz/{PINS['catch2']}", "catch2.tar.gz"),
        ("better-enums", f"https://codeload.github.com/aantron/better-enums/tar.gz/{PINS['better_enums_commit']}", "better-enums.tar.gz"),
    ]
    for name, url, filename in sources:
        extract(download(url, filename), LOCAL / name)
        print(f"Ready: {name}", flush=True)
    buildtools = LOCAL / "buildtools"
    if not (buildtools / "cmake" / "data" / "bin" / "cmake.exe").is_file():
        subprocess.run([sys.executable, "-m", "pip", "install", "--disable-pip-version-check", "--target", str(buildtools),
                        f"cmake=={PINS['cmake']}", f"ninja=={PINS['ninja']}"], check=True)
    # The SDK activation is embedded under .local; never use --global.
    sdk = LOCAL / "emsdk"
    env = os.environ.copy()
    env["EMSDK_PYTHON"] = sys.executable
    for action in ("install", "activate"):
        subprocess.run([sys.executable, str(sdk / "emsdk.py"), action, PINS["emscripten"]], cwd=sdk, env=env, check=True)
    print("Workspace-local WASM dependencies are ready.", flush=True)


if __name__ == "__main__":
    main()
