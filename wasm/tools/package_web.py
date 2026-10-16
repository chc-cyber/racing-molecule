"""Create a static site using the shared renderer and a browser-only transport."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import shutil
import re
import zipfile
import argparse
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--no-archive',action='store_true',help='Refresh local files without rebuilding the ZIP during development')
    args=parser.parse_args()
    target = ROOT / "wasm" / "dist"
    frontend = ROOT / "frontend"
    target.mkdir(parents=True, exist_ok=True)
    for required in ("browser/runtime/molecular-car.wasm", "browser/runtime/molecular-car.mjs", "browser/assets/model.json",
                     "browser/assets/fourier-atlas.json", "browser/assets/fourier-atlas-0.bin.gz"):
        if not (target / required).is_file():
            raise RuntimeError(f"Build and export browser assets before packaging: {required}")
    for filename in ("app.js", "style.css", "gold-surface.js", "transport.js", "start-pose.js", "start-pose.css", "race-course.mjs", "competition.mjs", "competition.css", "race-theme.css", "design-system.css", "app-state.mjs", "race-result.mjs", "race-result.css", "playback.mjs", "local-scores.mjs", "local-scores-ui.mjs", "local-scores.css"):
        shutil.copyfile(frontend / filename, target / filename)
    for source in (frontend / "browser").glob("*.mjs"):
        shutil.copyfile(source, target / "browser" / source.name)
    shutil.copytree(frontend / "vendor", target / "vendor", dirs_exist_ok=True)
    shutil.copytree(frontend / "browser/vendor", target / "browser/vendor", dirs_exist_ok=True)
    # Keep home, workshop and race under one origin and under relative paths,
    # including a GitHub Pages repository prefix.
    html = (frontend / "index.html").read_text(encoding="utf-8")
    html = html.replace('<html lang="zh-CN">', '<html lang="zh-CN" data-engine="wasm">')
    html = html.replace('href="/', 'href="./').replace('src="/app.js"', 'src="./app.js"')
    html = html.replace('"three":"/vendor/', '"three":"./vendor/')
    html = html.replace('href="./" aria-label="Racing Molecule 首页"', 'href="./#home" data-app-route="home" aria-label="Racing Molecule 首页"')
    html = html.replace('href="./editor/index.html"', 'href="./#garage" data-app-route="garage"')
    gear = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="m9.5 3-.6 2.4-1.5.9L5 5.6 2.5 9.9l1.8 1.7v1.8l-1.8 1.7L5 19.4l2.4-.7 1.5.9.6 2.4h5l.6-2.4 1.5-.9 2.4.7 2.5-4.3-1.8-1.7v-1.8l1.8-1.7L19 5.6l-2.4.7-1.5-.9L14.5 3Z"/><circle cx="12" cy="12.5" r="3"/></svg>'
    html = html.replace('<span class="engine-badge"', '<a class="text-button shell-link" target="_top" href="./#settings" data-app-route="settings">'+gear+'设置</a><span class="engine-badge"')
    html = html.replace('连接中…', '加载中…')
    (target / "race.html").write_text(html, encoding="utf-8")
    ui = ROOT / "web-ui" / "dist"
    (target / "ui").mkdir(exist_ok=True)
    for name in ("index.html", "app.js", "style.css", "minimal.css", "integration.css"):
        shutil.copyfile(ui / name, target / "ui" / name)
    shutil.copytree(ui / "assets", target / "ui" / "assets", dirs_exist_ok=True)
    home = (ui / "index.html").read_text(encoding="utf-8")
    home = home.replace('href="./', 'href="./ui/').replace('src="./app.js"', 'src="./ui/app.js"')
    home = home.replace('href="../design-system.css"', 'href="./design-system.css"').replace('href="../competition.css"', 'href="./competition.css"')
    (target / "index.html").write_text(home, encoding="utf-8")
    editor = target / "editor"
    editor.mkdir(exist_ok=True)
    for source in (frontend / "editor").iterdir():
        if source.is_file() and source.suffix in (".js", ".mjs", ".css", ".json", ".md"):
            shutil.copyfile(source, editor / source.name)
    for directory in ("jsmol", "ketcher-c60", "ketcher-icons"):
        shutil.copytree(frontend / "editor" / "vendor" / directory, editor / "vendor" / directory, dirs_exist_ok=True)
    # The pinned native editor was built for /editor/. Adapt only its asset
    # resolution in this generated copy; never assume a domain-root deployment.
    ketcher=editor/'vendor/ketcher-c60'
    entry=ketcher/'index.html'
    entry.write_text(entry.read_text(encoding='utf-8').replace('/editor/vendor/ketcher-c60/','./'),encoding='utf-8')
    for asset in (ketcher/'assets').glob('*.js'):
        text=asset.read_text(encoding='utf-8')
        if '/editor/vendor/ketcher-c60/' not in text:
            continue
        text=text.replace('return`/editor/vendor/ketcher-c60/`+e','return new URL(`../`+e,import.meta.url).href')
        text=re.sub(r'new Worker\(`/editor/vendor/ketcher-c60/assets/([^`]+)`',r'new Worker(new URL(`./\1`,import.meta.url)',text)
        if '/editor/vendor/ketcher-c60/' in text:
            raise RuntimeError(f'Unrecognized Ketcher asset base: {asset.name}')
        asset.write_text(text,encoding='utf-8')
    shutil.copytree(frontend / "editor" / "examples", editor / "examples", dirs_exist_ok=True)
    workshop = (frontend / "editor" / "index.html").read_text(encoding="utf-8")
    workshop = workshop.replace('<html lang="zh-CN">', '<html lang="zh-CN" data-engine="wasm">')
    workshop = workshop.replace('href="/" aria-label="Racing Molecule 首页"', 'href="../#home" data-app-route="home" aria-label="Racing Molecule 首页"')
    workshop = workshop.replace('href="/" title="返回赛道"', 'href="../#track" data-app-route="track" title="返回赛道"')
    (editor / "index.html").write_text(workshop, encoding="utf-8")
    (target / ".nojekyll").write_text("")
    shutil.copyfile(ROOT / "LICENSE", target / "LICENSE")
    shutil.copyfile(ROOT / "RELEASING.md", target / "SOURCE.md")
    notices = target / "licenses"
    notices.mkdir(exist_ok=True)
    dependencies = [(ROOT / ".local" / "wasm" / "lammps" / "LICENSE", "LAMMPS-GPL-2.0.txt"),
                    (ROOT / ".local" / "wasm" / "rdkit" / "license.txt", "RDKIT-BSD.txt"),
                    (ROOT / ".local" / "wasm" / "boost" / "LICENSE_1_0.txt", "BOOST.txt"),
                    (ROOT / ".local" / "wasm" / "emsdk" / "upstream" / "emscripten" / "LICENSE", "EMSCRIPTEN.txt"),
                    (ROOT / ".local" / "wasm" / "better-enums" / "LICENSE.md", "BETTER-ENUMS.txt")]
    for source, filename in dependencies:
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copyfile(source, notices / filename)
    zlib_source = ROOT / ".local" / "wasm" / "emsdk" / "upstream" / "emscripten" / "cache" / "ports" / "zlib"
    for license_file in zlib_source.rglob("LICENSE"):
        shutil.copyfile(license_file, notices / "ZLIB.txt")
    (target / "THIRD_PARTY.md").write_text((ROOT / "THIRD_PARTY.md").read_text(encoding="utf-8"), encoding="utf-8")
    tool_html = (editor / "vendor/ketcher-c60/index.html").read_text(encoding="utf-8")
    paths = sorted(set(re.findall(r'(?:src|href)="(?:\./)?(assets/[^"?]+\.(?:js|css))"', tool_html)))
    paths = sorted(set(paths) | {"assets/"+p.name for p in (editor/"vendor/ketcher-c60/assets").glob("*.js")})
    tools = {"2d": [{"path": "vendor/ketcher-c60/"+p, "bytes": (editor/"vendor/ketcher-c60"/p).stat().st_size} for p in paths],
             "3d": [{"path": p, "bytes": (editor/p).stat().st_size} for p in ["vendor/jsmol/JSmol.min.js", "vendor/jsmol/j2s/core/corejmol.z.js"]]}
    (editor/"tool-assets.json").write_text(json.dumps(tools), encoding="utf-8")
    # Independently verified compressed chunks survive interrupted visits.
    for metadata_name in ("model.json", "fourier-atlas.json"):
        metadata_path = target / "browser" / "assets" / metadata_name
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        descriptors = (metadata["surface"]["arrays"].values() if metadata_name == "model.json" else metadata["arrays"])
        for descriptor in descriptors:
            packed = (target / "browser" / descriptor["path"]).read_bytes()
            if len(packed) < 1024 * 1024:
                continue
            descriptor["chunks"] = []
            for offset in range(0, len(packed), 1024 * 1024):
                data = packed[offset:offset + 1024 * 1024]
                digest = hashlib.sha256(data).hexdigest()
                relative = "assets/chunks/" + digest + ".bin"
                chunk = target / "browser" / relative
                chunk.parent.mkdir(parents=True, exist_ok=True)
                chunk.write_bytes(data)
                descriptor["chunks"].append({"path": relative, "bytes": len(data), "download_bytes": len(data), "sha256": digest})
        if metadata_name == "fourier-atlas.json":
            metadata["base_model_sha256"] = hashlib.sha256((target / "browser/assets/model.json").read_bytes()).hexdigest()
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    subprocess.run(["node", str(ROOT / "wasm/tools/compress_tables.cjs")], check=True)
    atlas_path = target / "browser/assets/fourier-atlas.json"
    atlas = json.loads(atlas_path.read_text())
    atlas["base_model_sha256"] = hashlib.sha256((target / "browser/assets/model.json").read_bytes()).hexdigest()
    atlas_path.write_text(json.dumps(atlas, indent=2) + "\n", encoding="utf-8")
    legacy_proxy = target / "browser/assets/fourier.json"
    if legacy_proxy.exists():
        metadata = json.loads(legacy_proxy.read_text())
        metadata["base_model_sha256"] = hashlib.sha256((target / "browser/assets/model.json").read_bytes()).hexdigest()
        legacy_proxy.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    runtime_bytes = (target / "browser/runtime/molecular-car.wasm").read_bytes()
    (target / "browser/runtime/download.json").write_text(json.dumps({"path": "runtime/molecular-car.wasm", "bytes": len(runtime_bytes), "download_bytes": len(runtime_bytes), "sha256": hashlib.sha256(runtime_bytes).hexdigest()}), encoding="utf-8")
    manifest = {"schema": "molecular-car-static-site-v1", "execution": "browser-wasm",
        "backend_required": False, "files": []}
    for path in sorted(target.rglob("*")):
        if path.is_file() and path.name != "site-manifest.json":
            manifest["files"].append({"path": path.relative_to(target).as_posix(), "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    (target / "site-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if args.no_archive:
        print(json.dumps({'site':str(target),'archive':None,'files':len(manifest['files'])}))
        return
    archive = ROOT / "outputs" / "wasm-preparation" / "racing-molecule-web-preview.zip"
    archive.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as output:
        for path in sorted(target.rglob("*")):
            if path.is_file():output.write(path,path.relative_to(target))
    print(json.dumps({"site": str(target), "archive": str(archive),
                      "total_bytes": sum(f["bytes"] for f in manifest["files"])}, indent=2))


if __name__ == "__main__":
    main()
