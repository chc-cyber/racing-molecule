"""20 nm square game acceptance on a new isolated local service.

Never touches an existing race. Timing checks belong in the physics benchmark.
"""
from __future__ import annotations

import os
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import socket
import subprocess
import sys
import time
import traceback
from urllib.error import HTTPError, URLError
from urllib.request import ProxyHandler, Request, build_opener

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/20nm-game-check"
PROTECTED_PORTS = {8765, 8766, 8767, 8768, 8769, 8777}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keep-server", action="store_true")
    parser.add_argument("--port", type=int, default=8770)
    args = parser.parse_args()
    if not 1 <= args.port <= 65535 or args.port in PROTECTED_PORTS:
        parser.error("Use a new test port, excluding the existing 8765–8769 and 8777 scenes.")
    with socket.socket() as probe:
        if os.name == "nt":
            probe.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
        try:
            probe.bind(("127.0.0.1", args.port))
        except OSError as error:
            parser.error(f"Test port occupied; no existing service changed: {error}")
    base = f"http://127.0.0.1:{args.port}"
    OUT.mkdir(parents=True, exist_ok=True)
    opener = build_opener(ProxyHandler({}))
    report = {"passed": False, "created_utc": datetime.now(timezone.utc).isoformat(),
              "base_url": base, "checks": [], "purpose": "isolated API/MD acceptance, not performance benchmark"}

    def request(path, data=None, expected=200, raw=False):
        body = None if data is None else json.dumps(data).encode("utf-8")
        req = Request(base + path, data=body, headers={"Content-Type": "application/json"})
        try:
            response = opener.open(req, timeout=20)
        except HTTPError as error:
            response = error
        with response:
            payload = response.read()
            require(response.status == expected, f"HTTP {response.status}, expected {expected}: {payload[:400]!r}")
            return payload.decode("utf-8") if raw else json.loads(payload)

    def wait(predicate, timeout=240, allow_error=False):
        deadline = time.monotonic() + timeout
        latest = None
        while time.monotonic() < deadline:
            try:
                latest = request("/api/state")
            except (URLError, ConnectionError):
                time.sleep(.05)
                continue
            if not allow_error and latest["status"] == "error":
                raise AssertionError(latest["error"])
            if predicate(latest):
                return latest
            time.sleep(.05)
        raise AssertionError(f"Timeout waiting for isolated state: {latest}")

    def ready(mode):
        return wait(lambda s: s["status"] == "paused" and s["frame"] is not None
                    and s["settings"]["surface"] == mode)

    def record(name, **details):
        report["checks"].append({"name": name, "passed": True, **details})
        print(f"PASS {name}", flush=True)

    stdout = (OUT / "server.stdout.log").open("w", encoding="utf-8")
    stderr = (OUT / "server.stderr.log").open("w", encoding="utf-8")
    process = subprocess.Popen([sys.executable, "-m", "molecular_car.server", "--port", str(args.port)],
                              cwd=ROOT, stdout=stdout, stderr=stderr,
                              creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
    report["server_pid"] = process.pid
    (OUT / "api-report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    try:
        initial = ready("continuum")
        require(len(initial["molecule"]["atoms"]) == 52 and initial["frame"]["step"] == 0,
                "Default 52-atom car changed or auto-started")
        record("original continuum and default 52-atom car remain paused")
        kernel = initial["runtime"]["external_field_kernel"]
        require(kernel["accelerator"] == "numba" and kernel["fastmath"] is False
                and initial["runtime"]["timestep_fs"] == .5,
                "Compiled external forces or unchanged timestep not recorded")
        report["external_field_kernel"] = kernel
        record("compiled external kernel retains 0.5 fs and the existing dynamics")
        expected_track = {"start_x": 0.0, "finish_x": 200.0, "half_width": 100.0, "surface_z": 0.0}
        require(initial["track"] == expected_track and request("/api/replay")["track"] == expected_track,
                "20 nm square track missing from state/replay")
        report["track"] = initial["track"]
        record("API and replay use a 200 Angstrom finish and 100 Angstrom half-width")
        original_tip = initial["control"]["tip"]
        request("/api/control", {"tip": {"x": 206, "y": 100}})
        wait(lambda s: s["control"]["tip"]["x"] == 206 and s["control"]["tip"]["y"] == 100)
        for tip in ({"x": 206.001}, {"x": -5.001}, {"y": 100.001}, {"y": -100.001}):
            request("/api/control", {"tip": tip}, 400)
        request("/api/control", {"tip": {"x": original_tip["x"], "y": original_tip["y"]}})
        limits_restored = wait(lambda s: s["control"]["tip"]["x"] == original_tip["x"] and s["control"]["tip"]["y"] == original_tip["y"])
        require(limits_restored["generation"] == initial["generation"] and limits_restored["frame"]["positions"] == initial["frame"]["positions"],
                "Paused needle limit checks changed the molecular coordinates")
        record("Needle API accepts new rectangle limits and rejects targets outside them")
        generation = initial["generation"]
        for payload in ({"surface": "mobile-au"}, {"surface": True},
                        {"surface": "fourier-au", "temperature": 600}):
            request("/api/control", payload, 400)
        require(request("/api/state")["generation"] == generation, "Invalid request reset race")
        record("invalid or mixed surface changes do not mutate race")

        request("/api/control", {"temperature":600})
        for drive in ("gradient", "tip"):
            request("/api/control", {"drive":{"mode":drive, "strength":.5}})
            wait(lambda s: s["status"] == "paused" and s["control"]["drive"]["mode"] == drive)
            before = request("/api/state")["frame"]["step"]
            request("/api/action", {"action":"start"})
            running = wait(lambda s: s["frame"] is not None and s["frame"]["step"] >= before+300)
            request("/api/action", {"action":"pause"})
            stopped = wait(lambda s: s["status"] == "paused" and s["frame"]["step"] >= running["frame"]["step"])
            replay = request("/api/replay")
            require(stopped["settings"]["surface"] == "continuum" and
                    all(math.isfinite(v) for v in stopped["frame"]["positions"]) and
                    replay["frames"][-1]["positions"] == stopped["frame"]["positions"] and
                    replay["runtime"]["external_field_kernel"]["accelerator"] == "numba",
                    "Continuous MD or replay changed")
            record(f"optimized continuum runs real {drive} MD and pauses with exact replay")
        request("/api/control", {"drive":{"mode":"gradient", "strength":.5}})
        wait(lambda s: s["control"]["drive"]["mode"] == "gradient")
        before = request("/api/state")["frame"]["step"]
        request("/api/action", {"action":"pulse"})
        pulse = wait(lambda s: s["status"] == "paused" and s["frame"]["step"] >= before+800)
        require(pulse["frame"]["pulse_active"] is False and pulse["control"]["drive"]["strength"] == .5,
                "Pulse failed to expire or altered the base drive")
        record("compiled forces retain the temporary pulse and its exact expiry")

        report["surface_metadata"] = {}
        for mode in ("periodic-au", "fourier-au"):
            previous = request("/api/state")
            changed = request("/api/control", {"surface": mode})
            require(changed["generation"] > previous["generation"] and changed["frame"] is None,
                    f"{mode} retained previous race generation or frame")
            require(changed["molecule"] is None and request("/api/frames?after=-1")["frames"] == [],
                    f"{mode} retained previous molecule or frame queue")
            active = ready(mode)
            require(active["frame"]["step"] == 0 and active["frame"]["surface_mode"] == mode,
                    f"{mode} fresh frame invalid")
            runtime, metadata = active["runtime"], active["runtime"]["surface"]
            require(runtime["lammps_atoms"] == runtime["integrated_atoms"] == len(active["molecule"]["atoms"]) == 52,
                    "Gold entered molecular MD worker")
            require(runtime["mobile_surface_atoms"] == 0 and metadata["mode"] == mode,
                    "Actual surface metadata incorrect")
            require(runtime["model"] == f"uff-{mode}-gradient-soft-tip-v1", "Wrong actual model identity")
            require(metadata["algorithm_version"] and len(metadata["parameter_sha256"]) == 64,
                    "Algorithm or parameter hash missing")
            require(metadata["units"] == {"energy": "kcal/mol", "force": "kcal/mol/angstrom", "length": "angstrom"},
                    "Surface units incompatible with game")
            require(len(metadata["layer_heights_angstrom"]) == 6 and math.isclose(metadata["layer_heights_angstrom"][-1], 0),
                    "Centre-plane convention incorrect")
            require(metadata["accelerator"] == "numba", "Game did not use compiled surface")
            require(runtime["external_field_kernel"]["accelerator"] == "numba"
                    and runtime["external_field_kernel"]["surface_enabled"] is False,
                    "Periodic mode did not use the shared compiled drive and wall kernel")
            if mode == "fourier-au":
                require(metadata["proxy_height_support_angstrom"] and metadata["fallback"], "Proxy domain not recorded")
                require(set(active["frame"]["surface_evaluation"]) == {"full_proxy", "blended_proxy", "exact_direction_fallback"},
                        "Current proxy/fallback counts not recorded in frame")
            report["surface_metadata"][mode] = metadata
            request("/api/control", {"corrugation": .2}, 400)
            paused = request("/api/state")
            request("/api/control", {"surface": mode})
            same = request("/api/state")
            require(same["generation"] == paused["generation"] and same["frame"] == paused["frame"],
                    "Same mode reset paused race")
            record(f"{mode} prepares with actual metadata and only 52 racing atoms", generation=active["generation"])
            for drive in ("tip", "gradient"):
                request("/api/control", {"drive": {"mode": drive, "strength": .5}})
                wait(lambda s: s["status"] == "paused" and s["control"]["drive"]["mode"] == drive)
                request("/api/action", {"action": "start"})
                before = request("/api/state")["frame"]["step"]
                running = wait(lambda s: s["frame"] is not None and s["frame"]["step"] >= before + 100)
                request("/api/action", {"action": "pause"})
                stopped = wait(lambda s: s["status"] == "paused" and s["frame"]["step"] >= running["frame"]["step"])
                require(all(math.isfinite(v) for v in stopped["frame"]["positions"]), "Nonfinite MD coordinates")
                if mode == "fourier-au":
                    counts = stopped["frame"]["surface_evaluation"]
                    require(sum(counts.values()) >= 30 and all(isinstance(v, int) and v >= 0 for v in counts.values()),
                            "Current proxy count invalid")
                replay = request("/api/replay")
                require(replay["surface"] == metadata and replay["runtime"]["surface"] == metadata,
                        "Replay surface differs from active physical model")
                require(replay["frames"][-1]["positions"] == stopped["frame"]["positions"], "Replay coordinates differ")
                record(f"{mode} runs real {drive} MD and pauses with exact replay", step=stopped["frame"]["step"])

        paused = request("/api/state")
        imported = request("/api/editor/import", {"content": (ROOT / "examples/benzene.mol").read_text(), "filename": "benzene.mol"})
        require(imported["project"]["atoms"], "Editor import route failed")
        unchanged = request("/api/state")
        require(unchanged["generation"] == paused["generation"] and unchanged["frame"] == paused["frame"],
                "Independent editor route changed MD")
        record("independent molecule editor API preserves paused race")

        from rdkit import Chem
        from rdkit.Chem import AllChem
        water = Chem.AddHs(Chem.MolFromSmiles("O"))
        require(AllChem.EmbedMolecule(water, randomSeed=2026) == 0, "Water fixture failed")
        request("/api/load", {"content": Chem.MolToMolBlock(water), "filename": "unsupported-water.mol"})
        rejected = wait(lambda s: s["status"] == "error", allow_error=True)
        require("含有 O" in rejected["error"] and rejected["settings"]["surface"] == "fourier-au"
                and rejected["runtime"] is None, "Unsupported structure silently changed physical surface")
        request("/api/control", {"surface": "continuum"})
        restored = ready("continuum")
        require(len(restored["molecule"]["atoms"]) == 3 and request("/api/replay")["surface"]["mode"] == "continuum",
                "Surface error recovery failed")
        request("/api/load", {"example": "tetraphenylbenzene.mol"})
        ready("continuum")
        request("/api/control", {"temperature":10, "drive":{"mode":"gradient", "strength":.75}})
        request("/api/action", {"action":"reset"})
        final = ready("continuum")
        require(len(final["molecule"]["atoms"]) == 52 and final["frame"]["step"] == 0,
                "Default car not restored for browser acceptance")
        record("unsupported element reports error without fallback; explicit continuum recovery works")
        html = request("/", raw=True)
        require(all(f'value="{mode}"' in html for mode in ("continuum", "periodic-au", "fourier-au"))
                and 'id="space-filling-model"' in html and 'id="track-size-label">20 nm × 20 nm' in html, "UI choices or 20 nm track copy missing")
        record("same-origin UI serves all three surface choices and both representations")
        report["passed"] = True
    except BaseException as error:
        report["error"] = str(error)
        report["traceback"] = traceback.format_exc()
        traceback.print_exc()
    finally:
        if not args.keep_server or not report["passed"]:
            try:
                request("/api/shutdown", {})
                process.wait(timeout=15)
            except Exception:
                if process.poll() is None:
                    if os.name == "nt":
                        subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
                                       creationflags=subprocess.CREATE_NO_WINDOW)
                    else:
                        process.terminate()
                    process.wait(timeout=10)
        else:
            print(f"Isolated 20 nm square game: {base} (PID {process.pid})", flush=True)
        stdout.close(); stderr.close()
        (OUT / "api-report.json").write_text(json.dumps(report, ensure_ascii=False, allow_nan=False, indent=2) + "\n", encoding="utf-8")
        print(f"Report: {OUT / 'api-report.json'}", flush=True)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

