"""Smoke-check a new isolated local server, never an existing user scene.

This is a short API/MD integration test, not a performance benchmark.
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
OUT = ROOT / "outputs/surface-mode-check"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keep-server", action="store_true", help="Keep this isolated server for UI QA")
    parser.add_argument("--port", type=int, default=0, help="New test port; zero selects an empty port")
    args = parser.parse_args()
    if not 0 <= args.port <= 65535:
        parser.error("port must be zero or in 1..65535")
    with socket.socket() as probe:
        if os.name == "nt":
            probe.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
        try:
            probe.bind(("127.0.0.1", args.port))
        except OSError as error:
            parser.error(f"Test port is occupied; no existing service was changed: {error}")
        port = probe.getsockname()[1]
    base = f"http://127.0.0.1:{port}"
    OUT.mkdir(parents=True, exist_ok=True)
    opener = build_opener(ProxyHandler({}))
    report = {"passed": False, "created_utc": datetime.now(timezone.utc).isoformat(),
              "base_url": base, "checks": [], "purpose": "isolated API smoke, not performance benchmark"}

    def request(path, data=None, expected=200, raw=False):
        body = None if data is None else json.dumps(data).encode("utf-8")
        req = Request(base+path, data=body, headers={"Content-Type": "application/json"})
        try:
            response = opener.open(req, timeout=20)
        except HTTPError as error:
            response = error
        with response:
            payload = response.read()
            require(response.status == expected, f"HTTP {response.status}, expected {expected}: {payload[:400]!r}")
            return payload.decode("utf-8") if raw else json.loads(payload)

    def wait(predicate, timeout=120, allow_error=False):
        deadline = time.monotonic()+timeout
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

    def record(name, **details):
        report["checks"].append({"name": name, "passed": True, **details})
        print(f"PASS {name}", flush=True)

    stdout = (OUT/"server.stdout.log").open("w", encoding="utf-8")
    stderr = (OUT/"server.stderr.log").open("w", encoding="utf-8")
    process = subprocess.Popen([sys.executable, "-m", "molecular_car.server", "--port", str(port),
        "--molecule", str(ROOT/"examples/benzene.mol")], cwd=ROOT, stdout=stdout, stderr=stderr,
        creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
    report["server_pid"] = process.pid
    (OUT/"api-report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    try:
        initial = wait(lambda s: s["status"] == "paused" and s["frame"] is not None)
        require(initial["settings"]["surface"] == "continuum", "Default surface changed")
        require(initial["frame"]["step"] == 0, "Server auto-started race")
        record("continuum remains default and loads paused")
        generation = initial["generation"]
        for payload in ({"surface": "mobile-au"}, {"surface": True},
                        {"surface": "periodic-au", "temperature": 600}):
            request("/api/control", payload, 400)
        require(request("/api/state")["generation"] == generation, "Invalid request reset race")
        record("invalid or mixed surface requests are rejected without reset")
        response = request("/api/control", {"surface": "periodic-au"})
        require(response["generation"] > generation, "Surface switch did not create new generation")
        require(response["frame"] is None and response["molecule"] is None, "Old molecule or frames retained")
        periodic = wait(lambda s: s["status"] == "paused" and s["frame"] is not None)
        require(periodic["settings"]["surface"] == "periodic-au", "Wrong requested surface")
        require(periodic["runtime"]["lammps_atoms"] == len(periodic["molecule"]["atoms"]) == 12, "Au atoms entered MD worker")
        require(periodic["runtime"]["mobile_surface_atoms"] == 0, "Surface was not frozen")
        require(periodic["frame"]["step"] == 0 and periodic["frame"]["surface_mode"] == "periodic-au", "Fresh periodic state invalid")
        metadata = periodic["runtime"]["surface"]
        require(metadata["mode"] == "periodic-au" and len(metadata) > 2, "Missing actual surface metadata")
        require(metadata["algorithm_version"] and len(metadata["parameter_sha256"]) == 64,
                "Surface version or parameter hash missing")
        require(metadata["units"] == {"energy": "kcal/mol", "force": "kcal/mol/angstrom", "length": "angstrom"},
                "Surface units incompatible with game real units")
        require(metadata["cache_requested"] and metadata["cache_loaded"] and metadata["cache_grid_shape"],
                "Game did not load the validated shared cache")
        require(len(metadata["layer_heights_angstrom"]) == 6 and math.isclose(metadata["layer_heights_angstrom"][-1], 0),
                "Surface centre-plane convention incorrect")
        report["surface_metadata"] = metadata
        record("periodic Au resets safely, creates only racing atoms, and remains paused", generation=periodic["generation"])
        request("/api/control", {"corrugation": .2}, 400)
        paused = request("/api/state")
        request("/api/control", {"surface": "periodic-au"})
        same = request("/api/state")
        require(same["generation"] == paused["generation"] and same["frame"] == paused["frame"], "Same mode reset race")
        record("fixed Au corrugation rejected; same mode preserves paused state")
        request("/api/action", {"action": "start"})
        running = wait(lambda s: s["frame"] is not None and s["frame"]["step"] >= 100)
        request("/api/action", {"action": "pause"})
        stopped = wait(lambda s: s["status"] == "paused" and s["frame"]["step"] >= running["frame"]["step"])
        require(all(math.isfinite(v) for v in stopped["frame"]["positions"]), "Nonfinite periodic MD coordinates")
        time.sleep(.15)
        stable = request("/api/state")
        require(stable["frame"] == stopped["frame"], "Paused MD advanced")
        replay = request("/api/replay")
        require(replay["surface"] == metadata and replay["runtime"]["surface"] == metadata,
                "Replay metadata differs from active physical model")
        require(replay["frames"][-1]["positions"] == stopped["frame"]["positions"], "Replay used different coordinates")
        record("short real MD, pause, and replay preserve finite physical coordinates", steps=stopped["frame"]["step"])
        old_generation = stopped["generation"]
        request("/api/control", {"surface": "continuum"})
        restored = wait(lambda s: s["status"] == "paused" and s["frame"] is not None)
        require(restored["generation"] > old_generation and restored["frame"]["step"] == 0,
                "Switching back did not reset race")
        require(restored["settings"]["surface"] == "continuum" and request("/api/replay")["surface"]["mode"] == "continuum",
                "Replay retained prior surface")
        record("return to continuum clears prior race and metadata")
        html = request("/", raw=True)
        require('id="surface-mode"' in html and 'value="periodic-au"' in html and 'id="space-filling-model"' in html,
                "Surface or representation controls missing")
        record("same-origin UI serves surface choice and existing representations")
        report["passed"] = True
    except BaseException as error:
        report["error"] = str(error); report["traceback"] = traceback.format_exc()
        traceback.print_exc()
    finally:
        if not args.keep_server or not report["passed"]:
            try:
                request("/api/shutdown", {})
                process.wait(timeout=15)
            except Exception:
                if process.poll() is None:
                    if os.name == "nt":
                        # This PID is the Popen child created above. Include its
                        # daemon MD worker if loopback itself is unavailable.
                        subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
                                       creationflags=subprocess.CREATE_NO_WINDOW)
                    else:
                        process.terminate()
                    process.wait(timeout=10)
        else:
            print(f"Isolated UI server: {base} (PID {process.pid})", flush=True)
        stdout.close(); stderr.close()
        (OUT/"api-report.json").write_text(json.dumps(report, ensure_ascii=False, allow_nan=False, indent=2)+"\n", encoding="utf-8")
        print(f"Report: {OUT/'api-report.json'}", flush=True)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
