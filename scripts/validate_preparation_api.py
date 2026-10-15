"""Validate force-gated preparation and mixed-element races on an owned port."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import socket
import subprocess
import sys
import time
import traceback
from urllib.error import HTTPError, URLError
from urllib.request import ProxyHandler, Request, build_opener

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from molecular_car.relaxation import PROTOCOL_VERSION
OUT = ROOT/"outputs/md-preparation-elements/api"
PROTECTED = set(range(8765, 8771)) | {8777, 8778, 8781}
SOURCES = ["molecular_car/engine.py", "molecular_car/relaxation.py", "molecular_car/preparation_bfgs.py", "molecular_car/molecule.py",
           "molecular_car/element_models.py", "molecular_car/periodic_au.py", "molecular_car/periodic_au_kernel.py",
           "molecular_car/fourier_au.py", "molecular_car/fourier_au_kernel.py", "molecular_car/hetero_au.py",
           "molecular_car/hetero_au_kernel.py", "molecular_car/data/hetero_au_v1.json",
           "molecular_car/data/au111_reference.json", "molecular_car/server.py", "frontend/app.js", "frontend/index.html"]


def hashes():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8782)
    parser.add_argument("--keep-server", action="store_true")
    args = parser.parse_args()
    require(args.port not in PROTECTED and 1 <= args.port <= 65535, "Use a new isolated port.")
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", args.port))
    out = OUT/f"port-{args.port}"; out.mkdir(parents=True, exist_ok=True)
    before = hashes()
    report = {"passed": False, "created_utc": datetime.now(timezone.utc).isoformat(),
              "base_url": f"http://127.0.0.1:{args.port}", "source_sha256": before, "checks": []}
    opener = build_opener(ProxyHandler({}))
    process = None

    def request(path, data=None, raw=False):
        req = Request(report["base_url"]+path, data=None if data is None else json.dumps(data).encode(),
                      headers={"Content-Type": "application/json"})
        with opener.open(req, timeout=20) as response:
            body = response.read()
            return body.decode() if raw else json.loads(body)

    def wait(predicate, *, error_ok=False):
        deadline = time.monotonic()+90
        while time.monotonic() < deadline:
            try:
                state = request("/api/state")
            except (URLError, ConnectionError) as error:
                if getattr(getattr(error, "reason", None), "winerror", None) == 10013:
                    raise RuntimeError("Local HTTP is blocked; rerun this owned-service verifier with local-network permission.") from error
                time.sleep(.05)
                continue
            require(error_ok or state["status"] != "error", state.get("error"))
            if predicate(state):
                return state
            time.sleep(.05)
        raise AssertionError("Owned service did not reach the requested state.")

    def ready():
        return wait(lambda s: s["status"] == "paused" and s.get("frame") is not None)

    def gate(state):
        p = state["runtime"]["preparation"]
        require(p["protocol_version"] == PROTOCOL_VERSION, "Wrong preparation version.")
        require(p["accepted"] and p["force_tolerance_met"] and p["final_coordinate_force_verified"], "Preparation gate failed.")
        require(p["global_force_norm_kcal_mol_A"] <= 1e-4, "Residual force exceeds the release gate.")
        require(state["frame"]["step"] == 0 and state["runtime"]["timestep_fs"] == .5, "Race time step/start changed.")
        return p

    def record(name, **details):
        report["checks"].append(dict(name=name, passed=True, **details))
        print("PASS "+name, flush=True)

    try:
        with (out/"server.stdout.log").open("w", encoding="utf-8") as stdout, (out/"server.stderr.log").open("w", encoding="utf-8") as stderr:
            process = subprocess.Popen([sys.executable, "-E", "-B", "-m", "molecular_car.server", "--port", str(args.port)],
                cwd=ROOT, stdout=stdout, stderr=stderr,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        report["server_pid"] = process.pid
        initial = ready(); p = gate(initial)
        require(initial["settings"]["surface"] == "periodic-au" and initial["settings"]["temperature"] == 600.
                and initial["control"]["drive"]["strength"] == .75, "Defaults changed.")
        require(len(initial["molecule"]["atoms"]) == 52, "Default fixture changed.")
        record("default 52-atom car strictly converges before race velocities", preparation=p)
        html = request("/", raw=True); js = request("/app.js", raw=True)
        require('value="600"' in html and 'value="periodic-au" selected' in html and "尚未标定真实吸附" in js,
                "Frontend default/element description missing.")
        record("frontend shows current defaults and provisional interface qualification")
        fixtures = ROOT/"outputs/md-preparation-elements/elements/inputs"
        for surface, filename in (("periodic-au", "pyridine.mol"), ("fourier-au", "trimethylphosphine-oxide.mol"),
                                  ("continuum", "iodobenzene.mol")):
            if ready()["settings"]["surface"] != surface:
                request("/api/control", {"surface": surface}); ready()
            request("/api/load", {"content": (fixtures/filename).read_text(encoding="utf-8"), "filename": filename})
            loaded = ready(); prep = gate(loaded)
            if surface != "continuum":
                require(loaded["runtime"]["surface"]["hetero_au"]["model"] == "hetero-au-uff-vdw-game-v1", "Mixed surface provenance missing.")
            request("/api/action", {"action": "start"})
            wait(lambda s: s["frame"] is not None and s["frame"]["step"] >= 1000)
            request("/api/action", {"action": "pause"}); paused = ready()
            require(all(math.isfinite(v) for v in paused["frame"]["positions"]), "Nonfinite race positions.")
            replay = request("/api/replay")
            require(replay["runtime"]["preparation"]["accepted"] and replay["runtime"]["elements"]["used_elements"], "Replay lost preparation/element metadata.")
            record("mixed-element import, strict preparation, live MD and replay: "+surface,
                   filename=filename, preparation=prep, integrated_step=paused["frame"]["step"])
        request("/api/control", {"surface": "periodic-au"}); ready()
        for filename, explanation in (("methyl-sulfur-pentafluoride.mol", "UFF"), ("trimethylphosphine.mol", "三配位磷"), ("iodoxybenzene.mol", "UFF")):
            request("/api/load", {"content": (fixtures/filename).read_text(encoding="utf-8"), "filename": filename})
            failed = wait(lambda s: s["status"] == "error", error_ok=True)
            require(explanation in failed["error"] and failed.get("frame") is None, "Invalid environment was not rejected before racing.")
            record("unsupported/pathological environment blocked before racing", filename=filename, error=failed["error"])
        request("/api/load", {"example": "tetraphenylbenzene.mol"}); final = ready(); gate(final)
        require(final["control"]["drive"] == {"mode": "gradient", "strength": .75}, "Final drive changed.")
        record("failed preparation recovers to converged default car paused at start")
        report["preview"] = {"status": final["status"], "step": final["frame"]["step"],
                             "settings": final["settings"], "control": final["control"],
                             "preparation": final["runtime"]["preparation"]}
        require(hashes() == before, "Production source changed during API validation.")
        report["passed"] = True
    except Exception as error:
        report.update(error=str(error), traceback=traceback.format_exc()); traceback.print_exc()
    finally:
        if process is not None and (not args.keep_server or not report["passed"]):
            try:
                request("/api/shutdown", {"owned_validation": True}); process.wait(timeout=15)
            except Exception:
                process.terminate(); process.wait(timeout=15)
        report["server_kept"] = bool(process is not None and args.keep_server and report["passed"])
        (out/"report.json").write_text(json.dumps(report, ensure_ascii=False, allow_nan=False, indent=2)+"\n", encoding="utf-8")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
