"""Verify race masses and the new defaults on a newly owned local service."""
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

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from molecular_car.engine import Simulation, DEFAULT_SETTINGS, DEFAULT_DRIVE, TIMESTEP_FS

OUT = ROOT / "outputs/md-release-validation/runtime"
PROTECTED_PORTS = set(range(8765, 8771)) | {8777, 8778}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8780)
    parser.add_argument("--keep-server", action="store_true")
    args = parser.parse_args()
    require(1 <= args.port <= 65535 and args.port not in PROTECTED_PORTS, "Use an isolated test port")
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", args.port))
    out = OUT / f"port-{args.port}"
    out.mkdir(parents=True, exist_ok=True)
    report = {"passed": False, "created_utc": datetime.now(timezone.utc).isoformat(),
              "defaults": DEFAULT_SETTINGS, "drive": DEFAULT_DRIVE, "checks": [],
              "base_url": f"http://127.0.0.1:{args.port}"}

    def record(name, **details):
        report["checks"].append({"name": name, "passed": True, **details})
        print(f"PASS {name}", flush=True)

    process = None
    opener = build_opener(ProxyHandler({}))

    def request(path, data=None, expected=200, raw=False):
        req = Request(report["base_url"] + path,
                      data=None if data is None else json.dumps(data).encode("utf-8"),
                      headers={"Content-Type": "application/json"})
        try:
            response = opener.open(req, timeout=20)
        except HTTPError as error:
            response = error
        with response:
            body = response.read()
            require(response.status == expected, f"HTTP {response.status}: {body[:500]!r}")
            return body.decode("utf-8") if raw else json.loads(body)

    def wait(predicate, timeout=240, allow_error=False):
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            try:
                state = request("/api/state")
            except (URLError, ConnectionError) as error:
                if getattr(getattr(error, "reason", None), "winerror", None) == 10013:
                    raise RuntimeError("Execution environment blocked the local HTTP connection; rerun with local-network permission.") from error
                time.sleep(.05)
                continue
            require(allow_error or state["status"] != "error", state.get("error"))
            if predicate(state):
                return state
            time.sleep(.05)
        raise AssertionError("Isolated service timed out")

    def ready():
        return wait(lambda s: s["status"] == "paused" and s["frame"] is not None)

    try:
        content = (ROOT / "examples/benzene.mol").read_text(encoding="utf-8")
        mol = Chem.MolFromMolBlock(content, removeHs=False)
        mol.GetAtomWithIdx(0).SetIsotope(13)
        mol.GetAtomWithIdx(6).SetIsotope(2)
        labelled = Chem.MolToMolBlock(mol)
        simulations = []
        for name, text in (("standard", content), ("labelled", labelled)):
            sim = Simulation(text, f"{name}.mol")
            try:
                types = np.ctypeslib.as_array(sim.lmp.gather_atoms("type", 0, 1), shape=(sim.n,)).copy()
                mass_table = sim.lmp.extract_atom("mass")
                integrated = np.array([mass_table[int(t)] for t in types])
                require(np.array_equal(integrated, sim.masses), "LAMMPS and COM masses differ")
                velocities = np.ctypeslib.as_array(sim.lmp.gather_atoms("v", 1, 3), shape=(3*sim.n,)).copy()
                simulations.append({"name": name, "mass": integrated.tolist(),
                                    "positions": sim.coordinates(), "velocities": velocities})
            finally:
                sim.close()
        np.testing.assert_array_equal(simulations[0]["mass"], simulations[1]["mass"])
        # MOL export changes coordinate precision; allow that input round trip.
        np.testing.assert_allclose(simulations[0]["positions"], simulations[1]["positions"], atol=2e-6, rtol=0)
        np.testing.assert_array_equal(simulations[0]["velocities"], simulations[1]["velocities"])
        record("real LAMMPS standard masses match preparation and mixed-isotope input", mass_u=sum(simulations[0]["mass"]))

        with (out / "server.stdout.log").open("w", encoding="utf-8") as stdout, \
                (out / "server.stderr.log").open("w", encoding="utf-8") as stderr:
            process = subprocess.Popen([sys.executable, "-E", "-B", "-m", "molecular_car.server", "--port", str(args.port)],
                cwd=ROOT, stdout=stdout, stderr=stderr,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        report["server_pid"] = process.pid
        initial = ready()
        require(initial["settings"] == DEFAULT_SETTINGS, "Service defaults differ")
        require(initial["control"]["drive"] == DEFAULT_DRIVE, "Drive default differs")
        require(initial["frame"]["step"] == 0 and len(initial["molecule"]["atoms"]) == 52, "Initial race is not paused at start")
        require(initial["runtime"]["timestep_fs"] == TIMESTEP_FS == .5, "Game time step changed")
        require(initial["runtime"]["surface"]["normal_policy"] == "bonded-saip-game-v1", "Game does not use topology normals")
        record("fresh service defaults to periodic Au / 600 K / 0.75 and unchanged 0.5 fs", runtime=initial["runtime"])
        html = request("/", raw=True)
        require('value="periodic-au" selected' in html and 'value="600"' in html, "Frontend defaults differ")
        record("frontend initial controls agree with server defaults")

        for mode in ("gradient", "tip"):
            request("/api/control", {"drive": {"mode": mode, "strength": .75}})
            ready_state = wait(lambda s: s["status"] == "paused" and s["control"]["drive"]["mode"] == mode)
            before = ready_state["frame"]["step"]
            request("/api/action", {"action": "start"})
            advanced = wait(lambda s: s["frame"] is not None and s["frame"]["step"] >= before+500)
            request("/api/action", {"action": "pause"})
            paused = ready()
            require(all(math.isfinite(x) for x in paused["frame"]["positions"]), "Nonfinite MD coordinates")
            replay = request("/api/replay")
            require(replay["surface"]["normal_policy"] == "bonded-saip-game-v1", "Replay lost normal policy")
            record(f"real {mode} MD advances and pauses with model metadata", step=paused["frame"]["step"])

        request("/api/load", {"content": labelled, "filename": "mixed-isotopes.mol"})
        loaded = ready()
        require(any("同位素标签" in note for note in loaded["molecule"]["notes"]), "Missing isotope policy notice")
        require(np.isclose(sum(a["mass"] for a in loaded["molecule"]["atoms"]), 78.114), "Isotope mass survived race preparation")
        record("API import reports ignored isotope labels and standard race mass")

        long_file = ROOT / "outputs/md-release-validation/inputs/folded-eicosane.mol"
        if not long_file.is_file():
            candidates = sorted((ROOT / "outputs/md-release-validation/inputs").glob("*eicosane*.mol"))
            require(bool(candidates), "Long-molecule regression fixture missing")
            long_file = candidates[0]
        request("/api/load", {"content": long_file.read_text(encoding="utf-8"), "filename": long_file.name})
        long_state = ready()
        record("long four-coordinate CH molecule prepares through game API", atom_count=len(long_state["molecule"]["atoms"]),
               normal_metadata=long_state["runtime"]["surface"])

        request("/api/control", {"drive": {"mode": "gradient", "strength": .75}})
        request("/api/load", {"example": "tetraphenylbenzene.mol"})
        final = ready()
        require(final["settings"] == DEFAULT_SETTINGS and final["frame"]["step"] == 0, "Final preview defaults changed")
        require(final["control"]["drive"] == DEFAULT_DRIVE, "Final drive differs")
        report["preview"] = {"generation": final["generation"], "settings": final["settings"],
                             "control": final["control"], "status": final["status"], "step": final["frame"]["step"]}
        record("preview restored to default 52-atom car and paused at start")
        report["passed"] = True
    except Exception as error:
        report.update(error=str(error), traceback=traceback.format_exc())
        traceback.print_exc()
    finally:
        report["source_sha256"] = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in ("molecular_car/engine.py", "molecular_car/molecule.py", "molecular_car/periodic_au.py",
                         "molecular_car/periodic_au_kernel.py", "molecular_car/fourier_au.py",
                         "molecular_car/fourier_au_kernel.py", "molecular_car/data/au111_reference.json",
                         "frontend/app.js", "frontend/index.html")}
        if process is not None and (not args.keep_server or not report["passed"]):
            try:
                request("/api/shutdown", {"owned_validation": True})
                process.wait(timeout=15)
            except Exception:
                process.terminate()
                process.wait(timeout=15)
        report["server_kept"] = bool(process is not None and args.keep_server and report["passed"])
        (out / "report.json").write_text(json.dumps(report, ensure_ascii=False, allow_nan=False, indent=2)+"\n", encoding="utf-8")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
