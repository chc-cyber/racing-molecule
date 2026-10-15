"""Generate an isolated native MD reference; never touches running game services."""
from __future__ import annotations
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from molecular_car.engine import Simulation


def phase_space(sim):
    result = {}
    for field, key in (("x", "positions"), ("v", "velocities")):
        values = np.ctypeslib.as_array(sim.lmp.gather_atoms(field, 1, 3), shape=(sim.n * 3,))
        result[key] = values.copy().tolist()
    return result


def main():
    source = ROOT / "examples" / "tetraphenylbenzene.mol"
    started = time.perf_counter()
    sim = Simulation(source.read_text(), source.name)
    try:
        report = {"schema": "molecular-car-native-browser-reference-v1",
                  "runtime": sim.runtime, "prepare_wall_seconds": time.perf_counter() - started}
        started = time.perf_counter()
        for _ in range(100):
            frame = sim.advance(100)
        report["propagation"] = {"steps": frame["step"], "time_fs": frame["time_fs"],
                                 "wall_seconds": time.perf_counter() - started,
                                 "center": frame["center"], "energy": frame["energy"]}
        sim.lmp.command("unfix bath")
        sim.thermostat = False
        sim.drive_enabled = False
        sim.lmp.command("reset_timestep 0")
        sim.lmp.command("run 0")
        sim._check_callback()
        fixture = phase_space(sim)
        fixture["frames"] = [sim.frame()]
        started = time.perf_counter()
        for _ in range(100):
            fixture["frames"].append(sim.advance(100))
        fixture["wall_seconds"] = time.perf_counter() - started
        report["nve"] = fixture
    finally:
        sim.close()
    target = ROOT / "outputs" / "wasm-preparation" / "native-reference.json"
    target.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"file": str(target), "prepare_wall_seconds": report["prepare_wall_seconds"],
                      "propagation": report["propagation"], "nve_wall_seconds": report["nve"]["wall_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
