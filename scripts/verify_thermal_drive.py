"""Compare hotter baths and weaker drive using actual LAMMPS trajectories.

These short, three-seed samples describe trajectory spread. They do not fit
a converged diffusion coefficient, validate gold chemistry, or predict a
100 nm finish time for the weaker driver. Run without a web server.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from molecular_car.engine import Simulation, TIMESTEP_FS
from scripts.verify_gradient import free_com_rms

SEEDS = (98245165, 2718281, 1234567)


def trial(temperature, strength, seed, steps):
    content = (ROOT / "examples/tetraphenylbenzene.mol").read_text(encoding="utf-8")
    simulation = Simulation(content, "tetraphenylbenzene.mol",
                            {"temperature": temperature}, random_seed=seed,
                            drive={"mode": "gradient", "strength": strength})
    started = time.perf_counter()
    try:
        initial = simulation.frame()
        samples = []
        lo, hi = float("inf"), 0.0
        for offset in range(0, steps, 2000):
            frame = simulation.advance(min(2000, steps-offset))
            positions = np.array(frame["positions"]).reshape(-1, 3)
            lengths = [np.linalg.norm(positions[b["a"]-1]-positions[b["b"]-1])
                       for b in simulation.prepared.bonds]
            lo, hi = min(lo, min(lengths)), max(hi, max(lengths))
            if lo < .6 or hi > 2.5 or not np.isfinite(frame["energy"]):
                raise AssertionError("Nonfinite state or unstable molecular bond geometry")
            samples.append({key: frame[key] for key in
                            ("step", "time_fs", "center", "internal_temperature")})
        displacement = np.array(frame["center"])-initial["center"]
        delta_y = [row["center"][1]-initial["center"][1] for row in samples]
        return {"temperature_K": temperature, "strength": strength, "seed": seed,
                "runtime": simulation.runtime, "steps": steps,
                "time_ps": frame["time_fs"]/1000,
                "displacement_nm": (displacement/10).tolist(),
                "mean_forward_speed_m_s": displacement[0]/(frame["time_fs"])*1e5,
                "maximum_sampled_lateral_excursion_nm": max(map(abs, delta_y))/10,
                "sampled_height_range_angstrom": [float(positions[:, 2].min()), float(positions[:, 2].max())],
                "bond_length_range_angstrom": [lo, hi],
                "free_com_reference": free_com_rms(temperature, frame["time_fs"],
                                                    simulation.settings["damping_fs"], simulation.masses.sum()),
                "samples": samples, "wall_seconds": time.perf_counter()-started}
    finally:
        simulation.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duration-ps", type=float, default=50)
    args = parser.parse_args()
    if not 1 <= args.duration_ps <= 100:
        parser.error("duration must be 1..100 ps")
    steps = round(args.duration_ps*1000/TIMESTEP_FS)
    output = ROOT / "outputs/thermal-drive-check/report.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    report = {"passed": False, "created_utc": datetime.now(timezone.utc).isoformat(),
              "qualification": "Short three-seed trajectory spread, not converged diffusion or experimental Au validation",
              "groups": []}
    try:
        for temperature, strength in ((300, .75), (600, .75), (1000, .75), (600, .25)):
            rows = []
            for seed in SEEDS:
                row = trial(temperature, strength, seed, steps)
                rows.append(row)
                print(f"T={temperature}K S={strength} seed={seed}: "
                      f"dx={row['displacement_nm'][0]:.3f}nm dy={row['displacement_nm'][1]:.3f}nm "
                      f"v={row['mean_forward_speed_m_s']:.1f}m/s", flush=True)
            dx = np.array([row["displacement_nm"][0] for row in rows])
            dy = np.array([row["displacement_nm"][1] for row in rows])
            report["groups"].append({"temperature_K": temperature, "strength": strength,
                                     "duration_ps": args.duration_ps,
                                     "mean_forward_displacement_nm": float(dx.mean()),
                                     "rms_final_lateral_displacement_nm": float(np.sqrt(np.mean(dy*dy))),
                                     "final_lateral_range_nm": [float(dy.min()), float(dy.max())],
                                     "largest_sampled_lateral_excursion_nm": max(r["maximum_sampled_lateral_excursion_nm"] for r in rows),
                                     "trials": rows})
            output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        report["passed"] = True
    except Exception as error:
        report["error"] = str(error)
        raise
    finally:
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        print(f"Report: {output}", flush=True)


if __name__ == "__main__":
    main()
