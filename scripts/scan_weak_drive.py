"""Actual 20 nm game trajectories: weak fixed drive and the existing Au lattice.

Run with the game's .venv, without plotting packages in its import path.
Each checkpoint contains complete molecular coordinates; no rigid-body model,
force rescaling, larger timestep, defects or moving substrate are introduced.
"""
from __future__ import annotations

import os
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/weak-drive-20nm"
sys.path.insert(0, str(ROOT))
import numpy as np
from molecular_car.engine import Simulation, TIMESTEP_FS
from molecular_car.track import TRACK
from molecular_car.periodic_au import parameter_file

SEEDS = (98245165, 2718281, 1234567)
SOURCE_NAMES = (
    "molecular_car/engine.py", "molecular_car/track.py", "molecular_car/fields.py",
    "molecular_car/fields_kernel.py", "molecular_car/external_fields.py",
    "molecular_car/molecule.py", "molecular_car/periodic_au.py",
    "molecular_car/periodic_au_kernel.py", "scripts/lammps_runtime.py",
    "examples/tetraphenylbenzene.mol", "scripts/scan_weak_drive.py",
)


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def fingerprint():
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in SOURCE_NAMES}


def trajectory_metrics(time_ps, centers_nm, finish_x_nm=20., window_ps=1., threshold_nm=.02):
    """COM displacement categories in disjoint 1 ps windows, not atom vibration.

    Pause means small *longitudinal* progress, including sideways movement.
    Forward/backward use a symmetric threshold; the three fractions sum to one.
    A final partial window is omitted. Finish is the first sampled crossing.
    """
    t = np.asarray(time_ps, dtype=float)
    centers = np.asarray(centers_nm, dtype=float)
    relative = centers - centers[0]
    window_times = np.arange(int(np.floor((t[-1] + 1e-9) / window_ps)) + 1) * window_ps
    window_x = np.interp(window_times, t, centers[:, 0])
    dx = np.diff(window_x)
    finish_indices = np.flatnonzero(centers[:, 0] >= finish_x_nm)
    reached = bool(finish_indices.size)
    return {
        "time_ps": float(t[-1]), "displacement_nm": relative[-1].tolist(),
        "maximum_lateral_excursion_nm": float(np.max(np.abs(relative[:, 1]))),
        "time_rms_lateral_displacement_nm": float(np.sqrt(np.trapz(relative[:, 1]**2, t) / t[-1])),
        "sampled_lateral_path_nm": float(np.abs(np.diff(centers[:, 1])).sum()),
        "backward_fraction": float(np.mean(dx < -threshold_nm)) if len(dx) else None,
        "pause_fraction": float(np.mean(np.abs(dx) <= threshold_nm)) if len(dx) else None,
        "forward_fraction": float(np.mean(dx > threshold_nm)) if len(dx) else None,
        "coarse_window_count": len(dx), "coarse_window_ps": window_ps,
        "longitudinal_threshold_nm": threshold_nm,
        "finish_observed": reached,
        "finish_sample_time_ps": float(t[finish_indices[0]]) if reached else None,
        "right_censored": not reached,
    }


def trial_key(config):
    return (f"{config['phase']}_{config['surface']}_T{config['temperature_K']:g}"
            f"_S{config['strength']:g}_C{config['corrugation']:g}_seed{config['seed']}"
            f"_{config['duration_ps']:g}ps")


def run_trial(config, hashes):
    key = trial_key(config)
    folder = OUT / "trials"
    folder.mkdir(parents=True, exist_ok=True)
    checkpoint = folder / (key + ".json")
    coordinates_path = folder / (key + ".npz")
    if checkpoint.exists() and coordinates_path.exists():
        previous = json.loads(checkpoint.read_text(encoding="utf-8"))
        if previous.get("status") == "passed" and previous["source_sha256"] == hashes and previous["config"] == config:
            return previous
    simulation = None
    began = time.perf_counter()
    result = {"key": key, "config": config, "status": "running", "source_sha256": hashes}
    try:
        content = (ROOT / "examples/tetraphenylbenzene.mol").read_text(encoding="utf-8")
        simulation = Simulation(content, "tetraphenylbenzene.mol", {
            "surface": config["surface"], "temperature": config["temperature_K"],
            "corrugation": config["corrugation"], "damping_fs": 350.,
        }, random_seed=config["seed"], drive={"mode": "gradient", "strength": 0.})
        prepared_seconds = time.perf_counter() - began
        # Thermal preparation only: keep the same existing bath instance alive.
        # Constructor and subsequent run setups have identical call cadence.
        sample_steps = round(config["sample_ps"] * 1000 / TIMESTEP_FS)
        warmup_steps = round(config["warmup_ps"] * 1000 / TIMESTEP_FS)
        for offset in range(0, warmup_steps, sample_steps):
            simulation.advance(min(sample_steps, warmup_steps - offset))
        simulation.update({"drive": {"mode": "gradient", "strength": config["strength"]}})
        initial = simulation.frame()
        initial_step = initial["step"]
        samples, coordinates = [], []
        bond_ids = np.array([[b["a"] - 1, b["b"] - 1] for b in simulation.prepared.bonds])
        minimum_bond, maximum_bond = float("inf"), 0.
        minimum_z, maximum_z = float("inf"), -float("inf")

        def append(frame):
            nonlocal minimum_bond, maximum_bond, minimum_z, maximum_z
            positions = np.asarray(frame["positions"], dtype=float).reshape(-1, 3)
            lengths = np.linalg.norm(positions[bond_ids[:, 0]] - positions[bond_ids[:, 1]], axis=1)
            minimum_bond = min(minimum_bond, float(lengths.min()))
            maximum_bond = max(maximum_bond, float(lengths.max()))
            minimum_z = min(minimum_z, float(positions[:, 2].min()))
            maximum_z = max(maximum_z, float(positions[:, 2].max()))
            if not np.isfinite(positions).all() or not np.isfinite(frame["energy"]):
                raise AssertionError("Non-finite production state")
            if minimum_bond < .6 or maximum_bond > 2.5:
                raise AssertionError("Molecular bond geometry left the accepted diagnostic range")
            coordinates.append(positions.copy())
            samples.append({
                "time_ps": (frame["step"] - initial_step) * TIMESTEP_FS / 1000.,
                "absolute_step": frame["step"], "center_nm": (np.asarray(frame["center"]) / 10).tolist(),
                "energy_kcal_mol": frame["energy"], "internal_temperature_K": frame["internal_temperature"],
            })

        append(initial)
        production_steps = round(config["duration_ps"] * 1000 / TIMESTEP_FS)
        propagate_seconds = 0.
        actual_steps = 0
        for offset in range(0, production_steps, sample_steps):
            count = min(sample_steps, production_steps - offset)
            started = time.perf_counter()
            frame = simulation.advance(count)
            propagate_seconds += time.perf_counter() - started
            actual_steps += count
            append(frame)
            if frame["center"][0] >= TRACK["finish_x"]:
                break
        result.update(
            status="passed", runtime=simulation.runtime, atoms=simulation.n,
            mass_u=float(simulation.masses.sum()), preparation_seconds=prepared_seconds,
            total_wall_seconds=time.perf_counter() - began,
            production_wall_seconds=propagate_seconds, production_steps=actual_steps,
            production_us_per_step=propagate_seconds / actual_steps * 1e6,
            bond_length_range_angstrom=[minimum_bond, maximum_bond],
            atom_height_range_angstrom=[minimum_z, maximum_z],
            metrics=trajectory_metrics([s["time_ps"] for s in samples], [s["center_nm"] for s in samples]),
            samples=samples, coordinates_file=coordinates_path.relative_to(OUT).as_posix(),
        )
        np.savez_compressed(coordinates_path,
            time_ps=np.array([s["time_ps"] for s in samples]),
            positions_angstrom=np.stack(coordinates), masses_u=simulation.masses,
            bonds_zero_based=bond_ids)
        result["coordinates_sha256"] = hashlib.sha256(coordinates_path.read_bytes()).hexdigest()
    except BaseException as error:
        result.update(status="failed", error=str(error), traceback=traceback.format_exc())
    finally:
        if simulation is not None:
            simulation.close()
        write_json(checkpoint, result)
    return result


def summarize(rows):
    groups = {}
    for row in rows:
        if row["status"] != "passed":
            continue
        c = row["config"]
        group_key = (c["phase"], c["surface"], c["temperature_K"], c["strength"], c["corrugation"], c["duration_ps"])
        groups.setdefault(group_key, []).append(row)
    output = []
    for key, trials in sorted(groups.items()):
        phase, surface, temperature, strength, corrugation, duration = key
        dx = np.array([r["metrics"]["displacement_nm"][0] for r in trials])
        dy = np.array([r["metrics"]["displacement_nm"][1] for r in trials])
        output.append({
            "phase": phase, "surface": surface, "temperature_K": temperature,
            "strength": strength, "corrugation": corrugation, "duration_ps": duration,
            "n": len(trials), "mean_forward_nm": float(dx.mean()), "forward_range_nm": [float(dx.min()), float(dx.max())],
            "final_lateral_rms_nm": float(np.sqrt(np.mean(dy**2))),
            "mean_backward_fraction": float(np.mean([r["metrics"]["backward_fraction"] for r in trials])),
            "mean_pause_fraction": float(np.mean([r["metrics"]["pause_fraction"] for r in trials])),
            "maximum_lateral_excursion_nm": max(r["metrics"]["maximum_lateral_excursion_nm"] for r in trials),
            "finishers": sum(r["metrics"]["finish_observed"] for r in trials),
            "median_production_us_per_step": float(np.median([r["production_us_per_step"] for r in trials])),
            "trial_keys": [r["key"] for r in trials],
        })
    return output


def configurations(args):
    conditions = []
    if args.phase == "screen":
        conditions = [(surface, t, s, .08) for t in (150, 300, 600) for s in (.05, .10, .15, .25)
                      for surface in ("periodic-au", "continuum")]
        conditions += [(surface, 300, s, .08) for s in (0., .75) for surface in ("periodic-au", "continuum")]
    else:
        selected = json.loads(Path(args.selected).read_text(encoding="utf-8"))
        conditions = [(surface, c["temperature_K"], c["strength"], c.get("corrugation", .08))
                      for c in selected for surface in ("periodic-au", "continuum")]
    return [{"phase": args.phase, "surface": surface, "temperature_K": t, "strength": s,
             "corrugation": corr, "seed": seed, "duration_ps": args.duration_ps,
             "warmup_ps": 5., "sample_ps": .1}
            for surface, t, s, corr in conditions for seed in SEEDS]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("screen", "extended"), default="screen")
    parser.add_argument("--duration-ps", type=float, default=50.)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--selected", help="JSON list of temperature_K/strength conditions for an extended phase")
    args = parser.parse_args()
    if not 1 <= args.duration_ps <= 1000 or not 1 <= args.workers <= 3:
        parser.error("Use 1..1000 ps and 1..3 independent single-thread workers")
    if args.phase == "extended" and not args.selected:
        parser.error("Extended phase requires --selected")
    if TRACK["finish_x"] != 200. or TRACK["half_width"] != 100.:
        raise ValueError("This protocol requires the 20 nm square track")
    OUT.mkdir(parents=True, exist_ok=True)
    hashes = fingerprint()
    tasks = configurations(args)
    report = {
        "status": "running", "created_utc": datetime.now(timezone.utc).isoformat(), "phase": args.phase,
        "source_sha256": hashes, "parameter_file_sha256": hashlib.sha256(parameter_file().read_bytes()).hexdigest(),
        "track_angstrom": TRACK, "workers": args.workers, "planned_trials": len(tasks),
        "protocol": {"molecule": "52-atom C30H22 test car", "timestep_fs": TIMESTEP_FS,
            "damping_fs": 350., "warmup_ps": 5., "sample_ps": .1, "coarse_window_ps": 1.,
            "longitudinal_pause_threshold_nm": .02, "seeds": SEEDS,
            "surface": "Full existing frozen periodic Au SAIP; continuum corrugation .08 is a separate control",
            "drive": "Unchanged fixed broad gradient; only strength varied",
            "cost": "Concurrent-process end-to-end propagation and frame cost; not an isolated core benchmark",
            "qualification": "Finite trajectories of one molecule and three seeds; no converged diffusion or inferred finish times"},
        "trials": [], "summary": [],
    }
    report_path = OUT / (args.phase + "-report.json")
    write_json(report_path, report)
    began = time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_trial, config, hashes): config for config in tasks}
        for future in as_completed(futures):
            row = future.result()
            report["trials"].append(row)
            report["summary"] = summarize(report["trials"])
            report["elapsed_wall_seconds"] = time.perf_counter() - began
            write_json(report_path, report)
            metrics = row.get("metrics", {})
            print(f"{len(report['trials'])}/{len(tasks)} {row['key']} {row['status']} "
                  f"dx={metrics.get('displacement_nm', [None])[0]} "
                  f"back={metrics.get('backward_fraction')} error={row.get('error', '')}", flush=True)
    if fingerprint() != hashes:
        report["status"] = "failed"
        report["error"] = "Simulation source changed while the scan was running"
    else:
        report["status"] = "passed" if all(r["status"] == "passed" for r in report["trials"]) else "failed"
    report["trials"].sort(key=lambda r: r["key"])
    write_json(report_path, report)
    summary_path = OUT / (args.phase + "-summary.csv")
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        columns = [k for k in report["summary"][0] if k not in ("trial_keys", "forward_range_nm")]
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(report["summary"])
    print(f"{report['status']}: {report_path}", flush=True)
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
