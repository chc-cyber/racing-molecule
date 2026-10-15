"""Native SAIP accuracy, periodic safety and real-MD tests of cached Au.

All tests are isolated from the running game. Static errors compare the
same complete molecular coordinates to native baseline SAIP, not fitted targets.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import traceback

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "outputs" / "periodic-au-check"

import numpy as np
from molecular_car.surface_racing import benzene_input
from molecular_car.explicit_au import FrozenAuOracle
from molecular_car.periodic_au import create_periodic_surface
from molecular_car.periodic_au_racing import PeriodicAuRace
from molecular_car.molecule import prepare_molecule


def write_json(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def static_checks(*, use_cache=True):
    prepared, geometry, _ = benzene_input()
    spacing = 4.08 / np.sqrt(2)
    cell = np.array([spacing, np.sqrt(3) * spacing])
    anchor = np.array([5 * spacing, cell[1] / 3])
    points = []
    for tilt in (0., 30., 60., 85.):
        angle = np.deg2rad(tilt)
        rotation = np.array([[1., 0., 0.], [0., np.cos(angle), -np.sin(angle)],
                             [0., np.sin(angle), np.cos(angle)]])
        shape = geometry @ rotation.T
        for low_height in (2.8, 3.3, 5., 9.):
            for phase in (np.array([0., 0.]), .37 * cell):
                center_z = low_height - shape[:, 2].min()
                coords = shape + [*(anchor + phase), center_z]
                points.append((f"tilt{tilt:g}-minz{low_height:g}-phase{phase[0]:.3g}", coords,
                               {"tilt_deg": tilt, "min_z_A": low_height}))
    distorted = geometry.copy()
    distorted[0, 2] += .12
    distorted[2, 2] -= .08
    points.append(("nonplanar-carbon-chain-forces", distorted + [*anchor, 3.4],
                   {"tilt_deg": 0., "min_z_A": 3.32}))
    points.append(("corner-seam", geometry + [0., 0., 3.4],
                   {"tilt_deg": 0., "min_z_A": 3.4}))
    surface = create_periodic_surface(prepared.atoms, points[0][1], use_cache=use_cache)
    rows = []
    elements = [atom["element"] for atom in prepared.atoms]
    with FrozenAuOracle(points[0][1], elements, style="base", requested_xy=(40., 39.9)) as oracle:
        for name, coordinates, extra in points:
            reference = oracle.evaluate(coordinates)
            energy, forces = surface.energy_and_forces_eV(coordinates)
            row = {"case": name, **extra,
                   "atoms": len(coordinates),
                   "native_energy_eV": reference["energy_eV"], "fast_energy_eV": energy,
                   "energy_error_eV": energy - reference["energy_eV"],
                   "max_force_error_eV_A": float(np.max(np.abs(forces - reference["forces_eV_A"]))),
                   "native_max_abs_force_eV_A": float(np.max(np.abs(reference["forces_eV_A"])))}
            rows.append(row)
    # A larger, genuinely connected fixture spans several primitive cells.
    path = ROOT / "examples" / "tetraphenylbenzene.mol"
    larger = prepare_molecule(path.read_text(encoding="utf-8"), path.name)
    coords = larger.positions.copy()
    coords[:, :2] -= np.mean(coords[:, :2], axis=0)
    coords[:, :2] += anchor
    coords[:, 2] += 3.1 - coords[:, 2].min()
    larger_surface = create_periodic_surface(larger.atoms, coords, use_cache=use_cache)
    with FrozenAuOracle(coords, [a["element"] for a in larger.atoms], style="base",
                        requested_xy=(40., 39.9)) as oracle:
        reference = oracle.evaluate(coords)
        energy, forces = larger_surface.energy_and_forces_eV(coords)
        rows.append({"case": "connected-52-atom-multi-primitive-cell",
                     "atoms": len(coords),
                     "tilt_deg": None, "min_z_A": float(coords[:, 2].min()),
                     "native_energy_eV": reference["energy_eV"], "fast_energy_eV": energy,
                     "energy_error_eV": energy - reference["energy_eV"],
                     "max_force_error_eV_A": float(np.max(np.abs(forces - reference["forces_eV_A"]))),
                     "native_max_abs_force_eV_A": float(np.max(np.abs(reference["forces_eV_A"])))})
    translations = []
    e0, f0 = larger_surface.energy_and_forces_eV(coords)
    for multiples in ([1, 0], [-7, 11], [120, -50]):
        shifted = coords + [*(cell * multiples), 0.]
        energy, forces = larger_surface.energy_and_forces_eV(shifted)
        translations.append({"primitive_cell_translation": multiples,
                             "energy_error_eV": abs(energy - e0),
                             "force_error_eV_A": float(np.max(np.abs(forces - f0)))})
    if max(x["energy_error_eV"] for x in translations) > 1e-8 or max(x["force_error_eV_A"] for x in translations) > 1e-8:
        raise AssertionError("Translating the complete molecule over multiple surface cells changed its surface interaction.")
    # This single connected graph is longer than two nominal small Au boxes.
    # Native Au expands its domain; the environmental field never replicates it.
    from scripts.verify_au_optimization import long_polyphenyl
    long_positions, long_elements = long_polyphenyl()
    long_surface = create_periodic_surface(long_elements, long_positions, use_cache=use_cache)
    with FrozenAuOracle(long_positions, long_elements, style="base",
                        requested_xy=(40., 39.9)) as oracle:
        reference = oracle.evaluate(long_positions)
        energy, forces = long_surface.energy_and_forces_eV(long_positions)
        rows.append({"case": "connected-212-atom-spans-two-small-supercells",
                     "atoms": len(long_positions), "tilt_deg": 0.,
                     "min_z_A": float(long_positions[:, 2].min()),
                     "native_energy_eV": reference["energy_eV"], "fast_energy_eV": energy,
                     "energy_error_eV": energy - reference["energy_eV"],
                     "max_force_error_eV_A": float(np.max(np.abs(forces - reference["forces_eV_A"]))),
                     "native_max_abs_force_eV_A": float(np.max(np.abs(reference["forces_eV_A"])))})
        native_long_cell = reference["cell"]
    long_e0, long_f0 = long_surface.energy_and_forces_eV(long_positions)
    moved_energy, moved_force = long_surface.energy_and_forces_eV(long_positions + [*(cell * [-80, 123]), 0.])
    long_translation = {"primitive_cells": [-80, 123], "energy_error_eV": abs(moved_energy - long_e0),
                        "force_error_eV_A": float(np.max(np.abs(moved_force - long_f0)))}
    if long_translation["energy_error_eV"] > 1e-8 or long_translation["force_error_eV_A"] > 1e-8:
        raise AssertionError("The long complete molecule changed on a multiple-cell translation.")
    # Full Cartesian derivative includes the forces distributed to normal neighbours.
    coordinates = distorted + [*anchor, 3.4]
    energy, forces = surface.energy_and_forces_eV(coordinates)
    epsilon = 1e-5
    finite_difference = np.zeros_like(forces)
    for atom in range(len(coordinates)):
        for axis in range(3):
            plus, minus = coordinates.copy(), coordinates.copy()
            plus[atom, axis] += epsilon
            minus[atom, axis] -= epsilon
            ep = surface.energy_and_forces_eV(plus)[0]
            em = surface.energy_and_forces_eV(minus)[0]
            finite_difference[atom, axis] = -(ep - em) / (2 * epsilon)
    fd_error = float(np.max(np.abs(finite_difference - forces)))
    if fd_error > 5e-6:
        raise AssertionError(f"Conservative force check failed: {fd_error} eV/A")
    max_energy_error = max(abs(row["energy_error_eV"]) for row in rows)
    max_energy_error_per_atom = max(abs(row["energy_error_eV"]) / row["atoms"] for row in rows)
    max_force_error = max(row["max_force_error_eV_A"] for row in rows)
    accepted = max_energy_error_per_atom <= 1e-4 and max_force_error <= 2e-3
    prefix = "static" if use_cache else "direct-static"
    with (OUT / f"{prefix}-accuracy.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = {"cases": rows, "max_energy_error_eV": max_energy_error,
              "passed": accepted,
              "max_energy_error_per_atom_eV": max_energy_error_per_atom,
              "max_force_error_eV_A": max_force_error, "finite_difference_force_error_eV_A": fd_error,
              "translations": translations, "large_molecule_span_A": np.ptp(coords, axis=0).tolist(),
              "surface": surface.metadata(), "larger_surface": larger_surface.metadata(),
              "long_molecule": {"atoms": len(long_positions), "span_A": np.ptp(long_positions, axis=0).tolist(),
                                "native_autoexpanded_cell": native_long_cell, "translation": long_translation},
              "acceptance": {"energy_eV_per_atom": 1e-4, "force_eV_A": 2e-3, "force_derivative_eV_A": 5e-6}}
    write_json(OUT / f"{prefix}-report.json", result)
    if not accepted:
        worst = max(rows, key=lambda row: row["max_force_error_eV_A"])
        raise AssertionError(f"Native comparison error exceeds acceptance: E={max_energy_error}, F={max_force_error}, worst={worst['case']}")
    print(f"PASS static {len(rows)} cases; max E error {max_energy_error:.3g} eV, F error {max_force_error:.3g} eV/A", flush=True)
    return result


def run_md_checks(steps=2000):
    trials = []
    for seed in (98245165, 2718281, 1234567):
        race = PeriodicAuRace(strength=.75, temperature=600, seed=seed)
        try:
            samples = [race.sample()]
            total_wall, total_cpu = 0., 0.
            warmup = race.advance(1000)
            samples = [race.sample()]
            batches = []
            for batch in range(3):
                timing = race.advance(steps)
                total_wall += timing["wall_seconds"]
                total_cpu += timing["cpu_seconds"]
                batches.append({**timing, "batch": batch + 1, "steps": steps,
                                "seconds_per_step": timing["wall_seconds"] / steps})
                samples.append(race.sample())
            actual_steps = 3 * steps
            description = race.description()
            trial = {"description": description, "steps": actual_steps,
                     "wall_seconds": total_wall, "cpu_seconds": total_cpu,
                     "warmup_steps": 1000, "warmup_timing": warmup, "batches": batches,
                     "median_wall_seconds_per_step": float(np.median([b["seconds_per_step"] for b in batches])),
                     "microseconds_per_step": total_wall / actual_steps * 1e6,
                     "ps_per_wall_second": actual_steps * race.timestep_ps / total_wall,
                     "samples": samples, "monitor": race.monitor}
            trials.append(trial)
            write_json(OUT / f"md-{seed}.json", trial)
            print(f"PASS MD seed {seed}: {actual_steps} steps, {trial['microseconds_per_step']:.2f} us/step", flush=True)
        finally:
            race.close()
    race = PeriodicAuRace(strength=0, temperature=10, seed=98245165, thermostat=False)
    try:
        energies = [race.sample()["total_energy_eV"]]
        samples = [race.sample()]
        for _ in range(20):
            race.advance(100)
            samples.append(race.sample())
            energies.append(samples[-1]["total_energy_eV"])
        drift = float(np.max(np.abs(np.array(energies) - energies[0])))
        if drift > .005:
            raise AssertionError(f"Cached surface NVE energy drift {drift} eV exceeds .005 eV")
        nve = {"steps": 2000, "duration_ps": 1., "thermostat": False, "strength": 0.,
               "maximum_energy_drift_eV": drift, "acceptance_eV": .005,
               "description": race.description(), "samples": samples}
        write_json(OUT / "nve.json", nve)
        print(f"PASS NVE maximum drift {drift:.3g} eV", flush=True)
    finally:
        race.close()
    result = {"trials": trials, "nve": nve,
              "protocol": "600 K, strength 0.75, 0.5 fs, 350 fs Langevin; three seeds; propagation timing excludes sampling and table preparation"}
    write_json(OUT / "md-report.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("static", "md", "all"), default="all")
    parser.add_argument("--steps", type=int, default=2000)
    parser.add_argument("--direct", action="store_true", help="Check the exact lattice sum without interpolation (static only)")
    args = parser.parse_args()
    if args.steps < 100 or args.steps % 100:
        parser.error("--steps must be a positive multiple of 100")
    if args.direct and args.stage != "static":
        parser.error("--direct is supported with --stage static only")
    OUT.mkdir(parents=True, exist_ok=True)
    result = {"status": "running", "stage": args.stage}
    try:
        if args.stage in ("static", "all"):
            result["static"] = static_checks(use_cache=not args.direct)
        if args.stage in ("md", "all"):
            result["md"] = run_md_checks(args.steps)
        result["status"] = "passed"
        result["source_sha256"] = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                   for name in ("molecular_car/periodic_au.py", "molecular_car/periodic_au_kernel.py", "molecular_car/explicit_au.py",
                                                "molecular_car/periodic_au_racing.py", "scripts/verify_periodic_au.py",
                                                "molecular_car/surface_racing.py", "molecular_car/molecule.py",
                                                "molecular_car/fields.py", "molecular_car/au_surface.py",
                                                "scripts/lammps_runtime.py", "examples/benzene.mol",
                                                "outputs/au-surface-check/au111-relaxed.xyz")}
    except BaseException as error:
        result.update(status="failed", error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        suffix = f"{args.stage}-direct" if args.direct else args.stage
        write_json(OUT / f"report-{suffix}.json", result)


if __name__ == "__main__":
    main()
