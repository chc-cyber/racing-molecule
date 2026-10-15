"""Real LAMMPS acceptance of strict joint UFF/surface race preparation.

Compares the former 350-iteration, energy-or-force FIRE protocol with the new
force-gated protocol on identical input files. Legacy preparation is substituted
only inside this process for measurement. No service or production file changes.
"""
from __future__ import annotations
import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[key] = "1"
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback
import zipfile
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from molecular_car import engine
from molecular_car.relaxation import (relax_joint_geometry, JointRelaxationError,
                                      GLOBAL_FORCE_TOLERANCE, PROTOCOL_VERSION)
from scripts.lammps_runtime import open_lammps
FIXTURES = ("benzene-12", "test-car-52", "fullerene-c60", "ethynyl-rigid", "folded-eicosane")
SURFACES = ("periodic-au", "fourier-au", "continuum")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    paths = sorted((ROOT / "molecular_car").glob("*.py"))
    paths += sorted(path for path in (ROOT / "molecular_car/data").rglob("*") if path.is_file())
    paths += [Path(__file__).resolve(), ROOT / "scripts/lammps_runtime.py"]
    return {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path) for path in paths}


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def gather(lmp, name, count):
    raw = lmp.gather_atoms(name, 1, 3)
    return np.ctypeslib.as_array(raw, shape=(count * 3,)).reshape(count, 3).copy()


def force_statistics(force, energy):
    norm = float(np.linalg.norm(force))
    return {"global_force_norm_kcal_mol_A": norm,
            "maximum_atom_force_kcal_mol_A": float(np.linalg.norm(force, axis=1).max()),
            "rms_atom_force_kcal_mol_A": norm / len(force)**.5,
            "maximum_force_component_kcal_mol_A": float(np.abs(force).max()),
            "potential_energy_kcal_mol": float(energy)}


def former_joint_preparation(lmp, count, callback, check):
    """Exact former native minimize settings, with diagnostic force refresh."""
    start = time.perf_counter()
    lmp.commands_list(["min_style fire", "min_modify dmax 0.1", "minimize 1.0e-6 1.0e-4 350 3500"])
    check()
    iterations = int(lmp.get_thermo("step"))
    lmp.command("run 0"); check()
    values = force_statistics(gather(lmp, "f", count), lmp.get_thermo("pe"))
    return {"protocol_version": "former-fire-350-energy-or-force-diagnostic-only",
            "energy_tolerance": 1e-6, "global_force_tolerance_kcal_mol_A": 1e-4,
            "iteration_limit": 350, "evaluation_limit": 3500, "iterations": iterations,
            "wall_seconds": time.perf_counter() - start,
            "force_tolerance_met": values["global_force_norm_kcal_mol_A"] <= GLOBAL_FORCE_TOLERANCE,
            "qualification": "Former protocol may return due to energy tolerance or budget without meeting force tolerance; measured as the comparison, not accepted by current production.",
            **values}


def conservative_evaluation(simulation, coordinates):
    flat = coordinates.ravel().tolist()
    energy = float(simulation.prepared.forcefield.CalcEnergy(flat))
    force = -np.asarray(simulation.prepared.forcefield.CalcGrad(flat)).reshape(-1, 3)
    field_energy, field_force = simulation.field_evaluator.energy_and_forces(
        coordinates, 1, 0., 0., 0., 0., simulation.settings["corrugation"])
    energy += field_energy; force += field_force
    if simulation.surface is not None:
        surface_energy, surface_force = simulation.surface.energy_and_forces(coordinates)
        energy += surface_energy; force += surface_force
    return energy, force


def one_preparation(path, surface, legacy, output):
    start = time.perf_counter()
    # Both constructors use the real game bath/velocity/drive sequence. The
    # independent conservative force at their unchanged initial coordinates
    # excludes the Langevin/drive force installed only after geometry checks.
    protocol = former_joint_preparation if legacy else relax_joint_geometry
    with patch("molecular_car.engine.relax_joint_geometry", protocol):
        simulation = engine.Simulation(path.read_text(encoding="utf-8"), path.name,
                                       {"surface": surface, "temperature": 600., "damping_fs": 350.},
                                       thermostat=True, drive={"mode": "gradient", "strength": .75})
    try:
        coordinates = simulation.coordinates()
        velocities = gather(simulation.lmp, "v", simulation.n)
        energy, force = conservative_evaluation(simulation, coordinates)
        independent = force_statistics(force, energy)
        prep = simulation.runtime["preparation"]
        types = simulation.lmp.extract_atom("type")
        masses = simulation.lmp.extract_atom("mass")
        integrated_masses = np.array([float(masses[int(types[i])]) for i in range(simulation.n)])
        mass_error = float(np.abs(integrated_masses - simulation.masses).max())
        npz = output / f"{path.stem}-{surface}-{'former' if legacy else 'strict'}-initial.npz"
        np.savez_compressed(npz, prepared_positions_A=coordinates, velocities_A_fs=velocities,
                            masses_amu=simulation.masses,
                            bonds_1_based=np.array([(bond["a"], bond["b"]) for bond in simulation.prepared.bonds]))
        row = {"preparation": prep, "independent_conservative_force": independent,
               "independent_force_matches_gate": abs(independent["global_force_norm_kcal_mol_A"] - prep["global_force_norm_kcal_mol_A"]) < 1e-12,
               "independent_energy_matches_gate": abs(energy - prep["potential_energy_kcal_mol"]) < 1e-10,
               "timestep_fs": float(simulation.lmp.extract_global("dt")),
               "integrated_mass_max_error_amu": mass_error,
               "initial_temperature_K": float(simulation.lmp.get_thermo("temp")),
               "initial_COM_A": simulation.center().tolist(),
               "atom_coordinate_minimum_A": coordinates.min(axis=0).tolist(),
               "atom_coordinate_maximum_A": coordinates.max(axis=0).tolist(),
               "initial_state_file": str(npz), "initial_state_sha256": sha256(npz),
               "runtime": simulation.runtime}
        simulation.lmp.command("run 200 post no"); simulation._check_callback()
        row["smoke_integrated_steps"] = int(simulation.lmp.get_thermo("step"))
        row["smoke_positions_finite"] = bool(np.isfinite(simulation.coordinates()).all())
        row["timer_remaining_after_cleanup_seconds"] = float(simulation.lmp.get_thermo("timeremain"))
        row["wall_seconds_including_construction_and_0_1_ps_smoke"] = time.perf_counter() - start
        row["checks_passed"] = bool(row["independent_force_matches_gate"] and row["independent_energy_matches_gate"]
                                    and mass_error < 1e-12 and row["timestep_fs"] == .5
                                    and row["smoke_integrated_steps"] == 200 and row["smoke_positions_finite"]
                                    and row["timer_remaining_after_cleanup_seconds"] == 0.)
        if not legacy:
            row["checks_passed"] &= prep["accepted"] and independent["global_force_norm_kcal_mol_A"] <= GLOBAL_FORCE_TOLERANCE
        return row
    finally:
        simulation.close()


def native_failure_checks():
    """Genuine native time/evaluation soft stops and pathological admission."""
    rows = []
    for name, scale, options in (("wall-expired-before-first-stage", 1., {"wall_time_limit": 1e-9}),
                                 ("evaluation-soft-stop", 1., {"evaluation_limit": 2}),
                                 ("pathological-initial-force", 1e6, {})):
        with open_lammps() as (lmp, runtime):
            lmp.commands_list(["units real", "atom_style atomic", "boundary f f f", "atom_modify map array sort 0 0.0",
                               "region box block -20 20 -20 20 -20 20 units box", "create_box 1 box", "mass 1 12.011",
                               "pair_style zero 10", "pair_coeff * *", "fix molecular all external pf/callback 1 1",
                               "thermo_style custom step pe ke etotal", "thermo 1000000", "timestep 0.5"])
            lmp.create_atoms(2, [1, 2], [1, 1], [4., 0., 3., 0., 2., 5.])
            def force(caller, step, nlocal, tags, positions, output):
                output[:] = -scale * positions
                lmp.fix_external_set_energy_global("molecular", .5 * scale * float(np.sum(positions**2)))
            lmp.set_fix_external_callback("molecular", force)
            try:
                relax_joint_geometry(lmp, 2, force, lambda: None, **options)
                raise AssertionError("Numerical failure fixture unexpectedly accepted.")
            except JointRelaxationError as error:
                diagnostics = error.diagnostics
            coordinates = gather(lmp, "x", 2)
            actual_force = gather(lmp, "f", 2)
            force_error = float(np.abs(-scale * coordinates - actual_force).max())
            previous_step = int(lmp.get_thermo("step"))
            # Never integrate the rejected million-kcal spring: the numerical
            # fixture's same callback switches to its safe spring for cleanup.
            scale = 1.
            lmp.commands_list(["velocity all set 0 0 0", "fix integrate all nve", "run 5 post no"])
            steps = int(lmp.get_thermo("step")) - previous_step
            timer = float(lmp.get_thermo("timeremain"))
            rows.append({"case": name, "rejected": True, "diagnostics": diagnostics,
                         "independent_final_force_error_kcal_mol_A": force_error,
                         "post_failure_safe_surrogate_integrated_steps": steps, "timer_remaining_seconds": timer,
                         "passed": force_error < 1e-9 and steps == 5 and timer == 0.})
    return {"passed": all(row["passed"] for row in rows), "cases": rows,
            "qualification": "Numerical harmonic fixtures check native cleanup; invalid molecular structures are not integrated."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/md-preparation-elements/preparation")
    parser.add_argument("--fixtures-dir", type=Path, default=ROOT / "outputs/md-release-validation/inputs")
    parser.add_argument("--surfaces", choices=SURFACES, nargs="+", default=list(SURFACES))
    parser.add_argument("--fixtures", choices=FIXTURES, nargs="+", default=list(FIXTURES))
    args = parser.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    before = source_hashes()
    start = time.perf_counter()
    report = {"status": "running", "created_utc": datetime.now(timezone.utc).isoformat(),
              "protocol_version": PROTOCOL_VERSION, "source_sha256_before": before,
              "game_timestep_fs": engine.TIMESTEP_FS,
              "race_model_settings": {"temperature_K": 600., "damping_fs": 350., "gradient_strength": .75},
              "cases": [], "qualification": "Finite deterministic local geometry convergence, not proof of a global minimum, equilibrium ensembles or experimental surface/field calibration."}
    try:
        snapshot = out / "source-snapshot.zip"
        with zipfile.ZipFile(snapshot, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for name, expected in before.items():
                data = (ROOT / name).read_bytes()
                if hashlib.sha256(data).hexdigest() != expected:
                    raise RuntimeError("Source changed while preparing the source snapshot.")
                archive.writestr(name, data)
        report.update(source_snapshot_file=str(snapshot), source_snapshot_sha256=sha256(snapshot))
        for surface in args.surfaces:
            for index, name in enumerate(args.fixtures):
                path = args.fixtures_dir / (name + ".mol")
                row = {"fixture": name, "surface": surface, "input_sha256": sha256(path)}
                # Alternate order so neither protocol systematically gets all
                # cold constructions. Native geometry timings exclude cache
                # creation/preparation, as indicated in the plot/report.
                for legacy in ((True, False) if index % 2 == 0 else (False, True)):
                    row["former" if legacy else "strict"] = one_preparation(path, surface, legacy, out)
                row["passed"] = row["former"]["checks_passed"] and row["strict"]["checks_passed"]
                report["cases"].append(row)
                write_json(out / "convergence-report.json", report)
                print(f"{surface} {name}: old F={row['former']['preparation']['global_force_norm_kcal_mol_A']:.6g}; strict F={row['strict']['preparation']['global_force_norm_kcal_mol_A']:.6g}; passed={row['passed']}", flush=True)
        report["native_failure_checks"] = native_failure_checks()
        after = source_hashes()
        report["source_sha256_after"] = after
        report["sources_unchanged_during_run"] = before == after
        report["status"] = "passed" if (all(row["passed"] for row in report["cases"])
                                            and report["native_failure_checks"]["passed"] and before == after) else "failed"
    except BaseException as error:
        report.update(status="failed", error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        report["wall_seconds"] = time.perf_counter() - start
        write_json(out / "convergence-report.json", report)
    if report["status"] != "passed":
        raise SystemExit("Preparation convergence acceptance failed; inspect convergence-report.json")


if __name__ == "__main__":
    main()
