"""Reproduce workshop fullerene preparation without touching a live service."""
from __future__ import annotations

import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[key] = "1"

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from molecular_car import engine
from molecular_car.relaxation import JointRelaxationError, relax_joint_geometry


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "tests/data/connected-fullerene.mol")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/joint-fullerene-diagnosis")
    parser.add_argument("--surface", default="periodic-au")
    parser.add_argument("--finite-difference", action="store_true")
    parser.add_argument("--stages", help="JSON stage list for a process-local experiment")
    parser.add_argument("--legacy", action="store_true", help="Disable the rescue to reproduce the earlier rejection")
    parser.add_argument("--smoke-steps", type=int, default=200)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    captured = {"input": str(args.input), "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
                "surface": args.surface, "rescue_disabled": args.legacy,
                "source_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                  for name in ("molecular_car/engine.py", "molecular_car/molecule.py",
                                               "molecular_car/relaxation.py", "molecular_car/preparation_bfgs.py")}}
    started = time.perf_counter()

    def measured(lmp, count, callback, check):
        lmp.command('log "' + (args.output / "native.log").as_posix() + '"')
        try:
            options = {"stages": json.loads(args.stages)} if args.stages else {}
            if args.legacy:
                options["bfgs_rescue"] = False
            report = relax_joint_geometry(lmp, count, callback, check, **options)
        except JointRelaxationError as exc:
            report = exc.diagnostics
            captured["error"] = str(exc)
        coordinates = np.ctypeslib.as_array(lmp.gather_atoms("x", 1, 3), shape=(count * 3,)).reshape(count, 3).copy()
        forces = np.ctypeslib.as_array(lmp.gather_atoms("f", 1, 3), shape=(count * 3,)).reshape(count, 3).copy()
        simulation = callback.__self__

        def evaluate(x):
            flat = x.ravel().tolist()
            energy = simulation.prepared.forcefield.CalcEnergy(flat)
            force = -np.asarray(simulation.prepared.forcefield.CalcGrad(flat)).reshape(count, 3)
            fe, ff = simulation.field_evaluator.energy_and_forces(x, 1, 0., 0., 0., 0., simulation.settings["corrugation"])
            energy += fe
            force += ff
            if simulation.surface is not None:
                se, sf = simulation.surface.energy_and_forces(x)
                energy += se
                force += sf
            return float(energy), force

        energy, force = evaluate(coordinates)
        captured.update(preparation=report, atom_count=count,
                        independent_energy_kcal_mol=energy,
                        independent_force_norm=float(np.linalg.norm(force)),
                        independent_force_max_error=float(np.max(np.abs(force - forces))),
                        net_force=force.sum(axis=0).tolist(),
                        torque_about_centroid=np.cross(coordinates - coordinates.mean(axis=0), force).sum(axis=0).tolist(),
                        final_bounds=[coordinates.min(axis=0).tolist(), coordinates.max(axis=0).tolist()])
        np.savez_compressed(args.output / "final-state.npz", positions=coordinates, forces=forces,
                            initial_positions=simulation.prepared.positions)
        if args.finite_difference:
            numerical = np.empty_like(force)
            delta = 1e-5
            for atom in range(count):
                for axis in range(3):
                    plus, minus = coordinates.copy(), coordinates.copy()
                    plus[atom, axis] += delta
                    minus[atom, axis] -= delta
                    numerical[atom, axis] = -(evaluate(plus)[0] - evaluate(minus)[0]) / (2 * delta)
            differences = np.abs(force - numerical)
            worst = np.unravel_index(np.argmax(differences), differences.shape)
            captured["finite_difference"] = dict(maximum_component_error=float(differences.max()),
                                                 worst_atom_1_based=int(worst[0] + 1),
                                                 worst_axis="xyz"[worst[1]], delta_angstrom=delta,
                                                 numerical_force_norm=float(np.linalg.norm(numerical)))
        if not report["accepted"]:
            raise JointRelaxationError(report)
        return report

    simulation = None
    try:
        with patch("molecular_car.engine.relax_joint_geometry", measured):
            simulation = engine.Simulation(args.input.read_text(encoding="utf-8"), args.input.name,
                                           {"surface": args.surface})
        initial = simulation.coordinates().copy()
        captured["initial_race_step"] = int(simulation.lmp.get_thermo("step"))
        for start in range(0, args.smoke_steps, 100):
            simulation.lmp.command(f"run {min(100, args.smoke_steps - start)} post no")
            simulation._check_callback()
        final = simulation.coordinates().copy()
        simulation._check_callback()
        endpoints = np.array([(b["a"] - 1, b["b"] - 1) for b in simulation.prepared.bonds])
        initial_lengths = np.linalg.norm(initial[endpoints[:, 0]] - initial[endpoints[:, 1]], axis=1)
        final_lengths = np.linalg.norm(final[endpoints[:, 0]] - final[endpoints[:, 1]], axis=1)
        ratio = final_lengths / initial_lengths
        captured["smoke"] = dict(steps=int(simulation.lmp.get_thermo("step")),
                                 timestep_fs=float(simulation.lmp.extract_global("dt")),
                                 positions_finite=bool(np.isfinite(final).all()),
                                 final_center_A=simulation.center().tolist(),
                                 final_minimum_z_A=float(final[:, 2].min()),
                                 bond_length_ratio_range=[float(ratio.min()), float(ratio.max())])
    except Exception as exc:
        captured["error"] = str(exc)
    finally:
        if simulation is not None:
            simulation.close()
    captured["wall_seconds"] = time.perf_counter() - started
    (args.output / "report.json").write_text(json.dumps(captured, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: captured.get(key) for key in ("input", "surface", "atom_count",
                     "independent_force_norm", "finite_difference", "smoke", "error", "wall_seconds")},
                     ensure_ascii=True, indent=2))


if __name__ == "__main__":
    main()
