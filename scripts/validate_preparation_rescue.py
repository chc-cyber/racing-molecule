"""Real native preparation regression and rescue-budget cleanup checks."""
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
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from molecular_car.engine import Simulation
from molecular_car.relaxation import JointRelaxationError, relax_joint_geometry
from scripts.lammps_runtime import open_lammps
from scripts.validate_element_md import PRIMARY, INPUTS


def gather(lmp, name, count):
    raw = lmp.gather_atoms(name, 1, 3)
    return np.ctypeslib.as_array(raw, shape=(3 * count,)).reshape(count, 3).copy()


def budget_recovery():
    """Actually stop BFGS in native LAMMPS, verify checkpoint forces, then reuse it."""
    stiffness = np.array([1e-4, .01, 1., 10., 100., 1000.]).reshape(2, 3)
    target = np.array([[-2., 0., 0.], [2., 0., 0.]])
    with open_lammps() as (lmp, _):
        lmp.commands_list(["units real", "atom_style atomic", "boundary f f f",
                           "atom_modify map array sort 0 0.0",
                           "region world block -10 10 -10 10 -10 10", "create_box 1 world", "mass 1 12"])
        lmp.create_atoms(2, [1, 2], [1, 1], (target + .8).ravel().tolist())
        lmp.commands_list(["pair_style zero 10", "pair_coeff * *", "fix molecular all external pf/callback 1 1",
                           "thermo_style custom step pe", "thermo 1000000", "timestep 0.5"])

        def callback(caller, step, count, tags, x, force):
            ids = np.asarray(tags, dtype=int) - 1
            delta = x - target[ids]
            force[:] = -stiffness[ids] * delta
            lmp.fix_external_set_energy_global("molecular", float(.5 * np.sum(stiffness[ids] * delta * delta)))

        lmp.set_fix_external_callback("molecular", callback)
        try:
            relax_joint_geometry(lmp, 2, callback, lambda: None, stages=(), evaluation_limit=20)
        except JointRelaxationError as error:
            rejected = error.diagnostics
        else:
            raise AssertionError("BFGS bypassed its shared evaluation budget")
        assert rejected["termination_reason"] == "force-evaluation-budget"
        assert rejected["stages"][-1]["minimizer"] == "curvature-initialized-bfgs"
        x, f = gather(lmp, "x", 2), gather(lmp, "f", 2)
        force_error = float(np.max(np.abs(f + stiffness * (x - target))))
        assert force_error < 1e-10 and lmp.get_thermo("timeremain") == 0
        recovered = relax_joint_geometry(lmp, 2, callback, lambda: None, stages=())
        assert recovered["accepted"] and recovered["global_force_norm_kcal_mol_A"] <= 1e-4
        assert lmp.get_thermo("timeremain") == 0 and lmp.get_thermo("step") == 0
        return {"passed": True, "rejected_preparation": rejected,
                "checkpoint_force_max_error": force_error, "recovered_preparation": recovered}


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/joint-fullerene-diagnosis/native-regression.json")
    parser.add_argument("--fullerene-only", action="store_true")
    parser.add_argument("--surfaces", choices=("periodic-au", "fourier-au", "continuum"), nargs="+",
                        default=["periodic-au", "fourier-au", "continuum"])
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    report = {"passed": False, "element_surface_cases": [], "fullerene_surface_cases": []}
    started = time.perf_counter()
    sources = [ROOT / "molecular_car" / name for name in ("engine.py", "molecule.py", "relaxation.py", "preparation_bfgs.py")]
    report["source_sha256"] = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}

    def save():
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    def one(path, surface):
        simulation = Simulation(path.read_text(encoding="utf-8"), path.name, {"surface": surface})
        try:
            preparation = simulation.runtime["preparation"]
            assert preparation["accepted"] and preparation["global_force_norm_kcal_mol_A"] <= 1e-4
            assert preparation["optimization_force_callback_evaluations"] < 25000
            assert preparation["final_coordinate_force_verified"]
            assert simulation.lmp.get_thermo("step") == 0 and simulation.lmp.extract_global("dt") == .5
            if path.name == "connected-fullerene.mol" and surface == "periodic-au":
                x = simulation.coordinates()
                energy = float(simulation.prepared.forcefield.CalcEnergy(x.ravel().tolist()))
                force = -np.asarray(simulation.prepared.forcefield.CalcGrad(x.ravel().tolist())).reshape(-1, 3)
                fe, ff = simulation.field_evaluator.energy_and_forces(x, 1, 0., 0., 0., 0., simulation.settings["corrugation"])
                se, sf = simulation.surface.energy_and_forces(x)
                energy += fe + se; force += ff + sf
                assert abs(np.linalg.norm(force) - preparation["global_force_norm_kcal_mol_A"]) < 1e-10
                (args.output.parent / "native-static.json").write_text(json.dumps({"positions": x.ravel().tolist(),
                    "forces": force.ravel().tolist(), "energy": energy}, allow_nan=False) + "\n", encoding="utf-8")
            simulation.lmp.command("run 200 post no"); simulation._check_callback()
            assert simulation.lmp.get_thermo("step") == 200 and np.isfinite(simulation.coordinates()).all()
            return {"input": path.name, "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    "surface": surface, "atoms": simulation.n, "preparation": preparation,
                    "smoke_steps": 200, "passed": True}
        finally:
            simulation.close()

    try:
        for surface in args.surfaces:
            for _, name in (() if args.fullerene_only else PRIMARY):
                row = one(INPUTS / (name + ".mol"), surface)
                report["element_surface_cases"].append(row)
                save(); print(f"PASS {surface} {name}", flush=True)
            row = one(ROOT / "tests/data/connected-fullerene.mol", surface)
            report["fullerene_surface_cases"].append(row)
            save(); print(f"PASS {surface} connected-fullerene", flush=True)
        report["native_budget_recovery"] = budget_recovery()
        report["source_unchanged"] = all(hashlib.sha256(p.read_bytes()).hexdigest() == report["source_sha256"][p.relative_to(ROOT).as_posix()] for p in sources)
        assert report["source_unchanged"]
        report["passed"] = True
    except Exception as error:
        report["error"] = str(error)
        raise
    finally:
        report["wall_seconds"] = time.perf_counter() - started
        save()
    print(f"PASS {len(report['element_surface_cases'])} element/surface preparations, "
          f"{len(report['fullerene_surface_cases'])} connected-fullerene preparations and native budget recovery")


if __name__ == "__main__":
    main()
