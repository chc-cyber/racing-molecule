"""Strict-preparation and mixed-element MD acceptance on real game Simulation.

Run only after the production source is frozen. --plan-only does not import
Simulation or start MD. These finite numerical checks do not calibrate the
provisional heteroatom/Au LJ model or establish arbitrary-environment support.
"""
from __future__ import annotations

import os
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[variable] = "1"

import argparse
from contextlib import contextmanager
import ctypes
import hashlib
import json
from pathlib import Path
import re
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np

INPUTS = ROOT / "outputs/md-preparation-elements/elements/inputs"
OUT = ROOT / "outputs/md-preparation-elements/elements/md"
PRIMARY = (
    ("H", "methane"), ("C", "benzene"), ("B", "phenylboronic-acid"),
    ("N", "pyridine"), ("O", "ethanol"), ("F", "fluorobenzene"),
    ("Si", "tetramethylsilane"), ("P", "trimethylphosphine-oxide"),
    ("S", "dimethyl-sulfide"), ("Cl", "chlorobenzene"),
    ("Br", "bromobenzene"), ("I", "iodobenzene"),
)
FOURIER = ("pyridine", "trimethylphosphine-oxide", "iodobenzene")
NVE = ("pyridine", "tetramethylsilane", "dimethyl-sulfide", "iodobenzene")
SEED = 98245165
NVE_SEEDS = (SEED, 2718281)
THRESHOLDS = {
    "fd_absolute_kcal_mol_A": 1e-3, "fd_relative": 2e-6,
    "mass_error_u": 1e-8, "center_error_A": 1e-10,
    "conservative_callback_error_kcal_mol_A": 1e-7,
    "energy_accounting_error_kcal_mol": 1e-6,
    "bond_ratio_min": .65, "bond_ratio_max": 1.45,
    "nve_coarse_fraction_initial_ke": .01, "nve_fine_fraction_initial_ke": .004,
    "nve_fine_to_coarse_maximum": .8, "nve_reduction_absolute_floor_kcal_mol": 1e-6,
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def source_hashes():
    files = sorted((ROOT / "molecular_car").glob("*.py"))
    files += sorted(p for p in (ROOT / "molecular_car/data").rglob("*") if p.is_file())
    files += [ROOT / "scripts/lammps_runtime.py", ROOT / "scripts/validate_element_extension.py", Path(__file__).resolve()]
    files += [INPUTS / (name + ".mol") for _, name in PRIMARY]
    files += [INPUTS / (name + ".mol") for name in ("methyl-sulfur-pentafluoride", "iodoxybenzene", "iodosylbenzene", "trimethylphosphine")]
    return {str(path.relative_to(ROOT)).replace("\\", "/"): digest(path) for path in files}


def construct(name, model, log, *, bath=True, strength=.75, random_seed=SEED):
    # Native logging is passive and local to this validation process.
    from molecular_car import engine
    original = engine.open_lammps
    @contextmanager
    def logged_open():
        with original() as (lmp, runtime):
            lmp.command('log "' + str(log.resolve()).replace("\\", "/") + '"')
            yield lmp, runtime
    engine.open_lammps = logged_open
    try:
        path = INPUTS / (name + ".mol")
        return engine.Simulation(path.read_text(encoding="utf-8"), path.name,
            {"surface": model, "temperature": 600., "damping_fs": 350.},
            thermostat=bath, random_seed=random_seed, drive={"mode": "gradient", "strength": strength})
    finally:
        engine.open_lammps = original


def native_stops(log):
    text = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
    return re.findall(r"Stopping criterion\s*=\s*([^\r\n]+)", text)


def gather(sim, name):
    raw = sim.lmp.gather_atoms(name, 1, 3)
    return np.ctypeslib.as_array(raw, shape=(sim.n * 3,)).reshape(sim.n, 3).copy()


def actual_masses(sim):
    types = np.ctypeslib.as_array(sim.lmp.gather_atoms("type", 0, 1), shape=(sim.n,)).copy()
    mass_table = sim.lmp.extract_atom("mass")
    return np.array([float(mass_table[int(t)]) for t in types])


def evaluator(sim, strength):
    """Independent scalar energy and analytic conservative force assembly.

    The field uses its NumPy reference; SAIP/LJ uses the production surface.
    Finite differences differentiate the assembled scalar energy. Langevin's
    random/damping force is excluded, while drive and soft walls are included.
    """
    from molecular_car.external_fields import ExternalFieldEvaluator
    fields = ExternalFieldEvaluator(sim.charges, sim.alphas, surface_enabled=False, accelerator="numpy")
    def evaluate(x):
        flat = tuple(float(value) for value in x.ravel())
        energy = float(sim.prepared.forcefield.CalcEnergy(flat))
        force = -np.asarray(sim.prepared.forcefield.CalcGrad(flat)).reshape(sim.n, 3)
        surface_energy, surface_force = sim.surface.energy_and_forces(x)
        external_energy, external_force = fields.energy_and_forces(x, 1, strength, 0., 0., 0., 0.)
        energy += surface_energy + external_energy
        force = force + surface_force + external_force
        if not np.isfinite(energy) or not np.isfinite(force).all():
            raise AssertionError("Nonfinite independent conservative energy or forces")
        return float(energy), force
    return evaluate


def full_fd(evaluate, x, epsilon=1e-5):
    energy, force = evaluate(x)
    error_max = 0.
    normalized_max = 0.
    worst_component = None
    for atom in range(len(x)):
        for axis in range(3):
            plus, minus = x.copy(), x.copy()
            plus[atom, axis] += epsilon
            minus[atom, axis] -= epsilon
            fd = -(evaluate(plus)[0] - evaluate(minus)[0]) / (2 * epsilon)
            error = abs(float(force[atom, axis]) - fd)
            tolerance = THRESHOLDS["fd_absolute_kcal_mol_A"] + THRESHOLDS["fd_relative"] * abs(fd)
            normalized = error / tolerance
            if normalized >= normalized_max:
                worst_component = [atom, axis]
                normalized_max = normalized
            error_max = max(error_max, error)
    return {"potential_energy_kcal_mol": energy, "checked_components": x.size,
            "epsilon_A": epsilon, "maximum_force_error_kcal_mol_A": error_max,
            "maximum_error_over_tolerance": normalized_max,
            "worst_component_atom_axis_zero_based": worst_component,
            "passed": normalized_max <= 1.}


def preparation_check(sim):
    report = sim.runtime["preparation"]
    x = gather(sim, "x")
    energy, force = evaluator(sim, 0.)(x)
    norm = float(np.linalg.norm(force))
    tolerance = float(report["global_force_tolerance_kcal_mol_A"])
    accepted = bool(report.get("accepted") and report.get("force_tolerance_met") and
                    norm <= tolerance + 1e-10 and report.get("final_coordinate_force_verified"))
    independent = {"global_force_two_norm_kcal_mol_A": norm,
                   "maximum_atom_force_kcal_mol_A": float(np.linalg.norm(force, axis=1).max()),
                   "potential_energy_kcal_mol": energy,
                   "native_gate_norm_difference_kcal_mol_A": abs(norm - report["global_force_norm_kcal_mol_A"]),
                   "force_tolerance_met": norm <= tolerance + 1e-10,
                   "qualification": "Coordinates checked after actual zero-drive deterministic joint minimization and before any integration; constructor's velocity creation/run 0 does not move positions."}
    if not accepted:
        raise AssertionError("Actual constructor did not independently meet the joint conservative force gate")
    return {"game_metadata": report, "independent": independent, "passed": accepted}


def snapshot(sim, index, time_ps, initial_bonds, pairs, evaluate, *, check_fd=False):
    frame = sim.frame()  # Real game's finite/world bounds admission.
    x, v, native_f = gather(sim, "x"), gather(sim, "v"), gather(sim, "f")
    if not np.isfinite(x).all() or not np.isfinite(v).all() or not np.isfinite(native_f).all():
        raise AssertionError("Nonfinite actual LAMMPS x/v/f")
    bonds = np.linalg.norm(x[pairs[:, 0]] - x[pairs[:, 1]], axis=1)
    ratio = bonds / initial_bonds
    conservative_callback = sim._molecular_forces.copy()
    potential, conservative = evaluate(x)
    callback_error = float(np.abs(conservative_callback - conservative).max())
    center = np.average(x, axis=0, weights=sim.masses)
    center_error = float(np.abs(center - np.array(frame["center"])).max())
    kinetic = float(sim.lmp.get_thermo("ke"))
    total = float(sim.lmp.get_thermo("etotal"))
    energy_accounting = abs(potential + kinetic - total)
    if not np.isfinite([kinetic, total]).all():
        raise AssertionError("Nonfinite actual native energy")
    row = {"sample": index, "time_ps": time_ps, "native_step": int(sim.lmp.get_thermo("step")),
           "COM_A": center.tolist(), "center_error_A": center_error,
           "total_energy_kcal_mol": total, "kinetic_energy_kcal_mol": kinetic,
           "conservative_potential_energy_kcal_mol": potential,
           "energy_accounting_error_kcal_mol": energy_accounting,
           "conservative_callback_force_error_kcal_mol_A": callback_error,
           "native_total_force_max_atom_kcal_mol_A": float(np.linalg.norm(native_f, axis=1).max()),
           "conservative_force_max_atom_kcal_mol_A": float(np.linalg.norm(conservative, axis=1).max()),
           "min_atom_height_A": float(x[:, 2].min()), "max_atom_height_A": float(x[:, 2].max()),
           "bond_length_range_A": [float(bonds.min()), float(bonds.max())],
           "bond_ratio_range": [float(ratio.min()), float(ratio.max())]}
    if check_fd:
        row["full_conservative_force_fd"] = full_fd(evaluate, x)
    row["passed"] = (center_error < THRESHOLDS["center_error_A"] and
                     callback_error < THRESHOLDS["conservative_callback_error_kcal_mol_A"] and
                     energy_accounting < THRESHOLDS["energy_accounting_error_kcal_mol"] and
                     ratio.min() > THRESHOLDS["bond_ratio_min"] and ratio.max() < THRESHOLDS["bond_ratio_max"] and
                     row.get("full_conservative_force_fd", {"passed": True})["passed"])
    return row, x, v


def sample_sequence(sim, duration, sample_ps, dt, evaluate, *, fd_times=()):
    initial = gather(sim, "x")
    pairs = np.array([[b["a"] - 1, b["b"] - 1] for b in sim.prepared.bonds], dtype=int)
    lengths = np.linalg.norm(initial[pairs[:, 0]] - initial[pairs[:, 1]], axis=1)
    count = round(duration / sample_ps)
    steps = round(sample_ps * 1000 / dt)
    if not np.isclose(count * sample_ps, duration) or not np.isclose(steps * dt, sample_ps * 1000):
        raise ValueError("Duration/sampling must match an integer native step schedule")
    rows, xs, vs = [], [], []
    for index in range(count + 1):
        if index:
            sim.lmp.command(f"run {steps} post no")
            sim._check_callback()
        t = index * sample_ps
        check_fd = any(np.isclose(t, value) for value in fd_times)
        row, x, v = snapshot(sim, index, t, lengths, pairs, evaluate, check_fd=check_fd)
        rows.append(row); xs.append(x); vs.append(v)
    longest_off, current_off = 0., 0.
    for row in rows:
        if row["min_atom_height_A"] > 8.:
            current_off += sample_ps
            longest_off = max(longest_off, current_off)
        else:
            current_off = 0.
    summary = {
        "samples": rows,
        "bond_ratio_range": [min(r["bond_ratio_range"][0] for r in rows), max(r["bond_ratio_range"][1] for r in rows)],
        "atom_height_range_A": [min(r["min_atom_height_A"] for r in rows), max(r["max_atom_height_A"] for r in rows)],
        "max_consecutive_all_atoms_above_8A_ps": max(0., longest_off - sample_ps),
        "independent_energy_accounting_error_max_kcal_mol": max(r["energy_accounting_error_kcal_mol"] for r in rows),
        "conservative_callback_force_error_max_kcal_mol_A": max(r["conservative_callback_force_error_kcal_mol_A"] for r in rows),
        "center_error_max_A": max(r["center_error_A"] for r in rows),
        "passed": all(row["passed"] for row in rows),
    }
    return summary, np.asarray(xs), np.asarray(vs), pairs


def driven_case(element, name, model, out):
    key = model + "-" + name
    log = out / (key + ".log")
    sim = None
    began = time.perf_counter()
    result = {"key": key, "target_element": element, "fixture": name, "surface": model,
              "seed": SEED, "temperature_K": 600., "drive_strength": .75, "duration_ps": 5.,
              "native_log": str(log.relative_to(ROOT)).replace("\\", "/")}
    try:
        sim = construct(name, model, log)
        result.update(atoms=sim.n, runtime=sim.runtime, preparation=preparation_check(sim))
        native_mass = actual_masses(sim)
        error = float(np.abs(native_mass - sim.masses).max())
        result["native_mass_max_error_u"] = error
        data, xs, vs, pairs = sample_sequence(sim, 5., .05, .5, evaluator(sim, .75), fd_times=(0., 1., 5.))
        trajectory = out / (key + ".npz")
        np.savez_compressed(trajectory, positions_A=xs, velocities_A_fs=vs, times_ps=np.arange(len(xs)) * .05,
                            masses_u=native_mass, elements=np.array([a["element"] for a in sim.prepared.atoms]), bonds_zero_based=pairs)
        result.update(data, trajectory_file=str(trajectory.relative_to(ROOT)).replace("\\", "/"),
                      trajectory_sha256=digest(trajectory), displacement_A=(xs[-1].T @ native_mass / native_mass.sum() - xs[0].T @ native_mass / native_mass.sum()).tolist())
        result["passed"] = result["passed"] and error <= THRESHOLDS["mass_error_u"]
        result["status"] = "passed" if result["passed"] else "failed"
    except Exception as exc:
        result.update(status="failed", passed=False, exception_type=type(exc).__name__, error=str(exc), traceback=traceback.format_exc())
        if hasattr(exc, "diagnostics"):
            result["preparation_failure"] = exc.diagnostics
    finally:
        if sim is not None:
            sim.close()
        result["native_minimize_stopping_criteria"] = native_stops(log)
        result["wall_seconds"] = time.perf_counter() - began
        write_json(out / (key + ".json"), result)
    print(f"MIXED {key}: {result['status']}", flush=True)
    return result


def reset_state(sim, x, v, dt):
    sim.lmp.command("reset_timestep 0")
    sim.lmp.command(f"timestep {dt}")
    for name, values in (("x", x), ("v", v)):
        packed = (ctypes.c_double * values.size)(*values.ravel())
        sim.lmp.scatter_atoms(name, 1, 3, packed)
    sim.drive_enabled = False
    sim.lmp.command("run 0")
    sim._check_callback()
    if not np.array_equal(gather(sim, "x"), x) or not np.array_equal(gather(sim, "v"), v):
        raise AssertionError("Paired timestep checks do not begin from bit-identical x/v")


def nve_case(name, out, seed):
    key = "paired-nve-" + name + "-seed-" + str(seed)
    log = out / (key + ".log")
    sim = None
    result = {"key": key, "fixture": name, "surface": "periodic-au", "seed": seed,
              "duration_ps": 2., "initial_temperature_K": 600., "bath": False, "drive_strength": 0.,
              "native_log": str(log.relative_to(ROOT)).replace("\\", "/")}
    began = time.perf_counter()
    try:
        sim = construct(name, "periodic-au", log, bath=False, strength=0., random_seed=seed)
        result.update(runtime=sim.runtime, preparation=preparation_check(sim))
        initial_x, initial_v = gather(sim, "x"), gather(sim, "v")
        independent = evaluator(sim, 0.)
        methods, histories = [], {}
        for dt in (.5, .25):
            reset_state(sim, initial_x, initial_v, dt)
            method, xs, vs, pairs = sample_sequence(sim, 2., .05, dt, independent)
            method.update(dt_fs=dt, paired_initial_xv_bit_identical=True)
            energies = np.array([s["total_energy_kcal_mol"] for s in method["samples"]])
            ke = method["samples"][0]["kinetic_energy_kcal_mol"]
            maximum_error = float(np.abs(energies - energies[0]).max())
            method.update(initial_ke_kcal_mol=ke, maximum_energy_error_kcal_mol=maximum_error,
                          maximum_energy_error_fraction_initial_ke=maximum_error / max(ke, 1e-15))
            limit = THRESHOLDS["nve_coarse_fraction_initial_ke"] if dt == .5 else THRESHOLDS["nve_fine_fraction_initial_ke"]
            method["passed"] = method["passed"] and method["maximum_energy_error_fraction_initial_ke"] <= limit
            methods.append(method)
            histories["positions_dt_" + str(dt).replace(".", "_")] = xs
            histories["velocities_dt_" + str(dt).replace(".", "_")] = vs
        coarse, fine = methods
        reduction = fine["maximum_energy_error_kcal_mol"] <= (THRESHOLDS["nve_fine_to_coarse_maximum"] * coarse["maximum_energy_error_kcal_mol"] + THRESHOLDS["nve_reduction_absolute_floor_kcal_mol"])
        trajectory = out / (key + ".npz")
        np.savez_compressed(trajectory, **histories, initial_positions_A=initial_x, initial_velocities_A_fs=initial_v,
                            masses_u=actual_masses(sim), bonds_zero_based=pairs, times_ps=np.arange(len(methods[0]["samples"])) * .05)
        result.update(methods=methods, fine_step_reduces_maximum_energy_error=bool(reduction),
                      fine_to_coarse_energy_error_ratio=fine["maximum_energy_error_kcal_mol"] / max(coarse["maximum_energy_error_kcal_mol"], 1e-15),
                      trajectory_file=str(trajectory.relative_to(ROOT)).replace("\\", "/"), trajectory_sha256=digest(trajectory),
                      passed=bool(reduction and all(m["passed"] for m in methods)))
        result["status"] = "passed" if result["passed"] else "failed"
    except Exception as exc:
        result.update(status="failed", passed=False, error=str(exc), exception_type=type(exc).__name__, traceback=traceback.format_exc())
        if hasattr(exc, "diagnostics"):
            result["preparation_failure"] = exc.diagnostics
    finally:
        if sim is not None:
            sim.close()
        result["native_minimize_stopping_criteria"] = native_stops(log)
        result["wall_seconds"] = time.perf_counter() - began
        write_json(out / (key + ".json"), result)
    print(f"NVE {name}: {result['status']}", flush=True)
    return result


def rejection_case(name, out):
    log = out / ("rejection-" + name + ".log")
    sim = None
    result = {"fixture": name, "surface": "periodic-au", "expected": "explicit rejection before integration"}
    try:
        sim = construct(name, "periodic-au", log)
        result.update(passed=False, error="Unexpected constructor acceptance of the known problematic boundary fixture", runtime=sim.runtime)
    except Exception as exc:
        if name in ("methyl-sulfur-pentafluoride", "trimethylphosphine"):
            passed = type(exc).__name__ == "MoleculePreparationError" and "UFF" in str(exc)
        else:
            passed = ((type(exc).__name__ == "JointRelaxationError" and getattr(exc, "diagnostics", {}).get("termination_reason") == "pathological-initial-energy-or-force")
                      or (type(exc).__name__ == "MoleculePreparationError" and "UFF" in str(exc)))
        result.update(passed=bool(passed), exception_type=type(exc).__name__, rejection=str(exc))
        if hasattr(exc, "diagnostics"):
            result["preparation_failure"] = exc.diagnostics
    finally:
        if sim is not None:
            sim.close()
        result["native_minimize_stopping_criteria"] = native_stops(log)
    print(f"BOUNDARY {name}: passed={result['passed']}", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--stage", choices=("drive", "nve", "all"), default="all")
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    plan = {"periodic_5ps": [name for _, name in PRIMARY], "fourier_5ps": list(FOURIER),
            "periodic_paired_nve_2ps": list(NVE), "nve_seeds": list(NVE_SEEDS),
            "negative_boundaries": ["methyl-sulfur-pentafluoride", "iodoxybenzene", "iodosylbenzene", "trimethylphosphine"],
            "game_timestep_fs": .5, "diagnostic_fine_timestep_fs": .25, "thresholds": THRESHOLDS}
    if args.plan_only:
        print(json.dumps(plan, indent=2))
        return
    out = args.output.resolve()
    if not out.is_relative_to(ROOT / "outputs"):
        raise ValueError("Output must remain in the workspace outputs directory")
    out.mkdir(parents=True, exist_ok=True)
    before = source_hashes()
    report = {"status": "running", "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "source_sha256": before, "stage": args.stage, "plan": plan, "driven": [], "nve": [], "rejections": [],
              "protocol": {
                  "preparation": "Unmodified real Simulation constructor, deterministic strict joint UFF/surface/wall force convergence. No random thermal warming is substituted for minimization.",
                  "drive": "Immediately after preparation and velocity creation: 600 K Langevin, gradient .75, default .5 fs, 100-step batches/.05 ps for 5 ps; one seed per representative.",
                  "fd": "All 3N conservative force components at actual initial/1 ps/5 ps coordinates; independent UFF+surface+NumPy external scalar energy assembly, excludes Langevin forces and includes static drive/walls.",
                  "nve": "600 K initial kinetic condition only, no bath/drive; exact same prepared x/v reset for .5/.25 fs, 2 ps/.05 ps; finer timestep is validation-only.",
                  "random_number_consumption": "Constructor's exact production velocity/Langevin setup and run 0; then fixed 100-step MD batches. Independent energy/FD checks issue no native run and consume no thermostat RNG. NVE has no thermostat RNG; paired x/v are bit-identical.",
                  "mass": "Compare actual LAMMPS per-type integration masses with prepared masses; reconstruct full mass-weighted COM independently.",
              },
              "limitations": ["Heteroatom/Au LJ is a provisional nonreactive game approximation, distinct from CHAu SAIP. No claim of fitted adsorption, chemical bonding, image forces or charge transfer.",
                              "One seed and 5 ps are finite numerical/geometry acceptance, not converged race, diffusion or orientation statistics and not arbitrary environment coverage.",
                              "Bond/height diagnostics do not prove chemical validity or calibrated adsorption behaviour."]}
    began = time.perf_counter()
    try:
        if args.stage in ("drive", "all"):
            for element, name in PRIMARY:
                report["driven"].append(driven_case(element, name, "periodic-au", out))
                write_json(out / ("report-" + args.stage + ".json"), report)
            for name in FOURIER:
                element = next(e for e, n in PRIMARY if n == name)
                report["driven"].append(driven_case(element, name, "fourier-au", out))
                write_json(out / ("report-" + args.stage + ".json"), report)
            report["rejections"] = [rejection_case(name, out) for name in plan["negative_boundaries"]]
        if args.stage in ("nve", "all"):
            for name in NVE:
                for seed in NVE_SEEDS:
                    report["nve"].append(nve_case(name, out, seed))
                    write_json(out / ("report-" + args.stage + ".json"), report)
        after = source_hashes()
        report["source_sha256_after"] = after
        report["source_unchanged"] = before == after
        report["status"] = "passed" if before == after and all(r["passed"] for r in [*report["driven"], *report["nve"], *report["rejections"]]) else "failed"
    except BaseException as exc:
        report.update(status="failed", error=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        report["wall_seconds"] = time.perf_counter() - began
        write_json(out / ("report-" + args.stage + ".json"), report)
    print(json.dumps({"status": report["status"], "driven": len(report["driven"]), "nve_pairs": len(report["nve"]), "boundaries": len(report["rejections"]), "wall_seconds": report["wall_seconds"]}, indent=2))
    if report["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
