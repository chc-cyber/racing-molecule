"""Isolated pre/post verification of the game's external-field optimization.

The immutable archived Simulation is loaded beside the current Simulation.
Every result uses real single-thread LAMMPS and the same RDKit UFF object type;
no service is contacted and no integration model or time step is altered.
"""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
import traceback

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/continuous-optimization-check"
sys.path.insert(0, str(ROOT))

import numpy as np
from molecular_car import engine as optimized

SEEDS = (98245165, 1234567, 2718281)
FILES = ("benzene.mol", "tetraphenylbenzene.mol")
DT_FS = .5


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def source_hashes():
    names = ("molecular_car/engine.py", "molecular_car/fields.py", "molecular_car/external_fields.py",
             "molecular_car/fields_kernel.py", "molecular_car/molecule.py",
             "scripts/verify_continuum_optimization.py", "examples/benzene.mol", "examples/tetraphenylbenzene.mol",
             "outputs/continuous-optimization-check/reference/molecular_car/engine.py",
             "outputs/continuous-optimization-check/reference/molecular_car/fields.py")
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in names if (ROOT / name).exists()}


def load_reference():
    folder = OUT / "reference/molecular_car"
    required = {"engine.py": "08b1ec697ea8d73a2e4f7e7d58b07c5a3dc0c77cb24ee4ebe343594b28c5c062",
                "fields.py": "e21ba830461ce6f0ef8d472ded100d7f1d545583c2623949f1736078843d5d27"}
    for name, digest in required.items():
        if hashlib.sha256((folder / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Archived pre-optimization source hash differs: {name}")
    def load(name, path):
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        return module
    archived_fields = load("molecular_car._continuum_reference_fields", folder / "fields.py")
    old_binding = sys.modules["molecular_car.fields"]
    try:
        sys.modules["molecular_car.fields"] = archived_fields
        reference = load("molecular_car._continuum_reference_engine", folder / "engine.py")
    finally:
        sys.modules["molecular_car.fields"] = old_binding
    return reference


def create(module, filename, seed=SEEDS[0], *, temperature=600., thermostat=True, strength=.75, mode="gradient"):
    path = ROOT / "examples" / filename
    return module.Simulation(path.read_text(encoding="utf-8"), filename,
        settings={"surface": "continuum", "temperature": temperature, "corrugation": .08},
        thermostat=thermostat, random_seed=seed, drive={"mode": mode, "strength": strength})


def velocities(sim):
    raw = sim.lmp.gather_atoms("v", 1, 3)
    return np.ctypeslib.as_array(raw, shape=(sim.n * 3,)).copy().reshape(sim.n, 3)


def scatter(sim, field, values):
    flat = np.asarray(values, dtype=np.float64).ravel()
    raw = (ctypes.c_double * len(flat))(*flat)
    sim.lmp.scatter_atoms(field, 1, 3, raw)


def pair_start(old, new):
    before = float(np.max(np.abs(old.coordinates() - new.coordinates())))
    scatter(new, "x", old.coordinates())
    scatter(new, "v", velocities(old))
    # Langevin can draw random forces during a run-0 force refresh. Apply the
    # same extra setup to BOTH trajectories, otherwise their RNG streams are
    # offset even when x/v and the declared seed are identical.
    for sim in (old, new):
        sim.lmp.command("run 0")
        sim._check_callback()
    return before


def run_steps(sim, steps):
    sim.lmp.command(f"run {steps} pre no post no")
    sim._check_callback()


def call_forces(sim, positions, *, step=1, order=None):
    """Exercise the actual callback with the specified immutable coordinates."""
    order = np.arange(sim.n) if order is None else np.asarray(order)
    output = np.empty((sim.n, 3))
    class Capture:
        def fix_external_set_energy_global(self, name, energy):
            self.energy = float(energy)
    capture, native = Capture(), sim.lmp
    try:
        sim.lmp = capture
        sim.callback_error = None
        sim._forces(None, step, sim.n, order.astype(np.int64) + 1,
                    np.asarray(positions)[order].copy(), output)
        sim._check_callback()
        force = np.empty_like(output)
        force[order] = output
        return capture.energy, force
    finally:
        sim.lmp = native


def fixed_coordinate_checks(reference):
    rows, derivatives = [], []
    rng = np.random.default_rng(20261011)
    for filename in FILES:
        sims = []
        try:
            old, new = create(reference, filename), create(optimized, filename)
            sims.extend((old, new))
            geometry = old.coordinates()
            configs = [("gradient", .75, [0., 0., 0.], False),
                       ("gradient", 0., [0., 0., .7], False),
                       ("gradient", 1., [984., 243., .2], False),
                       ("tip", .75, [0., 0., .4], False),
                       ("tip", .35, [-12., -249., .1], False),
                       ("gradient", .75, [0., 0., 0.], True)]
            for case, (mode, strength, shift, pulse) in enumerate(configs):
                coords = geometry + np.asarray(shift)
                coords += rng.normal(0., .015, coords.shape)
                for sim in sims:
                    sim.update({"drive": {"mode": mode, "strength": strength}})
                    # The public prototype also permits dict controls. Use
                    # values different from defaults to catch stale caches.
                    sim.tip.update(x=11.5 if case == 3 else 6.,
                                   y=-3.25 if case == 3 else 0.,
                                   z=8.4 if case == 3 else 7.)
                    sim.pulse_until = 100 if pulse else -1
                expected_energy, expected_force = call_forces(old, coords)
                actual_energy, actual_force = call_forces(new, coords)
                permuted_energy, permuted_force = call_forces(new, coords, order=rng.permutation(new.n))
                identity_energy, identity_force = call_forces(new, coords)
                row = {"file": filename, "atoms": new.n, "case": case, "drive": mode,
                       "strength": strength, "pulse": pulse, "translation_A": shift,
                       "reference_energy_kcal_mol": expected_energy,
                       "energy_error_kcal_mol": actual_energy - expected_energy,
                       "max_force_error_kcal_mol_A": float(np.max(np.abs(actual_force - expected_force))),
                       "permuted_force_error_kcal_mol_A": float(np.max(np.abs(permuted_force - actual_force))),
                       "restored_identity_force_error_kcal_mol_A": float(np.max(np.abs(identity_force - actual_force))),
                       "permuted_energy_error_kcal_mol": permuted_energy - actual_energy,
                       "restored_identity_energy_error_kcal_mol": identity_energy - actual_energy}
                rows.append(row)
                # Full molecule Cartesian finite differences include native UFF
                # and both the external surface and driver. Twelve-atom gradient,
                # tip, boundary/pulse plus a 52-atom gradient fixture are checked.
                if (filename == FILES[0] and case in (0, 2, 3, 5)) or (filename == FILES[1] and case == 0):
                    for label, sim in zip(("reference", "optimized"), sims):
                        force = call_forces(sim, coords)[1]
                        numerical = np.zeros_like(force)
                        epsilon = 1e-5
                        for atom in range(sim.n):
                            for axis in range(3):
                                plus, minus = coords.copy(), coords.copy()
                                plus[atom, axis] += epsilon
                                minus[atom, axis] -= epsilon
                                numerical[atom, axis] = -(call_forces(sim, plus)[0] - call_forces(sim, minus)[0]) / (2 * epsilon)
                        derivatives.append({"file": filename, "case": case, "implementation": label,
                            "max_force_derivative_error_kcal_mol_A": float(np.max(np.abs(force - numerical)))})
        finally:
            for sim in reversed(sims):
                sim.close()
    maxima = {"energy_error_kcal_mol": max(abs(row["energy_error_kcal_mol"]) for row in rows),
              "force_error_kcal_mol_A": max(row["max_force_error_kcal_mol_A"] for row in rows),
              "permuted_force_error_kcal_mol_A": max(row["permuted_force_error_kcal_mol_A"] for row in rows),
              "restored_identity_force_error_kcal_mol_A": max(row["restored_identity_force_error_kcal_mol_A"] for row in rows),
              "force_derivative_error_kcal_mol_A": max(row["max_force_derivative_error_kcal_mol_A"] for row in derivatives)}
    result = {"cases": rows, "cartesian_derivative_checks": derivatives, "maxima": maxima,
              "acceptance": {"energy_kcal_mol": 1e-9, "force_kcal_mol_A": 1e-9, "derivative_kcal_mol_A": 1e-4},
              "passed": maxima["energy_error_kcal_mol"] < 1e-9 and maxima["force_error_kcal_mol_A"] < 1e-9
                and maxima["permuted_force_error_kcal_mol_A"] < 1e-9
                and maxima["restored_identity_force_error_kcal_mol_A"] < 1e-9
                and maxima["force_derivative_error_kcal_mol_A"] < 1e-4,
              "source_sha256": source_hashes()}
    write_json(OUT / "static-report.json", result)
    if not result["passed"]:
        raise AssertionError(f"Static callback check failed: {maxima}")
    print(f"PASS {len(rows)} fixed-coordinate cases and {len(derivatives)} full Cartesian derivative checks", flush=True)
    return result


def benchmark(reference, steps):
    trials = []
    for filename in FILES:
        for index, seed in enumerate(SEEDS):
            sims = []
            try:
                old, new = create(reference, filename, seed), create(optimized, filename, seed)
                sims.extend((old, new))
                natural_initial_difference = pair_start(old, new)
                order = [("reference", old), ("optimized", new)]
                if index % 2:
                    order.reverse()
                for _, sim in order:
                    run_steps(sim, 1000)
                batches = {label: [] for label, _ in order}
                for batch in range(3):
                    interval_order = order if batch % 2 == 0 else list(reversed(order))
                    for label, sim in interval_order:
                        started = time.perf_counter()
                        run_steps(sim, steps)
                        batches[label].append((time.perf_counter() - started) / steps * 1e6)
                costs = {label: float(np.median(values)) for label, values in batches.items()}
                trials.append({"file": filename, "atoms": old.n, "seed": seed,
                    "first_interval_order": [label for label, _ in order], "batches_us_per_step": batches,
                    "median_us_per_step": costs, "speedup": costs["reference"] / costs["optimized"],
                    "natural_initial_position_difference_A": natural_initial_difference,
                    "runtime": new.runtime})
                print(f"TIMING {filename} seed {seed}: {costs['reference']:.2f} -> {costs['optimized']:.2f} us/step", flush=True)
            finally:
                for sim in reversed(sims):
                    sim.close()
    summary = []
    for filename in FILES:
        subset = [row for row in trials if row["file"] == filename]
        costs = {label: float(np.median([row["median_us_per_step"][label] for row in subset]))
                 for label in ("reference", "optimized")}
        summary.append({"file": filename, "atoms": subset[0]["atoms"], "median_us_per_step": costs,
                        "speedup": costs["reference"] / costs["optimized"],
                        "cost_reduction_percent": 100 * (1 - costs["optimized"] / costs["reference"])})
    result = {"protocol": f"Actual game Simulation; paired identical initial state; one thread; 600 K; gradient .75; dt .5 fs; 1000-step warmup + 3 alternating {steps}-step intervals per seed; three seeds; compilation, setup, frames, rendering and sampling excluded",
              "trials": trials, "summary": summary, "passed": True, "source_sha256": source_hashes()}
    write_json(OUT / "benchmark-report.json", result)
    return result


def state_sample(sim):
    pos = sim.coordinates()
    carbon = pos[[i for i, atom in enumerate(sim.prepared.atoms) if atom["element"] == "C"]]
    _, _, axes = np.linalg.svd(carbon - carbon.mean(axis=0), full_matrices=False)
    return {"time_ps": float(sim.lmp.get_thermo("step")) * DT_FS / 1000.,
            "center_A": np.average(pos, axis=0, weights=sim.masses).tolist(),
            "minimum_z_A": float(pos[:, 2].min()),
            "carbon_plane_tilt_deg": float(np.degrees(np.arccos(np.clip(abs(axes[-1, 2]), 0., 1.)))),
            "total_energy_kcal_mol": float(sim.lmp.get_thermo("etotal")),
            "positions": pos}


def paired_motion(reference, total_steps):
    rows = []
    for filename in FILES:
        for seed in SEEDS:
            sims = []
            try:
                old, new = create(reference, filename, seed), create(optimized, filename, seed)
                sims.extend((old, new))
                natural_difference = pair_start(old, new)
                samples, comparison = {"reference": [], "optimized": []}, []
                for step in range(0, total_steps + 1, 100):
                    if step:
                        run_steps(old, 100)
                        run_steps(new, 100)
                    states = [state_sample(sim) for sim in sims]
                    delta = states[1]["positions"] - states[0]["positions"]
                    comparison.append({"time_ps": states[0]["time_ps"],
                        "max_position_difference_A": float(np.max(np.abs(delta))),
                        "rms_position_difference_A": float(np.sqrt(np.mean(delta**2))),
                        "center_difference_A": (np.array(states[1]["center_A"]) - states[0]["center_A"]).tolist(),
                        "total_energy_difference_kcal_mol": states[1]["total_energy_kcal_mol"] - states[0]["total_energy_kcal_mol"]})
                    for label, state in zip(samples, states):
                        state.pop("positions")
                        samples[label].append(state)
                maximum = max(row["max_position_difference_A"] for row in comparison)
                early = max(row["max_position_difference_A"] for row in comparison if row["time_ps"] <= 1.)
                row = {"file": filename, "atoms": old.n, "seed": seed,
                       "natural_initial_position_difference_A": natural_difference,
                       "samples": samples, "comparison": comparison,
                       "max_position_difference_A": maximum, "first_ps_max_position_difference_A": early,
                       "final_center_difference_A": comparison[-1]["center_difference_A"],
                       "passed": early < 1e-4 and maximum < .05}
                rows.append(row)
                write_json(OUT / "motion" / f"{Path(filename).stem}-{seed}.json", row)
                print(f"MOTION {filename} seed {seed}: max dx {maximum:.3g} A", flush=True)
            finally:
                for sim in reversed(sims):
                    sim.close()
    result = {"protocol": f"Paired identical post-minimization positions and velocities; same Langevin and velocity seed; 600 K, gradient .75, dt .5 fs, {total_steps * DT_FS / 1000:g} ps; 100-step sampling; all UFF and dynamical parameters unchanged",
              "acceptance": {"first_ps_max_position_difference_A": 1e-4, "whole_run_max_position_difference_A": .05},
              "limitations": "Identical equations need not preserve bitwise trajectories in chaotic MD. Thresholds test this short paired run; they do not establish universal equality at arbitrary time.",
              "trials": rows, "passed": all(row["passed"] for row in rows), "source_sha256": source_hashes()}
    write_json(OUT / "motion-report.json", result)
    if not result["passed"]:
        raise AssertionError("Paired short trajectories differ beyond the declared thresholds")
    return result


def nve_checks(reference):
    trials = []
    for filename in FILES:
        sims = []
        try:
            old = create(reference, filename, temperature=10., thermostat=False, strength=0.)
            new = create(optimized, filename, temperature=10., thermostat=False, strength=0.)
            sims.extend((old, new))
            pair_start(old, new)
            samples = {"reference": [], "optimized": []}
            for step in range(0, 2001, 100):
                if step:
                    for sim in sims:
                        run_steps(sim, 100)
                for label, sim in zip(samples, sims):
                    samples[label].append({"time_ps": step * DT_FS / 1000,
                                          "total_energy_kcal_mol": float(sim.lmp.get_thermo("etotal"))})
            drifts = {label: max(abs(value["total_energy_kcal_mol"] - values[0]["total_energy_kcal_mol"])
                                  for value in values) for label, values in samples.items()}
            delta = max(abs(a["total_energy_kcal_mol"] - b["total_energy_kcal_mol"])
                        for a, b in zip(samples["reference"], samples["optimized"]))
            trials.append({"file": filename, "atoms": old.n, "samples": samples,
                           "max_energy_drift_kcal_mol": drifts, "max_paired_energy_difference_kcal_mol": delta,
                           "passed": max(drifts.values()) < .2 and delta < 1e-5})
        finally:
            for sim in reversed(sims):
                sim.close()
    result = {"protocol": "10 K initial velocities, no thermostat, drive disabled, unchanged .5 fs step, 1 ps NVE, both molecules; drift is the existing velocity-Verlet integration error",
              "trials": trials, "passed": all(row["passed"] for row in trials), "source_sha256": source_hashes()}
    write_json(OUT / "nve-report.json", result)
    if not result["passed"]:
        raise AssertionError(f"NVE agreement check failed: {[(row['max_energy_drift_kcal_mol'], row['max_paired_energy_difference_kcal_mol']) for row in trials]}")
    print("PASS paired 1 ps NVE energy conservation", flush=True)
    return result


def profile_callback(module, filename, label):
    sim = create(module, filename)
    original_field = getattr(module, "energy_and_forces", None)
    totals = {"UFF_energy": 0., "UFF_gradient": 0., "surface_drive_walls": 0., "callback_total": 0., "callback_calls": 0}
    def timer(key, method):
        def wrapped(*args, **kwargs):
            began = time.perf_counter()
            value = method(*args, **kwargs)
            totals[key] += time.perf_counter() - began
            return value
        return wrapped
    forcefield = sim.prepared.forcefield
    class TimedUFF:
        CalcEnergy = staticmethod(timer("UFF_energy", forcefield.CalcEnergy))
        CalcGrad = staticmethod(timer("UFF_gradient", forcefield.CalcGrad))
    try:
        run_steps(sim, 1000)
        sim.prepared.forcefield = TimedUFF()
        if hasattr(sim, "field_evaluator"):
            evaluator = sim.field_evaluator
            class TimedField:
                energy_and_forces = staticmethod(timer("surface_drive_walls", evaluator.energy_and_forces))
                def __getattr__(self, name):
                    return getattr(evaluator, name)
            sim.field_evaluator = TimedField()
        else:
            module.energy_and_forces = timer("surface_drive_walls", original_field)
        callback = timer("callback_total", sim._forces)
        def timed_callback(*args):
            totals["callback_calls"] += 1
            callback(*args)
        sim.lmp.set_fix_external_callback("molecular", timed_callback)
        began = time.perf_counter()
        run_steps(sim, 2000)
        elapsed = time.perf_counter() - began
        if totals["callback_calls"] != 2000:
            raise AssertionError("Diagnostic callback count differs")
        components = {key: totals[key] / 2000 * 1e6 for key in ("UFF_energy", "UFF_gradient", "surface_drive_walls")}
        components["other_inside_callback"] = (totals["callback_total"] - sum(totals[key] for key in ("UFF_energy", "UFF_gradient", "surface_drive_walls"))) / 2000 * 1e6
        components["outside_callback"] = (elapsed - totals["callback_total"]) / 2000 * 1e6
        total = elapsed / 2000 * 1e6
        shares = {key: value / total for key, value in components.items()}
        p = shares["UFF_energy"] + shares["UFF_gradient"]
        callback_fraction = totals["callback_total"] / elapsed
        return {"file": filename, "atoms": sim.n, "implementation": label,
                "instrumented_us_per_step": total, "components_us_per_step": components,
                "component_fraction": shares, "callback_calls": totals["callback_calls"],
                "third_step_estimates": {"UFF_only_zero_cost_upper_bound_speedup": 1 / (1 - p),
                    "whole_callback_zero_cost_upper_bound_speedup": 1 / (1 - callback_fraction),
                    "UFF_only_scenarios": [{"assumed_UFF_speedup": speedup,
                        "whole_MD_speedup": 1 / ((1 - p) + p / speedup)} for speedup in (2., 4., 10.)]},
                "runtime": sim.runtime}
    finally:
        if original_field is not None:
            module.energy_and_forces = original_field
        sim.close()


def profile_checks(reference):
    rows = []
    for filename in FILES:
        for label, module in (("reference", reference), ("optimized", optimized)):
            rows.append(profile_callback(module, filename, label))
    result = {"protocol": "Separate instrumented 2000-step actual-game propagation after warmup, one seed; excluded from the three-seed performance benchmark",
              "limitations": "Timers and wrappers add overhead. Amdahl estimates are hypothetical component-removal scenarios, not a measured native UFF migration. Whole-callback zero cost is only a mathematical ceiling, not attainable speed.",
              "results": rows, "passed": True, "source_sha256": source_hashes()}
    write_json(OUT / "profile-report.json", result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("all", "static", "benchmark", "motion", "nve", "profile"), default="all")
    parser.add_argument("--steps", type=int, default=2000)
    parser.add_argument("--motion-steps", type=int, default=20000)
    args = parser.parse_args()
    if args.steps <= 0 or args.motion_steps <= 0 or args.motion_steps % 100:
        raise ValueError("Positive steps and a 100-step-divisible motion interval are required")
    OUT.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    initial_hashes = source_hashes()
    reference = load_reference()
    stages = {"static": lambda: fixed_coordinate_checks(reference),
              "benchmark": lambda: benchmark(reference, args.steps),
              "motion": lambda: paired_motion(reference, args.motion_steps),
              "nve": lambda: nve_checks(reference), "profile": lambda: profile_checks(reference)}
    try:
        results = {name: run() for name, run in stages.items() if args.stage in ("all", name)}
        if source_hashes() != initial_hashes:
            raise AssertionError("Source changed during this verification; repeat the affected stage")
        write_json(OUT / f"report-{args.stage}.json", {"status": "passed", "stage": args.stage,
            "reports": results, "wall_seconds": time.perf_counter() - start, "source_sha256": initial_hashes})
        print(f"PASS {args.stage} in {time.perf_counter() - start:.1f}s", flush=True)
    except Exception:
        write_json(OUT / f"report-{args.stage}.json", {"status": "failed", "stage": args.stage,
            "traceback": traceback.format_exc(), "wall_seconds": time.perf_counter() - start,
            "source_sha256": initial_hashes})
        raise


if __name__ == "__main__":
    main()
