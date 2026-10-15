"""Independent release checks of conservative forces and timestep convergence.

Only this process overrides the validation timestep. Game defaults, force-field
parameters, services and historical reports are untouched. NVE is deliberately
undriven and has no Langevin bath; its 600 K is an initial kinetic condition.
"""
from __future__ import annotations

import argparse
import ctypes
from datetime import datetime, timezone
import hashlib
import inspect
import json
import os
from pathlib import Path
import sys
import time
import traceback

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[name] = "1"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from molecular_car.engine import Simulation, TIMESTEP_FS
from molecular_car.external_fields import ExternalFieldEvaluator
from molecular_car.molecule import prepare_molecule
from molecular_car.periodic_au import create_periodic_surface, PeriodicAuPotential, _normal_geometry
from molecular_car.fourier_au import create_fourier_surface

GAME_POLICY = "bonded-saip-game-v1"
LEGACY_POLICY = "distance-saip-reference-v1"
SEEDS = (98245165, 2718281, 1234567)
THRESHOLDS = {"total_force_fd_absolute_kcal_mol_A": 1e-3,
              "total_force_fd_relative": 2e-6,
              "surface_period_energy_kcal_mol": 1e-6,
              "surface_period_force_kcal_mol_A": 1e-6,
              "seam_or_blend_force_jump_kcal_mol_A": 2e-3,
              "topology_cutoff_energy_jump_kcal_mol": 1e-3,
              "nve_coarse_drift_fraction_initial_ke": .01,
              "nve_fine_drift_fraction_initial_ke": .004,
              "fine_to_coarse_drift_maximum": .8,
              "drift_convergence_absolute_floor_kcal_mol": 1e-6,
              "bond_relative_minimum": .65, "bond_relative_maximum": 1.45,
              "full_energy_accounting_kcal_mol": 1e-6}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    names = sorted(ROOT.glob("molecular_car/*.py"))
    names += sorted(path for path in (ROOT / "molecular_car/data").rglob("*") if path.is_file())
    names += [Path(__file__).resolve()]
    return {str(path.relative_to(ROOT)).replace("\\", "/"): digest(path) for path in names}


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def load_fixture(path):
    path = Path(path).resolve()
    text = path.read_text(encoding="utf-8-sig")
    return {"path": path, "content": text, "prepared": prepare_molecule(text, path.name),
            "input_sha256": digest(path)}


def surface_factory(model, atoms, positions, bonds, policy):
    factory = create_periodic_surface if model == "periodic-au" else create_fourier_surface
    if "normal_policy" not in inspect.signature(PeriodicAuPotential).parameters:
        if policy != LEGACY_POLICY:
            raise RuntimeError("The bonded-normal repair is not available yet; wait for the production source to stabilize.")
        return factory(atoms, positions)
    return factory(atoms, positions, bonds=bonds, normal_policy=policy)


def total_evaluator(prepared, surface, *, strength=.75):
    q = np.array([atom["charge"] for atom in prepared.atoms])
    a = np.array([atom["polarizability"] for atom in prepared.atoms])
    external = ExternalFieldEvaluator(q, a, surface_enabled=False)
    def evaluate(positions):
        flat = tuple(float(value) for value in positions.ravel())
        uff_energy = float(prepared.forcefield.CalcEnergy(flat))
        uff_force = -np.asarray(prepared.forcefield.CalcGrad(flat)).reshape(-1, 3)
        field_energy, field_force = external.energy_and_forces(positions, 1, strength, 6., 0., 7., 0.)
        surface_energy, surface_force = surface.energy_and_forces(positions)
        return uff_energy + field_energy + surface_energy, uff_force + field_force + surface_force
    return evaluate


def reference_comparison(surface, prepared, snapshots, model):
    """Reference agreement is required only with unchanged, stable normal sets."""
    try:
        reference = surface_factory(model, prepared.atoms, snapshots[0], prepared.bonds, LEGACY_POLICY)
    except ValueError as error:
        return {"passed": True, "qualification": "Legacy distance reference rejects this fixture; adapted sp3/sp axis-normal models do not claim equivalence.",
                "legacy_rejection": str(error), "cases": []}
    rows = []
    for index, coords in enumerate(snapshots):
        try:
            old_geometry = _normal_geometry(coords, surface.normal_cutoffs)
            new_geometry = surface.normal_geometry(coords) if hasattr(surface, "normal_geometry") else old_geometry
            same = True
            for atom in range(len(coords)):
                old_count, new_count = int(old_geometry[2][atom]), int(new_geometry[2][atom])
                if old_count <= 1 and new_count <= 1:
                    continue
                if old_count != new_count or set(old_geometry[1][atom, :old_count]) != set(new_geometry[1][atom, :new_count]):
                    same = False; break
            normals_close = float(np.abs(old_geometry[0] - new_geometry[0]).max()) < 1e-8
            old_e, old_f = reference.energy_and_forces(coords)
            new_e, new_f = surface.energy_and_forces(coords)
            de, df = abs(old_e - new_e), float(np.abs(old_f - new_f).max())
            required = same and normals_close
            rows.append({"snapshot_index": index, "unchanged_normal_adjacency": same,
                         "stable_reference_normal_geometry": normals_close,
                         "reference_agreement_required": required,
                         "energy_difference_kcal_mol": de, "max_force_difference_kcal_mol_A": df,
                         "passed": (de < 1e-6 and df < 1e-6) if required else True})
        except ValueError as error:
            rows.append({"snapshot_index": index, "reference_agreement_required": False,
                         "legacy_rejection": str(error), "passed": True})
    return {"passed": all(row["passed"] for row in rows), "cases": rows,
            "qualification": "Native-distance-reference equivalence is checked only while the normal adjacency and well-conditioned normals coincide. Adapted sp3/sp/switch states are reported separately."}


def force_fd(evaluate, positions, *, selected=None, epsilon=1e-5):
    energy, force = evaluate(positions)
    axes = selected if selected is not None else [(i, a) for i in range(len(positions)) for a in range(3)]
    errors, normalized = [], []
    for atom, axis in axes:
        plus, minus = positions.copy(), positions.copy()
        plus[atom, axis] += epsilon
        minus[atom, axis] -= epsilon
        derivative = -(evaluate(plus)[0] - evaluate(minus)[0]) / (2 * epsilon)
        error = abs(float(force[atom, axis]) - derivative)
        tolerance = THRESHOLDS["total_force_fd_absolute_kcal_mol_A"] + THRESHOLDS["total_force_fd_relative"] * abs(derivative)
        errors.append(error)
        normalized.append(error / tolerance)
    return {"energy_kcal_mol": energy, "checked_components": len(axes), "epsilon_A": epsilon,
            "max_force_fd_error_kcal_mol_A": max(errors, default=0.),
            "max_error_over_tolerance": max(normalized, default=0.),
            "passed": max(normalized, default=0.) <= 1}


def translated_to_height(positions, height):
    coordinates = positions.copy()
    coordinates[:, 2] += height - float(coordinates[:, 2].min())
    return coordinates


def rigid_benzene(prepared):
    positions = prepared.positions.copy()
    carbon = [i for i, atom in enumerate(prepared.atoms) if atom["element"] == "C"]
    positions -= positions[carbon].mean(axis=0)
    _, _, axes = np.linalg.svd(positions[carbon])
    normal = axes[-1]
    x = positions[carbon[1]] - positions[carbon[0]]
    x -= normal * np.dot(x, normal)
    x /= np.linalg.norm(x)
    y = np.cross(normal, x)
    positions = positions @ np.stack([x, y, normal]).T
    positions[:, :2] += [.371, .719]
    return positions


def static_checks(fixtures, policy, out):
    reports = []
    for fixture in fixtures:
        prepared = fixture["prepared"]
        for model in ("periodic-au", "fourier-au"):
            positions = translated_to_height(prepared.positions, 3.4)
            surface = surface_factory(model, prepared.atoms, positions, prepared.bonds, policy)
            evaluate = total_evaluator(prepared, surface)
            cases = []
            for name, coords in (("prepared", positions),
                                ("local_distortion", positions + np.random.default_rng(6012).normal(0., .015, positions.shape))):
                cases.append({"case": name, **force_fd(evaluate, coords)})
            if fixture["path"].name == "benzene.mol":
                rigid = rigid_benzene(prepared)
                # The radial cache blends at 2..2.3 and 11.7..12 Å. The
                # Fourier height support additionally blends at 2.8..3.1
                # and 7.7..8 Å; its normal chart is bounded by |nx/nz|=2.
                heights = (1.999999, 2.000001, 2.15, 2.299999, 2.300001,
                           2.799999, 2.800001, 2.95, 3.099999, 3.100001,
                           7.699999, 7.700001, 7.85, 7.999999, 8.000001,
                           11.699999, 11.700001, 11.85, 11.999999, 12.000001)
                for height in heights:
                    coords = translated_to_height(rigid, height)
                    cases.append({"case": f"height_support_{height:.6f}",
                                  **force_fd(evaluate, coords, selected=[(0, 0), (0, 2), (6, 2)])})
                for angle in (57.9946, 58.0, 61., 63.43494, 63.43496, 70., 89.5):
                    radians = np.deg2rad(angle)
                    rotation = np.array([[1., 0., 0.], [0., np.cos(radians), -np.sin(radians)],
                                         [0., np.sin(radians), np.cos(radians)]])
                    coords = translated_to_height(rigid @ rotation.T, 3.4)
                    cases.append({"case": f"normal_chart_{angle:g}deg",
                                  **force_fd(evaluate, coords, selected=[(i, a) for i in (0, 1, 2, 6) for a in range(3)])})
            # Only the Au surface is periodic. UFF preserves whole unwrapped
            # coordinates, and the drive and track walls deliberately do not repeat.
            base_energy, base_force = surface.energy_and_forces(positions)
            translations = []
            for multiples in ((1, 0), (-7, 11), (120, -50)):
                moved = positions + [*(surface.period * multiples), 0.]
                energy, force = surface.energy_and_forces(moved)
                de, df = abs(energy - base_energy), float(np.abs(force - base_force).max())
                translations.append({"cell_multiples": list(multiples), "energy_error_kcal_mol": de,
                                     "max_force_error_kcal_mol_A": df,
                                     "passed": de < THRESHOLDS["surface_period_energy_kcal_mol"]
                                               and df < THRESHOLDS["surface_period_force_kcal_mol_A"]})
            seams = []
            for axis in (0, 1):
                coords = positions.copy()
                coords[:, axis] -= coords[0, axis]
                epsilon = 1e-6
                plus, minus = coords.copy(), coords.copy()
                plus[:, axis] += epsilon; minus[:, axis] -= epsilon
                e0, f0 = surface.energy_and_forces(coords)
                ep, fp = surface.energy_and_forces(plus)
                em, fm = surface.energy_and_forces(minus)
                jump = float(np.abs(fp - fm).max())
                derivative = -(ep - em) / (2 * epsilon)
                fd_error = abs(float(f0[:, axis].sum()) - derivative)
                seams.append({"axis": axis, "epsilon_A": epsilon, "force_jump_kcal_mol_A": jump,
                              "translation_fd_error_kcal_mol_A": fd_error,
                              "passed": jump < THRESHOLDS["seam_or_blend_force_jump_kcal_mol_A"]
                                        and fd_error < THRESHOLDS["total_force_fd_absolute_kcal_mol_A"]})
            row = {"fixture": fixture["path"].name, "input_sha256": fixture["input_sha256"],
                   "atoms": len(positions), "model": model, "normal_policy": policy,
                   "surface": surface.metadata(), "total_force_fd": cases,
                   "surface_period_translation": translations, "surface_seams": seams,
                   "passed": all(value["passed"] for value in [*cases, *translations, *seams])}
            row["legacy_reference_comparison"] = reference_comparison(surface, prepared, [positions], model)
            row["passed"] = row["passed"] and row["legacy_reference_comparison"]["passed"]
            reports.append(row)
            print(f"STATIC {model} {fixture['path'].name}: {len(cases)} FD cases passed={row['passed']}", flush=True)
    switching = normal_switch_checks(policy)
    report = {"passed": all(row["passed"] for row in reports) and switching["passed"],
              "normal_policy": policy, "cases": reports, "normal_switch": switching,
              "thresholds": THRESHOLDS, "source_sha256": source_hashes()}
    write_json(out / "static-report.json", report)
    return report


def normal_switch_checks(policy):
    """A numerical four-C graph deliberately crosses a nonbonded SAIP cutoff.

    This is a derivative fixture, not a chemical racing structure. The graph's
    bonded adjacency is unchanged on both sides; the distance reference changes
    local normals and remains available solely as the native-reference model.
    """
    base = np.array([[0., 0., 3.4], [1.4, 0., 3.55], [0., 1.4, 3.45], [-1.8, .7, 3.6]])
    atoms = [{"element": "C"}] * 4
    bonds = [{"a": 1, "b": 2, "order": 1}, {"a": 1, "b": 3, "order": 1}, {"a": 3, "b": 4, "order": 1}]
    rows = []
    for model in ("periodic-au", "fourier-au"):
        surface = surface_factory(model, atoms, base, bonds, policy)
        cutoff = float(surface.normal_cutoffs[0, 3])
        sweep = []
        for epsilon in (1e-3, 1e-4, 1e-5, 1e-6):
            pair = []
            for sign in (-1, 1):
                coords = base.copy()
                radius = cutoff + sign * epsilon
                coords[3, 0] = -np.sqrt(radius**2 - .7**2 - .2**2)
                e, f = surface.energy_and_forces(coords)
                pair.append((e, f))
            sweep.append({"epsilon_A": epsilon, "energy_jump_kcal_mol": abs(pair[1][0] - pair[0][0]),
                          "max_force_jump_kcal_mol_A": float(np.abs(pair[1][1] - pair[0][1]).max())})
        exact = base.copy(); exact[3, 0] = -np.sqrt(cutoff**2 - .7**2 - .2**2)
        fd = force_fd(surface.energy_and_forces, exact, selected=[(i, a) for i in range(4) for a in range(3)])
        latest = sweep[-1]
        passed = (fd["passed"] and latest["energy_jump_kcal_mol"] < THRESHOLDS["topology_cutoff_energy_jump_kcal_mol"]
                  and latest["max_force_jump_kcal_mol_A"] < THRESHOLDS["seam_or_blend_force_jump_kcal_mol_A"])
        # A distance-reference discontinuity is intentionally diagnostic and
        # is not accepted as evidence for the adapted game-normal policy.
        rows.append({"model": model, "normal_policy": policy, "cutoff_A": cutoff,
                     "sweep": sweep, "at_cutoff_fd": fd, "passed": passed})
    return {"passed": all(row["passed"] for row in rows), "cases": rows,
            "fixture_qualification": "Numerical derivative fixture; fixed connected graph with a nonbonded atom crossing the legacy spatial-normal cutoff."}


def gather(simulation, name):
    raw = simulation.lmp.gather_atoms(name, 1, 3)
    return np.ctypeslib.as_array(raw, shape=(simulation.n * 3,)).reshape(simulation.n, 3).copy()


def reset_initial_state(simulation, positions, velocities, dt):
    simulation.lmp.command("reset_timestep 0")
    simulation.lmp.command(f"timestep {dt:.12g}")
    for name, values in (("x", positions), ("v", velocities)):
        flat = (ctypes.c_double * values.size)(*values.ravel())
        simulation.lmp.scatter_atoms(name, 1, 3, flat)
    simulation.drive_enabled = False
    simulation.lmp.command("run 0")
    simulation._check_callback()


def nve_samples(simulation, duration_ps, sample_ps, dt, evaluate):
    steps_per_sample = round(sample_ps * 1000 / dt)
    if not np.isclose(steps_per_sample * dt, sample_ps * 1000):
        raise ValueError("Sampling interval must be an integer number of both timesteps.")
    intervals = round(duration_ps / sample_ps)
    bond_pairs = [(bond["a"] - 1, bond["b"] - 1) for bond in simulation.prepared.bonds]
    initial_x = gather(simulation, "x")
    initial_bonds = np.array([np.linalg.norm(initial_x[a] - initial_x[b]) for a, b in bond_pairs])
    samples, positions_history = [], []
    accounting = []
    started = time.perf_counter()
    for index in range(intervals + 1):
        if index:
            simulation.lmp.command(f"run {steps_per_sample} post no")
            simulation._check_callback()
        positions = gather(simulation, "x")
        bonds = np.array([np.linalg.norm(positions[a] - positions[b]) for a, b in bond_pairs])
        energy = float(simulation.lmp.get_thermo("etotal"))
        kinetic = float(simulation.lmp.get_thermo("ke"))
        if not np.isfinite(positions).all() or not np.isfinite([energy, kinetic]).all():
            raise AssertionError("Nonfinite NVE state")
        frame = simulation.frame()  # Retains the game's coordinate/world checks.
        samples.append({"time_ps": index * sample_ps, "step": int(simulation.lmp.get_thermo("step")),
                        "total_energy_kcal_mol": energy, "kinetic_energy_kcal_mol": kinetic,
                        "COM_A": frame["center"], "min_atom_height_A": float(positions[:, 2].min()),
                        "max_atom_height_A": float(positions[:, 2].max()),
                        "bond_min_A": float(bonds.min()), "bond_max_A": float(bonds.max()),
                        "bond_ratio_min": float((bonds / initial_bonds).min()),
                        "bond_ratio_max": float((bonds / initial_bonds).max())})
        positions_history.append(positions)
        if index in (0, intervals):
            potential, _ = evaluate(positions)
            accounting.append(abs(potential + kinetic - energy))
    energy = np.array([row["total_energy_kcal_mol"] for row in samples])
    delta = energy - energy[0]
    ke = samples[0]["kinetic_energy_kcal_mol"]
    drift = float(np.abs(delta).max())
    limit = THRESHOLDS["nve_coarse_drift_fraction_initial_ke" if dt == .5 else "nve_fine_drift_fraction_initial_ke"]
    report = {"dt_fs": dt, "duration_ps": duration_ps, "samples": samples,
              "max_energy_drift_kcal_mol": drift, "rms_energy_error_kcal_mol": float(np.sqrt(np.mean(delta**2))),
              "energy_error_slope_kcal_mol_ps": float(np.polyfit(np.arange(len(samples)) * sample_ps, delta, 1)[0]),
              "initial_ke_kcal_mol": ke, "max_drift_fraction_initial_ke": drift / max(ke, 1e-15),
              "maximum_full_energy_accounting_error_kcal_mol": max(accounting),
              "bond_ratio_range": [min(row["bond_ratio_min"] for row in samples), max(row["bond_ratio_max"] for row in samples)],
              "wall_seconds_including_sampling": time.perf_counter() - started}
    report["passed"] = (report["max_drift_fraction_initial_ke"] <= limit
                         and report["maximum_full_energy_accounting_error_kcal_mol"] < THRESHOLDS["full_energy_accounting_kcal_mol"]
                         and report["bond_ratio_range"][0] > THRESHOLDS["bond_relative_minimum"]
                         and report["bond_ratio_range"][1] < THRESHOLDS["bond_relative_maximum"])
    return report, np.asarray(positions_history)


def nve_checks(fixtures, policy, args, out):
    results = []
    for fixture in fixtures:
        for seed in args.seeds:
            simulation = Simulation(fixture["content"], fixture["path"].name,
                                    {"surface": "periodic-au", "temperature": args.temperature, "damping_fs": 350.},
                                    thermostat=False, random_seed=seed,
                                    drive={"mode": "gradient", "strength": 0.})
            try:
                actual_policy = simulation.runtime["surface"].get("normal_policy", LEGACY_POLICY)
                if actual_policy != policy:
                    raise RuntimeError(f"Game uses {actual_policy}, but requested validation policy is {policy}.")
                initial_x, initial_v = gather(simulation, "x"), gather(simulation, "v")
                evaluate = total_evaluator(simulation.prepared, simulation.surface, strength=0.)
                relaxed_energy, relaxed_force = evaluate(initial_x)
                force_norms = np.linalg.norm(relaxed_force, axis=1)
                preparation = simulation.runtime.get("preparation", {})
                residual = {"potential_energy_kcal_mol": relaxed_energy,
                            "max_atom_force_norm_kcal_mol_A": float(force_norms.max()),
                            "rms_atom_force_norm_kcal_mol_A": float(np.sqrt(np.mean(force_norms**2))),
                            "max_force_component_kcal_mol_A": float(np.abs(relaxed_force).max()),
                            "global_force_two_norm_kcal_mol_A": float(np.linalg.norm(relaxed_force)),
                            "configured_minimize": {"energy_tolerance": 1e-6, "force_tolerance_kcal_mol_A": 1e-4,
                                                    "maximum_iterations": 350, "maximum_force_evaluations": 3500},
                            "force_tolerance_met": bool(np.linalg.norm(relaxed_force) < 1e-4),
                            "actual_minimization_iterations": preparation.get("iterations"),
                            "actual_stop_reason": preparation.get("termination_reason"),
                            "game_preparation_diagnostics": preparation,
                            "independent_vs_preparation_force_norm_difference_kcal_mol_A": abs(float(np.linalg.norm(relaxed_force)) - preparation.get("global_force_norm_kcal_mol_A", float(np.linalg.norm(relaxed_force)))),
                            "stop_information_qualification": "Iterations and forces are captured by the game before timestep reset and thermal velocities. Its stopping reason is not exposed; native release logs provide that separately. Finite coordinates do not establish convergence."}
                methods, histories = [], {}
                for dt in (.5, .25):
                    reset_initial_state(simulation, initial_x, initial_v, dt)
                    if not np.array_equal(gather(simulation, "x"), initial_x) or not np.array_equal(gather(simulation, "v"), initial_v):
                        raise AssertionError("Paired NVE initial coordinates/velocities differ after reset.")
                    method, positions = nve_samples(simulation, args.duration_ps, args.sample_ps, dt, evaluate)
                    methods.append(method); histories[str(dt)] = positions
                    print(f"NVE {fixture['path'].name} seed={seed} dt={dt}: drift={method['max_energy_drift_kcal_mol']:.6g} kcal/mol passed={method['passed']}", flush=True)
                coarse, fine = methods
                convergence = fine["max_energy_drift_kcal_mol"] <= (THRESHOLDS["fine_to_coarse_drift_maximum"] * coarse["max_energy_drift_kcal_mol"]
                                                                              + THRESHOLDS["drift_convergence_absolute_floor_kcal_mol"])
                prefix = fixture["path"].stem + f"-seed-{seed}"
                trajectory = out / (prefix + "-paired.npz")
                np.savez_compressed(trajectory, initial_positions_A=initial_x, initial_velocities_A_fs=initial_v,
                                    positions_dt_0_5_A=histories["0.5"], positions_dt_0_25_A=histories["0.25"],
                                    time_ps=np.arange(len(histories["0.5"])) * args.sample_ps)
                delta = histories["0.5"] - histories["0.25"]
                comparison = reference_comparison(simulation.surface, simulation.prepared,
                                                  [histories["0.5"][i] for i in (0, len(histories["0.5"]) // 2, len(histories["0.5"]) - 1)],
                                                  "periodic-au")
                result = {"fixture": fixture["path"].name, "input_sha256": fixture["input_sha256"], "atoms": simulation.n,
                          "seed": seed, "runtime": simulation.runtime, "normal_policy": policy,
                          "initial_temperature_K": args.temperature, "thermostat": False, "drive_strength": 0.,
                          "joint_relaxation_residual": residual, "legacy_reference_comparison": comparison,
                          "same_initial_x_v": True, "initial_state_and_trajectories_npz": str(trajectory),
                          "trajectory_sha256": digest(trajectory), "dt_comparisons": methods,
                          "fine_step_reduces_energy_error": bool(convergence),
                          "fine_to_coarse_energy_drift_ratio": fine["max_energy_drift_kcal_mol"] / max(coarse["max_energy_drift_kcal_mol"], 1e-15),
                          "max_paired_coordinate_difference_A": float(np.abs(delta).max()),
                          "coordinate_difference_qualification": "Reported diagnostic only; chaotic finite-temperature trajectory separation is not the acceptance criterion.",
                          "passed": all(method["passed"] for method in methods) and bool(convergence) and comparison["passed"]}
                results.append(result); write_json(out / (prefix + "-report.json"), result)
            finally:
                simulation.close()
    report = {"passed": all(row["passed"] for row in results), "cases": results, "thresholds": THRESHOLDS,
              "all_joint_relaxations_met_configured_force_tolerance": all(row["joint_relaxation_residual"]["force_tolerance_met"] for row in results),
              "qualification": "Finite-duration conservative numerical tests, not calibration of UFF/SAIP, exact chaotic trajectory agreement, or converged race statistics.",
              "game_timestep_fs_unchanged": TIMESTEP_FS, "source_sha256": source_hashes()}
    write_json(out / "nve-report.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("static", "nve", "all"), default="all")
    parser.add_argument("--normal-policy", choices=(GAME_POLICY, LEGACY_POLICY), default=GAME_POLICY)
    parser.add_argument("--duration-ps", type=float, default=5.)
    parser.add_argument("--sample-ps", type=float, default=.05)
    parser.add_argument("--temperature", type=float, default=600.)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    parser.add_argument("--extra-mol", type=Path, action="append", default=[])
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/md-release-validation/numerics")
    args = parser.parse_args()
    if args.duration_ps <= 0 or args.sample_ps <= 0 or not np.isclose(args.duration_ps / args.sample_ps, round(args.duration_ps / args.sample_ps)):
        parser.error("Duration must be positive and an integer multiple of sample interval.")
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    fixtures = [load_fixture(ROOT / "examples" / name) for name in ("benzene.mol", "tetraphenylbenzene.mol")]
    fixtures += [load_fixture(path) for path in args.extra_mol]
    before = source_hashes(); started = time.perf_counter()
    result = {"status": "running", "created_utc": datetime.now(timezone.utc).isoformat(),
              "stage": args.stage, "normal_policy": args.normal_policy,
              "source_sha256": before, "reports": {}, "thresholds": THRESHOLDS}
    try:
        if args.stage in ("static", "all"):
            result["reports"]["static"] = static_checks(fixtures, args.normal_policy, out)
        if args.stage in ("nve", "all"):
            result["reports"]["nve"] = nve_checks(fixtures, args.normal_policy, args, out)
        if source_hashes() != before:
            raise RuntimeError("Production source changed while release numerical checks were running; rerun after it stabilizes.")
        result["status"] = "passed" if all(value["passed"] for value in result["reports"].values()) else "failed"
    except BaseException as error:
        result.update(status="failed", error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        result["wall_seconds"] = time.perf_counter() - started
        write_json(out / f"report-{args.stage}.json", result)
    if result["status"] != "passed":
        raise SystemExit("Numerical release checks did not satisfy the declared thresholds; inspect the report.")


if __name__ == "__main__":
    main()
