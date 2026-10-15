"""Verify the static broad gradient field with actual LAMMPS trajectories.

Full run: .venv/Scripts/python.exe scripts/verify_gradient.py
Short run: add --quick (does not claim completion of the 100 nm track).
Only outputs/gradient-check is written.  No server or browser is required.

The lateral samples are short, finite-time trajectories, not a converged
diffusion-coefficient measurement.  The comparison estimate is an unbounded,
flat, free center of mass under Langevin friction, initially at rest; the
actual simulation also has molecular degrees of freedom and corrugation.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time
import traceback

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import rdkit

from molecular_car import fields
from molecular_car.engine import DEFAULT_SETTINGS, MAX_STEPS, TIMESTEP_FS, TRACK, Simulation


FILENAME = "tetraphenylbenzene.mol"
CONTENT = (PROJECT_ROOT / "examples" / FILENAME).read_text(encoding="utf-8")
OUTPUT_ROOT = PROJECT_ROOT / "outputs" / "gradient-check"
SEEDS = (98245165, 2718281, 1234567)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def save_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def external(simulation: Simulation, *, drive=None, tip=None, positions=None):
    return fields.energy_and_forces(
        simulation.coordinates() if positions is None else positions,
        simulation.charges,
        simulation.alphas,
        dict(simulation.tip) if tip is None else tip,
        simulation.settings["corrugation"],
        dict(simulation.drive) if drive is None else drive,
    )


def refresh(simulation: Simulation) -> None:
    simulation.lmp.command("run 0")
    simulation._check_callback()


def sample(simulation: Simulation) -> dict:
    frame = simulation.frame()
    positions = np.asarray(frame["positions"]).reshape(simulation.n, 3)
    raw = simulation.lmp.gather_atoms("f", 1, 3)
    forces = np.ctypeslib.as_array(raw, shape=(simulation.n * 3,)).reshape(simulation.n, 3)
    require(np.isfinite(positions).all() and np.isfinite(forces).all(), "Nonfinite positions or forces")
    require(np.isfinite(frame["energy"]), "Nonfinite total energy")
    lengths = np.array([
        np.linalg.norm(positions[bond["a"] - 1] - positions[bond["b"] - 1])
        for bond in simulation.prepared.bonds
    ])
    require(lengths.min() > 0.6 and lengths.max() < 2.5, "Molecular bond lengths became unstable")
    return {
        "step": int(frame["step"]),
        "time_ps": float(frame["time_fs"] / 1000.0),
        "center_angstrom": list(frame["center"]),
        "energy_kcal_mol": float(frame["energy"]),
        "max_atom_force_kcal_mol_angstrom": float(np.linalg.norm(forces, axis=1).max()),
        "bond_length_range_angstrom": [float(lengths.min()), float(lengths.max())],
        "min_z_angstrom": float(positions[:, 2].min()),
        "max_z_angstrom": float(positions[:, 2].max()),
        "internal_temperature_K": float(frame["internal_temperature"]),
        "drive": dict(frame["drive"]),
        "tip": dict(frame["tip"]),
    }


def summarize(simulation: Simulation, samples: list[dict], wall: float, label: str) -> dict:
    initial = np.asarray(samples[0]["center_angstrom"])
    final = np.asarray(samples[-1]["center_angstrom"])
    return {
        "label": label,
        "filename": FILENAME,
        "atoms": simulation.n,
        "mass_u": float(simulation.masses.sum()),
        "random_seed": simulation.random_seed,
        "thermostat": simulation.thermostat,
        "settings": dict(simulation.settings),
        "runtime": dict(simulation.runtime),
        "initial": samples[0],
        "final": samples[-1],
        "displacement_angstrom": (final - initial).tolist(),
        "bond_length_range_angstrom": [
            min(row["bond_length_range_angstrom"][0] for row in samples),
            max(row["bond_length_range_angstrom"][1] for row in samples),
        ],
        "maximum_atom_force_kcal_mol_angstrom": max(row["max_atom_force_kcal_mol_angstrom"] for row in samples),
        "samples": samples,
        "wall_seconds": wall,
    }


def fixed_field_structure_check() -> dict:
    simulation = None
    try:
        simulation = Simulation(CONTENT, FILENAME)
        require(simulation.drive == {"mode": "gradient", "strength": 0.75}, "Default driver contract changed")
        positions = simulation.coordinates()
        energy, forces = external(simulation, positions=positions)
        alternative_tip = dict(simulation.tip, x=410.0, y=95.0, z=13.0, strength=0.0)
        other_energy, other_forces = external(simulation, positions=positions, tip=alternative_tip)
        require(abs(energy - other_energy) < 1e-10 and np.max(np.abs(forces-other_forces)) < 1e-10,
                "Static gradient still depends on the movable tip")
        zero_energy, zero_forces = external(simulation, positions=positions, drive={"mode": "gradient", "strength": 0.0})
        driver_forces = forces - zero_forces
        require(np.max(np.abs(driver_forces[:, 1])) < 1e-10, "Gradient driver adds a transverse restoring force")
        require(driver_forces[:, 0].sum() > 0.0, "Gradient driver does not exert forward force")
        return {
            "tip_independent": True,
            "max_driver_y_force_kcal_mol_angstrom": float(np.max(np.abs(driver_forces[:, 1]))),
            "initial_driver_net_force_kcal_mol_angstrom": driver_forces.sum(axis=0).tolist(),
            "initial_driver_potential_kcal_mol": float(energy-zero_energy),
            "track": dict(TRACK),
        }
    finally:
        if simulation is not None:
            simulation.close()


def gradient_trial(*, max_steps: int, undriven=False, finish=False) -> dict:
    simulation = None
    started = time.perf_counter()
    # A later quick run must not overwrite the complete race's saved tail.
    label = "undriven-control" if undriven else ("static-gradient" if finish else "static-gradient-short")
    samples = []
    try:
        simulation = Simulation(CONTENT, FILENAME)
        if undriven:
            simulation.drive_enabled = False
            simulation.update({"drive": {"mode": "gradient", "strength": 0.0}})
        # The callback has both drive modes; disabling only tip strength is wrong.
        refresh(simulation)
        original_tip = dict(simulation.tip)
        samples.append(sample(simulation))
        total = 0
        next_progress = 40_000
        while total < max_steps:
            batch = min(2000, max_steps-total)
            simulation.advance(batch)
            total += batch
            row = sample(simulation)
            samples.append(row)
            require(row["tip"] == original_tip, "The static-field trial moved the tip")
            if total >= next_progress:
                print(f"Gradient {total} steps: x={row['center_angstrom'][0]:.1f} A, y={row['center_angstrom'][1]:.2f} A", flush=True)
                next_progress += 40_000
            if finish and row["center_angstrom"][0] >= TRACK["finish_x"]:
                break
        result = summarize(simulation, samples, time.perf_counter()-started, label)
        result["finished"] = bool(samples[-1]["center_angstrom"][0] >= TRACK["finish_x"])
        result["tip_stationary"] = True
        coordinates = simulation.coordinates()
        save_json(OUTPUT_ROOT / f"{label}-tail.json", {
            "format": "molecular-car-gradient-tail-v1",
            "atoms": simulation.prepared.atoms,
            "bonds": simulation.prepared.bonds,
            "positions_angstrom": coordinates.tolist(),
            "frame": simulation.frame(),
            "runtime": simulation.runtime,
        })
        return result
    except Exception:
        # Partial trajectory is useful when a proposed field needs recalibration.
        if samples:
            save_json(OUTPUT_ROOT / f"{label}-partial.json", {"samples": samples})
        raise
    finally:
        if simulation is not None:
            simulation.close()


def nve_trial() -> dict:
    simulation = None
    started = time.perf_counter()
    try:
        simulation = Simulation(CONTENT, FILENAME, thermostat=False)
        simulation.autopilot = False
        simulation.update({"drive": {"mode": "gradient", "strength": 0.75}})
        refresh(simulation)
        samples = [sample(simulation)]
        initial_energy = samples[0]["energy_kcal_mol"]
        coordinates = simulation.coordinates()
        ext_energy, _ = external(simulation, positions=coordinates)
        molecular_energy = float(simulation.prepared.forcefield.CalcEnergy(coordinates.ravel().tolist()))
        kinetic_energy = float(simulation.lmp.get_thermo("ke"))
        accounting_error = abs(initial_energy-molecular_energy-ext_energy-kinetic_energy)
        require(accounting_error < 1e-7, "External potential is missing from total energy")
        for _ in range(50):
            simulation.advance(100)
            samples.append(sample(simulation))
        energies = np.array([row["energy_kcal_mol"] for row in samples])
        absolute = float(np.max(np.abs(energies-initial_energy)))
        scale = max(abs(initial_energy), 1.0)
        result = summarize(simulation, samples, time.perf_counter()-started, "static-gradient-nve")
        result.update({
            "initial_uff_energy_kcal_mol": molecular_energy,
            "initial_external_energy_kcal_mol": float(ext_energy),
            "initial_kinetic_energy_kcal_mol": kinetic_energy,
            "energy_accounting_error_kcal_mol": accounting_error,
            "max_absolute_energy_drift_kcal_mol": absolute,
            "max_relative_energy_drift": absolute/scale,
            "absolute_drift_limit_kcal_mol": 0.05,
            "relative_drift_limit": 1e-3,
        })
        # A large negative polarization energy must not hide a large absolute error.
        require(absolute < 0.05 and absolute/scale < 1e-3, "Static-gradient NVE drift exceeds the absolute or relative tolerance")
        return result
    finally:
        if simulation is not None:
            simulation.close()


def mode_switch_trial() -> dict:
    simulation = None
    started = time.perf_counter()
    try:
        simulation = Simulation(CONTENT, FILENAME)
        before = simulation.coordinates()
        gradient = dict(simulation.drive)
        simulation.update({"drive": {"mode": "tip", "strength": 0.75}})
        simulation.autopilot = True
        refresh(simulation)
        require(np.array_equal(before, simulation.coordinates()), "Switching modes changed atom positions")
        require(simulation.frame()["drive"]["mode"] == "tip", "Frame does not expose the selected tip mode")
        samples = [sample(simulation)]
        for _ in range(100):
            simulation.advance(100)
            samples.append(sample(simulation))
        result = summarize(simulation, samples, time.perf_counter()-started, "switched-to-tip")
        require(result["displacement_angstrom"][0] > 1.0, "The old movable-tip mode no longer drives the car")
        before_return = simulation.coordinates()
        simulation.update({"drive": gradient})
        refresh(simulation)
        require(np.array_equal(before_return, simulation.coordinates()), "Returning to gradient changed atom positions")
        require(simulation.frame()["drive"] == gradient, "Returning to gradient did not restore the driver")
        result["returned_to_gradient"] = True
        result["switch_preserves_coordinates"] = True
        return result
    finally:
        if simulation is not None:
            simulation.close()


def free_com_rms(temperature: float, time_fs: float, damping_fs: float, mass_u: float) -> dict:
    # real units: acceleration = F * 4.184e-4 / mass, with fs and angstrom.
    boltzmann_kcal_mol_K = 0.00198720425864083
    thermal_variance = boltzmann_kcal_mol_K * temperature * 4.184e-4 / mass_u
    diffusion = thermal_variance * damping_fs
    tau = damping_fs
    # Integrated Ornstein-Uhlenbeck process, v_COM(0)=0, not equilibrium v0.
    factor = time_fs - 2*tau*(-np.expm1(-time_fs/tau)) + 0.5*tau*(-np.expm1(-2*time_fs/tau))
    variance = max(0.0, 2*diffusion*factor)
    return {
        "model": "Free flat unbounded COM Langevin, initial COM velocity zero",
        "rms_dy_angstrom": float(np.sqrt(variance)),
        "asymptotic_free_D_angstrom_squared_per_ps": float(diffusion*1000),
        "mass_u": mass_u,
        "damping_fs": damping_fs,
        "time_ps": time_fs/1000,
        "formula": "Var[y]=2 D [t-2 tau (1-exp(-t/tau))+tau/2 (1-exp(-2t/tau))], D=kBT tau/M",
        "qualification": "Reference estimate only; no fitted diffusion coefficient from the short trajectories",
    }


def lateral_trials(duration_ps: float) -> list[dict]:
    grouped = []
    steps = int(round(duration_ps*1000/TIMESTEP_FS))
    for temperature in (10.0, 150.0):
        trials = []
        for seed in SEEDS:
            simulation = None
            started = time.perf_counter()
            try:
                simulation = Simulation(CONTENT, FILENAME, {"temperature": temperature}, random_seed=seed)
                refresh(simulation)
                original_tip = dict(simulation.tip)
                samples = [sample(simulation)]
                for offset in range(0, steps, 1000):
                    simulation.advance(min(1000, steps-offset))
                    samples.append(sample(simulation))
                    require(samples[-1]["tip"] == original_tip, "A lateral gradient trial moved the tip")
                result = summarize(simulation, samples, time.perf_counter()-started, f"gradient-{temperature:g}K-seed-{seed}")
                trials.append(result)
                print(f"Lateral {temperature:g} K seed {seed}: dx={result['displacement_angstrom'][0]:.2f} A, dy={result['displacement_angstrom'][1]:.3f} A", flush=True)
            finally:
                if simulation is not None:
                    simulation.close()
        displacements = np.array([row["displacement_angstrom"][1] for row in trials])
        path_excursions = [
            abs(sample_row["center_angstrom"][1]-row["initial"]["center_angstrom"][1])
            for row in trials for sample_row in row["samples"]
        ]
        grouped.append({
            "temperature_K": temperature,
            "duration_ps": duration_ps,
            "seeds": list(SEEDS),
            "final_dy_angstrom": displacements.tolist(),
            "mean_final_dy_angstrom": float(displacements.mean()),
            "rms_final_dy_angstrom": float(np.sqrt(np.mean(displacements**2))),
            "final_dy_range_angstrom": [float(displacements.min()), float(displacements.max())],
            "largest_sampled_lateral_excursion_angstrom": float(max(path_excursions)),
            "track_half_width_angstrom": TRACK["half_width"],
            "free_com_reference": free_com_rms(temperature, duration_ps*1000, DEFAULT_SETTINGS["damping_fs"], trials[0]["mass_u"]),
            "interpretation": "Three short trajectories; descriptive finite-time spread, not converged diffusion or a race failure probability",
            "trials": trials,
        })
    return grouped


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="10 ps lateral trials and short drive; skip track-completion claim")
    args = parser.parse_args()
    started = time.perf_counter()
    report = {
        "passed": False,
        "quick": args.quick,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model": "prototype-uff-static-gradient-v1",
        "versions": {"python": sys.version.split()[0], "numpy": np.__version__, "rdkit": rdkit.__version__},
        "parameters": {
            "timestep_fs": TIMESTEP_FS,
            "default_settings": dict(DEFAULT_SETTINGS),
            "track": dict(TRACK),
            "field_constants": {key: value for key, value in vars(fields).items() if key.isupper() and isinstance(value, (int, float))},
        },
    }
    path = OUTPUT_ROOT / ("quick-report.json" if args.quick else "report.json")
    try:
        require(TRACK["finish_x"] == 1000.0 and TRACK["half_width"] == 250.0, "Track must be 100 nm by 50 nm")
        report["fixed_field_structure"] = fixed_field_structure_check()
        print("Static driver ignores tip position and has no y driver force.", flush=True)
        control = gradient_trial(max_steps=10_000, undriven=True)
        report["undriven_control"] = control
        driven = gradient_trial(max_steps=10_000 if args.quick else min(MAX_STEPS, 2_000_000), finish=not args.quick)
        report["gradient_drive"] = driven
        early = next(row for row in driven["samples"] if row["step"] == 10_000)
        dx = early["center_angstrom"][0]-driven["initial"]["center_angstrom"][0]
        difference = dx-control["displacement_angstrom"][0]
        require(dx > 2.0 and difference > 2.0, "Static gradient fails the matched-seed motion causality control")
        report["drive_causality_10000_steps"] = {"driven_dx_angstrom": dx, "control_dx_angstrom": control["displacement_angstrom"][0], "difference_angstrom": difference}
        if not args.quick:
            require(driven["finished"], "Default molecule did not complete the 100 nm static-gradient track")
        print(f"Drive causality passed; finish={driven['finished']}, final x={driven['final']['center_angstrom'][0]:.2f} A", flush=True)
        report["static_field_nve"] = nve_trial()
        print(f"NVE absolute drift: {report['static_field_nve']['max_absolute_energy_drift_kcal_mol']:.6g} kcal/mol", flush=True)
        report["mode_switch"] = mode_switch_trial()
        print("Switch to the old movable tip and return to gradient passed.", flush=True)
        report["lateral_spread"] = lateral_trials(10.0 if args.quick else 20.0)
        report["passed"] = True
    except Exception as error:
        report["error"] = str(error)
        report["traceback"] = traceback.format_exc()
        print(f"Gradient verification failed: {error}", file=sys.stderr, flush=True)
    finally:
        report["wall_seconds"] = time.perf_counter()-started
        save_json(path, report)
        print(f"Report: {path}", flush=True)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
