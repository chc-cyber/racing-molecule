"""Reproducible checks of the actual LAMMPS prototype, without the web UI.

Run with ``.venv/Scripts/python.exe scripts/verify_prototype.py``.  The JSON
report records the exact engine and effective game-model parameters.  These
are integration/causality checks, not validation of experimental Au/STM data.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import io
import json
from pathlib import Path
import sys
import time
import traceback
import unittest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import rdkit

from molecular_car import fields
from molecular_car.engine import (
    DEFAULT_SETTINGS,
    DEFAULT_TIP,
    FRAME_STEPS,
    TIMESTEP_FS,
    Simulation,
)
from molecular_car.molecule import EFFECTIVE_POLARIZABILITY, prepare_molecule


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def molecular_checks() -> list[dict]:
    results = []
    for filename, expected in (("benzene.mol", 12), ("tetraphenylbenzene.mol", 52)):
        text = (PROJECT_ROOT / "examples" / filename).read_text(encoding="utf-8")
        molecule = prepare_molecule(text, filename)
        require(len(molecule.atoms) == expected, f"Unexpected atom count in {filename}")
        require({atom["element"] for atom in molecule.atoms} == {"C", "H"}, "Fixture composition changed")
        require(molecule.positions.shape == (expected, 3), "Coordinate order/shape mismatch")
        require(np.isfinite(molecule.positions).all(), "Nonfinite molecular positions")
        require(abs(sum(atom["charge"] for atom in molecule.atoms)) < 1e-8, "Neutral fixture charge not conserved")
        require(abs(molecule.positions[:, 2].min() - 2.9) < 1e-10, "Molecule starting height changed")
        # Check the callback's molecular-force sign and units against energy.
        coordinates = molecule.positions.ravel()
        direction = np.random.default_rng(17).normal(size=coordinates.size)
        direction /= np.linalg.norm(direction)
        gradient = np.asarray(molecule.forcefield.CalcGrad(tuple(coordinates)))
        epsilon = 1e-5
        plus = molecule.forcefield.CalcEnergy(tuple(coordinates + epsilon * direction))
        minus = molecule.forcefield.CalcEnergy(tuple(coordinates - epsilon * direction))
        difference = abs((plus - minus) / (2 * epsilon) - float(gradient @ direction))
        require(difference < 1e-5, "UFF analytic gradient disagrees with finite differences")
        require(np.linalg.norm(gradient.reshape(-1, 3).sum(axis=0)) < 1e-8, "Molecular internal forces have a net force")
        results.append({
            "filename": filename,
            "atoms": expected,
            "bonds": len(molecule.bonds),
            "charge_sum_e": float(sum(atom["charge"] for atom in molecule.atoms)),
            "uff_energy_kcal_mol": float(molecule.forcefield.CalcEnergy(tuple(coordinates))),
            "uff_directional_gradient_error": difference,
        })
    return results


def field_checks() -> dict:
    log = io.StringIO()
    suite = unittest.defaultTestLoader.discover(str(PROJECT_ROOT / "tests"), pattern="test_fields.py")
    result = unittest.TextTestRunner(stream=log, verbosity=1).run(suite)
    require(result.testsRun > 0, "No field tests discovered")
    require(result.wasSuccessful(), log.getvalue())
    return {"tests": result.testsRun, "passed": result.wasSuccessful()}


def _sample(simulation: Simulation) -> dict:
    frame = simulation.frame()
    coordinates = np.asarray(frame["positions"]).reshape(simulation.n, 3)
    bond_lengths = [
        np.linalg.norm(coordinates[bond["a"] - 1] - coordinates[bond["b"] - 1])
        for bond in simulation.prepared.bonds
    ]
    require(np.isfinite(frame["energy"]), "LAMMPS energy is not finite")
    require(np.isfinite(coordinates).all(), "LAMMPS coordinates are not finite")
    require(min(bond_lengths) > 0.6 and max(bond_lengths) < 2.5, "Fixture molecular bond geometry became unstable")
    return {
        "step": frame["step"],
        "time_fs": frame["time_fs"],
        "center_angstrom": frame["center"],
        "energy_kcal_mol": frame["energy"],
        "min_bond_angstrom": float(min(bond_lengths)),
        "max_bond_angstrom": float(max(bond_lengths)),
        "min_z_angstrom": float(coordinates[:, 2].min()),
        "max_z_angstrom": float(coordinates[:, 2].max()),
        "internal_temperature_K": frame["internal_temperature"],
    }


def trial(filename: str, steps: int, mode: str, *, thermostat: bool = True, tip_x: float | None = None) -> dict:
    """Each trial starts from identical fixture, velocity and bath seeds."""
    content = (PROJECT_ROOT / "examples" / filename).read_text(encoding="utf-8")
    simulation = None
    started = time.perf_counter()
    try:
        # This regression suite deliberately checks the retained original tip mode.
        simulation = Simulation(content, filename, thermostat=thermostat,
                                drive={"mode": "tip", "strength": DEFAULT_TIP["strength"]})
        simulation.autopilot = mode == "autopilot"
        if mode == "undriven":
            simulation.drive_enabled = False
            simulation.tip["strength"] = 0.0
        elif mode == "manual":
            simulation.tip.update({"x": tip_x, "y": 0.0, "z": DEFAULT_TIP["z"], "strength": DEFAULT_TIP["strength"]})
        simulation.lmp.command("run 0")
        samples = [_sample(simulation)]
        initial_energy = samples[0]["energy_kcal_mol"]

        # Ensure external potential energy really reaches LAMMPS thermodynamics.
        effective_tip = dict(simulation.tip)
        if not simulation.drive_enabled:
            effective_tip["strength"] = 0.0
        coordinates = simulation.coordinates()
        uff_energy = simulation.prepared.forcefield.CalcEnergy(coordinates.ravel().tolist())
        external_energy, _ = fields.energy_and_forces(
            coordinates, simulation.charges, simulation.alphas,
            effective_tip, simulation.settings["corrugation"],
        )
        potential = float(uff_energy + external_energy)
        reconstructed_energy = potential + float(simulation.lmp.get_thermo("ke"))
        energy_accounting_error = abs(initial_energy - reconstructed_energy)
        require(energy_accounting_error < 1e-7, "External potential missing from LAMMPS total energy")

        for offset in range(0, steps, FRAME_STEPS):
            simulation.advance(min(FRAME_STEPS, steps - offset))
            samples.append(_sample(simulation))
        initial_center = np.asarray(samples[0]["center_angstrom"])
        displacement = np.asarray(samples[-1]["center_angstrom"]) - initial_center
        energies = np.asarray([sample["energy_kcal_mol"] for sample in samples])
        maximum_drift = float(np.max(np.abs(energies - initial_energy)))
        energy_scale = max(abs(initial_energy), abs(potential), 1.0)
        return {
            "filename": filename,
            "mode": mode,
            "atoms": simulation.n,
            "steps": steps,
            "time_fs": steps * TIMESTEP_FS,
            "thermostat": thermostat,
            "settings": dict(simulation.settings),
            "initial_tip": effective_tip,
            "runtime": dict(simulation.runtime),
            "initial": samples[0],
            "final": samples[-1],
            "displacement_angstrom": displacement.tolist(),
            "max_absolute_energy_drift_kcal_mol": maximum_drift,
            "energy_scale_kcal_mol": energy_scale,
            "max_relative_energy_drift": maximum_drift / energy_scale,
            "initial_energy_accounting_error_kcal_mol": energy_accounting_error,
            "trajectory_samples": len(samples),
            "bond_length_range_angstrom": [
                min(sample["min_bond_angstrom"] for sample in samples),
                max(sample["max_bond_angstrom"] for sample in samples),
            ],
            "wall_seconds": time.perf_counter() - started,
        }
    finally:
        if simulation is not None:
            simulation.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "outputs" / "prototype-check" / "report.json")
    args = parser.parse_args()
    output = args.output.resolve()
    started = time.perf_counter()
    report = {
        "passed": False,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model": "prototype-uff-soft-tip-v1",
        "purpose": "Real LAMMPS integration, finite gradients, drive causality and undriven energy stability",
        "versions": {"python": sys.version.split()[0], "numpy": np.__version__, "rdkit": rdkit.__version__},
        "parameters": {
            "units": {"length": "angstrom", "time": "fs", "energy": "kcal/mol"},
            "timestep_fs": TIMESTEP_FS,
            "frame_steps": FRAME_STEPS,
            "default_settings": dict(DEFAULT_SETTINGS),
            "default_tip": dict(DEFAULT_TIP),
            "effective_polarizability": dict(EFFECTIVE_POLARIZABILITY),
            "field_constants": {key: value for key, value in vars(fields).items() if key.isupper() and isinstance(value, (int, float))},
            "velocity_seed": 17329,
            "langevin_seed": 98245165,
        },
        "trials": [],
    }
    try:
        report["molecule_checks"] = molecular_checks()
        report["field_checks"] = field_checks()
        print("Molecular preparation and analytic field checks passed.", flush=True)
        undriven = trial("tetraphenylbenzene.mol", 10_000, "undriven")
        report["trials"].append(undriven)
        print(f"Undriven 52-atom control: dx={undriven['displacement_angstrom'][0]:.3f} A", flush=True)
        driven = trial("tetraphenylbenzene.mol", 10_000, "autopilot")
        report["trials"].append(driven)
        difference = driven["displacement_angstrom"][0] - undriven["displacement_angstrom"][0]
        require(driven["displacement_angstrom"][0] > 10.0 and difference > 10.0, "Drive does not cause significant forward motion versus the matched control")
        report["drive_causality"] = {"driven_minus_undriven_dx_angstrom": difference, "minimum_required_difference_angstrom": 10.0}
        print(f"Automatic drive: dx={driven['displacement_angstrom'][0]:.3f} A; advantage={difference:.3f} A", flush=True)

        nve = trial("tetraphenylbenzene.mol", 2_000, "undriven", thermostat=False)
        report["trials"].append(nve)
        require(nve["max_relative_energy_drift"] < 1e-3, "Undriven NVE energy drift exceeds 0.1% of the initial energy scale")
        report["energy_conservation"] = {"max_relative_drift": nve["max_relative_energy_drift"], "maximum_allowed": 1e-3}
        print(f"Undriven NVE: maximum relative drift={nve['max_relative_energy_drift']:.3g}", flush=True)

        left = trial("benzene.mol", 2_000, "manual", tip_x=-6.0)
        report["trials"].append(left)
        right = trial("benzene.mol", 2_000, "manual", tip_x=6.0)
        report["trials"].append(right)
        left_dx, right_dx = left["displacement_angstrom"][0], right["displacement_angstrom"][0]
        require(left_dx < -0.1 and right_dx > 0.1 and right_dx-left_dx > 0.5, "Manual tip position does not steer the molecule in opposite directions")
        report["manual_tip_causality"] = {"left_dx_angstrom": left_dx, "right_dx_angstrom": right_dx, "right_minus_left_dx_angstrom": right_dx-left_dx}
        print(f"Manual tip: left dx={left_dx:.3f} A; right dx={right_dx:.3f} A", flush=True)
        report["passed"] = True
    except Exception as error:
        report["error"] = str(error)
        report["traceback"] = traceback.format_exc()
        print(f"Prototype check failed: {error}", file=sys.stderr, flush=True)
    finally:
        report["wall_seconds"] = time.perf_counter() - started
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
        print(f"Report: {output}", flush=True)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
