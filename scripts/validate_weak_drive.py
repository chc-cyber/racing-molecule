"""Independent coordinate, checkpoint and continuation audit of the weak scan."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/weak-drive-20nm"
sys.path.insert(0, str(ROOT))
from scripts.scan_weak_drive import trajectory_metrics, write_json, fingerprint
from molecular_car.molecule import prepare_molecule
from molecular_car.fields import energy_and_forces
from molecular_car.periodic_au import parameter_file


def signature(row):
    config = row["config"]
    return tuple(config[key] for key in ("surface", "temperature_K", "strength", "corrugation", "seed"))


def main():
    reports = [json.loads((OUT / "screen-report.json").read_text(encoding="utf-8"))]
    if (OUT / "extended-report.json").exists():
        reports.append(json.loads((OUT / "extended-report.json").read_text(encoding="utf-8")))
    audit = {"status": "running", "trials": [], "continuation_prefixes": [],
             "wall_zone_definition": "Any atom x < -9 A or x > 209 A, or |y| > 99 A; 1 A inside nominal soft wall",
             "initial_energy_sampling": "Production t=0 energy is the left-limit zero-drive thermal state. Driver controls are switched immediately after preparation; existing pre-yes run setup refreshes forces and energy before the first production MD step. No energy-conservation inference is made across this switch.",
             "qualification": "Checks sampled geometry and coordinates; callback already rejects nonfinite forces at every MD step"}
    coordinates = {}
    parameter_path = parameter_file()
    relaxed_path = ROOT / "outputs/au-surface-check/au111-relaxed.xyz"
    parameter_sha = hashlib.sha256(parameter_path.read_bytes()).hexdigest()
    gold = np.loadtxt(relaxed_path, skiprows=2, usecols=(1, 2, 3))
    heights = gold[:, 2].reshape(6, 480).mean(axis=1)
    heights -= heights[-1]
    audit["physical_input_sha256"] = {
        "CHAu.ILP": parameter_sha,
        "au111-relaxed.xyz": hashlib.sha256(relaxed_path.read_bytes()).hexdigest(),
    }
    prepared = prepare_molecule((ROOT / "examples/tetraphenylbenzene.mol").read_text(encoding="utf-8"), "tetraphenylbenzene.mol")
    charges = np.array([atom["charge"] for atom in prepared.atoms])
    alphas = np.array([atom["polarizability"] for atom in prepared.atoms])
    try:
        for report in reports:
            assert report["status"] == "passed" and len(report["trials"]) == report["planned_trials"]
            assert report["source_sha256"] == fingerprint()
            assert report["parameter_file_sha256"] == parameter_sha
            assert all(group["n"] == 3 for group in report["summary"])
            for row in report["trials"]:
                assert row["status"] == "passed"
                if row["config"]["surface"] == "periodic-au":
                    assert row["runtime"]["surface"]["parameter_sha256"] == parameter_sha
                    np.testing.assert_array_equal(row["runtime"]["surface"]["layer_heights_angstrom"], heights)
                path = OUT / row["coordinates_file"]
                assert hashlib.sha256(path.read_bytes()).hexdigest() == row["coordinates_sha256"]
                with np.load(path, allow_pickle=False) as data:
                    time = data["time_ps"]
                    positions = data["positions_angstrom"]
                    masses = data["masses_u"]
                    bonds = data["bonds_zero_based"]
                assert positions.shape == (len(time), 52, 3) and len(masses) == 52
                assert np.isfinite(positions).all() and np.all(np.diff(time) > 0)
                np.testing.assert_allclose(np.diff(time), .1, rtol=0, atol=3e-14)
                center = np.sum(positions * masses[None, :, None], axis=1) / masses.sum() / 10
                stored = np.array([sample["center_nm"] for sample in row["samples"]])
                np.testing.assert_allclose(center, stored, rtol=0, atol=1e-12)
                metrics = trajectory_metrics(time, center)
                np.testing.assert_allclose(metrics["displacement_nm"], row["metrics"]["displacement_nm"], rtol=0, atol=1e-12)
                for name in ("forward_fraction", "backward_fraction", "pause_fraction"):
                    assert metrics[name] == row["metrics"][name]
                lengths = np.linalg.norm(positions[:, bonds[:, 0], :] - positions[:, bonds[:, 1], :], axis=2)
                np.testing.assert_allclose([lengths.min(), lengths.max()], row["bond_length_range_angstrom"], rtol=0, atol=1e-12)
                assert .6 < lengths.min() and lengths.max() < 2.5
                left = np.any(positions[:, :, 0] < -9., axis=1)
                right = np.any(positions[:, :, 0] > 209., axis=1)
                side = np.any(np.abs(positions[:, :, 1]) > 99., axis=1)
                # Coordinates/velocities do not change at the instantaneous
                # field switch. Surface, walls and UFF therefore cancel in
                # this independent potential-energy correction of the t=0
                # left-limit raw sample. Preserve the original measurement.
                tip = {"x": 0., "y": 0., "z": 7., "strength": 0.}
                before_energy, _ = energy_and_forces(positions[0], charges, alphas, tip, row["config"]["corrugation"],
                    drive={"mode": "gradient", "strength": 0.})
                after_energy, _ = energy_and_forces(positions[0], charges, alphas, tip, row["config"]["corrugation"],
                    drive={"mode": "gradient", "strength": row["config"]["strength"]})
                audit["trials"].append({"key": row["key"], "sample_count": len(time),
                    "left_wall_zone_fraction": float(left.mean()), "right_wall_zone_fraction": float(right.mean()),
                    "side_wall_zone_fraction": float(side.mean()),
                    "COM_height_range_angstrom": [float(center[:, 2].min()*10), float(center[:, 2].max()*10)],
                    "min_atom_z_range_angstrom": [float(positions[:, :, 2].min(axis=1).min()), float(positions[:, :, 2].min(axis=1).max())],
                    "initial_driver_potential_change_kcal_mol": float(after_energy-before_energy),
                    "initial_energy_after_drive_switch_kcal_mol": float(row["samples"][0]["energy_kcal_mol"] + after_energy-before_energy),
                    "max_mass_weighted_COM_error_nm": float(np.max(np.abs(center-stored)))})
                coordinates[(report["phase"], signature(row))] = positions
        if len(reports) > 1:
            for extended in reports[1]["trials"]:
                key = signature(extended)
                original = coordinates[("screen", key)]
                longer = coordinates[("extended", key)]
                shared = min(len(original), len(longer))
                np.testing.assert_array_equal(original[:shared], longer[:shared])
                audit["continuation_prefixes"].append({"key": extended["key"],
                    "shared_samples": shared, "bitwise_equal": True})
        audit["status"] = "passed"
        inputs = OUT / "inputs"
        inputs.mkdir(exist_ok=True)
        for source in (parameter_path, relaxed_path):
            (inputs / source.name).write_bytes(source.read_bytes())
    except BaseException as error:
        audit["status"] = "failed"
        audit["error"] = repr(error)
        raise
    finally:
        write_json(OUT / "validation.json", audit)
    print(json.dumps({"status": audit["status"], "trials": len(audit["trials"]),
                     "bitwise_continuations": len(audit["continuation_prefixes"])}, ensure_ascii=False))


if __name__ == "__main__":
    main()
