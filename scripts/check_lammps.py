"""Exercise the official API and save a small, reproducible trajectory."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

from lammps_runtime import PROJECT_ROOT, open_lammps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lammps-home", help="Installation containing bin and Python folders")
    parser.add_argument("--expected-version", type=int, help="Required LAMMPS numeric version")
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "outputs" / "lammps-check")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    report_path = output / "report.json"
    trajectory_path = output / "trajectory.csv"
    # A failed rerun must not leave an earlier successful report looking current.
    report_path.write_text(json.dumps({"status": "running"}), encoding="utf-8")

    try:
        with open_lammps(args.lammps_home, args.expected_version) as (simulation, runtime):
            required_styles = {
                "atom": ["full"],
                "bond": ["harmonic"],
                "angle": ["harmonic"],
                "pair": ["lj/cut/coul/cut"],
                "fix": ["nve", "efield", "external"],
            }
            styles = {
                f"{kind}:{style}": simulation.has_style(kind, style)
                for kind, names in required_styles.items() for style in names
            }
            missing = [name for name, present in styles.items() if not present]
            if missing:
                raise RuntimeError(f"Missing required styles: {missing}")

            data_file = (PROJECT_ROOT / "examples" / "interface_check.data").as_posix()
            simulation.commands_string(f'''units real
atom_style full
boundary f f f
read_data "{data_file}"
pair_style lj/cut/coul/cut 8.0
pair_coeff 1 1 0.07 3.4
bond_style harmonic
bond_coeff 1 100.0 1.4
angle_style harmonic
angle_coeff 1 30.0 180.0
special_bonds lj/coul 0.0 0.0 0.0
velocity all set 0.001 0.0 0.0 units box
fix integration all nve
timestep 0.25
thermo 1000
thermo_style custom step atoms pe ke etotal
run 0''')
            atom_count = int(simulation.get_natoms())
            if atom_count != 4:
                raise RuntimeError(f"Expected four atoms, got {atom_count}.")
            energies = []
            centers = []
            with trajectory_path.open("w", newline="", encoding="utf-8") as trajectory:
                writer = csv.writer(trajectory)
                writer.writerow(["step", "time_fs", "atom_id", "x_angstrom", "y_angstrom", "z_angstrom"])
                for sample in range(51):
                    if sample:
                        simulation.command("run 20 pre no post no")
                    step = int(simulation.get_thermo("step"))
                    coordinates = list(simulation.gather_atoms("x", 1, 3)[:atom_count * 3])
                    energy = float(simulation.get_thermo("etotal"))
                    if not all(math.isfinite(value) for value in [*coordinates, energy]):
                        raise RuntimeError("Non-finite coordinates or energy.")
                    energies.append(energy)
                    centers.append(sum(coordinates[0::3]) / atom_count)
                    for index in range(atom_count):
                        writer.writerow([step, step * 0.25, index + 1, *coordinates[index * 3:index * 3 + 3]])

            drift = max(abs(energy - energies[0]) for energy in energies)
            drift_limit = max(1e-6, abs(energies[0]) * 1e-3)
            displacement = centers[-1] - centers[0]
            expected_displacement = 0.001 * 1000 * 0.25
            if int(simulation.get_thermo("step")) != 1000:
                raise RuntimeError("The simulation did not reach step 1000.")
            if drift > drift_limit:
                raise RuntimeError(f"Energy drift {drift} exceeds {drift_limit} kcal/mol.")
            if abs(displacement - expected_displacement) > 1e-6:
                raise RuntimeError(f"Unexpected center-of-mass displacement: {displacement} angstrom.")
            report = {
                "status": "passed",
                "model": "Four-particle bonded interface check; not a predictive nanocar model",
                "runtime": runtime,
                "available_styles": styles,
                "atoms": atom_count,
                "bonds": 3,
                "angles": 2,
                "steps": 1000,
                "timestep_fs": 0.25,
                "simulation_time_fs": 250.0,
                "samples": 51,
                "maximum_energy_drift_kcal_per_mol": drift,
                "energy_drift_limit_kcal_per_mol": drift_limit,
                "center_of_mass_displacement_angstrom": displacement,
                "trajectory": str(trajectory_path),
            }
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"PASS: LAMMPS {runtime['lammps_version']}, Python {runtime['python_version']} ({runtime['python_bits']}-bit)")
        print(f"4 atoms, 3 bonds, 1000 steps, 51 coordinate samples; COM displacement {displacement:.6f} angstrom")
        print(f"Report: {report_path}")
        print(f"Trajectory: {trajectory_path}")
        return 0
    except Exception as error:
        report_path.write_text(json.dumps({"status": "failed", "error": str(error)}, indent=2), encoding="utf-8")
        print(f"FAIL: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
