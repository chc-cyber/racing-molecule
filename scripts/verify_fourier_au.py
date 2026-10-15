"""Static native-SAIP acceptance and isolated real MD for the Fourier proxy.

Reuses the established geometry and MD protocol without altering historical
source files or game services. Timing methods run sequentially in one process.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import traceback

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "outputs/fourier-au-check"

import numpy as np
from scripts import verify_periodic_au as established
from molecular_car.periodic_au_racing import PeriodicAuRace
from molecular_car.fourier_au_racing import FourierAuRace, ContinuumComparisonRace
from molecular_car.surface_racing import SurfaceRace, benzene_input
from molecular_car.explicit_au import FrozenAuOracle
from molecular_car.fourier_au import create_fourier_surface


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def factory(atoms, positions=None, *, use_cache=True):
    if not use_cache:
        raise ValueError("Use the original periodic verification for its direct lattice sum.")
    return create_fourier_surface(atoms, positions)


def extra_static_checks():
    """Independent holdouts cover chart transitions and unrelated orientations."""
    prepared, geometry, _ = benzene_input()
    surface = create_fourier_surface(prepared.atoms, geometry + [0., 0., 3.4])
    random = np.random.default_rng(140516)
    d = 4.08 / np.sqrt(2)
    cell = np.array([d, np.sqrt(3) * d])
    rows, fd_rows = [], []
    with FrozenAuOracle(geometry, [a["element"] for a in prepared.atoms],
                        style="base", requested_xy=(40., 39.9)) as oracle:
        for index in range(32):
            # QR provides non-grid orientations, including normals close to
            # vertical, chart boundaries, and almost parallel to the surface.
            rotation, _ = np.linalg.qr(random.normal(size=(3, 3)))
            if np.linalg.det(rotation) < 0:
                rotation[:, 0] *= -1
            coords = geometry @ rotation.T
            coords[:, 2] += random.uniform(2.8, 9.) - coords[:, 2].min()
            coords[:, :2] += random.uniform(0, 1, 2) * cell
            expected = oracle.evaluate(coords)
            energy, force = surface.energy_and_forces_eV(coords)
            rows.append({"case": index, "min_z_A": float(coords[:, 2].min()),
                         "normal": rotation[:, 2].tolist(),
                         "energy_error_eV": energy - expected["energy_eV"],
                         "max_force_error_eV_A": float(np.max(np.abs(force - expected["forces_eV_A"]))),
                         "native_energy_eV": expected["energy_eV"]})
        # Cartesian FD tests the potential composition, not just interpolator
        # partial derivatives, including normal gradients on neighbouring atoms.
        for tilt, low_z in ((0., 2.8), (30., 3.4), (60., 3.4),
                            (62., 3.4), (65., 3.4), (85., 3.4), (25., 8.0)):
            angle = np.deg2rad(tilt)
            rotation = np.array([[1., 0., 0.], [0., np.cos(angle), -np.sin(angle)],
                                 [0., np.sin(angle), np.cos(angle)]])
            coords = geometry @ rotation.T
            coords[:, 2] += low_z - coords[:, 2].min()
            coords[:, :2] += .173 * cell
            energy, force = surface.energy_and_forces_eV(coords)
            numeric = np.zeros_like(force)
            epsilon = 1e-5
            for atom in range(len(coords)):
                for axis in range(3):
                    plus, minus = coords.copy(), coords.copy()
                    plus[atom, axis] += epsilon
                    minus[atom, axis] -= epsilon
                    numeric[atom, axis] = -(surface.energy_and_forces_eV(plus)[0]
                                           - surface.energy_and_forces_eV(minus)[0]) / (2 * epsilon)
            fd_rows.append({"tilt_deg": tilt, "min_z_A": low_z,
                            "maximum_force_derivative_error_eV_A": float(np.max(np.abs(numeric - force)))})
    # A numerical H-normal fixture (not a chemically proposed racing molecule)
    # exercises the full hydrogen direction correction and its chain forces.
    h_elements = ["H", "C", "C"]
    h_geometry = np.array([[0., 0., 0.], [1.1, 0., 0.], [0., 1.1, 0.]])
    h_surface = create_fourier_surface(h_elements, h_geometry + [0., 0., 3.7])
    with FrozenAuOracle(h_geometry, h_elements, style="base", requested_xy=(40., 39.9)) as oracle:
        for tilt in (0., 35., 60., 85.):
            angle = np.deg2rad(tilt)
            rotation = np.array([[1., 0., 0.], [0., np.cos(angle), -np.sin(angle)],
                                 [0., np.sin(angle), np.cos(angle)]])
            coords = h_geometry @ rotation.T + [*(.29 * cell), 3.6]
            expected = oracle.evaluate(coords)
            energy, force = h_surface.energy_and_forces_eV(coords)
            rows.append({"case": f"hydrogen-normal-{tilt:g}", "atoms": 3,
                         "min_z_A": float(coords[:, 2].min()), "normal": rotation[:, 2].tolist(),
                         "energy_error_eV": energy - expected["energy_eV"],
                         "max_force_error_eV_A": float(np.max(np.abs(force - expected["forces_eV_A"]))),
                         "native_energy_eV": expected["energy_eV"]})
            numeric = np.zeros_like(force)
            epsilon = 1e-5
            for atom in range(len(coords)):
                for axis in range(3):
                    plus, minus = coords.copy(), coords.copy()
                    plus[atom, axis] += epsilon
                    minus[atom, axis] -= epsilon
                    numeric[atom, axis] = -(h_surface.energy_and_forces_eV(plus)[0]
                                           - h_surface.energy_and_forces_eV(minus)[0]) / (2 * epsilon)
            fd_rows.append({"fixture": "hydrogen-two-neighbour", "tilt_deg": tilt,
                            "min_z_A": 3.6,
                            "maximum_force_derivative_error_eV_A": float(np.max(np.abs(numeric - force)))})
    max_force = max(row["max_force_error_eV_A"] for row in rows)
    max_energy_per_atom = max(abs(row["energy_error_eV"]) / row.get("atoms", len(geometry)) for row in rows)
    max_fd = max(row["maximum_force_derivative_error_eV_A"] for row in fd_rows)
    result = {"cases": rows, "cartesian_derivative_checks": fd_rows,
              "max_force_error_eV_A": max_force, "max_energy_error_per_atom_eV": max_energy_per_atom,
              "max_force_derivative_error_eV_A": max_fd, "surface": surface.metadata(),
              "passed": max_force <= .002 and max_energy_per_atom <= 1e-4 and max_fd <= 5e-6}
    write_json(OUT / "holdout-report.json", result)
    if not result["passed"]:
        raise AssertionError(f"Independent Fourier holdouts failed: F={max_force}, E/atom={max_energy_per_atom}, FD={max_fd}")
    print(f"PASS independent {len(rows)} native holdouts and {len(fd_rows)} Cartesian derivative checks", flush=True)
    return result


def paired_md_checks(steps):
    # Temporarily changing imported runner bindings is confined to this isolated
    # verification process, not the original source or any running game.
    reports = {}
    for label, constructor in (("continuum", ContinuumComparisonRace),
                               ("periodic", PeriodicAuRace), ("fourier", FourierAuRace)):
        folder = OUT / label
        folder.mkdir(parents=True, exist_ok=True)
        established.OUT = folder
        established.PeriodicAuRace = constructor
        reports[label] = established.run_md_checks(steps)
        reports[label]["source_sha256"] = source_hashes()
        write_json(folder / "md-report.json", reports[label])
    established.OUT = OUT
    cost = {label: float(np.median([trial["median_wall_seconds_per_step"]
                                  for trial in report["trials"]])) * 1e6
            for label, report in reports.items()}
    result = {"reports": reports, "median_us_per_step": cost,
              "fourier_speedup_vs_periodic": cost["periodic"] / cost["fourier"],
              "fourier_cost_vs_continuum": cost["fourier"] / cost["continuum"],
              "protocol": "single thread; 12-atom benzene; 600 K; 0.75 strength; 0.5 fs; 350 fs damping; three seeds; 1000 warmup + 3 timed batches; preparation/sampling excluded"}
    write_json(OUT / "paired-md-report.json", result)
    return result


def paired_motion_checks():
    """Short trajectories check dynamic contact without a claim of convergence."""
    reports = {}
    for label, constructor in (("continuum", ContinuumComparisonRace),
                               ("periodic", PeriodicAuRace), ("fourier", FourierAuRace)):
        trials = []
        for seed in (98245165, 2718281, 1234567):
            race = constructor(strength=.75, temperature=600, seed=seed)
            try:
                samples = [race.sample()]
                for _ in range(100):
                    race.advance(200)
                    sample = race.sample()
                    if label == "fourier":
                        sample["surface_evaluation"] = race.surface.last_evaluation_counts.copy()
                    samples.append(sample)
                trial = {"description": race.description(), "samples": samples,
                         "monitor": race.monitor, "source_sha256": source_hashes()}
                trials.append(trial)
                write_json(OUT / label / f"motion-{seed}.json", trial)
                print(f"PASS motion {label}, seed {seed}: 10 ps", flush=True)
            finally:
                race.close()
        reports[label] = trials
    result = {"reports": reports, "duration_ps": 10., "sample_interval_ps": .1,
              "protocol": "12-atom benzene; 600 K; strength 0.75; 0.5 fs; 350 fs Langevin; three paired initial velocity seeds; each model vertically relaxed; no pre-production warmup",
              "limitations": "Three short trajectories are a dynamic sanity check, not converged race times or diffusion statistics. Later stochastic trajectories may diverge; pointwise agreement is not required.",
              "passed": True}
    write_json(OUT / "paired-motion-report.json", result)
    return result


def source_hashes():
    names = ("molecular_car/fourier_au.py", "molecular_car/fourier_au_kernel.py",
             "molecular_car/fourier_au_racing.py", "molecular_car/periodic_au.py",
             "molecular_car/periodic_au_kernel.py", "molecular_car/periodic_au_racing.py",
             "molecular_car/explicit_au.py", "molecular_car/surface_racing.py",
             "molecular_car/fields.py", "molecular_car/molecule.py",
             "scripts/verify_fourier_au.py", "scripts/verify_periodic_au.py",
             "examples/benzene.mol", "outputs/au-surface-check/au111-relaxed.xyz")
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in names}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("static", "md", "motion", "all"), default="all")
    parser.add_argument("--steps", type=int, default=2000)
    args = parser.parse_args()
    if args.steps < 100 or args.steps % 100:
        parser.error("--steps must be a positive multiple of 100")
    OUT.mkdir(parents=True, exist_ok=True)
    result = {"status": "running", "stage": args.stage}
    start = time.perf_counter()
    try:
        established.OUT = OUT
        established.create_periodic_surface = factory
        if args.stage in ("static", "all"):
            result["static"] = established.static_checks()
            result["holdouts"] = extra_static_checks()
        if args.stage in ("md", "all"):
            result["md"] = paired_md_checks(args.steps)
        if args.stage in ("motion", "all"):
            result["motion"] = paired_motion_checks()
        result.update(status="passed", source_sha256=source_hashes())
    except BaseException as error:
        result.update(status="failed", error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        result["wall_seconds"] = time.perf_counter() - start
        write_json(OUT / f"report-{args.stage}.json", result)


if __name__ == "__main__":
    main()
