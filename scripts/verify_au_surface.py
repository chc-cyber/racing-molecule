"""Real LAMMPS Au(111)/benzene reference, isolated from the game model.

The bundled Foiles Au EAM describes only gold. Benzene/gold uses bundled
CHAu.ILP SAIP (Ouyang, Hod and Guerra, JCTC 17, 7215, 2021). Gold and benzene
have different molecule IDs; all gold layers share the same ID. UFF molecular
coordinates/topology are held fixed during static scans; its constant energy
is recorded separately and is not added to the adsorption energy.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import json
from pathlib import Path
import sys
import time
import traceback
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.lammps_runtime import open_lammps
from molecular_car.au_surface import (build_au111, layer_statistics,
    potential_metadata, write_atomic_data, write_xyz, write_height_svg,
    minimum_image_xy)
from molecular_car.molecule import prepare_molecule
from molecular_car.fields import energy_and_forces

OUTPUT = ROOT / "outputs" / "au-surface-check"
EV_TO_KCAL_MOL = 23.0605478306
EAM_DOI = "https://doi.org/10.1103/PhysRevB.33.7983"
SAIP_PAPER = "https://www.tau.ac.il/~odedhod/papers/paper84.pdf"
SAIP_DOI = "https://doi.org/10.1021/acs.jctc.1c00622"


def gather(lmp, name, count):
    return np.ctypeslib.as_array(lmp.gather_atoms(name, 1, 3),
                                shape=(count * 3,)).copy().reshape(count, 3)


def finite(values, label):
    if not np.isfinite(np.asarray(values)).all():
        raise AssertionError(f"Non-finite {label}")


def potentials(runtime):
    installation = Path(runtime["lammps_library"]).parents[1]
    return installation / "Potentials" / "Au_u3.eam", installation / "Potentials" / "CHAu.ILP"


def displacement_statistics(slab, positions, reference):
    """Actual atom excursions relative to the minimized state, without fitting."""
    displacement = minimum_image_xy(np.asarray(positions) - reference, slab.lengths_xy)
    rows = []
    for layer in range(slab.layer_count):
        delta = displacement[slab.layers == layer]
        rows.append({"layer_bottom_to_top": layer,
                     "rms_displacement_angstrom": float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
                     "rms_xyz_angstrom": np.sqrt(np.mean(delta**2, axis=0)).tolist(),
                     "maximum_atom_displacement_angstrom": float(np.linalg.norm(delta, axis=1).max())})
    return rows


def slab_sample(lmp, slab):
    positions = gather(lmp, "x", len(slab.positions))
    forces = gather(lmp, "f", len(slab.positions))
    values = {"step": int(lmp.get_thermo("step")),
              "potential_energy_eV": float(lmp.get_thermo("pe")),
              "kinetic_energy_eV": float(lmp.get_thermo("ke")),
              "temperature_K": float(lmp.get_thermo("temp")),
              "maximum_mobile_force_eV_per_angstrom": float(np.linalg.norm(
                  forces[slab.frozen_count:], axis=1).max()),
              "layer_statistics": layer_statistics(slab.with_positions(positions))}
    finite(positions, "Au coordinates")
    finite(forces, "Au forces")
    finite([values["potential_energy_eV"], values["kinetic_energy_eV"],
            values["temperature_K"]], "Au energies/temperature")
    frozen_delta = positions[:slab.frozen_count] - slab.positions[:slab.frozen_count]
    frozen_delta = minimum_image_xy(frozen_delta, slab.lengths_xy)
    values["maximum_frozen_displacement_angstrom"] = float(np.linalg.norm(frozen_delta, axis=1).max())
    if values["maximum_frozen_displacement_angstrom"] > 1e-10:
        raise AssertionError("The fixed gold layers moved.")
    return values, positions


def verify_slab(report, quick):
    slab = build_au111()
    initial_data = OUTPUT / "au111-initial.data"
    write_atomic_data(initial_data, slab)
    with open_lammps() as (lmp, runtime):
        eam, saip = potentials(runtime)
        styles = {style: lmp.has_style("pair", style)
                  for style in ("eam", "hybrid/overlay", "saip/metal")}
        if not all(styles.values()):
            raise RuntimeError(f"Required styles unavailable: {styles}")
        report.update(runtime=runtime, available_styles=styles,
                      relevant_packages=[p for p in lmp.installed_packages if p in
                                         ("MANYBODY", "INTERLAYER", "MOLECULE")],
                      potentials={"Au_EAM": potential_metadata(eam),
                                  "benzene_Au_SAIP": potential_metadata(saip)})
        setup = ["units metal", "atom_style atomic", "boundary p p f",
                 "atom_modify map array sort 0 0.0",
                 f'read_data "{initial_data.as_posix()}"',
                 "pair_style eam", f'pair_coeff * * "{eam.as_posix()}"',
                 "neighbor 2.0 bin", "neigh_modify every 1 delay 0 check yes",
                 f"group frozen id 1:{slab.frozen_count}", "group mobile subtract all frozen",
                 "fix hold frozen setforce 0 0 0", "velocity frozen set 0 0 0",
                 "compute mobile_temp mobile temp", "thermo_style custom step temp pe ke etotal fmax",
                 "thermo_modify temp mobile_temp", "thermo 1000",
                 "timestep 0.001", "run 0"]
        lmp.commands_list(setup)
        initial, _ = slab_sample(lmp, slab)
        minimization = ["min_style fire", "min_modify dmax 0.1",
                        "minimize 1e-12 1e-8 2000 20000", "run 0"]
        lmp.commands_list(minimization)
        optimized, coordinates = slab_sample(lmp, slab)
        relaxed = slab.with_positions(coordinates, top_zero=True)
        if optimized["potential_energy_eV"] > initial["potential_energy_eV"] + 1e-8:
            raise AssertionError("Gold minimization increased its potential energy.")
        if optimized["maximum_mobile_force_eV_per_angstrom"] > 1e-5:
            raise AssertionError("Gold minimization did not reach a small residual force.")
        write_atomic_data(OUTPUT / "au111-relaxed.data", relaxed)
        write_xyz(OUTPUT / "au111-relaxed.xyz", relaxed)
        write_height_svg(OUTPUT / "au111-height.svg", relaxed)
        thermal_steps = 1000 if quick else 5000
        thermal_setup = ["reset_timestep 0", "velocity mobile create 300 20261011 mom yes rot no dist gaussian",
                         "fix integration mobile nve", "fix bath mobile langevin 300 300 0.1 20261012 zero yes",
                         "run 0"]
        lmp.commands_list(thermal_setup)
        samples = []
        for _ in range(thermal_steps // 1000):
            lmp.command("run 1000 pre no post no")
            sample, thermal_coordinates = slab_sample(lmp, slab)
            sample["displacement_from_minimized_by_layer"] = displacement_statistics(slab, thermal_coordinates, coordinates)
            samples.append(sample)
            if max(row["rms_displacement_angstrom"] for row in sample["displacement_from_minimized_by_layer"]) > 1.0:
                raise AssertionError("The Au slab developed a large atomic excursion during the short 300 K check.")
            if any(not 1.5 < gap < 3.2 for gap in sample["layer_statistics"]["interlayer_gaps_angstrom"]):
                raise AssertionError("The Au slab lost its layered structure during the short 300 K check.")
        thermal_slab = slab.with_positions(thermal_coordinates, top_zero=True)
        write_xyz(OUTPUT / "au111-300K-tail.xyz", thermal_slab)
        lmp.commands_list(["unfix bath", "unfix integration", "velocity all set 0 0 0"] + minimization)
        quenched, quench_coordinates = slab_sample(lmp, slab)
        quenched["displacement_from_minimized_by_layer"] = displacement_statistics(slab, quench_coordinates, coordinates)
        write_xyz(OUTPUT / "au111-300K-quenched.xyz", slab.with_positions(quench_coordinates, top_zero=True))
        replay = setup[:-1] + minimization + thermal_setup
        replay.append(f"run {thermal_steps}")
        replay += ["unfix bath", "unfix integration", "velocity all set 0 0 0"] + minimization
        (OUTPUT / "au111-reference.in").write_text("\n".join(replay) + "\n", encoding="utf-8")
        report["slab"] = {"atom_count": len(slab.positions), "dimensions_xy_angstrom": slab.lengths_xy.tolist(),
                          "layer_count": slab.layer_count, "frozen_bottom_layers": slab.frozen_layers,
                          "frozen_atoms": slab.frozen_count, "lattice_constant_angstrom": slab.lattice_constant,
                          "nearest_neighbor_angstrom": slab.nearest_spacing,
                          "boundary": "p p f", "timestep_ps": 0.001,
                          "initial": initial, "minimized": optimized,
                          "thermal_check": {"temperature_target_K": 300, "steps": thermal_steps,
                                            "duration_ps": thermal_steps * 0.001,
                                            "samples": samples, "quenched": quenched},
                          "cross_scan_geometry": "zero-temperature minimized slab before thermal check"}
        print(f"PASS gold: {len(slab.positions)} Au; mobile fmax {optimized['maximum_mobile_force_eV_per_angstrom']:.3g} eV/A", flush=True)
        return relaxed, eam, saip


def benzene_geometry():
    path = ROOT / "examples" / "benzene.mol"
    prepared = prepare_molecule(path.read_text(encoding="utf-8"), path.name)
    if len(prepared.atoms) != 12 or {a["element"] for a in prepared.atoms} != {"C", "H"}:
        raise AssertionError("The benzene fixture must contain C6H6.")
    carbon = [i for i, atom in enumerate(prepared.atoms) if atom["element"] == "C"]
    center = np.mean(prepared.positions[carbon], axis=0)
    centered = prepared.positions - center
    _, _, vectors = np.linalg.svd(centered[carbon])
    normal = vectors[-1]
    direction = centered[carbon[1]] - centered[carbon[0]]
    ex = direction - normal * np.dot(direction, normal)
    ex /= np.linalg.norm(ex)
    ey = np.cross(normal, ex)
    planar = centered @ np.stack([ex, ey, normal]).T
    source_energy = float(prepared.forcefield.CalcEnergy(prepared.positions.ravel().tolist()))
    rotated_energy = float(prepared.forcefield.CalcEnergy(planar.ravel().tolist()))
    if abs(rotated_energy - source_energy) > 1e-7 or np.ptp(planar[:, 2]) > 1e-4:
        raise AssertionError("Rigid benzene orientation changed its UFF geometry/energy.")
    return prepared, planar, {"atoms": 12, "formula": "C6H6", "fixed_topology": prepared.bonds,
                             "UFF_internal_constant_kcal_per_mol": source_energy,
                             "input_height_span_angstrom": float(np.ptp(planar[:, 2])),
                             "source": str(path), "orientation": "carbon bond 1-2 along +x, plane parallel to Au(111)"}


class SAIPProbe:
    """One reused LAMMPS instance for rigid translations; no car/Au EAM mixing."""
    def __init__(self, slab, prepared, geometry, eam, saip, *, subsystem="both"):
        self.slab = slab
        self.prepared = prepared
        self.geometry = geometry if subsystem != "gold" else np.empty((0, 3))
        self.gold_positions = slab.positions if subsystem != "benzene" else np.empty((0, 3))
        self.eam_path = eam
        self.saip_path = saip
        self.n_gold = len(self.gold_positions)
        self.n = self.n_gold + len(self.geometry)
        self.context = open_lammps()
        self.lmp, _ = self.context.__enter__()
        try:
            lmp = self.lmp
            bounds = " ".join(f"{bound:.15g}" for axis in slab.box for bound in axis)
            lmp.commands_list(["units metal", "atom_style full", "boundary p p f",
                "newton on", "atom_modify map array sort 0 0.0",
                f"region world block {bounds} units box", "create_box 3 world",
                "mass 1 196.97", "mass 2 12.011", "mass 3 1.008"])
            positions = np.concatenate([self.gold_positions, self.geometry + [*slab.lengths_xy / 2, 3.3]])
            types = [1] * self.n_gold
            if len(self.geometry):
                types += [2 if a["element"] == "C" else 3 for a in prepared.atoms]
            count = lmp.create_atoms(self.n, list(range(1, self.n + 1)), types, positions.ravel().tolist())
            if count != self.n:
                raise AssertionError("LAMMPS did not create all gold/probe atoms.")
            lmp.commands_list(["group gold type 1", "group probe type 2 3", "set group gold mol 1",
                "set group probe mol 2", "pair_style hybrid/overlay zero 16.0 eam saip/metal 16.0 1",
                "pair_coeff * * zero",
                f'pair_coeff 1 1 eam "{eam.as_posix()}"',
                f'pair_coeff * * saip/metal "{saip.as_posix()}" Au C H',
                "neighbor 2.0 bin", "neigh_modify every 1 delay 0 check yes",
                "compute cross all pair saip/metal", "thermo_style custom step pe c_cross c_cross[1] c_cross[2]",
                "thermo 1000000", "run 0"])
        except BaseException:
            self.close()
            raise

    def close(self):
        if self.context is not None:
            self.context.__exit__(None, None, None)
            self.context = None

    def evaluate(self, center, angle_deg=0.0):
        angle = np.deg2rad(angle_deg)
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                             [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        positions = np.concatenate([self.gold_positions, self.geometry @ rotation.T + center])
        positions[:, :2] %= self.slab.lengths_xy
        raw = (ctypes.c_double * positions.size)(*positions.ravel())
        self.lmp.scatter_atoms("x", 1, 3, raw)
        self.lmp.command("run 0")
        energy = float(self.lmp.extract_compute("cross", 0, 0))
        components = self.lmp.extract_compute("cross", 0, 1)
        probe_forces = gather(self.lmp, "f", self.n)[self.n_gold:]
        finite([energy, components[0], components[1]], "SAIP energy")
        finite(probe_forces, "SAIP forces")
        if abs(energy - components[0] - components[1]) > 1e-8:
            raise AssertionError("SAIP component accounting is inconsistent.")
        return energy, np.sum(probe_forces, axis=0), probe_forces, [float(components[0]), float(components[1])]

    def grouping_checks(self, center):
        baseline = self.evaluate(center)[0]
        self.lmp.command("set group probe mol 1")
        same_id = self.evaluate(center)[0]
        self.lmp.command("set group probe mol 2")
        results = {"normal_cross_energy_eV": baseline, "same_molecule_ID_energy_eV": same_id}
        # ILP requires complete coefficient mappings even for unused atom types.
        # Separate systems therefore retain the complete Au/C/H mapping while
        # physically containing only the gold or only the benzene atoms.
        for label, subsystem in (("bare_gold", "gold"), ("isolated_benzene", "benzene")):
            isolated = SAIPProbe(self.slab, self.prepared, self.geometry,
                                 self.eam_path, self.saip_path, subsystem=subsystem)
            try:
                results[label + "_SAIP_energy_eV"] = isolated.evaluate(center)[0]
            finally:
                isolated.close()
        if any(abs(value) > 1e-10 for key, value in results.items() if key != "normal_cross_energy_eV"):
            raise AssertionError("SAIP contributed an intramolecular or intragold interaction.")
        return results


def golden_minimum(function, lower, upper, iterations=22):
    ratio = (np.sqrt(5) - 1) / 2
    left = upper - ratio * (upper - lower)
    right = lower + ratio * (upper - lower)
    e_left, e_right = function(left), function(right)
    for _ in range(iterations):
        if e_left <= e_right:
            upper, right, e_right = right, left, e_left
            left = upper - ratio * (upper - lower)
            e_left = function(left)
        else:
            lower, left, e_left = left, right, e_right
            right = lower + ratio * (upper - lower)
            e_right = function(right)
    coordinate = (lower + upper) / 2
    return float(coordinate), float(function(coordinate))


def write_rows(path, columns, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        writer.writerows(rows)


def verify_translation_gradient(probe, center):
    """Check a nonsymmetric translation, where lateral and vertical force matter."""
    center = np.asarray(center, dtype=float)
    step = 0.0001
    _, total_force, _, _ = probe.evaluate(center)
    differences = []
    for axis in range(3):
        delta = np.eye(3)[axis] * step
        energy_plus = probe.evaluate(center + delta)[0]
        energy_minus = probe.evaluate(center - delta)[0]
        differences.append(-(energy_plus - energy_minus) / (2 * step))
    gradient_error = float(np.max(np.abs(total_force - differences)))
    if gradient_error > 1e-5:
        raise AssertionError(f"SAIP force/gradient mismatch: {gradient_error} eV/A")
    return {"center_angstrom": center.tolist(), "orientation_angle_deg": 0,
            "finite_difference_step_angstrom": step, "force_eV_per_angstrom": total_force.tolist(),
            "negative_energy_gradient_eV_per_angstrom": differences,
            "maximum_absolute_error_eV_per_angstrom": gradient_error}


def verify_cross(report, slab, eam, saip, quick):
    prepared, geometry, molecule_info = benzene_geometry()
    probe = SAIPProbe(slab, prepared, geometry, eam, saip)
    try:
        top = slab.positions[slab.layers == slab.layer_count - 1]
        nearest_center = np.argmin(np.linalg.norm(top[:, :2] - slab.lengths_xy / 2, axis=1))
        anchor = top[nearest_center, :2]
        height_rows = []
        for height in np.arange(2.4, 5.40001, 0.1):
            energy, force, _, components = probe.evaluate([*anchor, height])
            height_rows.append([float(height), energy, *force, *components])
        write_rows(OUTPUT / "benzene-saip-height.csv", ["height_A", "energy_eV", "Fx_eV_A", "Fy_eV_A",
                   "Fz_eV_A", "vdW_eV", "repulsion_eV"], height_rows)
        z_atop, e_atop = golden_minimum(lambda z: probe.evaluate([*anchor, z])[0], 2.7, 4.2)
        reference = np.array([anchor[0] + 0.31 * slab.nearest_spacing,
                              anchor[1] + 0.23 * np.sqrt(3) * slab.nearest_spacing, z_atop + 0.15])
        gradient_check = verify_translation_gradient(probe, reference)
        grouping = probe.grouping_checks([*anchor, z_atop])
        d, h = slab.nearest_spacing, np.sqrt(3) * slab.nearest_spacing
        u, v = np.array([d, 0.0]), np.array([d / 2, h / 2])
        grid_count = 9 if quick else 17
        angles = [0.0] if quick else [0.0, 30.0]
        lateral = []
        all_rows = []
        for angle in angles:
            samples = []
            for fy in np.linspace(0, 1, grid_count):
                for fx in np.linspace(0, 1, grid_count):
                    xy = anchor + fx * u + fy * v
                    energy, force, _, _ = probe.evaluate([*xy, z_atop], angle)
                    samples.append((float(energy), float(fx), float(fy)))
                    all_rows.append([angle, float(fx), float(fy), *xy, z_atop, energy, *force])
            low, high = min(samples), max(samples)
            best_xy = anchor + low[1] * u + low[2] * v
            best_z, best_energy = golden_minimum(lambda z: probe.evaluate([*best_xy, z], angle)[0], 2.7, 4.2)
            edge_errors = []
            for vector in (u, v):
                edge_errors.append(abs(probe.evaluate([*(anchor + vector), z_atop], angle)[0]
                                       - probe.evaluate([*anchor, z_atop], angle)[0]))
            if max(edge_errors) > 1e-6:
                raise AssertionError("The Au(111) translational scan did not repeat across its primitive-cell edges.")
            lateral.append({"angle_deg": angle, "fixed_height_angstrom": z_atop,
                "grid_points_per_axis": grid_count, "sampled_energy_min_eV": low[0],
                "sampled_energy_max_eV": high[0], "sampled_corrugation_eV": high[0] - low[0],
                "sampled_corrugation_meV": (high[0] - low[0]) * 1000,
                "best_grid_registry_fractional": [low[1], low[2]],
                "best_registry_vertical_min_height_angstrom": best_z,
                "best_registry_vertical_min_energy_eV": best_energy,
                "cell_edge_repeat_error_eV": max(edge_errors)})
            print(f"PASS benzene {angle:g} deg: z {best_z:.4f} A; adsorption {-best_energy:.4f} eV; sampled corrugation {(high[0]-low[0])*1000:.3f} meV", flush=True)
        write_rows(OUTPUT / "benzene-saip-lateral.csv", ["angle_deg", "fx", "fy", "x_A", "y_A", "height_A",
                   "energy_eV", "Fx_eV_A", "Fy_eV_A", "Fz_eV_A"], all_rows)
        # Use the actual continuous game surface with no electric drive, at
        # central coordinates where soft boundary tails are negligible.
        charges = np.array([a["charge"] for a in prepared.atoms])
        alphas = np.array([a["polarizability"] for a in prepared.atoms])
        def continuous(center):
            energy, forces = energy_and_forces(geometry + center, charges, alphas,
                {"x": 0, "y": 0, "z": 7, "strength": 0}, 0.08,
                drive={"mode": "gradient", "strength": 0})
            return energy / EV_TO_KCAL_MOL, np.sum(forces, axis=0) / EV_TO_KCAL_MOL
        continuous_anchor = np.array([100.8, 0.0])
        c_u = np.array([2.88, 0.0])
        c_v = np.array([1.44, np.sqrt(3.0) * 1.44])
        continuous_rows = []
        for height in np.arange(2.4, 5.40001, 0.1):
            energy, force = continuous([*continuous_anchor, height])
            continuous_rows.append([float(height), energy, *force])
        cz, ce = golden_minimum(lambda z: continuous([*continuous_anchor, z])[0], 2.4, 4.2)
        continuous_lateral_rows = []
        c_energy_grids = []
        for height in (cz, z_atop):
            grid = []
            for fy in np.linspace(0, 1, grid_count):
                for fx in np.linspace(0, 1, grid_count):
                    x, y = continuous_anchor + fx * c_u + fy * c_v
                    energy, force = continuous([x, y, height])
                    grid.append(energy)
                    continuous_lateral_rows.append([float(fx), float(fy), float(x), float(y), height, energy, *force])
            c_energy_grids.append(grid)
        c_edge_errors = [abs(continuous([*(continuous_anchor + vector), z_atop])[0]
                            - continuous([*continuous_anchor, z_atop])[0]) for vector in (c_u, c_v)]
        if max(c_edge_errors) > 1e-8:
            raise AssertionError("The continuous surface did not repeat across its primitive-cell edges.")
        write_rows(OUTPUT / "benzene-continuous-lateral.csv", ["fx", "fy", "x_A", "y_A", "height_A",
                   "energy_eV", "Fx_eV_A", "Fy_eV_A", "Fz_eV_A"], continuous_lateral_rows)
        write_rows(OUTPUT / "benzene-continuous-height.csv", ["height_A", "energy_eV", "Fx_eV_A",
                   "Fy_eV_A", "Fz_eV_A"], continuous_rows)
        report["benzene"] = {**molecule_info, "gold_molecule_ID": 1, "benzene_molecule_ID": 2,
            "SAIP_cutoff_angstrom": 16.0, "taper": True, "newton_pair": True,
            "intramolecular_model": "rigid UFF-prepared coordinates; UFF energy constant excluded from adsorption",
            "registry_anchor_xy_angstrom": anchor.tolist(), "atop_vertical_minimum_height_angstrom": z_atop,
            "atop_vertical_minimum_energy_eV": e_atop, "translation_force_gradient_check": gradient_check,
            "grouping_checks": grouping, "lateral_scans": lateral,
            "continuous_game_surface_comparison": {"driver_strength": 0, "corrugation_kcal_per_mol": 0.08,
                "orientation_angle_deg": 0,
                "anchor_vertical_minimum_height_angstrom": cz, "anchor_vertical_minimum_energy_eV": ce,
                "at_anchor_vertical_minimum_height_sampled_corrugation_eV": max(c_energy_grids[0]) - min(c_energy_grids[0]),
                "at_SAIP_angle0_atop_reference_height_angstrom": z_atop,
                "at_SAIP_angle0_atop_reference_height_sampled_corrugation_eV": max(c_energy_grids[1]) - min(c_energy_grids[1]),
                "comparison_scope": "same angle 0 and same fixed reference height for like-height comparison; each model spans its own complete primitive cell",
                "registry_anchor_xy_angstrom": continuous_anchor.tolist(),
                "primitive_vectors_angstrom": [c_u.tolist(), c_v.tolist()],
                "cell_edge_repeat_error_eV": max(c_edge_errors)},
            "paper_reference_not_an_acceptance_target": {"benzene_Au111_height_angstrom": 3.3,
                "adsorption_depth_eV": 0.689, "corrugation_meV_approximately": 52,
                "geometry": "paper rigid three-layer Au versus this minimized six-layer slab",
                "source": SAIP_DOI}}
    finally:
        probe.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="1 ps Au thermal check and smaller translation grid")
    parser.add_argument("--slab-only", action="store_true", help="Only gold EAM; skip the independently parameterized cross test")
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    report = {"status": "running", "quick": args.quick, "slab_only": args.slab_only,
        "units": {"length": "angstrom", "energy": "eV", "time": "ps", "force": "eV/angstrom"},
        "citations": {"Au_EAM": EAM_DOI, "Au_EAM_LAMMPS": "https://docs.lammps.org/pair_eam.html",
                      "benzene_Au_SAIP_paper": SAIP_DOI, "benzene_Au_SAIP_author_PDF": SAIP_PAPER,
                      "SAIP_LAMMPS": "https://docs.lammps.org/pair_saip_metal.html"},
        "limitations": ["Independent reference only; the running game remains on its continuous effective surface.",
            "Small periodic unreconstructed Au(111) with fixed bottom layers; no herringbone, steps, defects or electronic surface states.",
            "A 1–5 ps 300 K check demonstrates short numerical stability, not equilibrium or validated thermal statistics.",
            "CHAu SAIP is fitted to graphene/benzene-Au PBE+D3 data; no arbitrary-element/car coverage is claimed.",
            "An isolated carbon probe has an assumed +z normal and is outside validated benzene/graphene use, so is omitted.",
            "The nonplanar 52-atom test car is outside the cited SAIP validation domain and is omitted from physical adsorption benchmarks.",
            "Benzene and gold are frozen during cross scans; geometry relaxation, reaction, charge transfer and STM dynamics are absent.",
            "Translation corrugation is sampled max-minus-min at fixed height/orientation, not a converged diffusion barrier or NEB path.",
            "SAIP parameter-file energies are meV with embedded scaling S; LAMMPS converts the reported energy to metal-unit eV."]}
    path = OUTPUT / "report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    try:
        suite = unittest.defaultTestLoader.discover(str(ROOT / "tests"), pattern="test_au_surface.py")
        result = unittest.TextTestRunner(verbosity=1).run(suite)
        report["geometry_tests"] = {"tests_run": result.testsRun, "passed": result.wasSuccessful()}
        if not result.wasSuccessful():
            raise AssertionError("Au(111) geometry unit tests failed.")
        slab, eam, saip = verify_slab(report, args.quick)
        if not args.slab_only:
            verify_cross(report, slab, eam, saip, args.quick)
        report["status"] = "passed"
    except Exception as error:
        report.update(status="failed", error=str(error), traceback=traceback.format_exc())
        traceback.print_exc()
    finally:
        report["elapsed_wall_seconds"] = time.perf_counter() - started
        path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"Report: {path}", flush=True)
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
