"""Independent explicit Au reference with OPT and safe periodic supercells.

This does not change the original racing benchmark.  The surface is the
same six-layer unreconstructed Au(111), with the layer heights taken from
the existing zero-temperature EAM relaxation.  One complete C/H molecule
has one molecule ID, and its unwrapped coordinates are never folded into
a single molecular fragment.  Native LAMMPS uses wrapped coordinates plus
image flags solely to find periodic Au and molecular normal neighbours.
"""

from __future__ import annotations

import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

from dataclasses import dataclass, asdict
from pathlib import Path
import ctypes
import time

import numpy as np

from molecular_car.au_surface import AuSlab, build_au111, potential_metadata
from molecular_car.surface_racing import (
    SurfaceRace, benzene_input, paired_velocities, decode_images, DT_PS, DAMPING_PS,
)
from scripts.lammps_runtime import open_lammps, PROJECT_ROOT

SAIP_CUTOFF = 16.0
NEIGHBOR_SKIN = 2.0
NORMAL_CUTOFF = 2.0
REFERENCE_XYZ = PROJECT_ROOT / "outputs" / "au-surface-check" / "au111-relaxed.xyz"
DEFAULT_CELL_XY = (40.0, 39.9)


@dataclass(frozen=True)
class CellPlan:
    nx: int
    ny: int
    lengths_xy: tuple[float, float]
    requested_xy: tuple[float, float]
    required_xy: tuple[float, float]
    molecular_span_xy: tuple[float, float]
    expanded: bool
    cutoff: float = SAIP_CUTOFF
    skin: float = NEIGHBOR_SKIN
    normal_cutoff: float = NORMAL_CUTOFF
    lattice_constant: float = 4.08

    def metadata(self):
        return {**asdict(self),
                "cutoff_margin_policy":"2*(cutoff+skin) is a conservative benchmark policy, not a LAMMPS hard limit",
                "normal_margin_policy":"complete unwrapped span + 2*(normal_cutoff+skin) keeps molecular copies out of normal neighbourhoods"}


def _coordinates(coordinates):
    result = np.asarray(coordinates, dtype=float)
    if result.ndim != 2 or result.shape[1:] != (3,) or len(result) == 0 or not np.isfinite(result).all():
        raise ValueError("A complete molecule needs finite N by 3 unwrapped coordinates.")
    return result


def required_cell_lengths(coordinates, *, cutoff=SAIP_CUTOFF,
                          skin=NEIGHBOR_SKIN, normal_cutoff=NORMAL_CUTOFF):
    """Conservative single-image and molecular-normal safety margins.

    L > 2*(cross cutoff+skin) avoids duplicate periodic images in the
    16 A cross neighbour range.  L > span+2*(normal cutoff+skin) separates
    the complete molecule from its own copies even when it crosses seams.
    This uses full unwrapped span, not an atom-by-atom wrapped bounding box.
    """
    coordinates = _coordinates(coordinates)
    parameters = np.array([cutoff, skin, normal_cutoff], dtype=float)
    if not np.isfinite(parameters).all() or cutoff <= 0 or skin < 0 or normal_cutoff <= 0:
        raise ValueError("Cutoff and normal cutoff must be positive; skin must be nonnegative.")
    span = np.ptp(coordinates[:, :2], axis=0)
    required = np.maximum(2*(cutoff+skin), span+2*(normal_cutoff+skin)) + 1e-6
    return span, required


def choose_safe_cell(coordinates, requested_xy=DEFAULT_CELL_XY, *, auto_expand=True,
                     cutoff=SAIP_CUTOFF, skin=NEIGHBOR_SKIN,
                     normal_cutoff=NORMAL_CUTOFF, lattice_constant=4.08):
    requested = np.asarray(requested_xy, dtype=float)
    if requested.shape != (2,) or not np.isfinite(requested).all() or np.any(requested <= 0):
        raise ValueError("Requested x/y lengths must be finite and positive.")
    if not np.isfinite(lattice_constant) or lattice_constant <= 0:
        raise ValueError("The lattice constant must be finite and positive.")
    span, required = required_cell_lengths(coordinates, cutoff=cutoff, skin=skin,
                                          normal_cutoff=normal_cutoff)
    spacings = np.array([lattice_constant/np.sqrt(2), lattice_constant*np.sqrt(1.5)])
    requested_counts = np.ceil(requested/spacings).astype(int)
    safe_counts = np.ceil(required/spacings).astype(int)
    expanded = bool(np.any(safe_counts > requested_counts))
    if expanded and not auto_expand:
        raise ValueError(f"Unsafe periodic cell: complete molecular span {span.tolist()} A needs at least {required.tolist()} A; enable auto_expand.")
    counts = np.maximum(requested_counts, safe_counts)
    lengths = counts*spacings
    return CellPlan(int(counts[0]), int(counts[1]), tuple(lengths), tuple(requested),
                    tuple(required), tuple(span), expanded, cutoff, skin,
                    normal_cutoff, lattice_constant)


def build_relaxed_supercell(plan: CellPlan, reference_path=REFERENCE_XYZ):
    """Reconstruct exact ABC phase using shared relaxed mean layer heights.

    Both comparison sizes use this reconstruction, not an arbitrary crop
    of a noncommensurate box.  Tiny numerical in-plane relaxation residuals
    of the original ideal flat slab are not tiled as a new superstructure.
    """
    reference = build_au111(lattice_constant=plan.lattice_constant)
    positions = np.loadtxt(Path(reference_path), skiprows=2, usecols=(1, 2, 3))
    if positions.shape != reference.positions.shape or not np.isfinite(positions).all():
        raise ValueError("The six-layer reference coordinates do not match the reference Au slab.")
    heights = np.array([positions[reference.layers == layer, 2].mean() for layer in range(6)])
    heights -= heights[-1]
    slab = build_au111(plan.nx, plan.ny, layers=6, frozen_layers=2,
                      lattice_constant=plan.lattice_constant)
    coordinates = slab.positions.copy()
    coordinates[:, 2] = heights[slab.layers]
    return slab.with_positions(coordinates)


def reference_geometry_metadata(reference_path=REFERENCE_XYZ):
    reference = build_au111()
    positions = np.loadtxt(Path(reference_path), skiprows=2, usecols=(1, 2, 3))
    residual = positions[:, :2]-reference.positions[:, :2]
    residual -= np.rint(residual/reference.lengths_xy)*reference.lengths_xy
    return {"source":str(Path(reference_path).resolve()),
            "maximum_xy_relaxation_residual_A":float(np.linalg.norm(residual, axis=1).max()),
            "maximum_layer_z_std_A":float(max(np.std(positions[reference.layers==layer,2]) for layer in range(6))),
            "reconstruction":"exact original ABC lattice phase with shared relaxed mean z for each of six layers"}


def _styles(style):
    if style == "eam_opt":
        return "saip/metal", "eam/opt"
    if style not in ("base", "opt"):
        raise ValueError("Explicit style must be 'base', 'opt' or 'eam_opt'.")
    suffix = "/opt" if style == "opt" else ""
    return "saip/metal"+suffix, "eam"+suffix


def _native_system(lmp, slab, coordinates, elements, potentials, style, *, mobile=False):
    saip_style, eam_style = _styles(style)
    if not lmp.has_style("pair", saip_style) or (mobile and not lmp.has_style("pair", eam_style)):
        raise RuntimeError(f"Installed LAMMPS does not support the requested {style} native pair styles.")
    coordinates = _coordinates(coordinates)
    n_gold = len(slab.positions)
    n = n_gold+len(coordinates)
    box = slab.box.copy()
    box[2] = [min(-40.0, float(coordinates[:, 2].min())-20),
              max(120.0, float(coordinates[:, 2].max())+20)]
    lengths = box[:, 1]-box[:, 0]
    bounds = " ".join(f"{value:.15g}" for row in box for value in row)
    lmp.commands_list(["units metal", "atom_style full", "boundary p p f", "newton on",
        "atom_modify map array sort 0 0.0", f"region world block {bounds} units box",
        "create_box 3 world", "mass 1 196.97", "mass 2 12.011", "mass 3 1.008"])
    original = np.concatenate([slab.positions, coordinates])
    wrapped = original.copy()
    images = np.zeros((n, 3), dtype=int)
    images[:, :2] = np.floor((original[:, :2]-box[:2, 0])/lengths[:2]).astype(int)
    wrapped[:, :2] = (original[:, :2]-box[:2, 0]) % lengths[:2]+box[:2, 0]
    packed = [lmp.encode_image_flags(*(int(v) for v in row)) for row in images]
    types = [1]*n_gold+[2 if element == "C" else 3 for element in elements]
    if lmp.create_atoms(n, list(range(1, n+1)), types, wrapped.ravel().tolist(), image=packed) != n:
        raise RuntimeError("LAMMPS failed to create the complete gold and molecular system.")
    pairs = f"hybrid/overlay zero {SAIP_CUTOFF} {saip_style} {SAIP_CUTOFF} 1"
    if mobile:
        pairs += " "+eam_style
    lmp.commands_list(["group car type 2 3", "group gold type 1", "set group car mol 2",
        "set group gold mol 1", f"pair_style {pairs}", "pair_coeff * * zero",
        f'pair_coeff * * {saip_style} "{potentials["SAIP"].as_posix()}" Au C H',
        "pair_coeff 1 1 none", f"neighbor {NEIGHBOR_SKIN} bin",
        "neigh_modify every 1 delay 0 check yes"])
    if mobile:
        lmp.command(f'pair_coeff 1 1 {eam_style} "{potentials["EAM"].as_posix()}"')
    lmp.commands_list([f"compute cross all pair {saip_style}",
                      "thermo_style custom step temp pe ke etotal c_cross", "thermo 1000000"])
    return lengths, saip_style, eam_style


class FrozenAuOracle:
    """Persistent native SAIP oracle for one complete frozen-surface molecule.

    evaluate() returns only C/H--Au cross energy/forces in eV and eV/A.
    No UFF, driver, thermostat or Au--Au EAM is included.  Coordinates may
    cross corners and many periodic cells.  An unsafe enlarged conformation
    expands the Au supercell, or raises when auto_expand is False.
    """
    def __init__(self, coordinates, elements, *, cell=None, requested_xy=DEFAULT_CELL_XY,
                 style="base", auto_expand=True, reference_path=REFERENCE_XYZ,
                 allow_unapproved_opt=False):
        self.elements = list(elements)
        coordinates = _coordinates(coordinates)
        if len(self.elements) != len(coordinates) or set(self.elements)-{"C", "H"}:
            raise ValueError("The validated native reference accepts one complete C/H molecule.")
        if style not in ("base", "opt"):
            raise ValueError("The frozen oracle accepts 'base' or experimental 'opt'.")
        if style == "opt" and not allow_unapproved_opt:
            raise ValueError("SAIP /opt failed per-atom force finite differences in this installation; use base. Diagnostic opt requires allow_unapproved_opt=True.")
        self.style, self.auto_expand, self.reference_path = style, auto_expand, Path(reference_path)
        self.plan = cell or choose_safe_cell(coordinates, requested_xy, auto_expand=auto_expand)
        if (self.plan.cutoff, self.plan.skin, self.plan.normal_cutoff) != (SAIP_CUTOFF, NEIGHBOR_SKIN, NORMAL_CUTOFF):
            raise ValueError("The native oracle preserves the fixed 16 A cutoff, 2 A skin and 2 A maximum normal cutoff.")
        self.context = None
        self._ensure_safe(coordinates)
        self._open(coordinates)

    def _ensure_safe(self, coordinates):
        _, required = required_cell_lengths(coordinates)
        if np.any(required > self.plan.lengths_xy):
            if not self.auto_expand:
                raise ValueError("Complete unwrapped molecular span is unsafe for this periodic cell.")
            self.plan = choose_safe_cell(coordinates, self.plan.lengths_xy)
            return True
        return False

    def _open(self, coordinates):
        self.slab = build_relaxed_supercell(self.plan, self.reference_path)
        self.n_gold, self.n_car = len(self.slab.positions), len(coordinates)
        self.n = self.n_gold+self.n_car
        self.context = open_lammps()
        self.lmp, self.runtime = self.context.__enter__()
        installation = Path(self.runtime["lammps_library"]).parents[1]
        self.potentials = {"EAM": installation/"Potentials"/"Au_u3.eam",
                           "SAIP": installation/"Potentials"/"CHAu.ILP"}
        try:
            self.box_lengths, self.saip_style, _ = _native_system(
                self.lmp, self.slab, coordinates, self.elements, self.potentials, self.style)
            self.image_settings = [self.lmp.extract_setting(key) for key in ("IMGMASK", "IMGMAX", "IMGBITS", "IMG2BITS")]
            self.lmp.command("run 0")
        except BaseException:
            self.close()
            raise

    def close(self):
        if self.context is not None:
            self.context.__exit__(None, None, None)
            self.context = None

    def __enter__(self):
        return self

    def __exit__(self, *arguments):
        self.close()

    def evaluate(self, coordinates, *, include_gold_forces=False):
        coordinates = _coordinates(coordinates)
        if len(coordinates) != self.n_car:
            raise ValueError("The oracle preserves molecular atom count and atom ordering.")
        if self._ensure_safe(coordinates):
            self.close()
            self._open(coordinates)
        complete = np.concatenate([self.slab.positions, coordinates])
        wrapped = complete.copy()
        images = np.zeros((self.n, 3), dtype=int)
        images[:, :2] = np.floor(complete[:, :2]/self.box_lengths[:2]).astype(int)
        wrapped[:, :2] %= self.box_lengths[:2]
        raw = (ctypes.c_double*wrapped.size)(*wrapped.ravel())
        self.lmp.scatter_atoms("x", 1, 3, raw)
        packed = np.array([self.lmp.encode_image_flags(*(int(v) for v in row)) for row in images])
        tags = self.lmp.numpy.extract_atom("id", nelem=self.n)
        self.lmp.numpy.extract_atom("image", nelem=self.n)[:] = packed[tags-1]
        self.lmp.command("run 0")
        energy = float(self.lmp.extract_compute("cross", 0, 0))
        forces = np.ctypeslib.as_array(self.lmp.gather_atoms("f", 1, 3), shape=(self.n*3,)).copy().reshape(self.n, 3)
        if not np.isfinite(energy) or not np.isfinite(forces).all():
            raise RuntimeError("Nonfinite native SAIP cross energy or forces.")
        tags = self.lmp.numpy.extract_atom("id", nelem=self.n)
        local = np.where(tags>self.n_gold)[0]
        native_unwrapped = np.empty_like(coordinates)
        native_unwrapped[tags[local]-self.n_gold-1] = (
            self.lmp.numpy.extract_atom("x", nelem=self.n)[local]+
            decode_images(self.lmp.numpy.extract_atom("image", nelem=self.n)[local],
                          *self.image_settings)*self.box_lengths)
        if not np.allclose(native_unwrapped, coordinates, rtol=0, atol=1e-8):
            raise RuntimeError("Native periodic image flags did not preserve the complete unwrapped molecule.")
        result = {"energy_eV": energy, "forces_eV_A": forces[self.n_gold:],
                  "unwrapped_positions": native_unwrapped, "cell": self.plan.metadata()}
        if include_gold_forces:
            result["all_forces_eV_A"] = forces
        return result


class ExplicitAuRace(SurfaceRace):
    """Short native performance experiment reusing the old callback/diagnostics."""
    def __init__(self, *, model="frozen_au", style="base", requested_xy=DEFAULT_CELL_XY,
                 strength=.75, temperature=600, seed=98245165, thermostat=True,
                 equilibration_steps=4000, timestep_ps=DT_PS, allow_unapproved_opt=False):
        self.preparation_start_wall, self.preparation_start_cpu = time.perf_counter(), time.process_time()
        if model not in ("frozen_au", "mobile_au"):
            raise ValueError("Explicit races require frozen_au or mobile_au.")
        if style == "opt" and not allow_unapproved_opt:
            raise ValueError("SAIP /opt is not approved by the per-atom force check; diagnostic timing needs allow_unapproved_opt=True.")
        _styles(style)
        self.model, self.style = model, style
        self.strength, self.temperature, self.seed = strength, temperature, seed
        self.thermostat, self.timestep_ps = thermostat, timestep_ps
        self.drive_enabled, self.origin_x, self.callback_error = False, 0.0, None
        self.prepared, self.geometry, self.carbons = benzene_input()
        self.n_car = len(self.geometry)
        self.masses = np.array([a["mass"] for a in self.prepared.atoms])
        self.charges = np.array([a["charge"] for a in self.prepared.atoms])
        self.alphas = np.array([a["polarizability"] for a in self.prepared.atoms])
        self.bond_a = np.array([b["a"]-1 for b in self.prepared.bonds])
        self.bond_b = np.array([b["b"]-1 for b in self.prepared.bonds])
        self.monitor = {"callback_calls": 0, "bond_min_angstrom": 1e20,
            "bond_max_angstrom": 0.0, "min_atom_z_angstrom": 1e20, "max_atom_z_angstrom": -1e20}
        self.plan = choose_safe_cell(self.geometry, requested_xy)
        self.slab = build_relaxed_supercell(self.plan)
        self.n_gold, self.n = len(self.slab.positions), len(self.slab.positions)+self.n_car
        self.gold_seed = 1+((seed^479001599)%900000000)
        self.neighbor_validation = {"method":"Au-Au SAIP removed; EAM retained only for mobile Au; independently verified by the optimization runner"}
        self.initial_orientation = "0 deg carbon bond along x; planar benzene"
        self.context = open_lammps()
        self.lmp, self.runtime = self.context.__enter__()
        try:
            self._initialize_explicit(equilibration_steps)
        except BaseException:
            self.close()
            raise

    def _initialize_explicit(self, equilibration_steps):
        lmp = self.lmp
        installation = Path(self.runtime["lammps_library"]).parents[1]
        self.eam, self.saip = installation/"Potentials"/"Au_u3.eam", installation/"Potentials"/"CHAu.ILP"
        # A fixed primitive registry independent of the chosen supercell size.
        d = self.slab.nearest_spacing
        anchor = np.array([5*d, np.sqrt(3)*d/3])
        car = self.geometry+[*anchor, 3.34256139116]
        self.box_lengths, self.saip_style, self.eam_style = _native_system(
            lmp, self.slab, car, [a["element"] for a in self.prepared.atoms],
            {"EAM":self.eam, "SAIP":self.saip}, self.style, mobile=self.model=="mobile_au")
        self.image_settings = [lmp.extract_setting(key) for key in ("IMGMASK", "IMGMAX", "IMGBITS", "IMG2BITS")]
        self.runtime.update(units="metal", timestep_ps=self.timestep_ps,
            omp_threads=lmp.extract_setting("nthreads"),
            potentials={"EAM":potential_metadata(self.eam), "SAIP":potential_metadata(self.saip)})
        lmp.commands_list(["fix molecular car external pf/callback 1 1", "fix_modify molecular energy yes",
            "compute car_temp car temp", "thermo_style custom step temp pe ke etotal c_car_temp",
            f"timestep {self.timestep_ps}", "fix gold_hold gold setforce 0 0 0"])
        lmp.set_fix_external_callback("molecular", self._callback)
        lmp.commands_list(["run 0", "fix xy_hold car setforce 0 0 NULL", "min_style fire",
            "min_modify dmax 0.02", "minimize 1e-10 1e-6 1000 10000",
            "unfix xy_hold", "velocity car set 0 0 0", "run 0"])
        self._check()
        self.vertical_relaxation = self.sample()
        if self.model=="mobile_au":
            lmp.commands_list(["unfix gold_hold", f"group frozen id 1:{self.slab.frozen_count}",
                "group au_mobile subtract gold frozen", "fix bottom_hold frozen setforce 0 0 0",
                "velocity frozen set 0 0 0", f"velocity au_mobile create {self.temperature} {self.gold_seed} mom yes rot no dist gaussian",
                "fix gold_integrate au_mobile nve", "compute au_temp au_mobile temp",
                "thermo_style custom step temp pe ke etotal c_car_temp c_au_temp"])
            if self.thermostat:
                lmp.command(f"fix gold_bath au_mobile langevin {self.temperature} {self.temperature} {DAMPING_PS} {self.gold_seed} zero yes")
            if self.thermostat and equilibration_steps:
                lmp.commands_list(["fix car_hold car setforce 0 0 0", f"run {equilibration_steps}", "unfix car_hold"])
                self._check()
        self.preparation_equilibration_steps = equilibration_steps if self.model=="mobile_au" and self.thermostat else 0
        velocity = paired_velocities(self.masses, self.temperature, self.seed)
        tags = lmp.numpy.extract_atom("id", nelem=self.n)
        local = np.where(tags>self.n_gold)[0]
        lmp.numpy.extract_atom("v", nelem=self.n)[local] = velocity[tags[local]-self.n_gold-1]
        lmp.command("fix car_integrate car nve")
        if self.thermostat:
            lmp.command(f"fix car_bath car langevin {self.temperature} {self.temperature} {DAMPING_PS} {self.seed} zero no")
        lmp.command("reset_timestep 0")
        self.origin_x = float(np.average(self.coordinates(), axis=0, weights=self.masses)[0])
        self.drive_enabled = True
        lmp.command("run 0")
        self._check()
        self.initial_sample = self.sample()
        self.preparation_wall_seconds = time.perf_counter()-self.preparation_start_wall
        self.preparation_cpu_seconds = time.process_time()-self.preparation_start_cpu
        self.monitor = {"callback_calls":0, "bond_min_angstrom":1e20, "bond_max_angstrom":0.0,
                        "min_atom_z_angstrom":1e20, "max_atom_z_angstrom":-1e20}

    def description(self):
        return {**super().description(), "style":self.style, "cell":self.plan.metadata(),
                "approved_native_force_mode":self.style != "opt",
                "reference_reconstruction":reference_geometry_metadata()}

