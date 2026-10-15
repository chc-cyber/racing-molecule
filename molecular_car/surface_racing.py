"""Independent dynamic surface benchmark; never imported by the game engine.

LAMMPS metal units: angstrom, ps, eV; RDKit UFF and the game's analytic
fields are evaluated in kcal/mol and divided by EV_TO_KCAL_MOL for both
energies and forces. Native SAIP/EAM already return eV. Periodic x/y image
flags unwrap the molecular coordinates, including UFF and the external
field. The driving x origin is the initial molecular COM, never a box wrap.
"""

from __future__ import annotations

import os
for _variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_variable] = "1"

import ctypes
from pathlib import Path
import time

import numpy as np

from scripts.lammps_runtime import open_lammps, PROJECT_ROOT
from molecular_car.au_surface import build_au111, potential_metadata
from molecular_car.molecule import prepare_molecule
from molecular_car import fields

EV_TO_KCAL_MOL = 23.0605478306
KB_EV_K = 8.617333262145e-5
MVV_TO_EV = 1.0364269652680505e-4
DT_PS = 0.0005
DAMPING_PS = 0.35
MODELS = ("continuous", "frozen_au", "mobile_au")


def gradient_energy_forces(coordinates, charges, alphas, strength, origin_x):
    """Exactly the existing broad gradient, isolated from surface/boundaries."""
    x, _, z = np.asarray(coordinates).T
    relative_x = x - origin_x
    base = fields.GRADIENT_BASE_FIELD * strength
    slope = fields.GRADIENT_FIELD_SLOPE * strength
    ex, ez = slope * z, base + slope * relative_x
    response = fields.GRADIENT_POLARIZATION_GAIN * alphas
    energy = np.sum(-charges * ez * z - 0.5 * response * (ex**2 + ez**2))
    forces = np.zeros_like(coordinates)
    forces[:, 0] = charges * ex + response * slope * ez
    forces[:, 2] = charges * ez + response * slope**2 * z
    return float(energy), forces


def continuous_surface_energy_forces(coordinates, corrugation=0.08):
    """Original Morse/triangular surface; remove only the race's soft walls."""
    x, y, z = np.asarray(coordinates).T
    height = z - fields.ADSORPTION_HEIGHT
    exponential = np.exp(-fields.ADSORPTION_ALPHA * height)
    energy = fields.ADSORPTION_DEPTH * (exponential**2 - 2 * exponential)
    forces = np.zeros_like(coordinates)
    forces[:, 2] = 2 * fields.ADSORPTION_DEPTH * fields.ADSORPTION_ALPHA * (exponential**2 - exponential)
    kx = 2 * np.pi / fields.SURFACE_SPACING
    ky = kx / np.sqrt(3)
    cx, sx, cy, sy = np.cos(kx*x), np.sin(kx*x), np.cos(ky*y), np.sin(ky*y)
    pattern = (2*cx*cy + np.cos(2*ky*y))/3
    envelope = corrugation * np.exp(-fields.CORRUGATION_DECAY * height)
    energy += envelope * pattern
    forces[:, 0] += envelope * 2*kx/3 * sx*cy
    forces[:, 1] += envelope * 2*ky/3 * (cx*sy + np.sin(2*ky*y))
    forces[:, 2] += fields.CORRUGATION_DECAY * envelope * pattern
    return float(np.sum(energy)), forces


def decode_images(packed, mask, maximum, bits, bits2):
    values = np.asarray(packed, dtype=np.int64)
    return np.stack([(values & mask) - maximum,
                     ((values >> bits) & mask) - maximum,
                     ((values >> bits2) & mask) - maximum], axis=-1)


def paired_velocities(masses, temperature, seed):
    rng = np.random.default_rng(seed)
    velocity = rng.standard_normal((len(masses), 3)) * np.sqrt(KB_EV_K * temperature / (masses[:, None] * MVV_TO_EV))
    velocity -= np.average(velocity, axis=0, weights=masses)
    if temperature > 0:
        energy = 0.5 * MVV_TO_EV * np.sum(masses[:, None] * velocity**2)
        velocity *= np.sqrt(0.5 * (3*len(masses)-3) * KB_EV_K * temperature / energy)
    return velocity


def benzene_input():
    path = PROJECT_ROOT / "examples" / "benzene.mol"
    prepared = prepare_molecule(path.read_text(encoding="utf-8"), path.name)
    carbons = np.array([i for i, atom in enumerate(prepared.atoms) if atom["element"] == "C"])
    positions = prepared.positions - np.mean(prepared.positions[carbons], axis=0)
    _, _, vectors = np.linalg.svd(positions[carbons])
    normal = vectors[-1]
    ex = positions[carbons[1]] - positions[carbons[0]]
    ex -= normal * np.dot(ex, normal)
    ex /= np.linalg.norm(ex)
    ey = np.cross(normal, ex)
    positions = positions @ np.stack([ex, ey, normal]).T
    return prepared, positions, carbons


class SurfaceRace:
    def __init__(self, model, strength, temperature, seed, *, thermostat=True,
                 equilibration_steps=4000, timestep_ps=DT_PS, optimize_neighbors=True):
        self.preparation_start_wall, self.preparation_start_cpu = time.perf_counter(), time.process_time()
        if model not in MODELS:
            raise ValueError("Unknown reference surface.")
        self.model, self.strength, self.temperature, self.seed = model, strength, temperature, seed
        self.thermostat, self.timestep_ps = thermostat, timestep_ps
        self.optimize_neighbors = optimize_neighbors
        self.drive_enabled = False
        self.origin_x = 0.0
        self.callback_error = None
        self.prepared, self.geometry, self.carbons = benzene_input()
        self.n_car = len(self.geometry)
        self.masses = np.array([a["mass"] for a in self.prepared.atoms])
        self.charges = np.array([a["charge"] for a in self.prepared.atoms])
        self.alphas = np.array([a["polarizability"] for a in self.prepared.atoms])
        self.bond_a = np.array([b["a"]-1 for b in self.prepared.bonds])
        self.bond_b = np.array([b["b"]-1 for b in self.prepared.bonds])
        self.monitor = {"callback_calls": 0, "bond_min_angstrom": 1e20, "bond_max_angstrom": 0.0,
                        "min_atom_z_angstrom": 1e20, "max_atom_z_angstrom": -1e20}
        self.slab = build_au111()
        reference_path = PROJECT_ROOT / "outputs" / "au-surface-check" / "au111-relaxed.xyz"
        if not reference_path.exists():
            raise FileNotFoundError("Run scripts/verify_au_surface.py first to create the independently minimized Au reference.")
        self.slab = self.slab.with_positions(np.loadtxt(reference_path, skiprows=2, usecols=(1,2,3)))
        self.n_gold = 0 if model == "continuous" else len(self.slab.positions)
        self.n = self.n_gold + self.n_car
        self.gold_seed = 1 + ((seed ^ 479001599) % 900000000)
        self.context = open_lammps()
        self.lmp, self.runtime = self.context.__enter__()
        try:
            self._initialize(equilibration_steps)
        except BaseException:
            self.close()
            raise

    def close(self):
        if self.context is not None:
            self.context.__exit__(None, None, None)
            self.context = None

    def _initialize(self, equilibration_steps):
        lmp = self.lmp
        installation = Path(self.runtime["lammps_library"]).parents[1]
        self.eam, self.saip = installation/"Potentials"/"Au_u3.eam", installation/"Potentials"/"CHAu.ILP"
        self.runtime.update(units="metal", timestep_ps=self.timestep_ps,
                            omp_threads=lmp.extract_setting("nthreads"),
                            potentials={"EAM":potential_metadata(self.eam), "SAIP":potential_metadata(self.saip)})
        box = self.slab.box.copy()
        box[2] = [-40, 120]
        self.box_lengths = box[:, 1] - box[:, 0]
        bounds = " ".join(f"{v:.15g}" for row in box for v in row)
        lmp.commands_list(["units metal", "atom_style full", "boundary p p f", "newton on",
            "atom_modify map array sort 0 0.0", f"region world block {bounds} units box",
            "create_box 3 world", "mass 1 196.97", "mass 2 12.011", "mass 3 1.008"])
        if self.n_gold:
            top = self.slab.positions[self.slab.layers == self.slab.layer_count-1]
            center_index = np.argmin(np.linalg.norm(top[:, :2] - self.slab.lengths_xy/2, axis=1))
            anchor = top[center_index, :2]
            height = 3.34256139116
        else:
            anchor, height = np.array([28.8, 0.0]), 2.91088660868
        self.initial_orientation = "0 deg carbon bond along x; planar benzene"
        car_positions = self.geometry + [*anchor, height]
        # The continuous box is bookkeeping only. Encode initial images as
        # necessary so even atoms initially crossing y=0 unwrap correctly.
        original = np.concatenate([self.slab.positions[:self.n_gold], car_positions])
        positions = original.copy()
        image_counts = np.zeros_like(positions, dtype=int)
        image_counts[:, :2] = np.floor(positions[:, :2] / self.box_lengths[:2]).astype(int)
        positions[:, :2] %= self.box_lengths[:2]
        packed = [lmp.encode_image_flags(*(int(value) for value in row)) for row in image_counts]
        types = [1]*self.n_gold + [2 if a["element"] == "C" else 3 for a in self.prepared.atoms]
        created = lmp.create_atoms(self.n, list(range(1,self.n+1)), types, positions.ravel().tolist(), image=packed)
        if created != self.n:
            raise RuntimeError("Not all reference atoms were created.")
        self.image_settings = [lmp.extract_setting(key) for key in ("IMGMASK", "IMGMAX", "IMGBITS", "IMG2BITS")]
        lmp.commands_list(["group car type 2 3", "group gold type 1", "set group car mol 2", "set group gold mol 1",
            "neighbor 2.0 bin", "neigh_modify every 1 delay 0 check yes"])
        if self.n_gold:
            styles = "hybrid/overlay zero 16.0 saip/metal 16.0 1"
            if self.model == "mobile_au": styles += " eam"
            lmp.commands_list([f"pair_style {styles}", "pair_coeff * * zero",
                f'pair_coeff * * saip/metal "{self.saip.as_posix()}" Au C H'])
            if self.model == "mobile_au": lmp.command(f'pair_coeff 1 1 eam "{self.eam.as_posix()}"')
        else:
            lmp.commands_list(["pair_style zero 16.0", "pair_coeff * *"])
        lmp.commands_list(["fix molecular car external pf/callback 1 1", "fix_modify molecular energy yes",
            "compute car_temp car temp", "thermo_style custom step temp pe ke etotal c_car_temp",
            "thermo 1000000", f"timestep {self.timestep_ps}"])
        lmp.set_fix_external_callback("molecular", self._callback)
        lmp.command("run 0")
        self._check()
        self.neighbor_validation = None
        if self.n_gold and self.optimize_neighbors:
            old_energy = float(lmp.get_thermo("pe"))
            old_forces = self._all_property("f")
            lmp.command("pair_coeff 1 1 none")
            if self.model == "mobile_au": lmp.command(f'pair_coeff 1 1 eam "{self.eam.as_posix()}"')
            lmp.command("run 0")
            self._check()
            new_energy, new_forces = float(lmp.get_thermo("pe")), self._all_property("f")
            e_error, f_error = abs(new_energy-old_energy), float(np.max(np.abs(new_forces-old_forces)))
            if e_error > 1e-8 or f_error > 1e-8:
                raise AssertionError("Removing zero Au-Au SAIP interactions changed energy or forces.")
            self.neighbor_validation = {"energy_absolute_difference_eV":e_error,
                "maximum_force_difference_eV_per_angstrom":f_error,
                "method":"Au-Au hybrid mapping removed; EAM restored only for mobile Au; car normal neighbours retained"}
        # All surfaces retain the same molecular x/y graph and orientation.
        # Only z coordinates relax against their own surface, with no drive.
        if self.n_gold: lmp.command("fix gold_hold gold setforce 0 0 0")
        lmp.commands_list(["fix xy_hold car setforce 0 0 NULL", "min_style fire", "min_modify dmax 0.02",
            "minimize 1e-10 1e-6 1000 10000", "unfix xy_hold", "velocity car set 0 0 0", "run 0"])
        self._check()
        self.vertical_relaxation = self.sample()
        if self.model == "mobile_au":
            lmp.command("unfix gold_hold")
            lmp.commands_list([f"group frozen id 1:{self.slab.frozen_count}", "group au_mobile subtract gold frozen",
                "fix bottom_hold frozen setforce 0 0 0", "velocity frozen set 0 0 0",
                f"velocity au_mobile create {self.temperature} {self.gold_seed} mom yes rot no dist gaussian",
                "fix gold_integrate au_mobile nve"])
            lmp.commands_list(["compute au_temp au_mobile temp", "thermo_style custom step temp pe ke etotal c_car_temp c_au_temp"])
            if self.thermostat:
                lmp.command(f"fix gold_bath au_mobile langevin {self.temperature} {self.temperature} {DAMPING_PS} {self.gold_seed} zero yes")
            if self.thermostat and equilibration_steps:
                lmp.commands_list(["fix car_hold car setforce 0 0 0", f"run {equilibration_steps}", "unfix car_hold"])
                self._check()
        self.preparation_equilibration_steps = equilibration_steps if self.model == "mobile_au" and self.thermostat else 0
        # Reuse the identical array across models and strengths, unlike native
        # velocity/create which can consume random draws for Au atoms/IDs.
        velocity = paired_velocities(self.masses, self.temperature, self.seed)
        tags = lmp.numpy.extract_atom("id", nelem=self.n)
        car_local = np.where(tags > self.n_gold)[0]
        order = tags[car_local] - self.n_gold - 1
        lmp.numpy.extract_atom("v", nelem=self.n)[car_local] = velocity[order]
        lmp.command("fix car_integrate car nve")
        if self.thermostat:
            lmp.command(f"fix car_bath car langevin {self.temperature} {self.temperature} {DAMPING_PS} {self.seed} zero no")
        lmp.command("reset_timestep 0")
        self.origin_x = float(np.average(self.coordinates(), axis=0, weights=self.masses)[0])
        self.drive_enabled = True
        lmp.command("run 0")
        self._check()
        self.initial_sample = self.sample()
        self.preparation_wall_seconds = time.perf_counter() - self.preparation_start_wall
        self.preparation_cpu_seconds = time.process_time() - self.preparation_start_cpu
        self.monitor = {"callback_calls": 0, "bond_min_angstrom": 1e20, "bond_max_angstrom": 0.0,
                        "min_atom_z_angstrom": 1e20, "max_atom_z_angstrom": -1e20}

    def _all_property(self, name):
        return np.ctypeslib.as_array(self.lmp.gather_atoms(name,1,3),shape=(self.n*3,)).copy().reshape(self.n,3)

    def coordinates(self):
        tags = self.lmp.numpy.extract_atom("id",nelem=self.n)
        local = np.where(tags > self.n_gold)[0]
        result = np.empty((self.n_car,3))
        images = decode_images(self.lmp.numpy.extract_atom("image",nelem=self.n)[local], *self.image_settings)
        result[tags[local]-self.n_gold-1] = self.lmp.numpy.extract_atom("x",nelem=self.n)[local] + images*self.box_lengths
        return result

    def external(self, positions):
        flat = positions.ravel().tolist()
        internal = float(self.prepared.forcefield.CalcEnergy(flat))
        force = -np.array(self.prepared.forcefield.CalcGrad(flat)).reshape(self.n_car,3)
        surface = 0.0
        if self.model == "continuous":
            surface, extra = continuous_surface_energy_forces(positions)
            force += extra
        driver, extra = gradient_energy_forces(positions,self.charges,self.alphas,
            self.strength if self.drive_enabled else 0.0,self.origin_x)
        force += extra
        return (internal+surface+driver)/EV_TO_KCAL_MOL, force/EV_TO_KCAL_MOL

    def _callback(self, caller, step, nlocal, tags, positions, output):
        try:
            if nlocal != self.n: raise RuntimeError("The benchmark requires single-process LAMMPS.")
            local = np.where(tags > self.n_gold)[0]
            order = tags[local]-self.n_gold-1
            images = decode_images(self.lmp.numpy.extract_atom("image",nelem=nlocal)[local],*self.image_settings)
            car = np.empty((self.n_car,3))
            car[order] = positions[local]+images*self.box_lengths
            energy, force = self.external(car)
            lengths = np.linalg.norm(car[self.bond_a]-car[self.bond_b],axis=1)
            if not np.isfinite(car).all() or not np.isfinite(energy) or not np.isfinite(force).all():
                raise RuntimeError("Nonfinite molecular coordinates/energy/forces.")
            if lengths.min()<0.6 or lengths.max()>2.5: raise RuntimeError("Unstable benzene bond lengths.")
            output[:]=0
            output[local]=force[order]
            self.lmp.fix_external_set_energy_global("molecular",energy)
            self.monitor["callback_calls"]+=1
            self.monitor["bond_min_angstrom"]=min(self.monitor["bond_min_angstrom"],float(lengths.min()))
            self.monitor["bond_max_angstrom"]=max(self.monitor["bond_max_angstrom"],float(lengths.max()))
            self.monitor["min_atom_z_angstrom"]=min(self.monitor["min_atom_z_angstrom"],float(car[:,2].min()))
            self.monitor["max_atom_z_angstrom"]=max(self.monitor["max_atom_z_angstrom"],float(car[:,2].max()))
        except BaseException as error:
            self.callback_error = str(error)
            output[:]=0

    def _check(self):
        if self.callback_error: raise RuntimeError(self.callback_error)

    def advance(self, steps):
        wall,cpu=time.perf_counter(),time.process_time()
        self.lmp.command(f"run {steps} pre no post no")
        timing={"wall_seconds":time.perf_counter()-wall,"cpu_seconds":time.process_time()-cpu}
        self._check()
        return timing

    def sample(self):
        positions=self.coordinates()
        center=np.average(positions,axis=0,weights=self.masses)
        carbons=positions[self.carbons]
        _,_,vectors=np.linalg.svd(carbons-np.mean(carbons,axis=0))
        normal=vectors[-1]
        plane_rms=float(np.sqrt(np.mean(((carbons-np.mean(carbons,axis=0))@normal)**2)))
        tags=self.lmp.numpy.extract_atom("id",nelem=self.n)
        local=np.where(tags>self.n_gold)[0]
        velocities=np.empty((self.n_car,3))
        velocities[tags[local]-self.n_gold-1]=self.lmp.numpy.extract_atom("v",nelem=self.n)[local]
        velocity=np.average(velocities,axis=0,weights=self.masses)
        energy=float(self.lmp.get_thermo("etotal"))
        if not np.isfinite(energy): raise RuntimeError("Nonfinite total energy.")
        lengths=np.linalg.norm(positions[self.bond_a]-positions[self.bond_b],axis=1)
        car_forces=self._all_property("f")[self.n_gold:]
        top_gold_z=0.0
        if self.n_gold:
            gold=self._all_property("x")[:self.n_gold]
            top_gold_z=float(np.mean(gold[self.slab.layers==self.slab.layer_count-1,2]))
        return {"step":int(self.lmp.get_thermo("step")),"time_ps":int(self.lmp.get_thermo("step"))*self.timestep_ps,
            "center_angstrom":center.tolist(),"velocity_angstrom_per_ps":velocity.tolist(),
            "total_energy_eV":energy,"potential_energy_eV":float(self.lmp.get_thermo("pe")),
            "kinetic_energy_eV":float(self.lmp.get_thermo("ke")),
            "car_temperature_K":float(self.lmp.extract_compute("car_temp",0,0)),
            "Au_temperature_K":float(self.lmp.extract_compute("au_temp",0,0)) if self.model=="mobile_au" and self.lmp.has_id("compute","au_temp") else 0.0,
            "bond_min_angstrom":float(lengths.min()),"bond_max_angstrom":float(lengths.max()),
            "max_car_force_eV_per_angstrom":float(np.linalg.norm(car_forces,axis=1).max()),
            "top_Au_mean_z_angstrom":top_gold_z,"COM_height_above_top_Au_angstrom":float(center[2]-top_gold_z),
            "carbon_plane_rms_angstrom":plane_rms,
            "carbon_plane_tilt_deg":float(np.degrees(np.arccos(np.clip(abs(normal[2]),0,1)))),
            "min_atom_z_angstrom":float(positions[:,2].min()),"max_atom_z_angstrom":float(positions[:,2].max()),
            "positions":positions.tolist()}

    def description(self):
        return {"surface":self.model,"strength":self.strength,"temperature_K":self.temperature,"seed":self.seed,
            "Au_seed":self.gold_seed,"thermostat":self.thermostat,"timestep_ps":self.timestep_ps,
            "damping_ps":DAMPING_PS,"car_atoms":self.n_car,"Au_atoms":self.n_gold,"total_atoms":self.n,
            "integrated_atoms":self.n_car+(self.n_gold-self.slab.frozen_count if self.model=="mobile_au" else 0),
            "runtime":self.runtime,"neighbor_validation":self.neighbor_validation,
            "preparation_wall_seconds":self.preparation_wall_seconds,"preparation_cpu_seconds":self.preparation_cpu_seconds,
            "Au_equilibration_steps":self.preparation_equilibration_steps,"Au_equilibration_ps":self.preparation_equilibration_steps*self.timestep_ps,
            "vertical_relaxation":self.vertical_relaxation,"drive_origin_initial_COM_x_angstrom":self.origin_x,
            "initial_sample":self.initial_sample,
            "paired_initial_velocities_angstrom_per_ps":paired_velocities(self.masses,self.temperature,self.seed).tolist()}
