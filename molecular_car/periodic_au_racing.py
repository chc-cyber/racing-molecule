"""Isolated frozen-Au potential benchmark, with no Au in the MD domain.

The original surface_racing.py and its historical 18 trajectories stay
unchanged. The inherited molecular callback and diagnostics use full
coordinates; only the environment potential has periodic x/y structure.
"""

from __future__ import annotations

import time
import numpy as np

from scripts.lammps_runtime import PROJECT_ROOT
from molecular_car.surface_racing import (
    SurfaceRace, benzene_input, gradient_energy_forces, paired_velocities,
    EV_TO_KCAL_MOL, DT_PS, DAMPING_PS,
)
from molecular_car.periodic_au import create_periodic_surface


class PeriodicAuRace(SurfaceRace):
    def __init__(self, *, strength=.75, temperature=600, seed=98245165,
                 thermostat=True, timestep_ps=DT_PS):
        constructor_wall, constructor_cpu = time.perf_counter(), time.process_time()
        prepared, geometry, _ = benzene_input()
        spacing = 4.08 / np.sqrt(2)
        self.anchor = np.array([5 * spacing, np.sqrt(3) * spacing / 3])
        started = time.perf_counter()
        self.surface = create_periodic_surface(
            prepared.atoms, geometry + [*self.anchor, 3.34256139116])
        self.environment_preparation_wall_seconds = time.perf_counter() - started
        super().__init__("continuous", strength, temperature, seed,
                         thermostat=thermostat, equilibration_steps=0,
                         timestep_ps=timestep_ps)
        self.preparation_wall_seconds = time.perf_counter() - constructor_wall
        self.preparation_cpu_seconds = time.process_time() - constructor_cpu

    def _initialize(self, equilibration_steps):
        lmp = self.lmp
        box = np.array([[-110., 1110.], [-350., 350.], [-40., 120.]])
        self.box_lengths = box[:, 1] - box[:, 0]
        bounds = " ".join(f"{value:.15g}" for row in box for value in row)
        lmp.commands_list([
            "units metal", "atom_style full", "boundary f f f", "newton on",
            "atom_modify map array sort 0 0.0", f"region world block {bounds} units box",
            "create_box 3 world", "mass 1 196.97", "mass 2 12.011", "mass 3 1.008",
        ])
        car = self.geometry + [*self.anchor, 3.34256139116]
        types = [2 if atom["element"] == "C" else 3 for atom in self.prepared.atoms]
        if lmp.create_atoms(self.n, list(range(1, self.n + 1)), types,
                            car.ravel().tolist()) != self.n:
            raise RuntimeError("Failed to create the complete benchmark molecule.")
        self.image_settings = [lmp.extract_setting(key) for key in
                               ("IMGMASK", "IMGMAX", "IMGBITS", "IMG2BITS")]
        lmp.commands_list([
            "group car type 2 3", "pair_style zero 10.0", "pair_coeff * *",
            "neighbor 2.0 bin", "neigh_modify every 1 delay 0 check yes",
            "fix molecular car external pf/callback 1 1", "fix_modify molecular energy yes",
            "compute car_temp car temp", "thermo_style custom step temp pe ke etotal c_car_temp",
            "thermo 1000000", f"timestep {self.timestep_ps}",
        ])
        self.runtime.update(units="metal", timestep_ps=self.timestep_ps,
                            omp_threads=lmp.extract_setting("nthreads"),
                            environment=self.surface.metadata())
        lmp.set_fix_external_callback("molecular", self._callback)
        lmp.commands_list([
            "run 0", "fix xy_hold car setforce 0 0 NULL", "min_style fire",
            "min_modify dmax 0.02", "minimize 1e-10 1e-6 1000 10000",
            "unfix xy_hold", "velocity car set 0 0 0", "run 0",
        ])
        self._check()
        self.vertical_relaxation = self.sample()
        self.preparation_equilibration_steps = 0
        self.neighbor_validation = {"method": "Only the complete molecule is in LAMMPS; no periodic car copies."}
        velocity = paired_velocities(self.masses, self.temperature, self.seed)
        tags = lmp.numpy.extract_atom("id", nelem=self.n)
        lmp.numpy.extract_atom("v", nelem=self.n)[:] = velocity[tags - 1]
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
        self.monitor = {"callback_calls": 0, "bond_min_angstrom": 1e20,
                        "bond_max_angstrom": 0., "min_atom_z_angstrom": 1e20,
                        "max_atom_z_angstrom": -1e20}

    def external(self, positions):
        flat = positions.ravel().tolist()
        energy = float(self.prepared.forcefield.CalcEnergy(flat))
        force = -np.array(self.prepared.forcefield.CalcGrad(flat)).reshape(self.n_car, 3)
        driver, extra = gradient_energy_forces(
            positions, self.charges, self.alphas,
            self.strength if self.drive_enabled else 0., self.origin_x)
        surface_energy, surface_force = self.surface.energy_and_forces_eV(positions)
        return ((energy + driver) / EV_TO_KCAL_MOL + surface_energy,
                (force + extra) / EV_TO_KCAL_MOL + surface_force)

    def description(self):
        return {**super().description(), "surface": "periodic_au",
                "environment": self.surface.metadata(),
                "environment_preparation_wall_seconds": self.environment_preparation_wall_seconds,
                "boundary": "f f f; surface repeats, molecule and driving field do not"}
