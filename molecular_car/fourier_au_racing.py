"""Isolated Fourier surface benchmark using the established MD protocol."""

from __future__ import annotations

import time
import numpy as np

from molecular_car.surface_racing import (SurfaceRace, benzene_input, DT_PS,
    continuous_surface_energy_forces, EV_TO_KCAL_MOL)
from molecular_car.periodic_au_racing import PeriodicAuRace
from molecular_car.fourier_au import create_fourier_surface


class FourierAuRace(PeriodicAuRace):
    def __init__(self, *, strength=.75, temperature=600, seed=98245165,
                 thermostat=True, timestep_ps=DT_PS, surface_options=None,
                 surface_factory=None):
        constructor_wall, constructor_cpu = time.perf_counter(), time.process_time()
        prepared, geometry, _ = benzene_input()
        spacing = 4.08 / np.sqrt(2)
        self.anchor = np.array([5 * spacing, np.sqrt(3) * spacing / 3])
        started = time.perf_counter()
        self.surface = (surface_factory or create_fourier_surface)(
            prepared.atoms, geometry + [*self.anchor, 3.34256139116],
            **(surface_options or {}))
        self.environment_preparation_wall_seconds = time.perf_counter() - started
        SurfaceRace.__init__(self, "continuous", strength, temperature, seed,
                             thermostat=thermostat, equilibration_steps=0,
                             timestep_ps=timestep_ps)
        self.preparation_wall_seconds = time.perf_counter() - constructor_wall
        self.preparation_cpu_seconds = time.process_time() - constructor_cpu

    def description(self):
        return {**super().description(), "surface": "fourier_au"}


class _ContinuousSurface:
    def energy_and_forces_eV(self, positions):
        energy, force = continuous_surface_energy_forces(positions)
        return energy / EV_TO_KCAL_MOL, force / EV_TO_KCAL_MOL

    def metadata(self):
        return {"model": "unchanged effective Morse and triangular corrugation",
                "mode": "continuum", "track_boundaries": False}


class ContinuumComparisonRace(FourierAuRace):
    """Same MD box/callback as both proxies; original continuous surface."""
    def __init__(self, **kwargs):
        super().__init__(surface_factory=lambda *args, **kw: _ContinuousSurface(), **kwargs)

    def description(self):
        return {**super().description(), "surface": "continuous"}
