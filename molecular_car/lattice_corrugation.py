"""Conservative crystal contrast control, independent of mean adhesion.

The six translated copies average out the leading Au(111) harmonics. This is
a finite phase filter, not an exactly flat continuum: higher harmonics remain.
Its continuous unit-cell average equals the original potential at every height
and molecular configuration. Complete molecular translations retain every
local-normal chain force and the separate hetero-Au interaction.
"""
from __future__ import annotations
import numpy as np


def validate_lattice_scale(value):
    if isinstance(value, bool):
        raise ValueError("晶格起伏倍率必须是 0–12 的有限数值。")
    value = float(value)
    if not np.isfinite(value) or not 0 <= value <= 12:
        raise ValueError("晶格起伏倍率必须是 0–12 的有限数值。")
    return value


class LatticeCorrugation:
    def __init__(self, base, scale=1.):
        self.base = base
        self.scale = validate_lattice_scale(scale)
        self.shifts = np.array([[i*base.period[0]/2, j*base.period[1]/3, 0.]
                               for i in range(2) for j in range(3)])
        if self.scale != 1. and hasattr(base, 'elements'):
            # Compile the weighted ABI before the joint-relaxation wall budget.
            warmup = np.column_stack([np.arange(len(base.elements))*3., np.zeros(len(base.elements)), np.full(len(base.elements),3.4)])
            base.energy_and_forces_eV(warmup, np.ones(len(base.elements)), np.zeros(len(base.elements)))

    def __getattr__(self, name):
        return getattr(self.base, name)

    def phase_mean_eV(self, positions, reference=None):
        energy, force = 0., np.zeros_like(positions, dtype=float)
        for k, shift in enumerate(self.shifts):
            e, f = reference if k == 0 and reference is not None else self.base.energy_and_forces_eV(positions+shift)
            energy += e/len(self.shifts)
            force += f/len(self.shifts)
        return float(energy), force

    def energy_and_forces_eV(self, positions):
        # Exact compatibility path: no filtering, summation, or extra calls.
        if self.scale == 1.:
            return self.base.energy_and_forces_eV(positions)
        if (getattr(self.base, "accelerator", None) == "numba"
                and hasattr(self.base, "energy_and_forces_corrugated_eV")):
            return self.base.energy_and_forces_corrugated_eV(positions, self.scale)
        # Every atomic interaction center has its own core gate. Weight its
        # normal gradient BEFORE transmitting forces to covalent neighbours.
        t = np.clip((np.asarray(positions)[:, 2]-2.)/.6, 0., 1.)
        weights = t**3*(10+t*(-15+6*t))
        derivative = 30*t*t*(t-1)**2/.6
        factor = 1.+(self.scale-1.)*weights
        dz = (self.scale-1.)*derivative
        energy, force = self.base.energy_and_forces_eV(positions, (1+5*factor)/6, 5*dz/6)
        for shift in self.shifts[1:]:
            e, f = self.base.energy_and_forces_eV(positions+shift, (1-factor)/6, -dz/6)
            energy += e; force += f
        return float(energy), force

    def energy_and_forces(self, positions):
        from .periodic_au import EV_TO_KCAL_MOL
        e, f = self.energy_and_forces_eV(positions)
        return e*EV_TO_KCAL_MOL, f*EV_TO_KCAL_MOL

    def metadata(self):
        return dict(self.base.metadata(), lattice_scale=self.scale,
                    corrugation_control={"version": "six-phase-atomic-crystal-contrast-v3",
                        "formula": "sum_i [U_i_reference + (lambda-1)*w(z_i)*(U_i_reference-U_i_phase)]",
                        "core_guard": "independent atomic quintic ramps, 0 at z<=2 A, 1 at z>=2.6 A; includes dw/dz and weighted normal-chain forces",
                        "phase_grid": [2, 3], "shifts_A": self.shifts.tolist(),
                        "mean_adhesion": "continuous lateral cell average unchanged at fixed geometry and height",
                        "zero_scale": "finite phase filter; higher lattice harmonics remain",
                        "literature_reference_scale": 1.,
                        "game_adjustment": self.scale != 1.,
                        "forces": "same scalar energy; complete normal-chain derivatives"})
