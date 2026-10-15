"""Reusable compiled external forces with the original NumPy reference.

Returned force memory belongs to the evaluator and is overwritten next call.
Only implementation changes; potentials, units and dynamics stay unchanged.
"""
from __future__ import annotations

import importlib.util
import time
import numpy as np
from . import fields


class ExternalFieldEvaluator:
    def __init__(self, charges, polarizabilities, *, surface_enabled=True,
                 accelerator="auto"):
        if accelerator not in ("auto", "numpy", "numba"):
            raise ValueError("External accelerator must be auto, numpy or numba.")
        self.charges = np.array(charges, dtype=np.float64, copy=True)
        self.polarizabilities = np.array(polarizabilities, dtype=np.float64, copy=True)
        if (self.charges.ndim != 1 or self.polarizabilities.shape != self.charges.shape
                or not np.isfinite(self.charges).all()
                or not np.isfinite(self.polarizabilities).all()
                or (self.polarizabilities < 0).any()):
            raise ValueError("External response coefficients must be finite matching vectors, with nonnegative polarizabilities.")
        self.surface_enabled = bool(surface_enabled)
        self.shape = (len(self.charges), 3)
        self.forces = np.empty(self.shape, dtype=np.float64)
        self.accelerator_requested = accelerator
        self.accelerator = "numpy"
        self._kernel = None
        self.warmup_seconds = 0.
        self.numba_version = None
        if accelerator != "numpy" and importlib.util.find_spec("numba"):
            start = time.perf_counter()
            from .fields_kernel import evaluate_external, NUMBA_VERSION
            self._kernel = evaluate_external
            self.numba_version = NUMBA_VERSION
            warmup = np.zeros(self.shape, dtype=np.float64)
            warmup[:, 2] = fields.ADSORPTION_HEIGHT
            self.energy_and_forces(warmup, 1, 0., 0., 0., 7., .08)
            self.warmup_seconds = time.perf_counter()-start
            self.accelerator = "numba"
        elif accelerator == "numba":
            raise ImportError("Compiled external forces require Numba.")

    def energy_and_forces(self, positions, mode, strength, tip_x, tip_y, tip_z,
                          corrugation, growth=1.):
        if positions.shape != self.shape:
            raise ValueError("External coordinates have an unexpected shape.")
        if self._kernel is not None:
            energy = self._kernel(positions, self.charges, self.polarizabilities,
                tip_x, tip_y, tip_z, mode, strength, corrugation,
                self.surface_enabled, self.forces, growth)
        else:
            tip = {"x":tip_x, "y":tip_y, "z":tip_z, "strength":strength}
            drive = {"mode":"gradient" if mode == 1 else "tip", "strength":strength, "growth":growth}
            energy, force = fields.energy_and_forces(positions, self.charges,
                self.polarizabilities, tip, corrugation if self.surface_enabled else 0., drive)
            if not self.surface_enabled:
                exp = np.exp(-fields.ADSORPTION_ALPHA*(positions[:, 2]-fields.ADSORPTION_HEIGHT))
                energy -= float(np.sum(fields.ADSORPTION_DEPTH*(exp**2-2*exp)))
                force[:, 2] -= 2*fields.ADSORPTION_DEPTH*fields.ADSORPTION_ALPHA*(exp**2-exp)
            np.copyto(self.forces, force)
        return float(energy), self.forces

    def metadata(self):
        return {"algorithm":"fused-external-fields-v1", "accelerator":self.accelerator,
                "accelerator_requested":self.accelerator_requested,
                "numba_version":self.numba_version, "warmup_seconds":self.warmup_seconds,
                "surface_enabled":self.surface_enabled, "threads":1, "fastmath":False,
                "energy_units":"kcal/mol", "force_units":"kcal/mol/angstrom",
                "time_step_or_dynamics_changed":False}
