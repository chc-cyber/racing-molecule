"""Compiled external potentials with the 20 nm square track's soft walls.

The caller prepares the immutable response coefficients and a force buffer once.
All values use the original ``fields`` units: Angstrom, kcal/mol, and elementary
charge.  This kernel does not integrate coordinates, alter the thermostat, or
change the timestep.  ``surface_enabled=False`` omits the Morse and triangular
surface terms so that a separate periodic Au surface can supply them directly.

There are no array allocations in the evaluation loop.  The supplied force
buffer is overwritten, and input coordinates and response arrays are read only.
Parameter and shape validation belongs to the public wrapper; dynamic finite
coordinate checks and finite output checks remain inside every kernel call.

The source revision also invalidates Numba's earlier disk cache, whose compiled
wall bounds were those of the previous 100 nm by 50 nm race geometry.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

os.environ.setdefault(
    "NUMBA_CACHE_DIR", str(Path(__file__).resolve().parents[1] / ".local/numba-cache")
)
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import numba

from .fields import (
    ADSORPTION_ALPHA,
    ADSORPTION_DEPTH,
    ADSORPTION_HEIGHT,
    BOUNDARY_STIFFNESS,
    BOUNDARY_WIDTH,
    CORRUGATION_DECAY,
    FIELD_SCALE,
    GRADIENT_BASE_FIELD,
    GRADIENT_FIELD_SLOPE,
    GRADIENT_POLARIZATION_GAIN,
    LANE_HALF_WIDTH,
    SURFACE_SPACING,
    TIP_SOFTENING,
    TRACK_X_MAX,
    TRACK_X_MIN,
)

NUMBA_VERSION = numba.__version__
TIP_MODE = 0
GRADIENT_MODE = 1
_KX = 2.0 * math.pi / SURFACE_SPACING
_KY = _KX / math.sqrt(3.0)
_TIP_SOFTENING_SQUARED = TIP_SOFTENING * TIP_SOFTENING


@numba.njit(cache=True, fastmath=False, inline="always")
def _softplus_sigmoid(value):
    """Return log(1 + exp(value)) and its derivative without overflow."""
    if value >= 0.0:
        exponential = math.exp(-value)
        return value + math.log1p(exponential), 1.0 / (1.0 + exponential)
    exponential = math.exp(value)
    return math.log1p(exponential), exponential / (1.0 + exponential)


@numba.njit(cache=True, fastmath=False, inline="always")
def _boundary(coordinate, lower, upper):
    low, low_slope = _softplus_sigmoid((lower - coordinate) / BOUNDARY_WIDTH)
    high, high_slope = _softplus_sigmoid((coordinate - upper) / BOUNDARY_WIDTH)
    energy = 0.5 * BOUNDARY_STIFFNESS * BOUNDARY_WIDTH**2 * (low**2 + high**2)
    force = BOUNDARY_STIFFNESS * BOUNDARY_WIDTH * (
        low * low_slope - high * high_slope
    )
    return energy, force


@numba.njit(cache=True, fastmath=False, parallel=False)
def evaluate_external(
    positions,
    charges,
    polarizabilities,
    tip_x,
    tip_y,
    tip_z,
    drive_mode,
    strength,
    corrugation,
    surface_enabled,
    forces_out,
    growth=1.0,
):
    """Write the external forces and return their total scalar energy.

    ``drive_mode`` is ``TIP_MODE`` or ``GRADIENT_MODE``.  The wrapper validates
    structured strengths in [0, 1]; the numeric kernel also supports the legacy
    unstructured tip interface, whose strength can be negative.  Tip coordinates
    are not read in gradient mode.  ``forces_out`` must have the same (N, 3)
    shape as ``positions`` and must not alias any input array.

    Numba uses one ordinary loop, no parallel worker pool, and no fast-math
    transformations.  Calling once with the prepared arrays before a timed MD
    run compiles or loads the cached signature; no simulation step is needed.
    """
    if drive_mode != TIP_MODE and drive_mode != GRADIENT_MODE:
        raise ValueError("drive mode must be gradient or tip")
    if not math.isfinite(strength) or not math.isfinite(corrugation):
        raise ValueError("field parameters must be finite")
    if drive_mode == TIP_MODE and not (
        math.isfinite(tip_x) and math.isfinite(tip_y) and math.isfinite(tip_z)
    ):
        raise ValueError("tip must be finite")

    amplitude = FIELD_SCALE * strength
    if not math.isfinite(growth) or growth < .1 or growth > 2.:
        raise ValueError("field growth must be in [0.1, 2]")
    root = math.sqrt(growth)
    base_field = GRADIENT_BASE_FIELD * strength/root
    slope = GRADIENT_FIELD_SLOPE * strength*root
    total_energy = 0.0
    for atom in range(positions.shape[0]):
        x, y, z = positions[atom, 0], positions[atom, 1], positions[atom, 2]
        charge, alpha = charges[atom], polarizabilities[atom]
        if not (
            math.isfinite(x)
            and math.isfinite(y)
            and math.isfinite(z)
            and math.isfinite(charge)
            and math.isfinite(alpha)
        ):
            raise ValueError("positions and response coefficients must be finite")
        if alpha < 0.0:
            raise ValueError("polarizabilities must be nonnegative")
        energy, fx, fy, fz = 0.0, 0.0, 0.0, 0.0
        if surface_enabled:
            height = z - ADSORPTION_HEIGHT
            exponential = math.exp(-ADSORPTION_ALPHA * height)
            squared = exponential * exponential
            energy = ADSORPTION_DEPTH * (squared - 2.0 * exponential)
            fz = 2.0 * ADSORPTION_DEPTH * ADSORPTION_ALPHA * (
                squared - exponential
            )
            if corrugation != 0.0:
                cx, sx = math.cos(_KX * x), math.sin(_KX * x)
                cy, sy = math.cos(_KY * y), math.sin(_KY * y)
                pattern = (2.0 * cx * cy + math.cos(2.0 * _KY * y)) / 3.0
                if CORRUGATION_DECAY == ADSORPTION_ALPHA:
                    envelope = corrugation * exponential
                else:
                    envelope = corrugation * math.exp(-CORRUGATION_DECAY * height)
                energy += envelope * pattern
                fx += envelope * (2.0 * _KX / 3.0) * sx * cy
                fy += envelope * (2.0 * _KY / 3.0) * (
                    cx * sy + math.sin(2.0 * _KY * y)
                )
                fz += CORRUGATION_DECAY * envelope * pattern

        boundary_energy, boundary_force = _boundary(x, TRACK_X_MIN, TRACK_X_MAX)
        energy += boundary_energy
        fx += boundary_force
        boundary_energy, boundary_force = _boundary(y, -LANE_HALF_WIDTH, LANE_HALF_WIDTH)
        energy += boundary_energy
        fy += boundary_force

        if drive_mode == GRADIENT_MODE:
            field_x = slope * z
            field_z = base_field + slope * x
            potential = -field_z * z
            response = GRADIENT_POLARIZATION_GAIN * alpha
            energy += charge * potential - 0.5 * response * (
                field_x**2 + field_z**2
            )
            fx += charge * field_x + response * slope * field_z
            fz += charge * field_z + response * slope**2 * z
        elif amplitude != 0.0:
            dx, dy, dz = x - tip_x, y - tip_y, z - tip_z
            distance_squared = dx * dx + dy * dy + dz * dz
            softened_squared = distance_squared + _TIP_SOFTENING_SQUARED
            inverse_radius = 1.0 / math.sqrt(softened_squared)
            inverse_squared = inverse_radius**2
            potential = amplitude * inverse_radius
            field_factor = amplitude * inverse_radius * inverse_squared
            field_squared = field_factor**2 * distance_squared
            energy += charge * potential - 0.5 * alpha * field_squared
            charge_force_factor = charge * field_factor
            induced_force_factor = (
                alpha
                * amplitude**2
                * (_TIP_SOFTENING_SQUARED - 2.0 * distance_squared)
                * inverse_squared**4
            )
            factor = charge_force_factor + induced_force_factor
            fx += factor * dx
            fy += factor * dy
            fz += factor * dz

        if not (
            math.isfinite(energy)
            and math.isfinite(fx)
            and math.isfinite(fy)
            and math.isfinite(fz)
        ):
            raise ValueError("coordinates exceed the supported field range")
        total_energy += energy
        forces_out[atom, 0] = fx
        forces_out[atom, 1] = fy
        forces_out[atom, 2] = fz
    if not math.isfinite(total_energy):
        raise ValueError("field calculation produced nonfinite energy or forces")
    return total_energy
