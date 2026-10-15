"""Continuous, deliberately effective surface and drive fields for the prototype.

Coordinates are in Angstrom, energies in kcal/mol, and forces in kcal/mol/Angstrom.
The surface parameters describe a game material, not a fitted gold force field.
Charges are expressed in elementary-charge units. ``polarizabilities`` are game
response coefficients compatible with the potential below, *not* polarizability
volumes in Angstrom cubed. In these effective units their dimensions are
charge**2 * Angstrom**2 / (kcal/mol).

The stationary broad-gradient mode uses a harmonic electrostatic potential,
phi = -(E0 + g*x)*z, zero on the plane z=0. Its polarization response is given
an explicit gameplay gain to provide useful speeds without making the field
acting on fixed partial charges excessively strong. This gain is not a measured
material property. Neither drive mode is a validated Au/STM model.

The tip potential is a softened point potential. Its induced contribution is
exactly -0.5 * alpha * |E|**2 for E = -grad(phi), rather than a separate attraction
that is labelled polarization. This remains a simplified linear-response model:
it does not include conductor images, charge transfer, or electronic excitation.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .track import LANE_HALF_WIDTH, TRACK_X_MIN, TRACK_X_MAX


ADSORPTION_DEPTH = 0.8
ADSORPTION_HEIGHT = 2.9
ADSORPTION_ALPHA = 1.0
SURFACE_SPACING = 2.88
CORRUGATION_DECAY = 1.0
TIP_SOFTENING = 1.5
# UI strength is nominally 0..1; this converts it to an effective phi amplitude.
FIELD_SCALE = 35.0
GRADIENT_BASE_FIELD = 0.5
GRADIENT_FIELD_SLOPE = 0.0015
GRADIENT_POLARIZATION_GAIN = 500.0
DEFAULT_GRADIENT_STRENGTH = 0.75
BOUNDARY_WIDTH = 1.0
BOUNDARY_STIFFNESS = 0.2


def gradient_coefficients(strength, growth=1.):
    """Keep E0*g at x=0 fixed; growth scales the relative along-track increase."""
    if isinstance(growth, bool):
        raise ValueError("电场沿程增长倍率必须是 0.1–2 的有限数值。")
    growth = float(growth)
    if not np.isfinite(growth) or not .1 <= growth <= 2.:
        raise ValueError("电场沿程增长倍率必须是 0.1–2 的有限数值。")
    root = np.sqrt(growth)
    return GRADIENT_BASE_FIELD*strength/root, GRADIENT_FIELD_SLOPE*strength*root


def _soft_boundary(
    coordinates: np.ndarray, lower: float, upper: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return smooth boundary energy and force without hard clipping positions."""
    low = (lower - coordinates) / BOUNDARY_WIDTH
    high = (coordinates - upper) / BOUNDARY_WIDTH
    low_soft = np.logaddexp(0.0, low)
    high_soft = np.logaddexp(0.0, high)
    # exp(-logaddexp(0, -x)) is a stable sigmoid, including large |x|.
    low_slope = np.exp(-np.logaddexp(0.0, -low))
    high_slope = np.exp(-np.logaddexp(0.0, -high))
    energy = (
        0.5
        * BOUNDARY_STIFFNESS
        * BOUNDARY_WIDTH**2
        * (low_soft**2 + high_soft**2)
    )
    force = BOUNDARY_STIFFNESS * BOUNDARY_WIDTH * (
        low_soft * low_slope - high_soft * high_slope
    )
    return energy, force


def energy_and_forces(
    positions: np.ndarray,
    charges: np.ndarray,
    polarizabilities: np.ndarray,
    tip: Mapping[str, float],
    corrugation: float,
    drive: Mapping[str, str | float] | None = None,
) -> tuple[float, np.ndarray]:
    """Return total external energy and a fresh N-by-3 force array.

    ``tip`` contains x, y, z in Angstrom and dimensionless ``strength`` (normally
    0..1). Its effective potential amplitude is FIELD_SCALE * strength, with
    units (kcal/mol)*Angstrom/charge. Negative strength reverses the charge force
    but leaves the induced force unchanged. ``corrugation`` is an energy
    amplitude in kcal/mol, with zero selecting a laterally smooth surface.

    ``drive=None`` preserves the original tip interface. A structured drive
    contains ``mode`` ("gradient" or "tip") and ``strength`` in [0, 1]. In tip
    mode its strength overrides tip["strength"]. Gradient mode does not read
    the tip position, a center of mass, or simulation time. It uses the fixed
    constants GRADIENT_BASE_FIELD and GRADIENT_FIELD_SLOPE, both scaled by
    strength, and GRADIENT_POLARIZATION_GAIN for the induced response only:
    U = q*phi - 0.5*gain*alpha*|E|**2. The base field has effective units of
    energy/(charge*length), and its slope has energy/(charge*length**2).

    All atoms retain independent surface and tip forces. Molecular bonds and
    atom-atom forces belong to the MD engine and are not calculated here.
    Boundaries are soft gameplay confinement at y=+/-250 and x=-10..1010 Angstrom.
    Their small, continuous tails avoid discontinuous contact forces.

    The softened tip has a field maximum at distance TIP_SOFTENING/sqrt(2).
    The normal operating tip height keeps the car outside this softened core.
    """
    positions = np.asarray(positions, dtype=np.float64)
    charges = np.asarray(charges, dtype=np.float64)
    polarizabilities = np.asarray(polarizabilities, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions must have shape (N, 3)")
    count = positions.shape[0]
    if charges.shape != (count,) or polarizabilities.shape != (count,):
        raise ValueError("charges and polarizabilities must each have shape (N,)")
    if not (
        np.isfinite(positions).all()
        and np.isfinite(charges).all()
        and np.isfinite(polarizabilities).all()
    ):
        raise ValueError("positions and response coefficients must be finite")
    if (polarizabilities < 0.0).any():
        raise ValueError("polarizabilities must be nonnegative")
    mode = "tip"
    drive_strength = None
    if drive is not None:
        if not isinstance(drive, Mapping) or set(drive) - {"mode", "strength", "growth"}:
            raise ValueError("drive accepts only mode, strength and growth")
        mode = drive.get("mode", "gradient")
        if not isinstance(mode, str) or mode not in {"gradient", "tip"}:
            raise ValueError("drive mode must be gradient or tip")
        value = drive.get("strength", DEFAULT_GRADIENT_STRENGTH)
        try:
            drive_strength = float(value)
        except (TypeError, ValueError) as error:
            raise ValueError("drive strength must be a finite number in [0, 1]") from error
        if isinstance(value, bool) or not np.isfinite(drive_strength) or not 0 <= drive_strength <= 1:
            raise ValueError("drive strength must be a finite number in [0, 1]")
    try:
        corrugation = float(corrugation)
    except (TypeError, ValueError) as error:
        raise ValueError("corrugation must be finite") from error
    if not np.isfinite(corrugation):
        raise ValueError("corrugation must be finite")
    if mode == "tip":
        try:
            tip_position = np.array([tip[axis] for axis in ("x", "y", "z")], dtype=float)
            strength = float(tip["strength"]) if drive_strength is None else drive_strength
            amplitude = FIELD_SCALE * strength
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("tip requires numeric x, y, z, and strength values") from error
        if not np.isfinite(tip_position).all() or not np.isfinite(amplitude):
            raise ValueError("tip must be finite")

    forces = np.zeros_like(positions)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            x, y, z = positions.T
            height = z - ADSORPTION_HEIGHT
            exp_height = np.exp(-ADSORPTION_ALPHA * height)
            energy = ADSORPTION_DEPTH * (exp_height**2 - 2.0 * exp_height)
            forces[:, 2] += 2.0 * ADSORPTION_DEPTH * ADSORPTION_ALPHA * (
                exp_height**2 - exp_height
            )

            if corrugation != 0.0:
                kx = 2.0 * np.pi / SURFACE_SPACING
                ky = kx / np.sqrt(3.0)
                cx, sx = np.cos(kx * x), np.sin(kx * x)
                cy, sy = np.cos(ky * y), np.sin(ky * y)
                # Three equal reciprocal-lattice cosines of a triangular plane.
                pattern = (2.0 * cx * cy + np.cos(2.0 * ky * y)) / 3.0
                envelope = corrugation * np.exp(-CORRUGATION_DECAY * height)
                energy += envelope * pattern
                forces[:, 0] += envelope * (2.0 * kx / 3.0) * sx * cy
                forces[:, 1] += envelope * (2.0 * ky / 3.0) * (
                    cx * sy + np.sin(2.0 * ky * y)
                )
                forces[:, 2] += CORRUGATION_DECAY * envelope * pattern

            for axis, lower, upper in (
                (0, TRACK_X_MIN, TRACK_X_MAX),
                (1, -LANE_HALF_WIDTH, LANE_HALF_WIDTH),
            ):
                boundary_energy, boundary_force = _soft_boundary(
                    positions[:, axis], lower, upper
                )
                energy += boundary_energy
                forces[:, axis] += boundary_force

            if mode == "gradient":
                base_field, slope = gradient_coefficients(drive_strength, drive.get("growth", 1.))
                field_x = slope * z
                field_z = base_field + slope * x
                potential = -field_z * z
                response = GRADIENT_POLARIZATION_GAIN * polarizabilities
                energy += charges * potential - 0.5 * response * (field_x**2 + field_z**2)
                forces[:, 0] += charges * field_x + response * slope * field_z
                forces[:, 2] += charges * field_z + response * slope**2 * z
            elif amplitude != 0.0:
                displacement = positions - tip_position
                distance_squared = np.einsum("ij,ij->i", displacement, displacement)
                softened_squared = distance_squared + TIP_SOFTENING**2
                inverse_radius = 1.0 / np.sqrt(softened_squared)
                inverse_squared = inverse_radius**2
                potential = amplitude * inverse_radius
                field_factor = amplitude * inverse_radius * inverse_squared
                field_squared = field_factor**2 * distance_squared
                energy += charges * potential - 0.5 * polarizabilities * field_squared
                charge_force_factor = charges * field_factor
                induced_force_factor = (
                    polarizabilities
                    * amplitude**2
                    * (TIP_SOFTENING**2 - 2.0 * distance_squared)
                    * inverse_squared**4
                )
                forces += (
                    charge_force_factor + induced_force_factor
                )[:, None] * displacement

            total_energy = float(np.sum(energy))
    except FloatingPointError as error:
        raise ValueError("coordinates exceed the supported field range") from error
    if not np.isfinite(total_energy) or not np.isfinite(forces).all():
        raise ValueError("field calculation produced nonfinite energy or forces")
    return total_energy, forces
