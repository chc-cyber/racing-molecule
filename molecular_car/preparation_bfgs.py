"""Budgeted, preparation-only BFGS rescue for disparate molecular stiffnesses.

The finite-difference Hessian only initializes an inverse preconditioner. All
energies and forces are supplied by the unchanged conservative race objective.
There is no auxiliary potential, MD integration, or force-tolerance relaxation.
"""
from __future__ import annotations

import numpy as np

BFGS_OPTIONS = {"maximum_iterations": 400, "hessian_step_angstrom": 1e-4,
                "maximum_atom_displacement_angstrom": .1,
                "maximum_line_search_trials": 30, "armijo_coefficient": 1e-4,
                "curvature_relative_threshold": 1e-12,
                "diagonal_shift_relative_start": 1e-8,
                "maximum_cholesky_attempts": 12}


class PreparationBudgetStop(Exception):
    """Abort an optimizer at its current accepted coordinate checkpoint."""


def dense_bfgs(initial_positions, evaluate, force_tolerance, *, options=None):
    """Return the last accepted geometry, including on a soft-budget stop.

    ``evaluate(x)`` returns total energy and *force*, and enforces the shared
    preparation budgets. The caller must independently refresh final-coordinate
    forces; neither this optimizer's return nor a line-search stop admits a race.
    """
    options = {**BFGS_OPTIONS, **(options or {})}
    current = np.array(initial_positions, dtype=np.float64, copy=True)
    dimension = current.size
    iterations, calls, shift, history = 0, 0, 0., []
    reason = "iteration-limit"

    def value(x):
        nonlocal calls
        energy, force = evaluate(x)
        calls += 1
        force = np.array(force, dtype=np.float64, copy=True)
        if force.shape != current.shape or not np.isfinite(force).all() or not np.isfinite(energy):
            raise RuntimeError("Preparation BFGS received nonfinite or inconsistent energy/forces.")
        return float(energy), force

    try:
        energy, force = value(current)
        if np.linalg.norm(force) <= force_tolerance:
            reason = "force-tolerance"
        else:
            hessian = np.empty((dimension, dimension))
            delta = options["hessian_step_angstrom"]
            for component in range(dimension):
                plus, minus = current.copy(), current.copy()
                plus.ravel()[component] += delta
                minus.ravel()[component] -= delta
                hessian[:, component] = -(value(plus)[1] - value(minus)[1]).ravel() / (2 * delta)
            hessian = (hessian + hessian.T) * .5
            scale = max(1., float(np.max(np.abs(np.diag(hessian)))))
            identity = np.eye(dimension)
            for attempt in range(options["maximum_cholesky_attempts"]):
                try:
                    lower = np.linalg.cholesky(hessian + shift * identity)
                    # The absolute pivot guard is also used by the browser solver.
                    if float(np.min(np.diag(lower) ** 2)) <= 1e-12 * scale:
                        raise np.linalg.LinAlgError("Near-singular Hessian pivot")
                    inverse_lower = np.linalg.solve(lower, identity)
                    inverse = inverse_lower.T @ inverse_lower
                    break
                except np.linalg.LinAlgError:
                    shift = (options["diagonal_shift_relative_start"] * scale if shift == 0.
                             else shift * 10.)
            else:
                return current, {"iterations": 0, "force_evaluations": calls,
                                 "termination_reason": "curvature-factorization-failed",
                                 "hessian_diagonal_shift": shift, "history": history}

            for _ in range(options["maximum_iterations"]):
                if np.linalg.norm(force) <= force_tolerance:
                    reason = "force-tolerance"
                    break
                direction = (inverse @ force.ravel()).reshape(current.shape)
                if not np.isfinite(direction).all():
                    reason = "nonfinite-search-direction"
                    break
                slope = -float(np.sum(force * direction))
                if slope >= 0:
                    reason = "non-descent-search-direction"
                    break
                size = float(np.linalg.norm(direction, axis=1).max())
                if size == 0:
                    reason = "zero-search-direction"
                    break
                factor = min(1., options["maximum_atom_displacement_angstrom"] / size)
                direction *= factor
                slope *= factor
                alpha = 1.
                for _ in range(options["maximum_line_search_trials"]):
                    trial = current + alpha * direction
                    trial_energy, trial_force = value(trial)
                    if trial_energy <= energy + options["armijo_coefficient"] * alpha * slope:
                        break
                    alpha *= .5
                else:
                    reason = "line-search-stalled"
                    break
                displacement = (trial - current).ravel()
                delta_gradient = (force - trial_force).ravel()
                curvature = float(displacement @ delta_gradient)
                if curvature > (options["curvature_relative_threshold"] * np.linalg.norm(displacement)
                                * np.linalg.norm(delta_gradient)):
                    product = inverse @ delta_gradient
                    inverse += ((1 + float(delta_gradient @ product) / curvature)
                                * np.outer(displacement, displacement) / curvature)
                    inverse -= ((np.outer(product, displacement) + np.outer(displacement, product)) / curvature)
                current, energy, force = trial, trial_energy, trial_force
                iterations += 1
                history.append({"iteration": iterations, "potential_energy_kcal_mol": energy,
                                "global_force_norm_kcal_mol_A": float(np.linalg.norm(force)),
                                "line_search_alpha": alpha})
            if np.linalg.norm(force) <= force_tolerance:
                reason = "force-tolerance"
    except PreparationBudgetStop as exc:
        reason = str(exc)
    return current, {"iterations": iterations, "force_evaluations": calls,
                     "termination_reason": reason, "hessian_diagonal_shift": shift,
                     "history": history}
