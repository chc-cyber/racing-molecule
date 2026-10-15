"""Preparation rescue must reduce the true objective and honor force/budget gates."""
import unittest
import numpy as np
from molecular_car.preparation_bfgs import PreparationBudgetStop, dense_bfgs


class PreparationBfgsTests(unittest.TestCase):
    def test_ill_conditioned_harmonic_coordinates_converge_by_actual_force(self):
        stiffness = np.array([1e-4, 1e-2, 1., 10., 100., 8000.]).reshape(2, 3)
        initial = np.full((2, 3), .4)
        observed = []

        def evaluate(x):
            observed.append(x.copy())
            return float(.5 * np.sum(stiffness * x * x)), -stiffness * x

        final, report = dense_bfgs(initial, evaluate, 1e-4)
        self.assertLessEqual(np.linalg.norm(stiffness * final), 1e-4)
        self.assertEqual(report["termination_reason"], "force-tolerance")
        self.assertTrue(np.array_equal(initial, np.full((2, 3), .4)))
        energies = [record["potential_energy_kcal_mol"] for record in report["history"]]
        self.assertTrue(all(a >= b for a, b in zip(energies, energies[1:])))
        self.assertEqual(report["force_evaluations"], len(observed))

    def test_budget_abort_during_curvature_probe_keeps_accepted_coordinates(self):
        initial = np.ones((2, 3))
        calls = 0

        def evaluate(x):
            nonlocal calls
            calls += 1
            if calls == 4:
                raise PreparationBudgetStop("force-evaluation-budget")
            return float(.5 * np.sum(x * x)), -x

        final, report = dense_bfgs(initial, evaluate, 1e-4)
        np.testing.assert_array_equal(final, initial)
        self.assertEqual(report["termination_reason"], "force-evaluation-budget")
        self.assertEqual(report["iterations"], 0)

    def test_budget_abort_during_line_search_keeps_latest_accepted_step(self):
        initial = np.ones((1, 3))
        calls, points = 0, []

        def evaluate(x):
            nonlocal calls
            calls += 1
            points.append(x.copy())
            if calls == 9:
                raise PreparationBudgetStop("wall-time-budget")
            return float(.5 * np.sum(x * x)), -x

        final, report = dense_bfgs(initial, evaluate, 1e-4)
        np.testing.assert_array_equal(final, points[7])
        self.assertEqual(report["iterations"], 1)
        self.assertEqual(report["termination_reason"], "wall-time-budget")

    def test_energy_flat_return_never_substitutes_for_force_convergence(self):
        initial = np.ones((1, 3))
        final, report = dense_bfgs(initial, lambda x: (-1., np.ones_like(x)), 1e-4)
        np.testing.assert_array_equal(final, initial)
        self.assertEqual(report["termination_reason"], "line-search-stalled")

    def test_nonfinite_force_is_an_explicit_failure(self):
        with self.assertRaisesRegex(RuntimeError, "nonfinite"):
            dense_bfgs(np.ones((1, 3)), lambda x: (0., np.full_like(x, np.nan)), 1e-4)


if __name__ == "__main__":
    unittest.main()
