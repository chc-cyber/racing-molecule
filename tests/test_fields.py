"""Check external field derivatives and the physical signs used by gameplay."""

import unittest

import numpy as np

from molecular_car.fields import (
    ADSORPTION_HEIGHT,
    GRADIENT_BASE_FIELD,
    GRADIENT_FIELD_SLOPE,
    GRADIENT_POLARIZATION_GAIN,
    LANE_HALF_WIDTH,
    TRACK_X_MAX,
    energy_and_forces,
)
from molecular_car.track import TRACK


class FieldTests(unittest.TestCase):
    def setUp(self):
        self.tip = {"x": 5.0, "y": -1.0, "z": 7.0, "strength": 0.75}

    def assert_gradient(self, positions, charges, polarizabilities, corrugation, drive=None):
        positions = np.array(positions, dtype=float)
        _, forces = energy_and_forces(
            positions, charges, polarizabilities, self.tip, corrugation, drive
        )
        step = 1.0e-5
        numerical = np.zeros_like(positions)
        for atom in range(len(positions)):
            for axis in range(3):
                plus = positions.copy()
                minus = positions.copy()
                plus[atom, axis] += step
                minus[atom, axis] -= step
                plus_energy, _ = energy_and_forces(
                    plus, charges, polarizabilities, self.tip, corrugation, drive
                )
                minus_energy, _ = energy_and_forces(
                    minus, charges, polarizabilities, self.tip, corrugation, drive
                )
                numerical[atom, axis] = -(plus_energy - minus_energy) / (2.0 * step)
        np.testing.assert_allclose(forces, numerical, rtol=2.0e-7, atol=2.0e-7)

    def test_combined_force_is_negative_energy_gradient(self):
        self.assert_gradient(
            [[0.37, 0.31, 2.95], [7.3, LANE_HALF_WIDTH + 0.2, 3.4], [-10.4, -LANE_HALF_WIDTH - 0.3, 2.55]],
            np.array([0.25, -0.4, 0.1]),
            np.array([1.0, 0.5, 1.2]),
            0.13,
        )

    def test_charge_field_gradient_without_induced_response(self):
        self.assert_gradient(
            [[3.2, -0.4, 3.1], [7.0, 2.2, 3.6]],
            np.array([0.6, -0.6]),
            np.zeros(2),
            0.0,
        )

    def test_induced_response_gradient_without_charges(self):
        self.assert_gradient(
            [[2.7, 0.8, 3.25], [6.1, -2.3, 4.2]],
            np.zeros(2),
            np.array([1.0, 0.7]),
            0.0,
        )

    def test_corrugation_has_consistent_vertical_force(self):
        positions = np.array([[0.37, 0.31, 2.95]])
        zero = np.zeros(1)
        _, flat_force = energy_and_forces(positions, zero, zero, self.tip, 0.0)
        _, corrugated_force = energy_and_forces(positions, zero, zero, self.tip, 0.13)
        step = 1.0e-5
        plus, minus = positions.copy(), positions.copy()
        plus[0, 2] += step
        minus[0, 2] -= step

        def corrugation_energy(point):
            curved, _ = energy_and_forces(point, zero, zero, self.tip, 0.13)
            flat, _ = energy_and_forces(point, zero, zero, self.tip, 0.0)
            return curved - flat

        derivative_force = -(corrugation_energy(plus) - corrugation_energy(minus)) / (
            2.0 * step
        )
        self.assertGreater(abs(derivative_force), 0.01)
        self.assertAlmostEqual(
            corrugated_force[0, 2] - flat_force[0, 2], derivative_force, places=8
        )

    def test_zero_strength_disables_both_tip_contributions(self):
        positions = np.array([[1.1, 0.2, 3.0], [2.0, -0.4, 3.2]])
        tip_off = {**self.tip, "strength": 0.0}
        loaded = energy_and_forces(
            positions, np.array([0.6, -0.3]), np.ones(2), tip_off, 0.08
        )
        neutral = energy_and_forces(positions, np.zeros(2), np.zeros(2), tip_off, 0.08)
        self.assertEqual(loaded[0], neutral[0])
        np.testing.assert_array_equal(loaded[1], neutral[1])

    def test_induced_tip_force_attracts_outside_softened_core(self):
        positions = np.array([[3.0, 0.0, 3.0]])
        _, loaded_force = energy_and_forces(
            positions, np.zeros(1), np.ones(1), self.tip, 0.0
        )
        _, flat_force = energy_and_forces(
            positions, np.zeros(1), np.zeros(1), self.tip, 0.0
        )
        displacement = positions[0] - np.array([self.tip[axis] for axis in "xyz"])
        self.assertLess(np.dot(loaded_force[0] - flat_force[0], displacement), 0.0)

    def test_induced_energy_fades_and_is_independent_of_strength_sign(self):
        near = np.array([[self.tip["x"], self.tip["y"], 3.0]])
        far = np.array([[self.tip["x"], self.tip["y"], 20.0]])

        def tip_energy(point):
            induced, _ = energy_and_forces(point, np.zeros(1), np.ones(1), self.tip, 0)
            baseline, _ = energy_and_forces(point, np.zeros(1), np.zeros(1), self.tip, 0)
            return induced - baseline

        self.assertLess(tip_energy(near), 0.0)
        self.assertLess(abs(tip_energy(far)), abs(tip_energy(near)) / 10.0)
        positive = energy_and_forces(near, np.zeros(1), np.ones(1), self.tip, 0)
        negative = energy_and_forces(
            near, np.zeros(1), np.ones(1), {**self.tip, "strength": -0.75}, 0
        )
        self.assertEqual(positive[0], negative[0])
        np.testing.assert_array_equal(positive[1], negative[1])

    def test_surface_binds_and_repels_at_the_plane(self):
        tip_off = {**self.tip, "strength": 0.0}
        energy, force = energy_and_forces(
            np.array([[20.0, 0.0, ADSORPTION_HEIGHT]]),
            np.zeros(1), np.zeros(1), tip_off, 0.0,
        )
        self.assertAlmostEqual(energy, -0.8, places=10)
        self.assertAlmostEqual(force[0, 2], 0.0, places=12)
        plane_energy, plane_force = energy_and_forces(
            np.array([[20.0, 0.0, 0.0]]), np.zeros(1), np.zeros(1), tip_off, 0.0
        )
        self.assertGreater(plane_energy, 100.0)
        self.assertGreater(plane_force[0, 2], 100.0)

    def test_soft_track_boundaries_push_inward(self):
        positions = np.array(
            [[20.0, LANE_HALF_WIDTH + 6.0, 3.0], [20.0, -LANE_HALF_WIDTH - 6.0, 3.0],
             [-20.0, 0.0, 3.0], [TRACK_X_MAX + 10.0, 0.0, 3.0]]
        )
        _, force = energy_and_forces(
            positions, np.zeros(4), np.zeros(4), {**self.tip, "strength": 0.0}, 0
        )
        self.assertLess(force[0, 1], 0)
        self.assertGreater(force[1, 1], 0)
        self.assertGreater(force[2, 0], 0)
        self.assertLess(force[3, 0], 0)

    def test_empty_input_and_nonmutation(self):
        energy, force = energy_and_forces(
            np.empty((0, 3)), np.empty(0), np.empty(0), self.tip, 0.08
        )
        self.assertEqual(energy, 0.0)
        self.assertEqual(force.shape, (0, 3))
        positions = np.array([[1.0, 0.0, 3.0]])
        original = positions.copy()
        energy_and_forces(positions, np.zeros(1), np.ones(1), self.tip, 0.08)
        np.testing.assert_array_equal(positions, original)

    def test_invalid_coordinates_fail_instead_of_clipping_force(self):
        for positions in (np.array([[0.0, 0.0, np.nan]]), np.array([[0.0, 0.0, -1000.0]])):
            with self.assertRaises(ValueError):
                energy_and_forces(positions, np.zeros(1), np.ones(1), self.tip, 0.08)
        with self.assertRaises(ValueError):
            energy_and_forces(np.zeros((2, 3)), np.zeros(1), np.ones(2), self.tip, 0.08)

    def test_broad_gradient_matches_energy_derivative_across_full_track(self):
        self.assert_gradient(
            [[0.37, 0.31, 2.95], [100.3, 20.2, 3.4], [197.4, -11.3, 2.55]],
            np.array([0.07, -0.04, 0.02]), np.array([1.0, 0.5, 1.2]), 0.08,
            {"mode": "gradient", "strength": 0.75},
        )

    def test_broad_gradient_charge_force_and_induction_separately(self):
        positions = [[11.3, -2.0, 3.1], [176.4, 7.0, 3.4]]
        drive = {"mode": "gradient", "strength": 0.6}
        self.assert_gradient(positions, np.array([0.3, -0.2]), np.zeros(2), 0.0, drive)
        self.assert_gradient(positions, np.zeros(2), np.array([1.0, 0.7]), 0.0, drive)

    def gradient_only(self, positions, charges, polarizabilities, strength=0.75, tip=None):
        drive = {"mode": "gradient", "strength": strength}
        arguments = (positions, charges, polarizabilities, self.tip if tip is None else tip, 0.0)
        loaded_energy, loaded_force = energy_and_forces(*arguments, drive)
        baseline_energy, baseline_force = energy_and_forces(
            *arguments, {"mode": "gradient", "strength": 0.0}
        )
        return loaded_energy - baseline_energy, loaded_force - baseline_force

    def test_broad_gradient_ignores_tip_position_and_has_no_lateral_steering(self):
        positions = np.array([[3.0, -15.0, 3.0], [140.0, 18.0, 3.2]])
        charges, alphas = np.array([0.1, -0.1]), np.ones(2)
        original = self.gradient_only(positions, charges, alphas)
        unrelated_tip = {"x": -1000.0, "y": 1500.0, "z": 100.0, "strength": 1.0}
        changed = self.gradient_only(positions, charges, alphas, tip=unrelated_tip)
        self.assertEqual(original[0], changed[0])
        np.testing.assert_array_equal(original[1], changed[1])
        np.testing.assert_array_equal(original[1][:, 1], np.zeros(2))
        # A stationary broad field has no dependency on any tip configuration.
        no_tip = self.gradient_only(positions, charges, alphas, tip={})
        self.assertEqual(original[0], no_tip[0])
        np.testing.assert_array_equal(original[1], no_tip[1])

    def test_translated_car_samples_fixed_gradient_without_recentering_field(self):
        positions = np.array([[1.0, 0.0, 3.0], [2.0, 1.0, 3.1]])
        translated = positions + np.array([100.0, 0.0, 0.0])
        alphas = np.array([1.0, 0.2])
        _, start_force = self.gradient_only(positions, np.zeros(2), alphas)
        _, shifted_force = self.gradient_only(translated, np.zeros(2), alphas)
        slope = GRADIENT_FIELD_SLOPE * 0.75
        expected_change = GRADIENT_POLARIZATION_GAIN * alphas * slope**2 * 100.0
        np.testing.assert_allclose(shifted_force[:, 0] - start_force[:, 0], expected_change)
        np.testing.assert_allclose(shifted_force[:, 2], start_force[:, 2], atol=1e-14)
        self.assertTrue((start_force[:, 0] > 0).all())

    def test_gradient_strength_zero_turns_off_force_and_negative_is_rejected(self):
        positions = np.array([[3.0, 0.0, 3.0]])
        energy, force = self.gradient_only(positions, np.array([0.2]), np.ones(1), strength=0)
        self.assertEqual(energy, 0.0)
        np.testing.assert_array_equal(force, np.zeros((1, 3)))
        for strength in (-0.1, 1.1, float("nan"), True):
            with self.assertRaises(ValueError):
                energy_and_forces(
                    positions, np.zeros(1), np.ones(1), self.tip, 0,
                    {"mode": "gradient", "strength": strength},
                )

    def test_gradient_strength_scales_charge_linearly_and_induction_quadratically(self):
        positions = np.array([[43.0, 0.0, 3.0]])
        for charges, alphas, power in ((np.array([0.2]), np.zeros(1), 1),
                                      (np.zeros(1), np.ones(1), 2)):
            full = self.gradient_only(positions, charges, alphas, strength=1.0)
            half = self.gradient_only(positions, charges, alphas, strength=0.5)
            self.assertAlmostEqual(half[0], full[0] * 0.5**power, places=11)
            np.testing.assert_allclose(half[1], full[1] * 0.5**power, atol=1e-14)

    def test_structured_tip_mode_preserves_original_driver(self):
        positions = np.array([[4.0, 0.0, 3.0]])
        charges, alphas = np.array([0.1]), np.ones(1)
        old = energy_and_forces(positions, charges, alphas, self.tip, 0.08)
        structured = energy_and_forces(
            positions, charges, alphas, self.tip, 0.08,
            {"mode": "tip", "strength": self.tip["strength"]},
        )
        self.assertEqual(old[0], structured[0])
        np.testing.assert_array_equal(old[1], structured[1])

    def test_gradient_parameters_produce_expected_start_and_finish_force(self):
        positions = np.array([[TRACK["start_x"], 0.0, 3.0], [TRACK["finish_x"], 0.0, 3.0]])
        _, force = self.gradient_only(positions, np.zeros(2), np.ones(2))
        slope = GRADIENT_FIELD_SLOPE * 0.75
        base = GRADIENT_BASE_FIELD * 0.75
        expected = GRADIENT_POLARIZATION_GAIN * slope * (base + slope * positions[:, 0])
        np.testing.assert_allclose(force[:, 0], expected)
        self.assertAlmostEqual(force[0, 0], 0.2109375)
        self.assertAlmostEqual(force[1, 0], 0.3375)


if __name__ == "__main__":
    unittest.main()
