"""Check the compiled external potentials against the preserved reference.

The reference is the pre-optimization NumPy implementation in ``fields.py``.
Cartesian energy derivatives are checked independently, including the softened
tip core, soft track walls, and the surface-disabled periodic-Au use case.
"""

import unittest

import numpy as np

from molecular_car.fields import (
    ADSORPTION_ALPHA,
    ADSORPTION_DEPTH,
    ADSORPTION_HEIGHT,
    LANE_HALF_WIDTH,
    TRACK_X_MIN,
    TRACK_X_MAX,
    energy_and_forces,
)
from molecular_car.fields_kernel import GRADIENT_MODE, TIP_MODE, evaluate_external


class CompiledFieldTests(unittest.TestCase):
    tip = {"x": 5.0, "y": -1.0, "z": 7.0, "strength": 0.75}

    def evaluate(self, positions, charges, alphas, corrugation, mode=TIP_MODE,
                 strength=0.75, surface=True, tip=None):
        positions = np.asarray(positions, dtype=np.float64)
        charges = np.asarray(charges, dtype=np.float64)
        alphas = np.asarray(alphas, dtype=np.float64)
        tip = self.tip if tip is None else tip
        output = np.full_like(positions, np.nan)
        energy = evaluate_external(
            positions, charges, alphas, tip.get("x", np.nan),
            tip.get("y", np.nan), tip.get("z", np.nan), mode, strength,
            corrugation, surface, output,
        )
        return energy, output

    def reference(self, positions, charges, alphas, corrugation, mode=TIP_MODE,
                  strength=0.75, surface=True, tip=None):
        tip = self.tip if tip is None else tip
        # The legacy tip interface intentionally also accepts a negative
        # strength, whereas a structured drive validates strengths in [0, 1].
        if mode == TIP_MODE:
            drive = None
            tip = {**tip, "strength": strength}
        else:
            drive = {"mode": "gradient", "strength": strength}
        energy, forces = energy_and_forces(
            positions, charges, alphas, tip, corrugation if surface else 0.0,
            drive,
        )
        if not surface:
            exponential = np.exp(-ADSORPTION_ALPHA * (
                positions[:, 2] - ADSORPTION_HEIGHT
            ))
            energy -= float(np.sum(ADSORPTION_DEPTH * (
                exponential**2 - 2.0 * exponential
            )))
            forces[:, 2] -= 2.0 * ADSORPTION_DEPTH * ADSORPTION_ALPHA * (
                exponential**2 - exponential
            )
        return energy, forces

    def assert_reference(self, positions, charges, alphas, corrugation, **kwargs):
        reference = self.reference(positions, charges, alphas, corrugation, **kwargs)
        compiled = self.evaluate(positions, charges, alphas, corrugation, **kwargs)
        np.testing.assert_allclose(compiled[0], reference[0], rtol=5e-14, atol=2e-10)
        np.testing.assert_allclose(compiled[1], reference[1], rtol=5e-14, atol=2e-11)

    def test_random_track_reference_both_drivers_and_surface_modes(self):
        rng = np.random.default_rng(761002)
        for count in (1, 12, 52, 212, 2049):
            positions = rng.uniform([TRACK_X_MIN - 20, -LANE_HALF_WIDTH - 30, 0.0],
                                    [TRACK_X_MAX + 30, LANE_HALF_WIDTH + 30, 30.0], (count, 3))
            charges = rng.uniform(-0.9, 0.9, count)
            alphas = rng.uniform(0.0, 2.0, count)
            for mode in (TIP_MODE, GRADIENT_MODE):
                for strength in (0.0, 0.25, 0.75, 1.0):
                    for surface in (False, True):
                        with self.subTest(count=count, mode=mode, strength=strength,
                                          surface=surface):
                            self.assert_reference(
                                positions, charges, alphas, 0.08, mode=mode,
                                strength=strength, surface=surface,
                            )

    def test_tip_core_boundaries_and_corrugation_sign_reference(self):
        tip = np.array([self.tip[axis] for axis in "xyz"])
        positions = np.vstack([
            tip, tip + [0.2, 0.3, 0.4], tip + [1.5 / np.sqrt(2), 0.0, 0.0],
            [[TRACK_X_MIN, -LANE_HALF_WIDTH, 2.9], [TRACK_X_MAX, LANE_HALF_WIDTH, 2.9],
             [TRACK_X_MIN - .001, -LANE_HALF_WIDTH - .001, 3.1],
             [TRACK_X_MAX + .001, LANE_HALF_WIDTH + .001, 3.1],
             [-1e6, -1e6, 20.0], [1e6, 1e6, 20.0],
             [40.0, 0.0, -10.0]],
        ])
        charges = np.linspace(-0.7, 0.7, len(positions))
        alphas = np.linspace(0.0, 1.1, len(positions))
        for corrugation in (-0.13, 0.0, 0.13):
            for strength in (-0.75, 0.0, 0.75):
                with self.subTest(corrugation=corrugation, strength=strength):
                    self.assert_reference(positions, charges, alphas, corrugation,
                                          strength=strength)

    def test_complete_cartesian_force_is_negative_energy_gradient(self):
        positions = np.array([
            [5.0, -1.0, 7.0], [5.3, -0.6, 7.2],
            [TRACK_X_MIN - .3, -LANE_HALF_WIDTH - .2, 2.75],
            [TRACK_X_MAX + .4, LANE_HALF_WIDTH + .3, 3.2],
            [100.2, 17.3, 3.4],
        ])
        charges = np.array([0.6, -0.4, 0.3, -0.2, 0.1])
        alphas = np.array([0.8, 1.1, 0.3, 0.9, 0.5])
        step = 1e-5
        for mode in (TIP_MODE, GRADIENT_MODE):
            for surface in (False, True):
                with self.subTest(mode=mode, surface=surface):
                    _, forces = self.evaluate(positions, charges, alphas, 0.13,
                                              mode=mode, surface=surface)
                    numerical = np.zeros_like(positions)
                    for atom in range(len(positions)):
                        for axis in range(3):
                            plus, minus = positions.copy(), positions.copy()
                            plus[atom, axis] += step
                            minus[atom, axis] -= step
                            ep, _ = self.evaluate(plus, charges, alphas, 0.13,
                                                  mode=mode, surface=surface)
                            em, _ = self.evaluate(minus, charges, alphas, 0.13,
                                                  mode=mode, surface=surface)
                            numerical[atom, axis] = -(ep - em) / (2.0 * step)
                    np.testing.assert_allclose(forces, numerical,
                                               rtol=2e-7, atol=4e-7)

    def test_gradient_mode_does_not_read_tip_position(self):
        positions = np.array([[0.1, -0.4, 2.9], [800.0, 3.0, 3.1]])
        charges, alphas = np.array([0.1, -0.1]), np.ones(2)
        real_tip = self.evaluate(positions, charges, alphas, 0.08, mode=GRADIENT_MODE)
        no_tip = self.evaluate(positions, charges, alphas, 0.08,
                               mode=GRADIENT_MODE, tip={})
        self.assertEqual(real_tip[0], no_tip[0])
        np.testing.assert_array_equal(real_tip[1], no_tip[1])

    def test_zero_tip_strength_has_no_charge_or_polarization_force(self):
        positions = np.array([[1.0, 0.0, 3.0], [5.0, -1.0, 7.0]])
        loaded = self.evaluate(positions, np.array([0.4, -0.3]), np.ones(2),
                               0.08, strength=0.0)
        neutral = self.evaluate(positions, np.zeros(2), np.zeros(2),
                                0.08, strength=0.0)
        self.assertEqual(loaded[0], neutral[0])
        np.testing.assert_array_equal(loaded[1], neutral[1])

    def test_empty_input_and_reused_buffer_do_not_mutate_inputs(self):
        energy, forces = self.evaluate(np.empty((0, 3)), np.empty(0),
                                       np.empty(0), 0.08)
        self.assertEqual(energy, 0.0)
        self.assertEqual(forces.shape, (0, 3))
        positions = np.array([[2.3, 4.1, 3.0], [12.0, -9.0, 4.3]])
        charges, alphas = np.array([0.25, -0.25]), np.ones(2)
        original = [array.copy() for array in (positions, charges, alphas)]
        output = np.full_like(positions, 1e8)
        arguments = (positions, charges, alphas, 5.0, -1.0, 7.0,
                     GRADIENT_MODE, 0.75, 0.08, True, output)
        first_energy = evaluate_external(*arguments)
        first_forces = output.copy()
        second_energy = evaluate_external(*arguments)
        self.assertEqual(first_energy, second_energy)
        np.testing.assert_array_equal(first_forces, output)
        for unchanged, source in zip((positions, charges, alphas), original):
            np.testing.assert_array_equal(unchanged, source)

    def test_dynamic_nonfinite_and_overflow_fail_instead_of_clipping(self):
        for positions in (np.array([[0.0, 0.0, np.nan]]),
                          np.array([[0.0, 0.0, np.inf]]),
                          np.array([[0.0, 0.0, -1000.0]]),
                          np.array([[1e308, 0.0, 3.0]])):
            with self.subTest(position=positions[0].tolist()):
                with self.assertRaises(ValueError):
                    self.evaluate(positions, np.zeros(1), np.ones(1), 0.08)
        for charges, alphas in (([np.nan], [1.0]), ([0.0], [np.inf]),
                               ([0.0], [-0.1])):
            with self.assertRaises(ValueError):
                self.evaluate([[0.0, 0.0, 3.0]], charges, alphas, 0.08)
        for kwargs in ({"mode": 3}, {"strength": np.inf},
                       {"tip": {"x": np.inf, "y": 0.0, "z": 7.0}}):
            with self.assertRaises(ValueError):
                self.evaluate([[0.0, 0.0, 3.0]], [0.0], [1.0], 0.08, **kwargs)
        with self.assertRaises(ValueError):
            self.evaluate([[0.0, 0.0, 3.0]], [0.0], [1.0], np.nan)


if __name__ == "__main__":
    unittest.main()
