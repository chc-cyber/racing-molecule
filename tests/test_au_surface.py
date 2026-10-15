"""Geometry checks for the independent explicit-gold reference experiment."""

import unittest
import numpy as np
from molecular_car.au_surface import build_au111, minimum_image_xy, layer_statistics


class AuGeometryTests(unittest.TestCase):
    def test_fcc_coordination_and_periodic_seams(self):
        slab = build_au111(nx=5, ny=4)
        # Interior FCC layers have 12 nearest neighbours, including across seams.
        for index in (80, 91, 119, 120, 159):
            distances = np.linalg.norm(minimum_image_xy(
                slab.positions - slab.positions[index], slab.lengths_xy), axis=1)
            nearest = np.count_nonzero(np.isclose(distances, slab.nearest_spacing, atol=1e-10))
            self.assertEqual(nearest, 12)
        self.assertEqual(len(slab.positions), 5 * 4 * 2 * 6)

    def test_layers_and_fixed_bottom(self):
        slab = build_au111()
        self.assertEqual(len(slab.positions), 2880)
        self.assertEqual(slab.frozen_count, 960)
        stats = layer_statistics(slab)
        np.testing.assert_allclose(stats["interlayer_gaps_angstrom"], 4.08 / np.sqrt(3))
        self.assertEqual(stats["layers"][-1]["mean_z_angstrom"], 0.0)
        self.assertTrue(np.all(slab.positions[:, :2] >= 0))
        self.assertTrue(np.all(slab.positions[:, :2] < slab.lengths_xy))

    def test_translation_to_top_zero_preserves_layer_gaps(self):
        slab = build_au111(nx=3, ny=3)
        moved = slab.with_positions(slab.positions + [0, 0, 11.3], top_zero=True)
        np.testing.assert_allclose(moved.positions, slab.positions, atol=1e-14)

    def test_invalid_dimensions(self):
        for options in ({"nx": True}, {"layers": 2}, {"frozen_layers": 7},
                        {"ny": 0}, {"lattice_constant": float("nan")}):
            with self.assertRaises(ValueError):
                build_au111(**options)


if __name__ == "__main__":
    unittest.main()
