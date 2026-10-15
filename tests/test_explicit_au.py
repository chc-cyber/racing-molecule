"""Guards for complete molecules on commensurate periodic Au supercells."""
import unittest
from unittest.mock import patch
import numpy as np
from molecular_car.explicit_au import (
    choose_safe_cell, required_cell_lengths, build_relaxed_supercell, FrozenAuOracle,
)
from molecular_car.au_surface import minimum_image_xy, build_au111


class ExplicitCellTests(unittest.TestCase):
    def test_full_span_is_translation_invariant(self):
        molecule = np.array([[0., 0., 3.], [87., 5., 3.]])
        a = choose_safe_cell(molecule)
        b = choose_safe_cell(molecule+[1473., -829., 0.])
        self.assertEqual(a.lengths_xy, b.lengths_xy)
        self.assertTrue(a.expanded)
        self.assertGreater(a.lengths_xy[0], 95.)

    def test_unsafe_long_molecule_rejected_without_expansion(self):
        with self.assertRaisesRegex(ValueError, "Unsafe periodic cell"):
            choose_safe_cell([[0,0,3],[90,0,3]], auto_expand=False)

    def test_cutoff_and_normal_margins_are_explicit(self):
        span, required = required_cell_lengths([[0,0,3],[90,4,3]])
        np.testing.assert_allclose(span, [90,4])
        np.testing.assert_allclose(required, [98.000001,36.000001])
        plan = choose_safe_cell([[0,0,3],[5,5,3]])
        self.assertEqual((plan.nx, plan.ny), (14,8))
        self.assertIn("not a LAMMPS hard limit", plan.metadata()["cutoff_margin_policy"])

    def test_wrapping_does_not_change_local_normal_neighbours(self):
        angle = np.arange(6)*np.pi/3
        geometry = np.stack([1.4*np.cos(angle),1.4*np.sin(angle),np.full(6,3.3)],axis=1)
        plan = choose_safe_cell(geometry)
        lengths = np.array(plan.lengths_xy)
        original = geometry+[lengths[0]-.2,lengths[1]-.3,0]
        wrapped = original.copy();wrapped[:,:2] %= lengths
        a = np.linalg.norm(original[:,None,:]-original[None,:,:],axis=-1)
        b = np.linalg.norm(minimum_image_xy(wrapped[:,None,:]-wrapped[None,:,:],lengths),axis=-1)
        np.testing.assert_array_equal(a<2.,b<2.)

    def test_six_layers_and_common_ABC_registry(self):
        # Unit tests do not require a previously generated integration output.
        reference = build_au111().positions
        with patch('molecular_car.explicit_au.np.loadtxt', return_value=reference):
            small = build_relaxed_supercell(choose_safe_cell([[0,0,3],[5,5,3]]))
            big = build_relaxed_supercell(choose_safe_cell([[0,0,3],[5,5,3]], (57.69,59.96)))
        self.assertEqual(len(small.positions),1344)
        self.assertEqual(len(big.positions),2880)
        self.assertEqual(small.layer_count,6)
        for layer in range(6):
            a=small.positions[small.layers==layer]
            b=big.positions[big.layers==layer]
            self.assertAlmostEqual(a[0,2],b[0,2],places=13)
            np.testing.assert_allclose(a[0],b[0],atol=1e-13)

    def test_invalid_geometry_and_parameters(self):
        for coordinates in ([], [[1,2]], [[1,2,float('nan')]]):
            with self.assertRaises(ValueError):choose_safe_cell(coordinates)
        with self.assertRaises(ValueError):choose_safe_cell([[1,2,3]], skin=-1)

    def test_unapproved_opt_requires_explicit_diagnostic_opt_in(self):
        with self.assertRaisesRegex(ValueError, "failed per-atom force"):
            FrozenAuOracle([[0,0,3]], ['C'], style='opt')


if __name__=='__main__':unittest.main()
