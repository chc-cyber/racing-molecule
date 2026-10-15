"""Conservative game-normal adaptation; original SAIP remains an explicit oracle.

Small, deterministic force checks include folds, sp/sp3, immutable initial
classes, the old distance cutoff, a degenerate plane, and Fourier fallback.
"""
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from molecular_car.periodic_au import (
    BONDED_NORMAL_POLICY, NORMAL_AREA_EPSILON, ROOT, PeriodicAuPotential,
    _bonded_normal_geometry, _normal_geometry, reference_layer_heights,
)
from molecular_car.fourier_au import FourierAuPotential


def benzene():
    phi = np.arange(6)*np.pi/3
    xy = np.stack([np.cos(phi), np.sin(phi)], axis=1)
    positions = np.column_stack([np.concatenate([1.4*xy, 2.48*xy]), np.zeros(12)])
    angle = .31
    rotation = np.array([[1., 0., 0.], [0., np.cos(angle), -np.sin(angle)],
                         [0., np.sin(angle), np.cos(angle)]])
    positions = positions@rotation.T+[.379, .713, 4.]
    bonds = ([dict(a=i+1, b=(i+1)%6+1, order=1.5) for i in range(6)]
             +[dict(a=i+1, b=i+7, order=1.) for i in range(6)])
    return positions, bonds


class BondedSurfaceNormalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="bonded-normal-", dir=ROOT/"outputs")
        cls.elements = ["C"]*6+["H"]*6
        cls.positions, cls.bonds = benzene()
        cls.options = dict(grid_nx=8, grid_ny=12, grid_dz=.25, cache_dir=cls.temporary.name)
        cls.fourier_options = dict(cls.options, fourier_cache_dir=cls.temporary.name,
            angular_step=.5, chart_bound=2., chart_blend=.4,
            height_min=2.8, height_max=4.6, height_blend=.3, height_step=.15)
        cls.game = PeriodicAuPotential(cls.elements, cls.positions, bonds=cls.bonds,
            normal_policy=BONDED_NORMAL_POLICY, accelerator="numpy", **cls.options)
        cls.compiled = PeriodicAuPotential(cls.elements, cls.positions, bonds=cls.bonds,
            normal_policy=BONDED_NORMAL_POLICY, accelerator="numba", **cls.options)
        cls.proxy = FourierAuPotential(cls.elements, cls.positions, bonds=cls.bonds,
            normal_policy=BONDED_NORMAL_POLICY, accelerator="numpy", **cls.fourier_options)
        cls.compiled_proxy = FourierAuPotential(cls.elements, cls.positions, bonds=cls.bonds,
            normal_policy=BONDED_NORMAL_POLICY, accelerator="numba", **cls.fourier_options)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def assert_fd(self, surface, positions, axes=None, step=1e-5, atol=3e-6, rtol=0.):
        _, force = surface.energy_and_forces_eV(positions)
        axes = axes or [(i, j) for i in range(len(positions)) for j in range(3)]
        for i, j in axes:
            plus, minus = positions.copy(), positions.copy()
            plus[i, j] += step
            minus[i, j] -= step
            derivative = (surface.energy_and_forces_eV(plus)[0]-surface.energy_and_forces_eV(minus)[0])/(2*step)
            np.testing.assert_allclose(force[i, j], -derivative, atol=atol, rtol=rtol,
                                       err_msg=f"{surface.__class__.__name__} {i}/{j}")

    def test_aromatic_reference_equivalence_numpy_and_numba(self):
        positions = self.positions.copy()
        positions[1, 2] += .08
        for cls, options, game in ((PeriodicAuPotential, self.options, self.game),
                                  (PeriodicAuPotential, self.options, self.compiled),
                                  (FourierAuPotential, self.fourier_options, self.proxy),
                                  (FourierAuPotential, self.fourier_options, self.compiled_proxy)):
            reference = cls(self.elements, self.positions, accelerator="numpy", **options)
            expected, actual = reference.energy_and_forces_eV(positions), game.energy_and_forces_eV(positions)
            self.assertAlmostEqual(expected[0], actual[0], delta=2e-11)
            np.testing.assert_allclose(expected[1], actual[1], atol=2e-11, rtol=0.)

    def test_full_chain_force_fd_for_both_engines_and_accelerators(self):
        positions = self.positions.copy()
        positions[1, 2] += .08
        for surface in (self.game, self.compiled, self.proxy, self.compiled_proxy):
            self.assert_fd(surface, positions)

    def test_nonbonded_close_fold_does_not_change_normal_graph(self):
        positions = self.positions.copy()
        positions[9] = positions[0]+[.1, .15, .2]  # unrelated H, not another bonded normal neighbour
        with self.assertRaisesRegex(ValueError, "three normal"):
            _normal_geometry(positions, self.game.normal_cutoffs)
        self.assertTrue(np.array_equal(self.game.normal_geometry(positions)[2], self.game._normal_counts))
        for surface in (self.game, self.compiled, self.proxy, self.compiled_proxy):
            self.assert_fd(surface, positions, [(0, 0), (0, 2), (3, 1), (9, 2)])

    def test_topology_remains_fixed_across_old_normal_distance_cutoff(self):
        for separation in (1.49999, 1.5, 1.50001):
            positions = self.positions.copy()
            positions[6] = positions[0]+[separation, 0., .2]
            for surface in (self.compiled, self.compiled_proxy):
                self.assertEqual(surface.normal_geometry(positions)[2][0], 3)
                self.assert_fd(surface, positions, [(0, 0), (0, 2), (6, 0), (6, 2)])

    def test_sp3_and_linear_sp_retain_direction_repulsion_without_truncation(self):
        methane = np.array([[0., 0., 4.], [1., 1., 5.], [1., -1., 3.],
                            [-1., 1., 3.], [-1., -1., 5.]])
        cases = [(["C"]+["H"]*4, methane,
                  [dict(a=1, b=j, order=1.) for j in range(2, 6)], "four-coordinate"),
                 (["H", "C", "C", "H"], np.array([[-1.8, 0., 4.], [-.6, 0., 4.],
                   [.6, 0., 4.], [1.8, 0., 4.]]),
                  [dict(a=1,b=2,order=1.), dict(a=2,b=3,order=3.), dict(a=3,b=4,order=1.)], "linear-sp")]
        for elements, positions, bonds, reason in cases:
            for cls, options in ((PeriodicAuPotential, self.options), (FourierAuPotential, self.fourier_options)):
                surface = cls(elements, positions, bonds=bonds, normal_policy=BONDED_NORMAL_POLICY,
                              accelerator="numba", **options)
                self.assertTrue(any(reason in mode for mode in surface.normal_modes))
                np.testing.assert_array_equal(surface.normal_geometry(positions)[0],
                                              np.tile([0.,0.,1.], (len(elements),1)))
                self.assert_fd(surface, positions)
                # Fixed-axis carbon retains the nonzero SAIP angular term;
                # its fitted coefficient can have either sign.
                if cls is PeriodicAuPotential:
                    from molecular_car.periodic_au import _anisotropic
                    _, delta, radius = surface._pairs(positions, np.array([elements.index("C")]))
                    normals = np.tile([0.,0.,1.], (len(delta),1))
                    angular = _anisotropic(delta,radius,normals,surface.parameters["C","Au"])[0].sum()
                    self.assertGreater(abs(angular), 1e-12)

    def test_initial_collinear_class_is_immutable(self):
        positions = np.array([[0.,0.,4.], [1.,0.,4.], [-1.,0.,4.]])
        surface = PeriodicAuPotential(["C"]*3, positions, use_cache=False, accelerator="numba",
            bonds=[dict(a=1,b=2,order=1.),dict(a=1,b=3,order=1.)], normal_policy=BONDED_NORMAL_POLICY)
        self.assertIn("initial-collinear", surface.normal_modes[0])
        positions[2] += [.2,.8,.3]
        self.assertEqual(surface.normal_geometry(positions)[2][0], 0)
        self.assert_fd(surface, positions)

    def test_dynamic_plane_degeneracy_is_finite_conservative_and_continuous(self):
        reference = np.array([[0.,0.,4.], [1.,0.,4.], [-1.,1.,4.3]])
        bonds = [dict(a=1,b=2,order=1.),dict(a=1,b=3,order=1.)]
        for cls, options in ((PeriodicAuPotential, self.options), (FourierAuPotential, self.fourier_options)):
            surfaces = [cls(["C"]*3, reference, bonds=bonds, normal_policy=BONDED_NORMAL_POLICY,
                accelerator=accelerator, **options) for accelerator in ("numpy", "numba")]
            # q in the Fourier norm-weight blend, the near-zero extension, and
            # exact collinearity. Each uses the same immutable class/graph.
            for area in (0., 2e-7, 1e-6, 3.4e-5, 6e-5, 1.1e-4):
                positions = reference.copy()
                positions[2] = [-1.,area,4.+.3*area]
                outcomes = [surface.energy_and_forces_eV(positions) for surface in surfaces]
                self.assertTrue(np.isfinite(outcomes[0][1]).all())
                self.assertAlmostEqual(outcomes[0][0],outcomes[1][0],delta=2e-11)
                np.testing.assert_allclose(outcomes[0][1],outcomes[1][1],atol=1e-7,rtol=1e-9)
                for surface in surfaces:
                    self.assert_fd(surface,positions,[(2,1),(2,2)],step=1e-10,atol=2e-4,rtol=3e-5)
            center = reference.copy(); center[2] = [-1.,0.,4.]
            minus, plus = center.copy(), center.copy()
            minus[2,1] = -1e-10; plus[2,1] = 1e-10
            em,fm = surfaces[1].energy_and_forces_eV(minus)
            ep,fp = surfaces[1].energy_and_forces_eV(plus)
            self.assertAlmostEqual(em,ep,delta=2e-9)
            self.assertTrue(np.isfinite(fm).all() and np.isfinite(fp).all())

    def test_metadata_and_graph_validation(self):
        metadata = self.compiled_proxy.metadata()
        json.dumps(metadata, allow_nan=False)
        self.assertTrue(metadata["normal_model_is_game_adaptation"])
        self.assertEqual(metadata["normal_policy"],BONDED_NORMAL_POLICY)
        self.assertEqual(metadata["normal_area_regularization_angstrom2"],NORMAL_AREA_EPSILON)
        with self.assertRaisesRegex(ValueError, "explicit bonds"):
            PeriodicAuPotential(self.elements,self.positions,normal_policy=BONDED_NORMAL_POLICY,use_cache=False)
        for bonds in ([dict(a=0,b=1,order=1.)], [dict(a=1,b=2,order=0.)],
                      [dict(a=1,b=2,order=1.),dict(a=2,b=1,order=1.)]):
            with self.assertRaises(ValueError):
                PeriodicAuPotential(self.elements,self.positions,bonds=bonds,
                                    normal_policy=BONDED_NORMAL_POLICY,use_cache=False)

    def test_bundled_layer_reference_is_independent_of_ignored_outputs(self):
        actual, source = reference_layer_heights()
        asset = Path(source)
        self.assertEqual(asset.parent.name,"data")
        data = json.loads(asset.read_text(encoding="utf-8"))
        np.testing.assert_array_equal(actual,data["layer_heights_angstrom"])
        with mock.patch.object(Path,"is_file",side_effect=AssertionError("output presence must not be probed")):
            np.testing.assert_array_equal(reference_layer_heights()[0],actual)
        ideal, source = reference_layer_heights(4.1)
        self.assertIn("ideal",source)
        np.testing.assert_allclose(ideal,np.arange(-5,1)*4.1/np.sqrt(3),rtol=0,atol=0)


if __name__ == "__main__":
    unittest.main()
