"""Independent v2 reference comparisons and conservative fused-phase checks."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from molecular_car.fourier_atlas import FourierAtlasPotential
from molecular_car.lattice_corrugation import LatticeCorrugation
from molecular_car.molecule import prepare_molecule
from molecular_car.periodic_au import ROOT, REFERENCE_NORMAL_POLICY


class AtlasKernelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cache = ROOT / 'outputs/fourier-au-check/unit-cache'
        cache.mkdir(parents=True, exist_ok=True)
        cls.temp = tempfile.TemporaryDirectory(prefix='atlas-', dir=cache)
        cls.molecule = prepare_molecule((ROOT/'examples/benzene.mol').read_text(), 'benzene.mol')
        cls.options = dict(angular_step=.6, height_min=2.4, height_max=5.4,
                           height_step=.15, grid_nx=8, grid_ny=12, grid_dz=.25,
                           cache_dir=cls.temp.name, fourier_cache_dir=cls.temp.name,
                           bonds=cls.molecule.bonds, normal_policy='bonded-saip-game-v1')
        cls.reference = FourierAtlasPotential(cls.molecule.atoms, cls.molecule.positions,
                                               accelerator='numpy', **cls.options)
        cls.compiled = FourierAtlasPotential(cls.molecule.atoms, cls.molecule.positions,
                                              accelerator='numba', **cls.options)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def compare(self, x, scale=1.5, reference=None, compiled=None, force_atol=2e-9):
        ref = LatticeCorrugation(reference or self.reference, scale)
        fast = LatticeCorrugation(compiled or self.compiled, scale)
        e, f = ref.energy_and_forces_eV(x)
        ec, fc = fast.energy_and_forces_eV(x)
        np.testing.assert_allclose(ec, e, atol=2e-10, rtol=1e-11)
        np.testing.assert_allclose(fc, f, atol=force_atol, rtol=1e-11)

    def rotated(self, seed, height=3.6):
        rotation, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
        x = self.molecule.positions.copy()
        x = (x-x.mean(axis=0))@rotation
        x[:, 2] += height-x[:, 2].min()
        return x

    def test_full_rotations_and_phase_fusion(self):
        for seed in range(8):
            for scale in (0., 1., 1.5, 2.):
                self.compare(self.rotated(seed), scale)

    def test_height_core_and_fallback_switches(self):
        for height in (1.99, 2., 2.1, 2.39, 2.4, 2.55, 2.6, 2.7, 3., 5.1, 5.25, 5.4, 5.41, 12.):
            x = self.molecule.positions.copy()
            x[:, 2] += height-x[:, 2].min()
            self.compare(x)

    def test_complete_force_is_energy_derivative(self):
        surface = LatticeCorrugation(self.compiled, 1.5)
        for x in (self.rotated(5, 2.5), self.rotated(9, 3.3)):
            _, f = surface.energy_and_forces_eV(x)
            for atom in range(len(x)):
                for axis in range(3):
                    a, b = x.copy(), x.copy()
                    a[atom, axis] += 1e-5; b[atom, axis] -= 1e-5
                    fd = (surface.energy_and_forces_eV(a)[0]-surface.energy_and_forces_eV(b)[0])/2e-5
                    self.assertAlmostEqual(f[atom, axis], -fd, delta=2e-6)

    def test_multicell_unwrapped_coordinates(self):
        x = self.rotated(7)
        shift = np.array([997*self.compiled.period[0], -613*self.compiled.period[1], 0.])
        self.compare(x+shift)
        s = LatticeCorrugation(self.compiled, 1.5)
        e, f = s.energy_and_forces_eV(x)
        es, fs = s.energy_and_forces_eV(x+shift)
        np.testing.assert_allclose(es, e, atol=1e-9, rtol=0)
        np.testing.assert_allclose(fs, f, atol=1e-8, rtol=0)

    def test_weighted_centers_and_regularized_normals(self):
        # A bonded C center with a nearly collinear plane exercises the
        # smooth normal-magnitude fallback and its neighbour chain forces.
        x = np.array([[0., 0., 3.4], [1.4, 0., 3.5], [2.8, .8, 3.6]])
        options = dict(self.options, bonds=[{'a': 1, 'b': 2, 'order': 1},
                                           {'a': 1, 'b': 3, 'order': 1}])
        ref = FourierAtlasPotential(['C']*3, x, accelerator='numpy', **options)
        fast = FourierAtlasPotential(['C']*3, x, accelerator='numba', **options)
        self.assertEqual(fast.normal_geometry(x)[2][0], 2)
        for y in (0., 1e-6, 2.4e-5, 4.3e-5, 8e-5):
            x[2, 1] = y
            # Dividing by a ~1e-5 A² plane area amplifies summation roundoff.
            self.compare(x, reference=ref, compiled=fast, force_atol=2e-8)
        x = self.rotated(2)
        w = np.linspace(.2, 1.4, len(x)); dz = np.linspace(-.1, .1, len(x))
        a = self.reference.energy_and_forces_eV(x, w, dz)
        b = self.compiled.energy_and_forces_eV(x, w, dz)
        np.testing.assert_allclose(a[0], b[0], atol=2e-10, rtol=1e-11)
        np.testing.assert_allclose(a[1], b[1], atol=2e-9, rtol=1e-11)

    def test_hydrogen_with_direction_dependent_normal(self):
        x = np.array([[0., 0., 3.4], [1.1, .1, 3.5], [.2, 1.1, 3.6]])
        options = dict(self.options, normal_policy=REFERENCE_NORMAL_POLICY)
        ref = FourierAtlasPotential(['H', 'C', 'C'], x, accelerator='numpy', **options)
        fast = FourierAtlasPotential(['H', 'C', 'C'], x, accelerator='numba', **options)
        self.compare(x, reference=ref, compiled=fast)

    def test_hetero_atoms_keep_same_weighted_potential(self):
        from rdkit import Chem
        content = Chem.MolToMolBlock(Chem.MolFromSmiles('c1ccncc1'))
        m = prepare_molecule(content, 'pyridine.mol')
        options = dict(self.options, bonds=m.bonds)
        ref = FourierAtlasPotential(m.atoms, m.positions, accelerator='numpy', **options)
        fast = FourierAtlasPotential(m.atoms, m.positions, accelerator='numba', **options)
        for height in (2.5, 3.4, 5.3):
            x = m.positions.copy(); x[:, 2] += height-x[:, 2].min()
            self.compare(x, reference=ref, compiled=fast)


if __name__ == '__main__':
    unittest.main()
