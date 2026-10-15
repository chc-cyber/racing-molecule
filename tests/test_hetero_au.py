"""Complete-force checks for the explicit provisional X-Au extension."""
from pathlib import Path
import tempfile
import unittest
import numpy as np
from rdkit import Chem
from molecular_car.molecule import prepare_molecule
from molecular_car.periodic_au import PeriodicAuPotential, BONDED_NORMAL_POLICY, ROOT
from molecular_car.fourier_au import FourierAuPotential


class HeteroAuTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix="hetero-au-", dir=ROOT/"outputs")
        mol = Chem.MolFromSmiles("c1ncccc1")
        cls.p = prepare_molecule(Chem.MolToMolBlock(mol), "pyridine.mol")
        cls.options = dict(bonds=cls.p.bonds, normal_policy=BONDED_NORMAL_POLICY,
                           grid_nx=8, grid_ny=12, grid_dz=.25, cache_dir=cls.tmp.name)
        cls.fourier_options = dict(cls.options, fourier_cache_dir=cls.tmp.name,
             angular_step=.5, height_min=2.8, height_max=4.6, height_step=.15)
        cls.surfaces = [kind(cls.p.atoms, cls.p.positions, accelerator=accelerator, **opts)
                       for kind, opts in ((PeriodicAuPotential, cls.options), (FourierAuPotential, cls.fourier_options))
                       for accelerator in ("numpy", "numba")]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_numpy_numba_and_full_cartesian_gradient(self):
        p = self.p.positions.copy(); p[0, 2] += .07; p[1, 0] += .025
        for a, b in ((self.surfaces[0], self.surfaces[1]), (self.surfaces[2], self.surfaces[3])):
            ea, fa = a.energy_and_forces_eV(p); eb, fb = b.energy_and_forces_eV(p)
            self.assertAlmostEqual(ea, eb, delta=2e-11)
            np.testing.assert_allclose(fa, fb, atol=2e-11, rtol=0.)
            for i in range(len(p)):
                for axis in range(3):
                    plus, minus = p.copy(), p.copy()
                    plus[i, axis] += 1e-5; minus[i, axis] -= 1e-5
                    fd = -(a.energy_and_forces_eV(plus)[0]-a.energy_and_forces_eV(minus)[0])/2e-5
                    self.assertAlmostEqual(fa[i, axis], fd, delta=3e-6)

    def test_hetero_coordinates_contribute_to_carbon_normal_chain(self):
        s = self.surfaces[0]; p = self.p.positions.copy(); p[0, 2] += .15
        n = next(i for i, a in enumerate(self.p.atoms) if a["element"] == "N")
        self.assertTrue(any(n in row for row in s._normal_neighbours))
        _, full = s.energy_and_forces_eV(p)
        _, own = s.hetero_surface.energy_and_forces_eV(p)
        self.assertGreater(np.linalg.norm(full[n]-own[n]), 1e-8)
        self.assertEqual(s._element_codes[n], -1)
        self.assertEqual(s.normal_modes[n], "isotropic-hetero-au-game")

    def test_hetero_cache_height_fallback_is_conservative(self):
        for s in self.surfaces[:2]:
            n = int(s.hetero_surface.indices[0])
            for z in (1.999, 2.001, 2.15, 2.299, 2.301, 11.70, 11.85, 11.999, 12.001):
                p = self.p.positions.copy(); p[n, 2] = z
                _, f = s.hetero_surface.energy_and_forces_eV(p)
                for axis in range(3):
                    plus, minus = p.copy(), p.copy()
                    plus[n, axis] += 1e-5; minus[n, axis] -= 1e-5
                    fd = -(s.hetero_surface.energy_and_forces_eV(plus)[0]-s.hetero_surface.energy_and_forces_eV(minus)[0])/2e-5
                    self.assertAlmostEqual(f[n, axis], fd, delta=3e-6)

    def test_periodic_query_does_not_fold_the_molecular_graph(self):
        s = self.surfaces[0]; p = self.p.positions.copy(); p[:, 0] += 20*s.period[0]
        e0, f0 = s.energy_and_forces_eV(self.p.positions); e1, f1 = s.energy_and_forces_eV(p)
        self.assertAlmostEqual(e0, e1, delta=2e-11)
        np.testing.assert_allclose(f0, f1, atol=2e-11, rtol=0.)

    def test_reference_rejects_hetero_and_game_requires_valid_vdw(self):
        with self.assertRaisesRegex(ValueError, "C/H"):
            PeriodicAuPotential(self.p.atoms, self.p.positions, use_cache=False)
        bad = [dict(a) for a in self.p.atoms]
        bad[1].pop("uff_vdw_rmin_angstrom")
        with self.assertRaisesRegex(ValueError, "vdW"):
            PeriodicAuPotential(bad, self.p.positions, use_cache=False, bonds=self.p.bonds,
                                normal_policy=BONDED_NORMAL_POLICY)

    def test_parameter_units_mixing_and_provisional_metadata(self):
        data = self.surfaces[0].metadata()
        self.assertEqual(data["hetero_au"]["model"], "hetero-au-uff-vdw-game-v1")
        pair = data["hetero_au"]["used_pairs"][0]
        self.assertAlmostEqual(pair["rmin_angstrom"]**2, pair["self_rmin_angstrom"]*3.293)
        self.assertAlmostEqual(pair["well_depth_kcal_mol"]**2, pair["self_well_depth_kcal_mol"]*.039)
        self.assertIn("not fitted", data["hetero_au"]["status"])
        self.assertEqual(len(data["hetero_au"]["parameter_sha256"]), 64)
        self.assertEqual(len(data["hetero_au"]["used_atom_parameters_sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
