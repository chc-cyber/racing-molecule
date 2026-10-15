"""Chemical environments are admitted by actual force-field consistency."""
import unittest
from unittest.mock import Mock
import numpy as np
from rdkit import Chem
from molecular_car.molecule import prepare_molecule, MoleculePreparationError, _check_uff_cartesian_gradient


class ElementAdmissionTests(unittest.TestCase):
    def test_three_coordinate_phosphorus_is_rejected_before_md(self):
        mol = Chem.MolFromSmiles("CP(C)C")
        with self.assertRaisesRegex(MoleculePreparationError, "三配位磷"):
            prepare_molecule(Chem.MolToMolBlock(mol), "phosphine.mol")

    def test_four_coordinate_phosphorus_and_hetero_atom_type_vdw_are_available(self):
        mol = Chem.MolFromSmiles("CP(=O)(C)C")
        p = prepare_molecule(Chem.MolToMolBlock(mol), "phosphine-oxide.mol")
        self.assertTrue(p.forcefield_validation["passed"])
        self.assertEqual(p.forcefield_validation["coordinates_checked"], 3*len(p.atoms))
        self.assertTrue(all(a["uff_vdw_rmin_angstrom"] > 0 and a["uff_vdw_well_depth_kcal_mol"] > 0 for a in p.atoms))

    def test_finite_but_wrong_gradient_is_not_admitted(self):
        forcefield = Mock()
        forcefield.CalcEnergy.side_effect = lambda p: float(np.dot(p, p))
        forcefield.CalcGrad.side_effect = lambda p: tuple(3*np.array(p))
        with self.assertRaisesRegex(MoleculePreparationError, "能量与受力不一致"):
            _check_uff_cartesian_gradient(forcefield, np.array([[.2, .3, .4]]))

    def test_unsupported_metal_has_an_explicit_element_boundary(self):
        mol = Chem.MolFromSmiles("C[Si](C)(C)C")
        p = prepare_molecule(Chem.MolToMolBlock(mol), "silane.mol")
        self.assertTrue(p.forcefield_validation["passed"])
        mol = Chem.MolFromSmiles("[Fe]")
        with self.assertRaisesRegex(MoleculePreparationError, "Fe"):
            prepare_molecule(Chem.MolToMolBlock(mol), "iron.mol")


if __name__ == "__main__":
    unittest.main()
