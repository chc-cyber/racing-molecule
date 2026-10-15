"""Race preparation deliberately removes isotope-specific masses on its copy."""
from pathlib import Path
import unittest

import numpy as np
from rdkit import Chem

from molecular_car.molecule import prepare_molecule

ROOT = Path(__file__).resolve().parents[1]


class RaceMassTests(unittest.TestCase):
    def test_mixed_carbon_and_hydrogen_isotopes_use_same_standard_race(self):
        content = (ROOT / "examples/benzene.mol").read_text(encoding="utf-8")
        source = Chem.MolFromMolBlock(content, removeHs=False)
        source.GetAtomWithIdx(0).SetIsotope(13)
        source.GetAtomWithIdx(1).SetIsotope(14)
        source.GetAtomWithIdx(6).SetIsotope(2)
        labelled = Chem.MolToMolBlock(source)
        prepared = prepare_molecule(labelled, "labelled.mol")
        standard = prepare_molecule(Chem.MolToMolBlock(
            Chem.MolFromMolBlock(content, removeHs=False)), "standard.mol")
        np.testing.assert_array_equal([a["mass"] for a in prepared.atoms],
                                      [a["mass"] for a in standard.atoms])
        np.testing.assert_allclose(prepared.positions, standard.positions, atol=1e-12, rtol=0)
        self.assertEqual(prepared.bonds, standard.bonds)
        self.assertTrue(any("3 个同位素标签" in note for note in prepared.notes))
        self.assertTrue(all(a.GetIsotope() == 0 for a in prepared._rdkit_molecule.GetAtoms()))
        # The editable input and its labels were not overwritten.
        restored = Chem.MolFromMolBlock(labelled, removeHs=False)
        self.assertEqual([restored.GetAtomWithIdx(i).GetIsotope() for i in (0, 1, 6)], [13, 14, 2])

    def test_unlabelled_mass_and_input_geometry_remain_unchanged(self):
        content = (ROOT / "examples/benzene.mol").read_text(encoding="utf-8")
        prepared = prepare_molecule(content, "benzene.mol")
        for atom in prepared.atoms:
            expected = Chem.GetPeriodicTable().GetAtomicWeight(atom["element"])
            self.assertEqual(atom["mass"], expected)
        self.assertFalse(any("同位素标签" in note for note in prepared.notes))


if __name__ == "__main__":
    unittest.main()
