"""Regression: private 2D depiction recognizes sp geometry without rejecting drafts."""
import copy
import math
import unittest

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor

from molecular_car.editor.service import geometry, import_file, to_mol


def drawing(smiles):
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    project = import_file({"filename": "alkyne.mol", "content": Chem.MolToMolBlock(mol)})["project"]
    project["has3d"] = True
    # A distinct manually designed conformer must survive a 2D layout exactly.
    for i, atom in enumerate(project["atoms"]):
        atom["pos2"] = [math.sin(i * .8) * 2, math.cos(i * .9) * 2]
        atom["pos3"] = [i * .71, math.sin(i), .123456789 + i * .34]
    return project


def sp_angles(project):
    mol = to_mol(project, 2)
    xyz = mol.GetConformer().GetPositions()
    values = []
    for atom in mol.GetAtoms():
        if atom.GetDegree() != 2 or not any(b.GetBondType() == Chem.BondType.TRIPLE for b in atom.GetBonds()):
            continue
        first, second = [n.GetIdx() for n in atom.GetNeighbors()]
        u, v = xyz[first] - xyz[atom.GetIdx()], xyz[second] - xyz[atom.GetIdx()]
        values.append(math.degrees(math.acos(np.clip(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)), -1, 1))))
    return values


class LinearLayoutTests(unittest.TestCase):
    def test_alkynes_and_diynes_are_straight_and_do_not_replace_3d(self):
        for smiles in ("CC#CC", "CC#CC#CC", "c1ccccc1C#Cc1ccccc1", "[H]C#CC"):
            with self.subTest(smiles=smiles):
                before = drawing(smiles)
                original = copy.deepcopy(before)
                result = geometry({"project": before, "operation": "layout2d", "revision": 28})
                after = result["project"]
                angles = sp_angles(after)
                self.assertTrue(angles)
                for angle in angles:
                    self.assertAlmostEqual(angle, 180, places=5)
                self.assertEqual(result["revision"], 28)
                self.assertEqual(before, original)
                self.assertEqual(after["bonds"], before["bonds"])
                for old, new in zip(before["atoms"], after["atoms"]):
                    self.assertEqual({k: v for k, v in old.items() if k != "pos2"}, {k: v for k, v in new.items() if k != "pos2"})

    def test_invalid_valence_draft_can_still_be_laid_out(self):
        before = drawing("CC#CC")
        anchor = before["atoms"][0]["id"]
        for i in range(4):
            atom_id = "draft-f-" + str(i)
            before["atoms"].append({"id": atom_id, "element": "F", "charge": 0, "pos2": [0, i], "pos3": [0, 0, i]})
            before["bonds"].append({"id": "draft-edge-" + str(i), "a": anchor, "b": atom_id, "order": 1})
        after = geometry({"project": before, "operation": "layout2d"})["project"]
        self.assertEqual([a["id"] for a in after["atoms"]], [a["id"] for a in before["atoms"]])
        self.assertEqual([a["pos3"] for a in after["atoms"]], [a["pos3"] for a in before["atoms"]])
        self.assertEqual([(b["id"], b["order"]) for b in after["bonds"]], [(b["id"], b["order"]) for b in before["bonds"]])
        for angle in sp_angles(after):
            self.assertAlmostEqual(angle, 180, places=5)


if __name__ == "__main__":
    unittest.main()
