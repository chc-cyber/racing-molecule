import unittest
from rdkit import Chem
from molecular_car.editor.superatoms import read_groups, write_groups


class SuperatomTests(unittest.TestCase):
    def test_real_atoms_bonds_attachment_and_display_roundtrip(self):
        for expanded in (False, True):
            mol = Chem.MolFromSmiles("CCC")
            ids = ["stable-a", "stable-b", "stable-c"]
            project = {"atoms": [{"id": i} for i in ids], "metadata": {"ketcherGroups": [{
                "type": "SUP", "name": "Et", "atomIds": ids[:2], "expanded": expanded,
                "attachmentPoints": [{"atomId": ids[1], "attachmentId": "1"}],
            }]}}
            write_groups(mol, project)
            # V3000 avoids the pinned RDKit V2000 short SAP-ID parser bug.
            block = Chem.MolToMolBlock(mol, forceV3000=True)
            parsed = Chem.MolFromMolBlock(block, removeHs=False)
            self.assertIsNotNone(parsed)
            self.assertEqual(parsed.GetNumAtoms(), 3)
            self.assertEqual(parsed.GetNumBonds(), 2)
            groups = read_groups(parsed, ids)
            self.assertEqual(groups[0]["atomIds"], ids[:2])
            self.assertEqual(groups[0]["name"], "Et")
            self.assertEqual(groups[0]["expanded"], expanded)
            self.assertEqual(groups[0]["attachmentPoints"], [{"atomId": ids[1], "attachmentId": "1"}])
            self.assertEqual(list(Chem.GetMolSubstanceGroups(parsed)[0].GetBonds()), [1])

    def test_other_group_semantics_are_rejected(self):
        mol = Chem.MolFromSmiles("CC")
        Chem.CreateMolSubstanceGroup(mol, "SRU").AddAtomWithIdx(0)
        with self.assertRaisesRegex(ValueError, "其他 S-group"):
            read_groups(mol, ["a", "b"])

    def test_mapping_cannot_invent_or_drop_a_member(self):
        mol = Chem.MolFromSmiles("CC")
        project = {"atoms": [{"id": "a"}, {"id": "b"}], "metadata": {"ketcherGroups": [{"type": "SUP", "atomIds": ["missing"]}]}}
        with self.assertRaisesRegex(ValueError, "编号映射"):
            write_groups(mol, project)


if __name__ == "__main__":
    unittest.main()
