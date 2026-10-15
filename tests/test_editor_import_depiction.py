"""Import depiction regressions: chemical 2D wedges never replace input 3D."""
import copy
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor

from molecular_car.editor import dispatch
from molecular_car.editor.service import _stereo_signature

ROOT = Path(__file__).resolve().parents[1]


class EditorImportDepictionTests(unittest.TestCase):
    def test_3d_chiral_import_adds_indexed_2d_wedges_without_moving_input(self):
        content = (ROOT / "frontend/editor/examples/chiral.mol").read_text()
        source = Chem.MolFromMolBlock(content, removeHs=False)
        expected = _stereo_signature(source)
        original_positions = np.array(source.GetConformer().GetPositions(), copy=True)
        result = dispatch("/api/editor/import", {"filename": "chiral.mol", "content": content})
        project = result["project"]
        self.assertTrue(project["has3d"])
        np.testing.assert_array_equal([a["pos3"] for a in project["atoms"]], original_positions)
        self.assertEqual(result["mapping"]["indexToAtomId"], [a["id"] for a in project["atoms"]])
        self.assertEqual(project["metadata"]["indexToAtomId"], result["mapping"]["indexToAtomId"])
        self.assertEqual([a["chiralTag"] for a in project["atoms"]], [str(a.GetChiralTag()) for a in source.GetAtoms()])
        wedges = [b for b in project["bonds"] if b["direction"] in ("BEGINWEDGE", "BEGINDASH")]
        self.assertTrue(wedges, "A tag-only 3D input must get a visible 2D stereo declaration")
        by_id = {a["id"]: a for a in project["atoms"]}
        self.assertTrue(all(by_id[b["a"]]["chiralTag"] != "CHI_UNSPECIFIED" for b in wedges))
        atom_indices = {a["id"]: i for i, a in enumerate(project["atoms"])}
        self.assertEqual(
            [(frozenset((atom_indices[b["a"]], atom_indices[b["b"]])), b["order"]) for b in project["bonds"]],
            [(frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())), b.GetBondTypeAsDouble()) for b in source.GetBonds()],
        )
        before_export = copy.deepcopy(project)
        for dimension in (2, 3):
            with self.subTest(dimension=dimension):
                exported = dispatch("/api/editor/export", {"project": project, "dimension": dimension})
                observed = Chem.MolFromMolBlock(exported["content"], removeHs=False)
                self.assertEqual(_stereo_signature(observed), expected)
                self.assertEqual(exported["mapping"]["indexToAtomId"], result["mapping"]["indexToAtomId"])
                if dimension == 3:
                    np.testing.assert_allclose(observed.GetConformer().GetPositions(), original_positions, atol=0.000051)
        self.assertEqual(project, before_export)

    def test_2d_import_keeps_the_supplied_layout_including_translation_and_orientation(self):
        molecule = Chem.MolFromSmiles("F[C@](Cl)(Br)I")
        rdDepictor.Compute2DCoords(molecule)
        conformer = molecule.GetConformer()
        for index in range(molecule.GetNumAtoms()):
            point = conformer.GetAtomPosition(index)
            conformer.SetAtomPosition(index, (11 + 0.8 * point.x - 0.6 * point.y,
                                             -7 + 0.6 * point.x + 0.8 * point.y, 0))
        content = Chem.MolToMolBlock(molecule)
        source = Chem.MolFromMolBlock(content, removeHs=False)
        project = dispatch("/api/editor/import", {"filename": "2d-stereo.mol", "content": content})["project"]
        self.assertFalse(project["has3d"])
        np.testing.assert_array_equal([a["pos2"] for a in project["atoms"]], source.GetConformer().GetPositions()[:, :2])
        np.testing.assert_array_equal([a["pos3"] for a in project["atoms"]], source.GetConformer().GetPositions())
        exported = dispatch("/api/editor/export", {"project": project, "dimension": 2})
        self.assertEqual(_stereo_signature(Chem.MolFromMolBlock(exported["content"], removeHs=False)), _stereo_signature(source))

    def test_depiction_failure_keeps_the_import_editable_and_reports_the_limit(self):
        content = (ROOT / "frontend/editor/examples/chiral.mol").read_text()
        source = Chem.MolFromMolBlock(content, removeHs=False)
        with patch("molecular_car.editor.service._apply_stereo_layout", side_effect=ValueError("test depiction unavailable")):
            result = dispatch("/api/editor/import", {"filename": "chiral.mol", "content": content})
        np.testing.assert_array_equal([a["pos3"] for a in result["project"]["atoms"]], source.GetConformer().GetPositions())
        self.assertTrue(any("未能生成完整的二维立体描画" in warning for warning in result["warnings"]))
        self.assertIn("stereoNotice", result["project"]["metadata"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
