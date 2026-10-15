"""Independent geometry actions must not secretly perform another action."""
import copy
import pathlib
import unittest
from unittest.mock import patch

from molecular_car.editor.service import geometry, import_file, EditorError
from molecular_car.editor.fragments import benzene

ROOT = pathlib.Path(__file__).resolve().parents[1]


def example(name):
    return import_file({"filename": name + ".mol", "content": (ROOT / "frontend/editor/examples" / (name + ".mol")).read_text()})["project"]


def bare_benzene():
    project = benzene()
    project["has3d"] = False
    return project


class SeparateOperationsTests(unittest.TestCase):
    def test_rebuild_is_coordinates_only_and_no_forcefield(self):
        before = bare_benzene()
        with patch("molecular_car.editor.service.AllChem.UFFHasAllMoleculeParams", side_effect=AssertionError("no UFF")), patch("molecular_car.editor.service.AllChem.UFFGetMoleculeForceField", side_effect=AssertionError("no UFF")):
            result = geometry({"project": before, "operation": "rebuild", "revision": 123})
        after = result["project"]
        self.assertTrue(after["has3d"])
        self.assertFalse(result["construction"]["optimized"])
        self.assertNotIn("energy", result)
        self.assertEqual(result["hydrogensAdded"], 0)
        self.assertEqual(result["revision"], 123)
        self.assertEqual(before["bonds"], after["bonds"])
        for old, new in zip(before["atoms"], after["atoms"]):
            self.assertEqual({k: v for k, v in old.items() if k != "pos3"}, {k: v for k, v in new.items() if k != "pos3"})
        self.assertEqual(len(after["atoms"]), len(before["atoms"]))

    def test_two_dimensional_hydrogenation_keeps_layout_and_old_coordinates(self):
        before = bare_benzene()
        self.assertFalse(before["has3d"])
        with patch("molecular_car.editor.service.AllChem.EmbedMolecule", side_effect=AssertionError("no embedding")), patch("molecular_car.editor.service.AllChem.UFFGetMoleculeForceField", side_effect=AssertionError("no UFF")):
            result = geometry({"project": before, "operation": "hydrogenate"})
        after = result["project"]
        self.assertFalse(after["has3d"])
        self.assertEqual(after["atoms"][:len(before["atoms"])], before["atoms"])
        self.assertEqual(after["bonds"][:len(before["bonds"])], before["bonds"])
        self.assertEqual(result["hydrogensAdded"], 6)
        self.assertEqual(len(set(a["id"] for a in after["atoms"])), 12)
        self.assertNotIn("energy", result)
        repeated = geometry({"project": after, "operation": "hydrogenate"})
        self.assertEqual(repeated["project"], after)

    def test_three_dimensional_hydrogenation_preserves_hand_coordinates(self):
        before = geometry({"project": bare_benzene(), "operation": "rebuild"})["project"]
        before["atoms"][0]["pos3"][2] += .123456789
        result = geometry({"project": before, "operation": "hydrogenate"})
        self.assertEqual(result["project"]["atoms"][:len(before["atoms"])], before["atoms"])
        self.assertEqual(result["hydrogensAdded"], 6)

    def test_uff_only_optimizes_and_honors_fixed_atom(self):
        before = geometry({"project": bare_benzene(), "operation": "rebuild"})["project"]
        atom = before["atoms"][0]
        before["atoms"][1]["pos3"][2] += .5
        with patch("molecular_car.editor.service.AllChem.EmbedMolecule", side_effect=AssertionError("no regeneration")):
            result = geometry({"project": before, "operation": "relax", "fixedIds": [atom["id"]]})
        after = result["project"]
        self.assertEqual(len(before["atoms"]), len(after["atoms"]))
        self.assertEqual(before["bonds"], after["bonds"])
        self.assertEqual(after["atoms"][0]["pos3"], atom["pos3"])
        self.assertEqual(result["hydrogensAdded"], 0)
        self.assertLessEqual(result["energy"]["after"], result["energy"]["before"])

    def test_layout_does_not_replace_three_dimensional_design(self):
        before = example("tetraphenylbenzene")
        before["atoms"][0]["pos3"][2] += .123456789
        result = geometry({"project": before, "operation": "layout2d"})
        self.assertEqual([a["id"] for a in result["project"]["atoms"]], [a["id"] for a in before["atoms"]])
        self.assertEqual([a["pos3"] for a in result["project"]["atoms"]], [a["pos3"] for a in before["atoms"]])

    def test_rebuild_failure_leaves_input_intact(self):
        before = bare_benzene()
        original = copy.deepcopy(before)
        with patch("molecular_car.editor.service.embed_initial", side_effect=RuntimeError("rejected")):
            with self.assertRaises(RuntimeError):
                geometry({"project": before, "operation": "rebuild"})
        self.assertEqual(before, original)


if __name__ == "__main__":
    unittest.main()
