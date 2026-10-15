"""Behavioral coverage for fast, non-optimizing automatic 2D -> 3D."""
import copy
from contextlib import ExitStack
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from molecular_car.editor import dispatch, EditorError
from molecular_car.editor import construction
from molecular_car.editor.service import FORMAT, to_mol, _stereo_signature
from tests.test_editor_service import from_smiles_fixture, mark_c60_seed


ROOT = Path(__file__).resolve().parents[1]


def native_cage(three_d=False):
    source = json.loads((ROOT / "frontend/editor/ketcher-source/c60-source.json").read_text(encoding="utf-8"))
    project = {"format": FORMAT, "name": "C60", "atoms": source["atoms"], "bonds": source["bonds"],
               "has3d": three_d, "fragments": [], "metadata": {}}
    ids = [a["id"] for a in project["atoms"]]
    mark_c60_seed(project, ids)
    project["metadata"]["geometryPending"] = {"atomIds": [] if three_d else list(ids)}
    return project


def attach_methyl(project):
    project = copy.deepcopy(project)
    anchor = project["atoms"][0]
    changed = next(b for b in project["bonds"] if b["order"] == 2 and anchor["id"] in (b["a"], b["b"]))
    changed["order"] = 1
    project["atoms"].append({"id": "new-methyl", "element": "C", "charge": 0,
                              "chiralTag": "CHI_UNSPECIFIED", "pos2": [anchor["pos2"][0] + 1, anchor["pos2"][1]],
                              "pos3": [anchor["pos3"][0] + 1, anchor["pos3"][1], anchor["pos3"][2]]})
    project["bonds"].append({"id": "new-methyl-attachment", "a": anchor["id"], "b": "new-methyl", "order": 1})
    pending = project["metadata"]["geometryPending"]
    pending["atomIds"].append("new-methyl")
    pending["changedAtomIds"] = [changed["a"], changed["b"]]
    pending["editedAtomIds"] = [changed["a"], changed["b"]]
    return project


def attach_new_fragment(project, smiles):
    """Explicit addition at the native cage, with exact fragment index mapping."""
    project = copy.deepcopy(project)
    molecule = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(molecule)
    fragment = dispatch("/api/editor/import", {"filename": "new.mol", "content": Chem.MolToMolBlock(molecule)})["project"]
    anchor = project["atoms"][0]
    changed = next(b for b in project["bonds"] if b["order"] == 2 and anchor["id"] in (b["a"], b["b"]))
    changed["order"] = 1
    for atom in fragment["atoms"]:
        atom["pos3"] = (np.array(atom["pos3"]) + anchor["pos3"] + np.array([1, 0, 0])).tolist()
    project["atoms"].extend(fragment["atoms"])
    project["bonds"].extend(fragment["bonds"])
    project["bonds"].append({"id": "new-fragment-attachment", "a": anchor["id"], "b": fragment["atoms"][0]["id"], "order": 1})
    project["metadata"]["geometryPending"] = {"atomIds": [a["id"] for a in fragment["atoms"]],
                                                   "changedAtomIds": [changed["a"], changed["b"]]}
    return project, fragment


def no_force_field():
    stack = ExitStack()
    for name in ("UFFHasAllMoleculeParams", "UFFGetMoleculeForceField", "UFFOptimizeMolecule",
                 "UFFOptimizeMoleculeConfs", "MMFFHasAllMoleculeParams", "MMFFGetMoleculeForceField"):
        stack.enter_context(patch.object(AllChem, name, side_effect=AssertionError(f"automatic prepare called {name}")))
    return stack


class InitialConstructionTests(unittest.TestCase):
    def assert_unchanged_cage(self, before, after):
        by_id = {a["id"]: a for a in after["atoms"]}
        for atom in before["atoms"][:60]:
            self.assertEqual(by_id[atom["id"]], atom)

    def test_native_c60_methyl_first_conversion_needs_neither_embedding_nor_force_field(self):
        project = attach_methyl(native_cage(False))
        original = copy.deepcopy(project)
        with no_force_field(), patch.object(AllChem, "EmbedMolecule", side_effect=AssertionError("whole cage re-embedded")):
            result = dispatch("/api/editor/geometry", {"project": project, "operation": "prepare", "revision": 101})
        after = result["project"]
        self.assertEqual(len(after["atoms"]), 65)
        self.assertEqual(result["hydrogensAdded"], 4)
        self.assertEqual(result["revision"], 101)
        self.assert_unchanged_cage(project, after)
        self.assertNotIn("energy", result)
        self.assertFalse(result["construction"]["optimized"])
        self.assertEqual(after["metadata"]["ketcherGroups"], project["metadata"]["ketcherGroups"])
        self.assertEqual(result["mapping"]["indexToAtomId"], [a["id"] for a in after["atoms"]])
        self.assertEqual(project, original)
        self.assert_valid_methyl_geometry(project, after)

    def assert_valid_methyl_geometry(self, before, after):
        by_id = {a["id"]: a for a in after["atoms"]}
        origin = np.array(by_id[before["atoms"][0]["id"]]["pos3"])
        methyl = np.array(by_id["new-methyl"]["pos3"])
        center = np.array([by_id[a["id"]]["pos3"] for a in before["atoms"][:60]]).mean(axis=0)
        self.assertAlmostEqual(np.linalg.norm(methyl - origin), 1.52, places=10)
        self.assertGreater(np.dot(methyl - origin, origin - center), 5)
        methyl_h = [by_id[b["b"] if b["a"] == "new-methyl" else b["a"]] for b in after["bonds"]
                    if "new-methyl" in (b["a"], b["b"])
                    and by_id[b["b"] if b["a"] == "new-methyl" else b["a"]]["element"] == "H"]
        self.assertEqual(len(methyl_h), 3)
        for hydrogen in methyl_h:
            vector = np.array(hydrogen["pos3"]) - methyl
            self.assertAlmostEqual(np.linalg.norm(vector), 1.1, places=10)
            cosine = np.dot(vector, origin - methyl) / np.linalg.norm(vector) / np.linalg.norm(origin - methyl)
            self.assertAlmostEqual(cosine, -1 / 3, places=4)

    def test_existing_c60_methyl_addition_preserves_even_changed_cage_atoms_and_is_idempotent(self):
        project = attach_methyl(native_cage(True))
        original = copy.deepcopy(project)
        with no_force_field(), patch.object(AllChem, "EmbedMolecule", side_effect=AssertionError("whole cage re-embedded")):
            first = dispatch("/api/editor/geometry", {"project": project, "operation": "prepare"})
            second = dispatch("/api/editor/geometry", {"project": first["project"], "operation": "prepare"})
        self.assert_unchanged_cage(project, first["project"])
        self.assertEqual(first["project"], second["project"])
        self.assertEqual(second["hydrogensAdded"], 0)
        self.assert_valid_methyl_geometry(project, first["project"])
        self.assertEqual(project, original)

    def test_first_ordinary_structure_embeds_but_never_uses_uff_parameters(self):
        molecule = Chem.MolFromSmiles("C[C@H](F)/C=C/Cl")
        AllChem.Compute2DCoords(molecule)
        project = dispatch("/api/editor/import", {"filename": "drawn.mol", "content": Chem.MolToMolBlock(molecule)})["project"]
        expected = _stereo_signature(to_mol(project, 2, sanitize=True))
        with no_force_field():
            result = dispatch("/api/editor/geometry", {"project": project, "operation": "prepare"})
        self.assertEqual(_stereo_signature(to_mol(result["project"], 3, sanitize=True)), expected)
        self.assertNotIn("energy", result)
        self.assertEqual([a["pos2"] for a in result["project"]["atoms"][:len(project["atoms"])]], [a["pos2"] for a in project["atoms"]])
        dispatch("/api/editor/export", {"project": result["project"], "dimension": 3})

    def test_elements_without_real_uff_parameters_can_still_get_initial_geometry(self):
        for smiles in ("[Fe+2]", "[He]"):
            with self.subTest(smiles=smiles):
                molecule = Chem.MolFromSmiles(smiles)
                with rdBase.BlockLogs():
                    self.assertFalse(AllChem.UFFHasAllMoleculeParams(molecule))
                AllChem.Compute2DCoords(molecule)
                project = dispatch("/api/editor/import", {"filename": "ion.mol", "content": Chem.MolToMolBlock(molecule)})["project"]
                with no_force_field():
                    result = dispatch("/api/editor/geometry", {"project": project, "operation": "prepare"})
                self.assertTrue(result["project"]["has3d"])
                self.assertEqual(result["project"]["atoms"][0]["id"], project["atoms"][0]["id"])
                self.assertEqual(result["project"]["atoms"][0]["charge"], project["atoms"][0]["charge"])
                self.assertNotIn("energy", result)

    def test_local_phenyl_and_chain_embed_only_new_fragment_and_one_cap(self):
        for smiles in ("c1ccccc1", "CCC"):
            with self.subTest(smiles=smiles):
                project, fragment = attach_new_fragment(native_cage(True), smiles)
                calls = []
                original_embed = AllChem.EmbedMolecule
                def local_only(mol, *args, **kwargs):
                    heavy = sum(a.GetAtomicNum() != 1 for a in mol.GetAtoms())
                    self.assertEqual(heavy, len(fragment["atoms"]) + 1)
                    calls.append(heavy)
                    return original_embed(mol, *args, **kwargs)
                with no_force_field(), patch.object(AllChem, "EmbedMolecule", side_effect=local_only):
                    result = dispatch("/api/editor/geometry", {"project": project, "operation": "prepare"})
                self.assertTrue(calls)
                self.assert_unchanged_cage(project, result["project"])
                by_id = {a["id"]: a for a in result["project"]["atoms"]}
                for bond in fragment["bonds"]:
                    length = np.linalg.norm(np.array(by_id[bond["a"]]["pos3"]) - by_id[bond["b"]]["pos3"])
                    self.assertTrue(1.25 < length < 1.7, length)
                self.assertAlmostEqual(np.linalg.norm(np.array(by_id[project["atoms"][0]["id"]]["pos3"])
                                                     - by_id[fragment["atoms"][0]["id"]]["pos3"]), 1.52, places=10)

    def test_new_tetrahedral_fragment_uses_proper_alignment_and_keeps_declared_stereo(self):
        for smiles in ("C[C@H](O)C(=O)O", "C[C@@H](O)C(=O)O"):
            with self.subTest(smiles=smiles):
                project, fragment = attach_new_fragment(native_cage(True), smiles)
                expected = _stereo_signature(to_mol(project, 3, sanitize=True))
                old_positions = [a["pos3"] for a in project["atoms"][:60]]
                with no_force_field():
                    after = dispatch("/api/editor/geometry", {"project": project, "operation": "prepare"})["project"]
                self.assertEqual(_stereo_signature(to_mol(after, 3, sanitize=True)), expected)
                self.assertEqual([a["pos3"] for a in after["atoms"][:60]], old_positions)
                for dimension in (2, 3):
                    output = dispatch("/api/editor/export", {"project": after, "dimension": dimension})
                    restored = dispatch("/api/editor/import", output)["project"]
                    # MOL import allocates new IDs, but atom row indices are exact.
                    observed = _stereo_signature(to_mol(restored, dimension, sanitize=True))
                    if dimension == 2:
                        self.assertEqual(observed, expected)
                    else:
                        # External readers may additionally infer previously
                        # unspecified cage stereo from 3D coordinates.
                        for index, signature in expected[0].items():
                            self.assertEqual(observed[0].get(index), signature)
                        self.assertEqual(observed[1], expected[1])

    def test_edited_old_centers_and_all_existing_h_stay_literal(self):
        project = from_smiles_fixture("C=C")
        carbons = [a["id"] for a in project["atoms"] if a["element"] == "C"]
        next(b for b in project["bonds"] if b["order"] == 2)["order"] = 1
        for i, atom in enumerate(project["atoms"]):
            atom["pos3"][2] += (i + 1) * 0.000000123456789
            if atom["element"] == "H":
                atom["generatedHydrogen"] = i % 2 == 0
        project["metadata"]["geometryPending"] = {"atomIds": [], "changedAtomIds": carbons, "editedAtomIds": carbons}
        original = copy.deepcopy(project)
        with no_force_field(), patch.object(AllChem, "EmbedMolecule", side_effect=AssertionError("old molecule re-embedded")):
            after = dispatch("/api/editor/geometry", {"project": project, "operation": "prepare"})["project"]
        self.assertEqual(after["atoms"][:len(project["atoms"])], project["atoms"])
        self.assertEqual(project, original)

    def test_new_rigid_cage_preserves_all_pair_distances_and_old_manual_h(self):
        parent = from_smiles_fixture("CO")
        anchor = parent["atoms"][0]["id"]
        old_by_id = {a["id"]: a for a in parent["atoms"]}
        replaced_h = next(b["b"] if b["a"] == anchor else b["a"] for b in parent["bonds"]
                          if anchor in (b["a"], b["b"])
                          and old_by_id[b["b"] if b["a"] == anchor else b["a"]]["element"] == "H")
        parent["atoms"] = [a for a in parent["atoms"] if a["id"] != replaced_h]
        parent["bonds"] = [b for b in parent["bonds"] if replaced_h not in (b["a"], b["b"])]
        cage = native_cage(False)
        cage_anchor = cage["atoms"][0]["id"]
        next(b for b in cage["bonds"] if b["order"] == 2 and cage_anchor in (b["a"], b["b"]))["order"] = 1
        old = copy.deepcopy(parent)
        parent["atoms"].extend(cage["atoms"])
        parent["bonds"].extend(cage["bonds"])
        parent["bonds"].append({"id": "cage-parent", "a": anchor, "b": cage_anchor, "order": 1})
        parent["metadata"]["ketcherGroups"] = cage["metadata"]["ketcherGroups"]
        parent["metadata"]["geometryPending"] = {"atomIds": [a["id"] for a in cage["atoms"]], "changedAtomIds": [anchor]}
        with no_force_field(), patch.object(AllChem, "EmbedMolecule", side_effect=AssertionError("known sphere re-embedded")):
            after = dispatch("/api/editor/geometry", {"project": parent, "operation": "prepare"})["project"]
        self.assertEqual(after["atoms"][:len(old["atoms"])], old["atoms"])
        before_points = np.array([a["pos3"] for a in cage["atoms"]])
        after_points = np.array([a["pos3"] for a in after["atoms"][len(old["atoms"]):len(old["atoms"]) + 60]])
        np.testing.assert_allclose(np.linalg.norm(before_points[:, None] - before_points[None, :], axis=2),
                                   np.linalg.norm(after_points[:, None] - after_points[None, :], axis=2), atol=1e-12, rtol=0)

    def test_complex_multi_anchor_region_fails_clearly_without_changing_draft(self):
        project, fragment = attach_new_fragment(native_cage(True), "CCC")
        second_cage = project["atoms"][10]["id"]
        next(b for b in project["bonds"] if b["order"] == 2 and second_cage in (b["a"], b["b"]))["order"] = 1
        project["bonds"].append({"id": "second-bridge", "a": second_cage, "b": fragment["atoms"][-1]["id"], "order": 1})
        original = copy.deepcopy(project)
        with no_force_field(), patch.object(AllChem, "EmbedMolecule", side_effect=AssertionError("invalid bridge embedded")):
            with self.assertRaisesRegex(EditorError, "多个原有固定连接点"):
                dispatch("/api/editor/geometry", {"project": project, "operation": "prepare"})
        self.assertEqual(project, original)

    def test_old_coordinate_stereo_conflict_is_rejected_without_silent_reflection(self):
        project = from_smiles_fixture("C[C@H](O)C(=O)O")
        center = next(a for a in project["atoms"] if a["chiralTag"] != "CHI_UNSPECIFIED")
        center["chiralTag"] = "CHI_TETRAHEDRAL_CW" if center["chiralTag"] == "CHI_TETRAHEDRAL_CCW" else "CHI_TETRAHEDRAL_CCW"
        for bond in project["bonds"]:
            bond["direction"] = "NONE"
        project["metadata"]["geometryPending"] = {"atomIds": [], "editedAtomIds": [center["id"]]}
        original = copy.deepcopy(project)
        with no_force_field(), patch.object(AllChem, "EmbedMolecule", side_effect=AssertionError("old chiral conformer changed")):
            with self.assertRaisesRegex(EditorError, "手性标记与坐标不一致"):
                dispatch("/api/editor/geometry", {"project": project, "operation": "prepare"})
        self.assertEqual(project, original)

    def test_proper_rotations_have_positive_determinant(self):
        for source, target in (([1, 0, 0], [-1, 0, 0]), ([1, 2, 3], [-3, 4, 1]), ([0, 1, 0], [0, 1, 0])):
            rotation = construction._rotation(source, target)
            self.assertAlmostEqual(np.linalg.det(rotation), 1, places=10)
            np.testing.assert_allclose(construction._unit(source) @ rotation.T, construction._unit(target), atol=1e-12, rtol=0)


if __name__ == "__main__":
    unittest.main()
