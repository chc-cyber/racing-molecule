"""Behavioral checks for independent editor data and RDKit geometry services."""
import copy
from pathlib import Path
import unittest

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from molecular_car.editor import dispatch, EditorError
from molecular_car.editor.fragments import benzene, ethynyl, fullerene
from molecular_car.editor.service import to_mol, _stereo_signature

ROOT = Path(__file__).resolve().parents[1]


def load_fixture(name="tetraphenylbenzene.mol"):
    return dispatch("/api/editor/import", {"filename": name, "content": (ROOT / "examples" / name).read_text()})["project"]


def from_smiles_fixture(smiles, three_d=True):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if three_d:
        AllChem.EmbedMolecule(mol, randomSeed=20261011)
        AllChem.UFFOptimizeMolecule(mol)
    else:
        AllChem.Compute2DCoords(mol)
    return dispatch("/api/editor/import", {"filename": "stereo.mol", "content": Chem.MolToMolBlock(mol)})["project"]


def append_unit_phenyl(project, anchor_id):
    """A real native-scale 2D ring plus deliberately provisional planar pos3."""
    project = copy.deepcopy(project)
    by_id = {atom["id"]: atom for atom in project["atoms"]}
    replaced_h = next((bond["b"] if bond["a"] == anchor_id else bond["a"] for bond in project["bonds"]
                       if anchor_id in (bond["a"], bond["b"])
                       and by_id[bond["b"] if bond["a"] == anchor_id else bond["a"]]["element"] == "H"), None)
    if replaced_h:
        project["atoms"] = [atom for atom in project["atoms"] if atom["id"] != replaced_h]
        project["bonds"] = [bond for bond in project["bonds"] if replaced_h not in (bond["a"], bond["b"])]
    fragment = benzene()
    mapping = {atom["id"]: "new-" + atom["id"] for atom in fragment["atoms"]}
    origin = np.array(by_id[anchor_id]["pos3"])
    for atom in fragment["atoms"]:
        atom["id"] = mapping[atom["id"]]
        atom["pos2"] = [value / 1.4 for value in atom["pos2"]]
        atom["pos3"] = (origin + np.array([*atom["pos2"], 0]) + [2, 0, 0]).tolist()
    for bond in fragment["bonds"]:
        bond.update(id="new-" + bond["id"], a=mapping[bond["a"]], b=mapping[bond["b"]])
    project["atoms"].extend(fragment["atoms"])
    project["bonds"].extend(fragment["bonds"])
    project["bonds"].append({"id": "new-attachment", "a": anchor_id, "b": fragment["atoms"][3]["id"], "order": 1})
    project.setdefault("metadata", {})["geometryPending"] = {
        "atomIds": [atom["id"] for atom in fragment["atoms"]], "changedAtomIds": [anchor_id]}
    return project, fragment


def mark_c60_seed(project, atom_ids):
    project.setdefault("metadata", {})["ketcherGroups"] = [{
        "id": "native-c60", "type": "SUP", "kind": "c60", "name": "C60", "expanded": False,
        "atomIds": atom_ids, "attachmentPoints": [{"atomId": atom_ids[0], "attachmentId": "1"}],
        "template": {"id": "c60", "version": 2, "indexToAtomId": atom_ids}}]


class EditorServiceTests(unittest.TestCase):
    def test_import_preserves_input_3d(self):
        p = load_fixture()
        original = Chem.MolFromMolFile(str(ROOT / "examples" / "tetraphenylbenzene.mol"), removeHs=False)
        np.testing.assert_array_equal([a["pos3"] for a in p["atoms"]], original.GetConformer().GetPositions())
        self.assertEqual((len(p["atoms"]), len(p["bonds"])), (52, 56))

    def test_generate_benzene_materializes_hydrogens_and_preserves_old_graph_ids(self):
        p = benzene()
        p["has3d"] = False
        original = copy.deepcopy(p)
        result = dispatch("/api/editor/geometry", {"project": p, "operation": "generate", "revision": 15})
        q = result["project"]
        self.assertEqual(result["revision"], 15)
        self.assertTrue(q["has3d"])
        self.assertEqual(len(q["atoms"]), 12)
        self.assertEqual(q["bonds"][:len(p["bonds"])], p["bonds"])
        self.assertEqual([a["id"] for a in q["atoms"][:6]], [a["id"] for a in p["atoms"]])
        self.assertEqual([a["pos2"] for a in q["atoms"][:6]], [a["pos2"] for a in p["atoms"]])
        self.assertEqual(result["hydrogensAdded"], 6)
        self.assertEqual(len(result["addedBondIds"]), 6)
        self.assertEqual([a["element"] for a in q["atoms"][6:]], ["H"] * 6)
        self.assertEqual(p, original)

    def test_dihedral_survives_layout_and_export(self):
        p = load_fixture()
        mol = to_mol(p, sanitize=True)
        torsion = None
        for bond in mol.GetBonds():
            if bond.GetBondTypeAsDouble() == 1 and not bond.IsInRing():
                b, c = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                aa = [a.GetIdx() for a in mol.GetAtomWithIdx(b).GetNeighbors() if a.GetIdx() != c and a.GetAtomicNum() != 1]
                dd = [a.GetIdx() for a in mol.GetAtomWithIdx(c).GetNeighbors() if a.GetIdx() != b and a.GetAtomicNum() != 1]
                if aa and dd:
                    torsion = (aa[0], b, c, dd[0])
                    break
        self.assertIsNotNone(torsion)
        rdMolTransforms.SetDihedralDeg(mol.GetConformer(), *torsion, 73)
        for atom, point in zip(p["atoms"], mol.GetConformer().GetPositions()):
            atom["pos3"] = point.tolist()
        modified = copy.deepcopy(p)
        q = dispatch("/api/editor/geometry", {"project": p, "operation": "layout2d", "revision": 7})["project"]
        self.assertEqual([a["pos3"] for a in q["atoms"]], [a["pos3"] for a in modified["atoms"]])
        self.assertAlmostEqual(rdMolTransforms.GetDihedralDeg(to_mol(q).GetConformer(), *torsion), 73, places=5)
        exported = dispatch("/api/editor/export", {"project": q, "dimension": 3})
        self.assertIn("3D", exported["content"].splitlines()[1])
        r = dispatch("/api/editor/import", exported)["project"]
        np.testing.assert_allclose([a["pos3"] for a in q["atoms"]], [a["pos3"] for a in r["atoms"]], atol=0.000051)
        self.assertEqual([b["order"] for b in q["bonds"]], [b["order"] for b in r["bonds"]])

    def test_fixed_relaxation_exact_and_nonmutating(self):
        p = load_fixture("benzene.mol")
        p["atoms"][2]["pos3"][2] += 0.45
        original = copy.deepcopy(p)
        fixed = [a["id"] for a in p["atoms"][:2]]
        result = dispatch("/api/editor/geometry", {"project": p, "operation": "relax", "fixedIds": fixed, "revision": 123})
        q = result["project"]
        self.assertEqual(result["revision"], 123)
        self.assertEqual([a["pos3"] for a in p["atoms"][:2]], [a["pos3"] for a in q["atoms"][:2]])
        self.assertNotEqual(p["atoms"][2]["pos3"], q["atoms"][2]["pos3"])
        self.assertEqual(q["bonds"], p["bonds"])
        self.assertEqual(p, original)

    def test_c60_real_graph_valence_and_parameters(self):
        p = fullerene()
        self.assertEqual((len(p["atoms"]), len(p["bonds"])), (60, 90))
        self.assertTrue(p["rigid"])
        self.assertTrue(all(port["maxOrder"] == 1 for port in p["ports"]))
        mol = to_mol(p, sanitize=True)
        self.assertEqual(len(Chem.GetMolFrags(mol)), 1)
        self.assertTrue(all(a.GetDegree() == 3 and a.GetTotalValence() == 4 for a in mol.GetAtoms()))
        self.assertTrue(AllChem.UFFHasAllMoleculeParams(mol))
        valid = dispatch("/api/editor/validate", {"project": p})
        self.assertTrue(valid["eligible"], valid["errors"])
        extra = copy.deepcopy(ethynyl()["atoms"][0])
        p["atoms"].append(extra)
        p["bonds"].append({"id": "bad-port-bond", "a": p["atoms"][0]["id"], "b": extra["id"], "order": 1})
        invalid = dispatch("/api/editor/validate", {"project": p})
        self.assertTrue(invalid["editable"])
        self.assertFalse(invalid["eligible"])

    def test_stereo_charge_isotope_roundtrip(self):
        for smiles in ("F[C@](Cl)(Br)I", "F/C=C/Cl", "F/C=C\\Cl", "[13CH3][NH3+]"):
            with self.subTest(smiles=smiles):
                p = from_smiles_fixture(smiles)
                q = dispatch("/api/editor/import", dispatch("/api/editor/export", {"project": p, "dimension": 3, "format": "sdf"}))["project"]
                self.assertEqual([a["charge"] for a in p["atoms"]], [a["charge"] for a in q["atoms"]])
                self.assertEqual([a["isotope"] for a in p["atoms"]], [a["isotope"] for a in q["atoms"]])
                self.assertEqual([a["chiralTag"] for a in p["atoms"]], [a["chiralTag"] for a in q["atoms"]])
                self.assertEqual([b["stereo"] for b in p["bonds"]], [b["stereo"] for b in q["bonds"]])
                np.testing.assert_allclose([a["pos3"] for a in p["atoms"]], [a["pos3"] for a in q["atoms"]], atol=0.000051)

    def test_conflicting_3d_chirality_rejected(self):
        p = from_smiles_fixture("F[C@](Cl)(Br)I")
        for atom in p["atoms"]:
            atom["pos3"][0] *= -1
        with self.assertRaisesRegex(EditorError, "手性"):
            dispatch("/api/editor/export", {"project": p, "dimension": 3})

    def test_generate_3d_retains_defined_stereo(self):
        for smiles in ("F[C@](Cl)(Br)I", "F/C=C/Cl", "F/C=C\\Cl"):
            with self.subTest(smiles=smiles):
                p = from_smiles_fixture(smiles, three_d=False)
                q = dispatch("/api/editor/geometry", {"project": p, "operation": "generate"})["project"]
                self.assertEqual([a["chiralTag"] for a in p["atoms"]], [a["chiralTag"] for a in q["atoms"]])
                dispatch("/api/editor/export", {"project": q, "dimension": 3})

    def test_multiple_sdf_records_rejected(self):
        text = (ROOT / "examples" / "benzene.mol").read_text()
        with self.assertRaisesRegex(EditorError, "2"):
            dispatch("/api/editor/import", {"filename": "multi.sdf", "content": text + "$$$$\n" + text + "$$$$\n"})

    def test_disconnected_and_empty_designs_are_editable(self):
        p = benzene()
        p["bonds"] = []
        result = dispatch("/api/editor/validate", {"project": p})
        self.assertTrue(result["editable"])
        self.assertFalse(result["eligible"])
        self.assertTrue(result["errors"])
        p["atoms"] = []
        result = dispatch("/api/editor/validate", {"project": p})
        self.assertFalse(result["eligible"])

    def test_race_validation_on_private_copy(self):
        p = load_fixture()
        original = copy.deepcopy(p)
        result = dispatch("/api/editor/validate", {"project": p, "revision": 42})
        self.assertTrue(result["eligible"], result["errors"])
        self.assertEqual(result["preparedAtomCount"], 52)
        self.assertEqual(result["revision"], 42)
        self.assertEqual(p, original)

    def test_malformed_requests_have_readable_value_errors(self):
        for path, payload in (
            ("/api/editor/import", {"filename": "x.mol", "content": "junk"}),
            ("/api/editor/export", {"project": None}),
            ("/api/editor/geometry", {"project": benzene(), "operation": "explode"}),
            ("/api/editor/geometry", {"project": {**benzene(), "fragments": [None]}, "operation": "relax"}),
            ("/api/editor/unknown", {}),
        ):
            with self.subTest(path=path, payload=payload):
                with self.assertRaises(EditorError) as caught:
                    dispatch(path, payload)
                self.assertTrue(str(caught.exception))

    def test_invalid_chirality_marker_is_not_silently_dropped(self):
        p = benzene()
        p["atoms"][0]["chiralTag"] = "CHI_TETRAHEDRAL_CW"
        with self.assertRaisesRegex(EditorError, "手性"):
            dispatch("/api/editor/export", {"project": p, "dimension": 3})

    def test_user_drawn_wedge_and_dash_generate_opposite_correct_stereo(self):
        for has_3d in (False, True):
            observed = []
            for direction in ("BEGINWEDGE", "BEGINDASH"):
                with self.subTest(has_3d=has_3d, direction=direction):
                    p = from_smiles_fixture("FC(Cl)(Br)I", three_d=False)
                    p["has3d"] = has_3d
                    for atom in p["atoms"]:
                        atom["chiralTag"] = "CHI_UNSPECIFIED"
                    for bond in p["bonds"]:
                        bond["direction"] = "NONE"
                    p["bonds"][0].update(a=p["atoms"][1]["id"], b=p["atoms"][0]["id"], direction=direction)
                    original = copy.deepcopy(p)
                    drawn_mol = dispatch("/api/editor/export", {"project": p, "dimension": 2})["content"]
                    expected = _stereo_signature(Chem.MolFromMolBlock(drawn_mol, removeHs=False))[0]
                    q = dispatch("/api/editor/geometry", {"project": p, "operation": "generate"})["project"]
                    generated = dispatch("/api/editor/export", {"project": q, "dimension": 3})["content"]
                    actual = _stereo_signature(Chem.MolFromMolBlock(generated, removeHs=False))[0]
                    self.assertEqual(actual, expected)
                    self.assertNotEqual(q["atoms"][1]["chiralTag"], "CHI_UNSPECIFIED")
                    self.assertEqual(p, original)
                    observed.append(actual[1])
                    # Explicit 2D layout may choose a different wedge bond;
                    # it must regenerate depiction consistently with the tag.
                    laid_out = dispatch("/api/editor/geometry", {"project": q, "operation": "layout2d"})["project"]
                    after_layout = dispatch("/api/editor/export", {"project": laid_out, "dimension": 2})["content"]
                    self.assertEqual(_stereo_signature(Chem.MolFromMolBlock(after_layout, removeHs=False))[0], expected)
            self.assertNotEqual(*observed)

    def test_mirrored_2d_wedge_layout_conflict_is_rejected(self):
        p = from_smiles_fixture("FC(Cl)(Br)I", three_d=False)
        p["bonds"][0].update(a=p["atoms"][1]["id"], b=p["atoms"][0]["id"], direction="BEGINWEDGE")
        p = dispatch("/api/editor/geometry", {"project": p, "operation": "generate"})["project"]
        for atom in p["atoms"]:
            atom["pos2"][0] *= -1
        for path, arguments in (("/api/editor/export", {"dimension": 2}),
                                ("/api/editor/export", {"dimension": 3}),
                                ("/api/editor/geometry", {"operation": "generate"})):
            with self.subTest(path=path, arguments=arguments):
                with self.assertRaisesRegex(EditorError, "楔"):
                    dispatch(path, {"project": p, **arguments})

    def test_hydrogenate_preserves_existing_coordinates_ids_and_is_idempotent(self):
        p = load_fixture()
        keep = {a["id"] for a in p["atoms"] if a["element"] != "H"}
        # Remove only some real H from the input fixture: existing H must be
        # retained and only missing neighbors appended with new stable IDs.
        keep.add(next(a["id"] for a in p["atoms"] if a["element"] == "H"))
        p["atoms"] = [a for a in p["atoms"] if a["id"] in keep]
        p["bonds"] = [b for b in p["bonds"] if b["a"] in keep and b["b"] in keep]
        p["atoms"][0]["pos3"] = [coordinate + 0.123456789012 for coordinate in p["atoms"][0]["pos3"]]
        original = copy.deepcopy(p)
        result = dispatch("/api/editor/geometry", {"project": p, "operation": "hydrogenate", "revision": 24})
        q = result["project"]
        self.assertEqual(result["revision"], 24)
        self.assertEqual(result["hydrogensAdded"], 21)
        self.assertEqual(q["atoms"][:len(p["atoms"])], p["atoms"])
        self.assertEqual(q["bonds"][:len(p["bonds"])], p["bonds"])
        self.assertEqual(len(q["atoms"]), 52)
        self.assertEqual(p, original)
        again = dispatch("/api/editor/geometry", {"project": q, "operation": "hydrogenate"})
        self.assertEqual(again["project"], q)
        self.assertEqual(again["hydrogensAdded"], 0)
        self.assertEqual(again["addedAtomIds"], [])
        self.assertEqual(again["addedBondIds"], [])

    def test_new_hydrogens_have_real_single_bonds_and_follow_rigid_fragments(self):
        p = benzene()
        p["fragments"] = [{"id": "ring-fragment", "name": "rigid benzene", "atomIds": [a["id"] for a in p["atoms"]],
                           "ports": copy.deepcopy(p["ports"]), "rigid": True}]
        q = dispatch("/api/editor/geometry", {"project": p, "operation": "generate"})["project"]
        old_ids = {a["id"] for a in p["atoms"]}
        self.assertEqual(len(q["fragments"][0]["atomIds"]), 12)
        self.assertEqual(q["fragments"][0]["ports"], p["fragments"][0]["ports"])
        self.assertEqual(q["fragments"][0]["id"], "ring-fragment")
        self.assertEqual(len({a["id"] for a in q["atoms"]}), 12)
        for atom in q["atoms"][6:]:
            bonds = [b for b in q["bonds"] if atom["id"] in (b["a"], b["b"])]
            self.assertEqual(len(bonds), 1)
            self.assertEqual(bonds[0]["order"], 1)
            parent_id = bonds[0]["b"] if bonds[0]["a"] == atom["id"] else bonds[0]["a"]
            self.assertIn(parent_id, old_ids)
            parent = next(a for a in q["atoms"] if a["id"] == parent_id)
            self.assertGreater(np.linalg.norm(np.array(atom["pos3"]) - parent["pos3"]), 0.7)
            self.assertLess(np.linalg.norm(np.array(atom["pos3"]) - parent["pos3"]), 1.5)
        self.assertEqual(len(q["bonds"]), 12)

    def test_persistent_hydrogens_preserve_tetrahedral_and_double_bond_stereo(self):
        for smiles in ("C[C@H](O)C(=O)O", "C[C@@H](F)/C=C/Cl", "C/C=C\\C", "[13CH3][NH3+]"):
            with self.subTest(smiles=smiles):
                mol = Chem.MolFromSmiles(smiles)
                AllChem.Compute2DCoords(mol)
                p = dispatch("/api/editor/import", {"filename": "implicit-H.mol", "content": Chem.MolToMolBlock(mol)})["project"]
                before_stereo = _stereo_signature(to_mol(p, 2, sanitize=True))
                q = dispatch("/api/editor/geometry", {"project": p, "operation": "generate"})["project"]
                self.assertGreater(len(q["atoms"]), len(p["atoms"]))
                self.assertEqual(_stereo_signature(to_mol(q, 3, sanitize=True)), before_stereo)
                self.assertEqual([a["pos2"] for a in q["atoms"][:len(p["atoms"])]], [a["pos2"] for a in p["atoms"]])
                self.assertEqual([a["id"] for a in q["atoms"][:len(p["atoms"])]], [a["id"] for a in p["atoms"]])
                self.assertEqual([b["id"] for b in q["bonds"][:len(p["bonds"])]], [b["id"] for b in p["bonds"]])
                self.assertTrue(all(a.get("explicitHs", 0) == 0 for a in q["atoms"]))
                for dimension in (2, 3):
                    out = dispatch("/api/editor/export", {"project": q, "dimension": dimension})
                    restored = dispatch("/api/editor/import", out)["project"]
                    self.assertEqual(_stereo_signature(to_mol(restored, dimension, sanitize=True)), before_stereo)
                second = dispatch("/api/editor/geometry", {"project": q, "operation": "hydrogenate"})
                self.assertEqual(second["project"], q)
                self.assertEqual(second["hydrogensAdded"], 0)

    def test_existing_3d_chirality_and_old_hydrogen_coordinates_survive_hydrogenate(self):
        for smiles in ("C[C@H](O)C(=O)O", "C[C@@H](F)/C=C/Cl"):
            with self.subTest(smiles=smiles):
                mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                AllChem.EmbedMolecule(mol, randomSeed=20261011)
                AllChem.UFFOptimizeMolecule(mol)
                reduced = Chem.RemoveHs(mol)
                p = dispatch("/api/editor/import", {"filename": "heavy-3D.mol", "content": Chem.MolToMolBlock(reduced)})["project"]
                expected = _stereo_signature(to_mol(p, 3, sanitize=True))
                result = dispatch("/api/editor/geometry", {"project": p, "operation": "hydrogenate"})
                q = result["project"]
                self.assertGreater(result["hydrogensAdded"], 0)
                self.assertEqual([a["pos3"] for a in q["atoms"][:len(p["atoms"])]], [a["pos3"] for a in p["atoms"]])
                self.assertEqual([a["pos2"] for a in q["atoms"][:len(p["atoms"])]], [a["pos2"] for a in p["atoms"]])
                self.assertEqual(_stereo_signature(to_mol(q, 3, sanitize=True)), expected)
                dispatch("/api/editor/export", {"project": q, "dimension": 3})

    def test_hydrogen_operations_reject_invalid_valence_without_mutating_editable_drafts(self):
        p = ethynyl()
        p["atoms"][0]["explicitHs"] = 3
        original = copy.deepcopy(p)
        for operation in ("generate", "hydrogenate"):
            with self.subTest(operation=operation):
                with self.assertRaisesRegex(EditorError, "价态"):
                    dispatch("/api/editor/geometry", {"project": p, "operation": operation})
                self.assertEqual(p, original)
        p = benzene()
        p["has3d"] = False
        original = copy.deepcopy(p)
        result = dispatch("/api/editor/geometry", {"project": p, "operation": "hydrogenate"})
        self.assertFalse(result["project"]["has3d"])
        self.assertEqual(result["hydrogensAdded"], 6)
        self.assertEqual(p, original)

    def test_relax_and_layout_do_not_persist_private_hydrogens(self):
        p = benzene()
        for operation in ("relax", "layout2d"):
            with self.subTest(operation=operation):
                q = dispatch("/api/editor/geometry", {"project": p, "operation": operation})["project"]
                self.assertEqual([a["id"] for a in q["atoms"]], [a["id"] for a in p["atoms"]])
                self.assertEqual([b["id"] for b in q["bonds"]], [b["id"] for b in p["bonds"]])
                self.assertEqual(len(q["atoms"]), 6)

    def test_prepare_first_generation_builds_initial_geometry_and_marks_real_h(self):
        p = benzene()
        p["has3d"] = False
        p["metadata"] = {"geometryPending": {"atomIds": [a["id"] for a in p["atoms"]]},
                         "geometryNotice": "provisional coordinates"}
        original = copy.deepcopy(p)
        prepared = dispatch("/api/editor/geometry", {"project": p, "operation": "prepare", "revision": 91})
        q = prepared["project"]
        self.assertEqual(prepared["revision"], 91)
        self.assertNotIn("energy", prepared)
        self.assertFalse(prepared["construction"]["optimized"])
        for original_atom, constructed_atom in zip(p["atoms"], q["atoms"]):
            self.assertEqual({k: v for k, v in original_atom.items() if k != "pos3"},
                             {k: v for k, v in constructed_atom.items() if k != "pos3"})
        by_id = {a["id"]: a for a in q["atoms"]}
        lengths = [np.linalg.norm(np.array(by_id[b["a"]]["pos3"]) - by_id[b["b"]]["pos3"])
                   for b in p["bonds"]]
        self.assertTrue(all(1.3 < length < 1.55 for length in lengths), lengths)
        self.assertEqual(q["bonds"][:6], p["bonds"])
        self.assertEqual(prepared["hydrogensAdded"], 6)
        self.assertTrue(all(a["generatedHydrogen"] for a in q["atoms"][6:]))
        self.assertNotIn("geometryPending", q["metadata"])
        self.assertNotIn("geometryNotice", q["metadata"])
        self.assertEqual(p, original)

    def test_prepare_existing_without_topology_edits_never_relaxes_manual_conformer(self):
        p = load_fixture()
        # Deliberately unoptimized geometry must survive mere view switching.
        p["atoms"][0]["pos3"][2] += 0.123456789012345
        original = copy.deepcopy(p)
        result = dispatch("/api/editor/geometry", {"project": p, "operation": "prepare"})
        self.assertEqual(result["project"], original)
        self.assertEqual(result["hydrogensAdded"], 0)
        self.assertNotIn("energy", result)
        self.assertEqual(p, original)

    def test_prepare_new_unit_scale_ring_has_real_bond_lengths_and_all_old_coordinates_fixed(self):
        old = dispatch("/api/editor/geometry", {"project": benzene(), "operation": "generate"})["project"]
        p, fragment = append_unit_phenyl(old, old["atoms"][0]["id"])
        original = copy.deepcopy(p)
        new_ids = {atom["id"] for atom in fragment["atoms"]}
        result = dispatch("/api/editor/geometry", {"project": p, "operation": "prepare", "revision": 93})
        q = result["project"]
        by_id = {atom["id"]: atom for atom in q["atoms"]}
        for atom in p["atoms"]:
            self.assertEqual(by_id[atom["id"]]["pos2"], atom["pos2"])
            if atom["id"] not in new_ids:
                self.assertEqual(by_id[atom["id"]], atom)
        self.assertEqual(q["bonds"][:len(p["bonds"])], p["bonds"])
        lengths = [np.linalg.norm(np.array(by_id[bond["a"]]["pos3"]) - by_id[bond["b"]]["pos3"])
                   for bond in fragment["bonds"]]
        self.assertTrue(all(1.3 < length < 1.55 for length in lengths), lengths)
        join = p["bonds"][-1]
        self.assertGreater(np.linalg.norm(np.array(by_id[join["a"]]["pos3"]) - by_id[join["b"]]["pos3"]), 1.4)
        self.assertLess(np.linalg.norm(np.array(by_id[join["a"]]["pos3"]) - by_id[join["b"]]["pos3"]), 1.7)
        self.assertEqual(result["hydrogensAdded"], 5)
        self.assertNotIn("energy", result)
        self.assertFalse(result["construction"]["optimized"])
        self.assertEqual(result["revision"], 93)
        self.assertEqual(p, original)
        again = dispatch("/api/editor/geometry", {"project": q, "operation": "prepare"})
        self.assertEqual(again["project"], q)
        self.assertNotIn("energy", again)

    def test_prepare_new_ring_preserves_an_existing_tpb_manual_torsion_exactly(self):
        p = load_fixture()
        old_ids = [a["id"] for a in p["atoms"]]
        mol = to_mol(p, sanitize=True)
        torsion = next((a.GetIdx(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), d.GetIdx())
                       for bond in mol.GetBonds() if not bond.IsInRing() and bond.GetBondTypeAsDouble() == 1
                       for a in bond.GetBeginAtom().GetNeighbors()
                       for d in bond.GetEndAtom().GetNeighbors()
                       if a.GetIdx() != bond.GetEndAtomIdx() and d.GetIdx() != bond.GetBeginAtomIdx()
                       and a.GetAtomicNum() != 1 and d.GetAtomicNum() != 1)
        rdMolTransforms.SetDihedralDeg(mol.GetConformer(), *torsion, 71.23456789)
        for atom, point in zip(p["atoms"], mol.GetConformer().GetPositions()):
            atom["pos3"] = point.tolist()
        anchor = next(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6
                      and any(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors()))
        p, fragment = append_unit_phenyl(p, p["atoms"][anchor]["id"])
        new_ids = {a["id"] for a in fragment["atoms"]}
        q = dispatch("/api/editor/geometry", {"project": p, "operation": "prepare"})["project"]
        by_id = {a["id"]: a for a in q["atoms"]}
        for atom in p["atoms"]:
            if atom["id"] not in new_ids:
                self.assertEqual(by_id[atom["id"]], atom)
        indices = {a["id"]: i for i, a in enumerate(q["atoms"])}
        observed = rdMolTransforms.GetDihedralDeg(to_mol(q).GetConformer(), *(indices[old_ids[i]] for i in torsion))
        self.assertAlmostEqual(observed, 71.23456789, places=8)

    def test_prepare_c60_first_generation_uses_rigid_native_sphere_not_2d_projection(self):
        p = fullerene()
        p["has3d"] = False
        atom_ids = [a["id"] for a in p["atoms"]]
        mark_c60_seed(p, atom_ids)
        p["metadata"]["geometryPending"] = {"atomIds": atom_ids}
        original = copy.deepcopy(p)
        q = dispatch("/api/editor/geometry", {"project": p, "operation": "prepare"})["project"]
        self.assertEqual(len(q["atoms"]), 60)
        self.assertEqual(len(q["bonds"]), 90)
        self.assertEqual(q["bonds"], p["bonds"])
        self.assertEqual([a["pos2"] for a in q["atoms"]], [a["pos2"] for a in p["atoms"]])
        self.assertEqual([a["id"] for a in q["atoms"]], atom_ids)
        before, after = np.array([a["pos3"] for a in p["atoms"]]), np.array([a["pos3"] for a in q["atoms"]])
        np.testing.assert_allclose(np.linalg.norm(before[:, None] - before[None, :], axis=2),
                                   np.linalg.norm(after[:, None] - after[None, :], axis=2), atol=1e-10, rtol=0)
        radius = np.linalg.norm(after - after.mean(axis=0), axis=1)
        self.assertTrue(np.all((radius > 3.4) & (radius < 3.7)))
        self.assertGreater(np.linalg.svd(after - after.mean(axis=0), compute_uv=False)[-1], 10)
        self.assertEqual(q["metadata"]["ketcherGroups"], p["metadata"]["ketcherGroups"])
        self.assertEqual(p, original)

    def test_prepare_appended_c60_keeps_old_manual_coordinates_and_template_pair_distances(self):
        p = from_smiles_fixture("C[C@H](O)C(=O)O")
        p["atoms"][0]["pos3"][2] += 0.000123456789
        old = copy.deepcopy(p)
        fragment = fullerene()
        for atom in fragment["atoms"]:
            atom["pos3"] = (np.array(atom["pos3"]) + [12.3456789, 0, 2]).tolist()
        p["atoms"].extend(fragment["atoms"])
        p["bonds"].extend(fragment["bonds"])
        atom_ids = [a["id"] for a in fragment["atoms"]]
        mark_c60_seed(p, atom_ids)
        p["metadata"]["geometryPending"] = {"atomIds": atom_ids, "changedAtomIds": []}
        expected_stereo = _stereo_signature(to_mol(p, sanitize=True))
        result = dispatch("/api/editor/geometry", {"project": p, "operation": "prepare"})
        q = result["project"]
        self.assertEqual(q["atoms"][:len(old["atoms"])], old["atoms"])
        self.assertEqual(q["bonds"], p["bonds"])
        self.assertEqual(_stereo_signature(to_mol(q, sanitize=True)), expected_stereo)
        before = np.array([a["pos3"] for a in fragment["atoms"]])
        after = np.array([a["pos3"] for a in q["atoms"][len(old["atoms"]):]])
        np.testing.assert_allclose(np.linalg.norm(before[:, None] - before[None, :], axis=2),
                                   np.linalg.norm(after[:, None] - after[None, :], axis=2), atol=1e-10, rtol=0)
        self.assertEqual(result["hydrogensAdded"], 0)

    def test_prepare_connected_c60_preserves_rigid_sphere_and_reasonable_attachment(self):
        p = dispatch("/api/editor/geometry", {"project": benzene(), "operation": "generate"})["project"]
        anchor = p["atoms"][0]["id"]
        by_id = {a["id"]: a for a in p["atoms"]}
        replaced_h = next(bond["b"] if bond["a"] == anchor else bond["a"] for bond in p["bonds"]
                          if anchor in (bond["a"], bond["b"])
                          and by_id[bond["b"] if bond["a"] == anchor else bond["a"]]["element"] == "H")
        p["atoms"] = [a for a in p["atoms"] if a["id"] != replaced_h]
        p["bonds"] = [b for b in p["bonds"] if replaced_h not in (b["a"], b["b"])]
        old = copy.deepcopy(p)
        fragment = fullerene()
        cage_anchor = fragment["atoms"][0]["id"]
        # This is an explicit topology edit, not permission for the geometry
        # service to silently change a saturated C60 attachment's bond orders.
        next(b for b in fragment["bonds"] if b["order"] == 2 and cage_anchor in (b["a"], b["b"]))["order"] = 1
        vertex = np.array(fragment["atoms"][0]["pos3"])
        translation = np.array(by_id[anchor]["pos3"]) - 1.52 * vertex / np.linalg.norm(vertex) - vertex
        for atom in fragment["atoms"]:
            atom["pos3"] = (np.array(atom["pos3"]) + translation).tolist()
        p["atoms"].extend(fragment["atoms"])
        p["bonds"].extend(fragment["bonds"])
        p["bonds"].append({"id": "cage-attachment", "a": anchor, "b": cage_anchor, "order": 1})
        atom_ids = [a["id"] for a in fragment["atoms"]]
        mark_c60_seed(p, atom_ids)
        p["metadata"]["geometryPending"] = {"atomIds": atom_ids, "changedAtomIds": [anchor]}
        result = dispatch("/api/editor/geometry", {"project": p, "operation": "prepare"})
        q = result["project"]
        by_id = {a["id"]: a for a in q["atoms"]}
        self.assertEqual(q["atoms"][:len(old["atoms"])], old["atoms"])
        self.assertEqual(q["bonds"][:len(p["bonds"])], p["bonds"])
        length = np.linalg.norm(np.array(by_id[anchor]["pos3"]) - by_id[cage_anchor]["pos3"])
        self.assertGreater(length, 1.4)
        self.assertLess(length, 1.7)
        self.assertEqual(result["hydrogensAdded"], 1)
        seed = np.array([a["pos3"] for a in fragment["atoms"]])
        after = np.array([by_id[atom_id]["pos3"] for atom_id in atom_ids])
        np.testing.assert_allclose(np.linalg.norm(seed[:, None] - seed[None, :], axis=2),
                                   np.linalg.norm(after[:, None] - after[None, :], axis=2), atol=1e-10, rtol=0)

    def test_prepare_explicit_old_bond_order_edit_keeps_every_existing_atom(self):
        p = from_smiles_fixture("C=C")
        carbon_ids = [a["id"] for a in p["atoms"] if a["element"] == "C"]
        next(b for b in p["bonds"] if b["order"] == 2)["order"] = 1
        p["metadata"]["geometryPending"] = {"atomIds": [], "changedAtomIds": carbon_ids, "editedAtomIds": carbon_ids}
        original = copy.deepcopy(p)
        result = dispatch("/api/editor/geometry", {"project": p, "operation": "prepare"})
        q = result["project"]
        by_id = {a["id"]: a for a in q["atoms"]}
        for atom in p["atoms"]:
            self.assertEqual(by_id[atom["id"]]["pos2"], atom["pos2"])
            self.assertEqual(by_id[atom["id"]], atom)
        length = np.linalg.norm(np.array(by_id[carbon_ids[0]]["pos3"]) - by_id[carbon_ids[1]]["pos3"])
        old_length = np.linalg.norm(np.array(p["atoms"][0]["pos3"]) - p["atoms"][1]["pos3"])
        self.assertEqual(length, old_length)
        self.assertFalse(result["construction"]["optimized"])
        self.assertNotIn("energy", result)
        self.assertEqual(result["hydrogensAdded"], 2)
        self.assertEqual(q["bonds"][:len(p["bonds"])], p["bonds"])
        self.assertEqual(p, original)

    def test_prepare_reports_strained_fixed_boundary_without_rewriting_the_conformer(self):
        p = {"format": "molecular-car-workshop-v1", "name": "fixed boundary",
             "has3d": True, "fragments": [],
             "atoms": [{"id": "old-c1", "element": "C", "charge": 0,
                        "pos2": [0, 0], "pos3": [0, 0, 0]},
                       {"id": "old-c2", "element": "C", "charge": 0,
                        "pos2": [1, 0], "pos3": [10, 0, 0]}],
             "bonds": [{"id": "new-connection", "a": "old-c1", "b": "old-c2", "order": 1}],
             "metadata": {"geometryPending": {"atomIds": [], "changedAtomIds": ["old-c1", "old-c2"]}}}
        original = copy.deepcopy(p)
        result = dispatch("/api/editor/geometry", {"project": p, "operation": "prepare"})
        q = result["project"]
        self.assertEqual(q["atoms"][:2], p["atoms"])
        self.assertEqual(q["bonds"][:1], p["bonds"])
        self.assertEqual(result["hydrogensAdded"], 6)
        # Automatic construction reports the fixed impossible C–C connection
        # without implying an energy minimum or moving either endpoint.
        self.assertNotIn("energy", result)
        self.assertTrue(any("C–C 10.00 Å" in warning and "连接键长异常" in warning
                            and "请移动片段或运行 UFF 优化" in warning for warning in result["warnings"]))
        self.assertEqual(p, original)

    def test_prepare_fixed_connection_warning_is_limited_to_this_pending_region(self):
        base = {"format": "molecular-car-workshop-v1", "name": "boundary scope",
                "has3d": True, "fragments": [],
                "atoms": [{"id": "c1", "element": "C", "charge": 0,
                           "pos2": [0, 0], "pos3": [0, 0, 0]},
                          {"id": "c2", "element": "C", "charge": 0,
                           "pos2": [1, 0], "pos3": [10, 0, 0]}],
                "bonds": [{"id": "cc", "a": "c1", "b": "c2", "order": 1}], "metadata": {}}
        for case in ("unchanged manual stretch", "unrelated pending region", "normal length", "movable endpoint"):
            with self.subTest(case=case):
                p = copy.deepcopy(base)
                if case == "unrelated pending region":
                    p["atoms"].append({"id": "unrelated-c", "element": "C", "charge": 0,
                                       "pos2": [5, 0], "pos3": [20, 0, 0]})
                    p["metadata"]["geometryPending"] = {"atomIds": [], "changedAtomIds": ["unrelated-c"]}
                elif case == "normal length":
                    p["atoms"][1]["pos3"] = [1.5, 0, 0]
                    p["metadata"]["geometryPending"] = {"atomIds": [], "changedAtomIds": ["c1", "c2"]}
                elif case == "movable endpoint":
                    p["metadata"]["geometryPending"] = {"atomIds": [], "changedAtomIds": ["c1", "c2"],
                                                         "editedAtomIds": ["c1"]}
                original = copy.deepcopy(p)
                result = dispatch("/api/editor/geometry", {"project": p, "operation": "prepare"})
                self.assertEqual(any("连接键长异常" in warning for warning in result["warnings"]),
                                 case == "movable endpoint")
                self.assertEqual(p, original)
                self.assertEqual(result["project"]["atoms"][1], p["atoms"][1])
                self.assertEqual(result["project"]["atoms"][:len(p["atoms"])], p["atoms"])

    def test_prepare_does_not_silently_delete_even_generated_h_after_invalid_new_bond(self):
        p = dispatch("/api/editor/geometry", {"project": benzene(), "operation": "generate"})["project"]
        fragment = benzene()
        mapping = {a["id"]: "new-" + a["id"] for a in fragment["atoms"]}
        for atom in fragment["atoms"]:
            atom["id"] = mapping[atom["id"]]
        for bond in fragment["bonds"]:
            bond.update(id="new-" + bond["id"], a=mapping[bond["a"]], b=mapping[bond["b"]])
        p["atoms"].extend(fragment["atoms"])
        p["bonds"].extend(fragment["bonds"])
        p["bonds"].append({"id": "bad-attachment", "a": p["atoms"][0]["id"], "b": fragment["atoms"][0]["id"], "order": 1})
        p["metadata"] = {"geometryPending": {"atomIds": list(mapping.values())}}
        original = copy.deepcopy(p)
        with self.assertRaisesRegex(EditorError, "价态"):
            dispatch("/api/editor/geometry", {"project": p, "operation": "prepare"})
        self.assertEqual(p, original)

    def test_prepare_rejects_invalid_pending_or_flat_c60_seed_without_mutating_draft(self):
        p = load_fixture()
        p["metadata"]["geometryPending"] = {"atomIds": ["missing-editor-UUID"]}
        original = copy.deepcopy(p)
        with self.assertRaisesRegex(EditorError, "原子 ID"):
            dispatch("/api/editor/geometry", {"project": p, "operation": "prepare"})
        self.assertEqual(p, original)
        p = fullerene()
        atom_ids = [a["id"] for a in p["atoms"]]
        mark_c60_seed(p, atom_ids)
        p["metadata"]["geometryPending"] = {"atomIds": atom_ids}
        for atom in p["atoms"]:
            atom["pos3"] = [*atom["pos2"], 0]
        original = copy.deepcopy(p)
        with self.assertRaisesRegex(EditorError, "球壳"):
            dispatch("/api/editor/geometry", {"project": p, "operation": "prepare"})
        self.assertEqual(p, original)

    def test_sup_mol_roundtrip_keeps_actual_chemistry_and_native_display_group(self):
        p = fullerene()
        mark_c60_seed(p, [a["id"] for a in p["atoms"]])
        for dimension in (2, 3):
            with self.subTest(dimension=dimension):
                out = dispatch("/api/editor/export", {"project": p, "dimension": dimension})
                self.assertIn("V3000", out["content"])
                q = dispatch("/api/editor/import", out)["project"]
                self.assertEqual(len(q["atoms"]), 60)
                self.assertEqual(len(q["bonds"]), 90)
                self.assertEqual([a["element"] for a in q["atoms"]], ["C"] * 60)
                self.assertEqual([b["order"] for b in q["bonds"]], [b["order"] for b in p["bonds"]])
                group = q["metadata"]["ketcherGroups"][0]
                self.assertEqual(group["type"], "SUP")
                self.assertEqual(group["name"], "C60")
                self.assertFalse(group["expanded"])
                self.assertEqual(group["atomIds"], [a["id"] for a in q["atoms"]])
                self.assertEqual(group["attachmentPoints"], [{"atomId": q["atoms"][0]["id"], "attachmentId": "1"}])
                self.assertNotIn("kind", group)
                if dimension == 3:
                    np.testing.assert_allclose([a["pos3"] for a in q["atoms"]], [a["pos3"] for a in p["atoms"]], atol=0.000051)


if __name__ == "__main__":
    unittest.main(verbosity=2)
