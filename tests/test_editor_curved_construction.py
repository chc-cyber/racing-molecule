"""Actual RDKit construction of strained fused aromatics, without UFF."""
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch
import numpy as np
from molecular_car.editor import dispatch


class CurvedConstructionTests(TestCase):
    def test_curved_fused_aromatics(self):
        for name in ("hexahelicene", "spiral"):
            with self.subTest(name=name):
                content = (Path(__file__).parent / "data/workshop" / (name + ".mol")).read_text()
                p = dispatch("/api/editor/import", {"content": content, "filename": name + ".mol"})["project"]
                with patch("rdkit.Chem.AllChem.UFFOptimizeMolecule", side_effect=AssertionError("automatic UFF")), \
                     patch("rdkit.Chem.AllChem.UFFGetMoleculeForceField", side_effect=AssertionError("automatic UFF")):
                    result = dispatch("/api/editor/geometry", {"project": p, "operation": "prepare", "revision": 7})
                q = result["project"]
                self.assertEqual(result["revision"], 7)
                self.assertTrue(q["has3d"])
                self.assertFalse(result["construction"]["optimized"])
                self.assertTrue(any("弯曲" in warning for warning in result["warnings"]))
                self.assertEqual([(a["id"], a["pos2"]) for a in p["atoms"]],
                                 [(a["id"], a["pos2"]) for a in q["atoms"][:len(p["atoms"])]])
                xyz = np.array([a["pos3"] for a in q["atoms"][:len(p["atoms"])]])
                self.assertGreater(np.linalg.svd(xyz - xyz.mean(axis=0), compute_uv=False)[-1], .5)
                again = dispatch("/api/editor/geometry", {"project": q, "operation": "prepare"})["project"]
                self.assertEqual(q, again)
