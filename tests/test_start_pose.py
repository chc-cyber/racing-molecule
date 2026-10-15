"""Rigid placement invariants, ambiguous axes, and preparation request lifecycle."""
import copy
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
from molecular_car.start_pose import (place_start_pose, validate_start_pose, PRESERVE_START_POSE,
                                      _principal_axes)
from molecular_car.molecule import prepare_molecule
from molecular_car.server import GameController


class RigidStartPoseTests(unittest.TestCase):
    def assert_rigid(self, original, placed, report):
        old = np.linalg.norm(original[:, None] - original[None, :], axis=2)
        new = np.linalg.norm(placed[:, None] - placed[None, :], axis=2)
        np.testing.assert_allclose(new, old, atol=2e-12, rtol=0)
        rotation = np.array(report["rotation_matrix"])
        np.testing.assert_allclose(rotation @ rotation.T, np.eye(3), atol=2e-12)
        self.assertAlmostEqual(np.linalg.det(rotation), 1., places=12)
        self.assertAlmostEqual(placed[:, 2].min(), 2.9, places=12)

    def test_manual_three_axis_rotation_keeps_distances_and_chirality(self):
        x = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
        result, report = place_start_pose(x, [12, 12, 16, 14], ["C", "C", "O", "N"],
                                         {"mode": "manual", "yaw_deg": 35, "pitch_deg": -70, "roll_deg": 90})
        self.assert_rigid(x, result, report)
        np.testing.assert_allclose(np.average(result[:, :2], axis=0, weights=[12, 12, 16, 14]), [0, 0], atol=1e-12)

    def test_line_along_surface_normal_becomes_horizontal(self):
        x = np.array([[0., 0., z] for z in (-3., -1., 1., 3.)])
        result, report = place_start_pose(x, [12]*4, ["C"]*4, {"mode": "auto"})
        self.assert_rigid(x, result, report)
        self.assertLess(np.ptp(result[:, 2]), 1e-12)
        self.assertAlmostEqual(np.ptp(result[:, 0]), 6.)
        self.assertEqual(report["strategy"], "automatic-heading-preserve-roll")

    def test_already_forward_axis_is_not_reversed_due_to_atom_order(self):
        x = np.array([[v,0.,0.] for v in (-3.,-1.,1.,3.)])
        result, report = place_start_pose(x,[12]*4,["C"]*4,{"mode":"auto"})
        np.testing.assert_array_equal(report["rotation_matrix"],np.eye(3))
        np.testing.assert_array_equal(result[:,0],x[:,0])

    def test_planar_symmetric_shape_gets_a_plane_without_invented_front(self):
        angles = np.arange(6) * np.pi / 3
        x = np.column_stack([np.cos(angles), np.zeros(6), np.sin(angles)])
        result, report = place_start_pose(x, [12]*6, ["C"]*6, {"mode": "auto"})
        self.assert_rigid(x, result, report)
        self.assertLess(np.ptp(result[:, 2]), 1e-12)
        self.assertFalse(report["heading_confident"])
        self.assertEqual(report["heading_atom_ids"], [])

    def test_near_spherical_shape_keeps_file_rotation(self):
        x = np.array([[1., 1., 1.], [1., -1., -1.], [-1., 1., -1.], [-1., -1., 1.]])
        result, report = place_start_pose(x, [12]*4, ["C"]*4, {"mode": "auto"})
        self.assert_rigid(x, result, report)
        np.testing.assert_array_equal(report["rotation_matrix"], np.eye(3))
        self.assertEqual(report["strategy"], "ambiguous-geometry-input-preserved")

    def test_explicit_heading_overrides_spherical_ambiguity(self):
        x = np.array([[1., 1., 1.], [1., -1., -1.], [-1., 1., -1.], [-1., -1., 1.]])
        result, report = place_start_pose(x, [12]*4, ["C"]*4, {"mode": "auto", "axis_atoms": [1, 2]})
        self.assert_rigid(x, result, report)
        direction = result[1] - result[0]
        self.assertGreater(direction[0], 0)
        np.testing.assert_allclose(direction[1:], [0, 0], atol=1e-12)

    def test_principal_axes_match_independent_symmetric_eigensolver(self):
        x = np.random.default_rng(42).normal(size=(80, 3)) * [5, 2, .3]
        values, axes = _principal_axes(x)
        covariance = (x-x.mean(axis=0)).T @ (x-x.mean(axis=0)) / len(x)
        np.testing.assert_allclose(values, np.linalg.eigvalsh(covariance)[::-1], atol=1e-12)
        for value, axis in zip(values, axes):
            np.testing.assert_allclose(covariance @ axis, value * axis, atol=1e-12)

    def test_preserve_mode_matches_default_molecule_preparation_exactly(self):
        path = Path(__file__).resolve().parents[1] / "examples/tetraphenylbenzene.mol"
        original = prepare_molecule(path.read_text(), path.name)
        preserved = prepare_molecule(path.read_text(), path.name, start_pose=PRESERVE_START_POSE)
        np.testing.assert_array_equal(original.positions, preserved.positions)
        self.assertEqual(original.atoms, preserved.atoms)
        self.assertEqual(original.bonds, preserved.bonds)

    def test_invalid_pose_and_heading_are_rejected(self):
        for value in ({"mode": "other"}, {"yaw_deg": float("nan")}, {"pitch_deg": True},
                      {"roll_deg": 181}, {"axis_atoms": [1, 1]}, {"axis_atoms": [1]},
                      {"axis_atoms": [True, 2]}, {"axis_atoms": [1, 257]}, {"unknown": 1}):
            with self.subTest(value=value), self.assertRaises(ValueError):
                validate_start_pose(value)
        with self.assertRaises(ValueError):
            place_start_pose(np.zeros((2, 3)), [12, 12], ["C", "C"], {"mode": "auto", "axis_atoms": [1, 2]})


class PoseLifecycleTests(unittest.TestCase):
    def controller(self):
        controller = GameController()
        controller.competition["mode"] = "custom"
        controller.input = ("valid text", "valid.mol")
        controller.state.update(status="paused", molecule={"atoms": [{"id": 1}, {"id": 2}]})
        return controller

    def test_invalid_pose_never_replaces_active_generation(self):
        controller = self.controller();before=controller.snapshot()
        with patch.object(controller, "_load_locked") as reload:
            with self.assertRaises(ValueError):
                controller.start_pose({"mode": "auto", "axis_atoms": [1, 3]})
            reload.assert_not_called()
        self.assertEqual(before, controller.snapshot())

    def test_pose_changes_require_pausing_and_reset_preserves_pose(self):
        controller = self.controller();controller.state["status"]="running"
        with patch.object(controller, "_load_locked") as reload:
            with self.assertRaises(ValueError):
                controller.start_pose({"mode": "manual", "yaw_deg": 90})
            reload.assert_not_called()
            controller.state["status"]="paused"
            controller.start_pose({"mode": "manual", "yaw_deg": 90})
            requested=copy.deepcopy(controller.requested_start_pose)
            controller.action("reset")
            self.assertEqual(controller.requested_start_pose, requested)
            self.assertEqual(reload.call_count, 2)


if __name__ == "__main__":
    unittest.main()
