"""Gameplay uses the shortened track for controls, replay and finish progress."""

import ctypes
import unittest
from unittest.mock import Mock

import numpy as np

from molecular_car.engine import DEFAULT_SETTINGS, DEFAULT_TIP, Simulation
from molecular_car.server import GameController


class TrackGameplayTests(unittest.TestCase):
    def test_removed_tip_and_twenty_nm_replay(self):
        controller = GameController()
        controller.state["status"] = "paused"
        controller.commands = Mock()
        with self.assertRaises(ValueError):
            controller.control({"tip": {"x": 206., "y": 100.}})
        with self.assertRaises(ValueError):
            controller.action("pulse")
        controller.commands.put.assert_not_called()
        track = controller.replay()["track"]
        self.assertEqual((track["finish_x"] - track["start_x"]) / 10, 20.)
        self.assertEqual(2 * track["half_width"] / 10, 20.)

    def test_finish_requires_twenty_nm_net_forward_displacement(self):
        simulation = Simulation.__new__(Simulation)
        simulation.n = 1
        simulation.masses = np.ones(1)
        simulation.lmp = Mock()
        simulation.lmp.gather_atoms.return_value = (ctypes.c_double * 3)(0., 0., 0.)
        simulation.lmp.get_thermo.side_effect = lambda name: {"step": 100, "etotal": 0.}[name]
        simulation.lmp.extract_compute.return_value = 300.
        simulation.sequence = 1
        simulation.settings = dict(DEFAULT_SETTINGS)
        simulation.tip = dict(DEFAULT_TIP)
        simulation.drive = {"mode": "gradient", "strength": .15}
        simulation.drive_enabled = True
        simulation.pulse_until = -1
        simulation.surface = None
        for start in (0., .7, -.7, 5.):
            simulation.initial_center = np.array([start, 0., 3.])
            # Lateral motion is not forward distance, and reversing subtracts
            # from progress rather than accumulating travelled path length.
            for distance, y in ((0., 0.), (0., 80.), (199.3, 0.), (200., 0.),
                                (201., 0.), (30., 60.), (-2., 0.)):
                with self.subTest(start=start, distance=distance, y=y):
                    simulation.coordinates = lambda: np.array([[start+distance, y, 3.]])
                    frame = simulation.frame()
                    self.assertAlmostEqual(frame["progress"], min(1., max(0., distance/200.)))
                    self.assertAlmostEqual(frame["distance_angstrom"], distance)
                    self.assertEqual(frame["race_start_x"], start)
                    self.assertEqual(frame["race_length_angstrom"], 200.)
                    self.assertEqual(frame["progress"] >= 1, distance >= 200.)

    def test_truncated_replay_preserves_prepared_origin(self):
        controller = GameController()
        controller.state['runtime'] = {'race': {'rule': 'longitudinal-displacement-v1',
            'start_x': .7, 'finish_x': 200.7, 'length_angstrom': 200.}}
        controller.frames.append({'seq': 800, 'step': 40000, 'race_start_x': .7,
                                  'race_length_angstrom': 200., 'distance_angstrom': 70.})
        replay = controller.replay()
        self.assertTrue(replay['truncated'])
        self.assertEqual(replay['runtime']['race']['finish_x'], 200.7)
        self.assertEqual(replay['frames'][0]['race_start_x'], .7)


if __name__ == "__main__":
    unittest.main()
