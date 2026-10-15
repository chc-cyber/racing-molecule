"""Surface replacement and lifecycle checks; no LAMMPS process is started."""

import unittest
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np

from molecular_car import fields
from molecular_car.engine import (DEFAULT_SETTINGS, DEFAULT_TIP, Simulation,
                                  drive_and_walls_energy_and_forces, validate_surface_mode,
                                  periodic_surface_error_message)
from molecular_car.server import GameController


class SurfaceCompositionTests(unittest.TestCase):
    def setUp(self):
        self.positions = np.array([[.3, .8, 3.1], [-13, 103, 3.4], [213, -102, 2.7]])
        self.charges = np.array([.07, -.1, .03])
        self.alphas = np.array([1., .2, .9])

    def test_default_and_invalid_modes(self):
        self.assertEqual(DEFAULT_SETTINGS["surface"], "fourier-au")
        self.assertEqual(DEFAULT_SETTINGS["temperature"], 150)
        for mode in ("continuum", "periodic-au", "fourier-au"):
            self.assertEqual(validate_surface_mode(mode), mode)
        for mode in (None, True, {}, "mobile-au", "PERIODIC-AU"):
            with self.assertRaises(ValueError):
                validate_surface_mode(mode)

    def test_replacement_keeps_both_drivers_and_walls_without_morse(self):
        zeros = np.zeros(len(self.positions))
        for mode in ("gradient", "tip"):
            drive = {"mode": mode, "strength": .75}
            total, force = fields.energy_and_forces(
                self.positions, self.charges, self.alphas, DEFAULT_TIP, 0, drive)
            baseline, baseline_force = fields.energy_and_forces(
                self.positions, zeros, zeros, DEFAULT_TIP, 0, {"mode": mode, "strength": 0})
            walls = 0.; wall_force = np.zeros_like(self.positions)
            for axis, low, high in ((0, fields.TRACK_X_MIN, fields.TRACK_X_MAX),
                                    (1, -fields.LANE_HALF_WIDTH, fields.LANE_HALF_WIDTH)):
                e, f = fields._soft_boundary(self.positions[:, axis], low, high)
                walls += float(e.sum()); wall_force[:, axis] = f
            actual, actual_force = drive_and_walls_energy_and_forces(
                self.positions, self.charges, self.alphas, DEFAULT_TIP, drive)
            self.assertAlmostEqual(actual, total - baseline + walls, places=11)
            np.testing.assert_allclose(actual_force, force - baseline_force + wall_force, atol=1e-12)

    def test_replacement_force_matches_energy_gradient(self):
        step = 1e-5
        for mode in ("gradient", "tip"):
            drive = {"mode": mode, "strength": .75}
            _, actual = drive_and_walls_energy_and_forces(
                self.positions, self.charges, self.alphas, DEFAULT_TIP, drive)
            numerical = np.zeros_like(actual)
            for atom in range(len(actual)):
                for axis in range(3):
                    plus = self.positions.copy(); minus = self.positions.copy()
                    plus[atom, axis] += step; minus[atom, axis] -= step
                    ep, _ = drive_and_walls_energy_and_forces(plus, self.charges, self.alphas, DEFAULT_TIP, drive)
                    em, _ = drive_and_walls_energy_and_forces(minus, self.charges, self.alphas, DEFAULT_TIP, drive)
                    numerical[atom, axis] = -(ep - em) / (2*step)
            np.testing.assert_allclose(actual, numerical, atol=2e-7, rtol=2e-7)

    def test_absolute_driver_does_not_repeat_with_surface_cell(self):
        original = np.array([[1., 2., 3.]])
        moved = original + [100, 0, 0]
        drive = {"mode": "gradient", "strength": .75}
        _, f0 = drive_and_walls_energy_and_forces(original, np.zeros(1), np.ones(1), DEFAULT_TIP, drive)
        _, f1 = drive_and_walls_energy_and_forces(moved, np.zeros(1), np.ones(1), DEFAULT_TIP, drive)
        self.assertGreater(f1[0, 0], f0[0, 0])
        self.assertAlmostEqual(f1[0, 0]-f0[0, 0], fields.GRADIENT_POLARIZATION_GAIN*(fields.GRADIENT_FIELD_SLOPE*.75)**2*100)

    def test_unsupported_element_rejected_before_cache_and_lammps(self):
        atoms = [{"element": e, "mass": 12, "charge": 0, "polarizability": 1} for e in ("C", "Fe")]
        prepared = SimpleNamespace(atoms=atoms, positions=np.array([[0., 0., 3.], [1., 0., 3.]]))
        for mode in ("periodic-au", "fourier-au"):
            with self.subTest(mode=mode), patch("molecular_car.engine.prepare_molecule", return_value=prepared), \
                    patch("molecular_car.engine.open_lammps") as open_lammps:
                with self.assertRaisesRegex(ValueError, "Fe"):
                    Simulation("unused", "test.mol", {"surface": mode})
                open_lammps.assert_not_called()

    def test_fourier_factory_records_actual_model_without_substitution(self):
        atoms = [{"element": "C", "mass": 12, "charge": 0, "polarizability": 1}]
        prepared = SimpleNamespace(atoms=atoms, positions=np.array([[0., 0., 3.]]), bonds=[], start_pose={})
        surface = Mock()
        surface.elements = ['C']
        surface.period = np.array([2.88, 4.99])
        surface.last_evaluation_counts = {"full_proxy": 1, "blended_proxy": 0, "exact_direction_fallback": 0}
        surface.metadata.return_value = {"algorithm_version": "test-fourier", "actual_kind": "harmonic-proxy",
                                         "last_evaluation_counts": dict(surface.last_evaluation_counts)}
        factory = Mock(return_value=surface)
        context = Mock(); context.__enter__ = Mock(return_value=(Mock(), {})); context.__exit__ = Mock()
        with patch("molecular_car.engine.prepare_molecule", return_value=prepared), \
                patch("molecular_car.engine.open_lammps", return_value=context), \
                patch.dict("sys.modules", {"molecular_car.fourier_atlas": SimpleNamespace(FourierAtlasPotential=factory)}), \
                patch.object(Simulation, "_initialize"):
            simulation = Simulation("unused", "test.mol", {"surface": "fourier-au"})
        factory.assert_called_once_with(prepared.atoms, prepared.positions, bonds=[],
                                       normal_policy="bonded-saip-game-v1")
        self.assertIs(simulation.surface.base, surface)
        self.assertEqual(simulation.runtime["model"], "uff-fourier-au-bonded-normal-gradient-v2")
        self.assertEqual(simulation.runtime["surface"]["mode"], "fourier-au")
        self.assertEqual(simulation.runtime["surface"]["actual_kind"], "harmonic-proxy")
        self.assertNotIn("last_evaluation_counts", simulation.runtime["surface"])
        self.assertEqual(simulation.runtime["surface"]["preparation_evaluation_counts"], surface.last_evaluation_counts)
        simulation.close()

    def test_worker_rejects_hot_surface_change(self):
        simulation = Simulation.__new__(Simulation)
        simulation.settings = {**DEFAULT_SETTINGS, "surface": "continuum"}
        for mode in ("periodic-au", "fourier-au"):
            with self.assertRaisesRegex(ValueError, "重新准备"):
                simulation.update({"surface": mode})

    def test_normal_rejection_explains_reason_without_hiding_original_error(self):
        for error, chinese in (("SAIP atom 1 has more than three normal neighbours; this geometry is unsupported.", "超过 3"),
                               ("A SAIP local plane has zero normal magnitude.", "共线")):
            actual = periodic_surface_error_message(ValueError(error))
            self.assertIn(chinese, actual)
            self.assertIn(error, actual)

    def test_callback_adds_one_periodic_surface_and_keeps_atom_identity(self):
        simulation = Simulation.__new__(Simulation)
        simulation.n = len(self.positions)
        simulation.charges, simulation.alphas = self.charges, self.alphas
        simulation.settings = {**DEFAULT_SETTINGS, "surface": "periodic-au", "corrugation": .6}
        simulation.tip = dict(DEFAULT_TIP)
        simulation.drive = {"mode": "gradient", "strength": .75}
        simulation.drive_enabled = True; simulation.pulse_until = -1
        simulation.callback_error = None
        uff = Mock()
        uff.CalcEnergy.side_effect = lambda p: .5*float(np.dot(p, p))
        uff.CalcGrad.side_effect = lambda p: p
        simulation.prepared = SimpleNamespace(forcefield=uff)
        surface = Mock()
        surface.energy_and_forces.side_effect = lambda p: (float(np.sum(p[:, 2]**2)),
            np.column_stack((np.zeros(len(p)), np.zeros(len(p)), -2*p[:, 2])))
        simulation.surface = surface; simulation.lmp = Mock()
        order = np.array([2, 0, 1]); output = np.empty_like(self.positions)
        simulation._forces(None, 0, 3, order+1, self.positions[order], output)
        energy, force = drive_and_walls_energy_and_forces(
            self.positions, self.charges, self.alphas, DEFAULT_TIP, simulation.drive)
        expected = force - self.positions
        expected[:, 2] -= 2*self.positions[:, 2]
        np.testing.assert_allclose(output, expected[order], atol=1e-12)
        surface.energy_and_forces.assert_called_once()
        np.testing.assert_array_equal(surface.energy_and_forces.call_args.args[0], self.positions)
        actual_energy = simulation.lmp.fix_external_set_energy_global.call_args.args[1]
        self.assertAlmostEqual(actual_energy, energy + .5*float(np.sum(self.positions**2)) + float(np.sum(self.positions[:, 2]**2)))
        self.assertIsNone(simulation.callback_error)

    def test_fourier_callback_failure_is_reported_without_switching_surface(self):
        simulation = Simulation.__new__(Simulation)
        simulation.n = len(self.positions)
        simulation.charges, simulation.alphas = self.charges, self.alphas
        simulation.settings = {**DEFAULT_SETTINGS, "surface": "fourier-au"}
        simulation.tip = dict(DEFAULT_TIP)
        simulation.drive = {"mode": "gradient", "strength": .75}
        simulation.drive_enabled = True; simulation.pulse_until = -1
        simulation.callback_error = None
        uff = Mock(); uff.CalcEnergy.return_value = 0.; uff.CalcGrad.return_value = np.zeros(self.positions.size)
        simulation.prepared = SimpleNamespace(forcefield=uff)
        surface = Mock(); surface.energy_and_forces.side_effect = ValueError("Unsupported Fourier geometry.")
        simulation.surface = surface; simulation.lmp = Mock()
        output = np.ones_like(self.positions)
        simulation._forces(None, 0, simulation.n, np.arange(simulation.n) + 1, self.positions, output)
        self.assertIs(simulation.surface, surface)
        self.assertEqual(simulation.settings["surface"], "fourier-au")
        np.testing.assert_array_equal(output, np.zeros_like(output))
        with self.assertRaisesRegex(RuntimeError, "Unsupported Fourier geometry"):
            simulation._check_callback()


class SurfaceLifecycleTests(unittest.TestCase):
    def setUp(self):
        self.game = GameController()
        # These lifecycle cases start from an explicitly selected continuum race.
        self.game.requested_settings["surface"] = "continuum"
        self.game.state["settings"]["surface"] = "continuum"
        self.game.input = ("content", "car.mol")
        self.game.state.update(status="running", frame={"step": 100, "seq": 1},
                               molecule={"name": "old"}, runtime={"surface": {"mode": "continuum"}})
        self.game.frames.append({"step": 100, "seq": 1})

    def test_surface_switch_restarts_worker_and_clears_old_generation(self):
        self.assert_mode_restarts_worker("periodic-au")

    def test_fourier_surface_switch_restarts_worker_and_clears_old_generation(self):
        self.assert_mode_restarts_worker("fourier-au")

    def assert_mode_restarts_worker(self, mode):
        context = Mock()
        with patch.object(self.game, "_stop") as stop, \
                patch("molecular_car.server.mp.get_context", return_value=context), \
                patch("molecular_car.server.threading.Thread"):
            self.game.load("content", "car.mol", competition={"mode":"custom", "surface":mode})
        stop.assert_called_once()
        context.Process.return_value.start.assert_called_once()
        self.assertEqual(self.game.generation, 1)
        self.assertEqual(self.game.state["status"], "preparing")
        self.assertEqual(self.game.state["settings"]["surface"], mode)
        self.assertEqual(len(self.game.frames), 0)
        for key in ("frame", "molecule", "runtime", "stats"):
            self.assertIsNone(self.game.state[key])
        worker_args = context.Process.call_args.kwargs["args"]
        self.assertEqual(worker_args[2]["surface"], mode)

    def test_same_surface_preserves_running_race(self):
        with patch.object(self.game, "_load_locked") as load:
            with self.assertRaises(ValueError):
                self.game.control({"surface": "continuum"})
            load.assert_not_called()
        self.assertEqual(self.game.state["frame"]["step"], 100)

    def test_invalid_or_mixed_changes_do_not_mutate_settings(self):
        for changes in ({"surface": True}, {"surface": "mobile-au"},
                        {"surface": "periodic-au", "temperature": 600},
                        {"surface": "fourier-au", "temperature": 600}):
            with self.assertRaises(ValueError):
                self.game.control(changes)
        self.assertEqual(self.game.requested_settings["surface"], "continuum")

    def test_error_can_recover_with_new_preparation(self):
        self.game.state["status"] = "error"
        with patch.object(self.game, "_load_locked") as load:
            self.game.load("content", "car.mol", competition={"mode":"custom","surface":"continuum"})
            load.assert_called_once_with("content", "car.mol")

    def test_periodic_corrugation_rejected_and_replay_records_actual_metadata(self):
        for mode in ("periodic-au", "fourier-au"):
            self.game.requested_settings["surface"] = mode
            with self.assertRaisesRegex(ValueError, "赛事参数"):
                self.game.control({"corrugation": .2})
            self.game.state["runtime"] = {"model": f"uff-{mode}-gradient-soft-tip-v1",
                                          "surface": {"mode": mode, "version": "test-cache-v1"}}
            replay = self.game.replay()
            self.assertEqual(replay["surface"]["version"], "test-cache-v1")
            self.assertEqual(replay["model"], f"uff-{mode}-gradient-soft-tip-v1")


if __name__ == "__main__":
    unittest.main()
