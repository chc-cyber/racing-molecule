"""A native minimizer's early return must not admit an unrelaxed race."""
import ctypes
import unittest
from unittest.mock import patch
import numpy as np

from molecular_car.relaxation import (relax_joint_geometry, JointRelaxationError,
                                      INITIAL_MAX_ATOM_FORCE)


class NativeMinimizerDouble:
    """Emulate native termination states, not the minimization algorithm."""
    def __init__(self, states, initial=.5):
        self.states = iter(states)
        self.force = initial
        self.step = 0
        self.commands = []
        self.callback = None
        self.soft_stopped = False
        self.clock = 0.
        self.stage_seconds = 0.
        self.positions = np.zeros((1, 3))

    def set_fix_external_callback(self, fix, callback):
        self.callback = callback

    def command(self, command):
        self.commands.append(command)
        if command == "run 0":
            if not self.soft_stopped:
                self.callback(None, self.step, 1, np.array([1]), np.zeros((1, 3)), np.zeros((1, 3)))
        elif command.startswith("minimize"):
            self.step += 1
            self.force = next(self.states)
            self.callback(None, self.step, 1, np.array([1]), np.zeros((1, 3)), np.zeros((1, 3)))
            self.clock += self.stage_seconds
        elif command == "timer timeout off every 10":
            self.soft_stopped = False

    def force_timeout(self):
        self.soft_stopped = True

    def gather_atoms(self, name, kind, count):
        if name == "x":
            return (ctypes.c_double * 3)(*self.positions.ravel())
        return (ctypes.c_double * 3)(self.force, 0., 0.)

    def scatter_atoms(self, name, kind, count, values):
        self.positions[:] = np.array(list(values)).reshape(1, 3)

    def get_thermo(self, name):
        # Identical energies deliberately emulate an energy-flat early stop.
        return self.step if name == "step" else -1.


class JointForceGateTests(unittest.TestCase):
    def setUp(self):
        self.callback = lambda *args: None

    def assert_cleanup(self, native):
        self.assertIs(native.callback, self.callback)
        self.assertIn("timer timeout off every 10", native.commands)
        self.assertFalse(native.soft_stopped)
        self.assertFalse(any(command.startswith("timestep") for command in native.commands))

    def test_energy_flat_fire_return_needs_another_force_converged_stage(self):
        native = NativeMinimizerDouble([.2, 3e-5])
        report = relax_joint_geometry(native, 1, self.callback, lambda: None)
        self.assertTrue(report["force_tolerance_met"])
        self.assertEqual([stage["minimizer"] for stage in report["stages"]], ["fire", "hftn"])
        self.assertFalse(report["stages"][0]["force_tolerance_met"])
        self.assertEqual(report["termination_reason"], "independently-verified-force-tolerance")
        self.assertTrue(all(command.split()[1] == "0.0" for command in native.commands if command.startswith("minimize")))
        self.assert_cleanup(native)

    def test_finite_energy_and_native_returns_cannot_start_unrelaxed_race(self):
        native = NativeMinimizerDouble([.4, .3, .2, .1])
        with self.assertRaises(JointRelaxationError) as caught:
            relax_joint_geometry(native, 1, self.callback, lambda: None)
        diagnostics = caught.exception.diagnostics
        self.assertFalse(diagnostics["force_tolerance_met"])
        self.assertEqual(diagnostics["global_force_norm_kcal_mol_A"], .1)
        self.assertEqual([stage["minimizer"] for stage in diagnostics["stages"][:4]],
                         ["fire", "hftn", "cg", "hftn"])
        self.assertEqual(diagnostics["stages"][-1]["minimizer"], "curvature-initialized-bfgs")
        self.assertFalse(diagnostics["stages"][-1]["accepted"])
        self.assertEqual(diagnostics["stages"][-1]["termination_reason"], "line-search-stalled")
        self.assert_cleanup(native)

    def test_total_evaluation_limit_soft_stops_without_zeroing_real_force(self):
        native = NativeMinimizerDouble([.2, .1])
        with self.assertRaises(JointRelaxationError) as caught:
            relax_joint_geometry(native, 1, self.callback, lambda: None, evaluation_limit=2)
        self.assertEqual(caught.exception.diagnostics["termination_reason"], "force-evaluation-budget")
        self.assertEqual(caught.exception.diagnostics["global_force_norm_kcal_mol_A"], .2)
        self.assert_cleanup(native)

    def test_wall_budget_fails_before_trying_another_stage(self):
        native = NativeMinimizerDouble([.2])
        native.stage_seconds = 35.
        with patch("molecular_car.relaxation.time.perf_counter", side_effect=lambda: native.clock):
            with self.assertRaises(JointRelaxationError) as caught:
                relax_joint_geometry(native, 1, self.callback, lambda: None, wall_time_limit=30.)
        self.assertEqual(caught.exception.diagnostics["termination_reason"], "wall-time-budget")
        self.assertEqual(len(caught.exception.diagnostics["stages"]), 1)
        self.assert_cleanup(native)

    def test_budget_cannot_admit_a_race_even_when_last_force_is_small(self):
        native = NativeMinimizerDouble([3e-5])
        with self.assertRaises(JointRelaxationError) as caught:
            relax_joint_geometry(native, 1, self.callback, lambda: None, evaluation_limit=2)
        report = caught.exception.diagnostics
        self.assertTrue(report["force_tolerance_met"])
        self.assertTrue(report["final_coordinate_force_verified"])
        self.assertEqual(report["final_verification_force_callback_evaluations"], 1)
        self.assertFalse(report["accepted"])
        self.assert_cleanup(native)

    def test_timer_setup_failure_still_restores_callback(self):
        native = NativeMinimizerDouble([])
        original_command = native.command
        calls = 0
        def fail_first_timer(command):
            nonlocal calls
            calls += 1
            if calls == 1:
                raise RuntimeError("timer setup failed")
            original_command(command)
        native.command = fail_first_timer
        with self.assertRaisesRegex(RuntimeError, "timer setup failed"):
            relax_joint_geometry(native, 1, self.callback, lambda: None)
        self.assert_cleanup(native)

    def test_pathological_initial_force_is_rejected_before_minimization(self):
        native = NativeMinimizerDouble([], initial=2 * INITIAL_MAX_ATOM_FORCE)
        with self.assertRaises(JointRelaxationError) as caught:
            relax_joint_geometry(native, 1, self.callback, lambda: None)
        self.assertEqual(caught.exception.diagnostics["termination_reason"], "pathological-initial-energy-or-force")
        self.assertEqual(caught.exception.diagnostics["iterations"], 0)
        self.assertFalse(any(command.startswith("minimize") for command in native.commands))
        self.assert_cleanup(native)

    def test_callback_validation_error_also_restores_native_state(self):
        native = NativeMinimizerDouble([])
        def check():
            raise RuntimeError("actual force evaluator failed")
        with self.assertRaisesRegex(RuntimeError, "actual force evaluator failed"):
            relax_joint_geometry(native, 1, self.callback, check)
        self.assert_cleanup(native)

    def test_initial_force_gate_does_not_override_an_expired_empty_stage_budget(self):
        native = NativeMinimizerDouble([], initial=3e-5)
        def slow_force(*args):
            native.clock = 35.
        with patch("molecular_car.relaxation.time.perf_counter", side_effect=lambda: native.clock):
            with self.assertRaises(JointRelaxationError) as caught:
                relax_joint_geometry(native, 1, slow_force, lambda: None, stages=())
        self.assertTrue(caught.exception.diagnostics["force_tolerance_met"])
        self.assertFalse(caught.exception.diagnostics["accepted"])
        self.assertEqual(caught.exception.diagnostics["termination_reason"], "wall-time-budget")
        self.assertIs(native.callback, slow_force)
        self.assertFalse(native.soft_stopped)


if __name__ == "__main__":
    unittest.main()
