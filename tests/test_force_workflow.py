"""Regression of mutable controls, atom ordering and reference force semantics."""
import unittest
from types import SimpleNamespace
from unittest.mock import Mock
import numpy as np
from molecular_car import fields
from molecular_car.engine import Simulation, DEFAULT_SETTINGS, DEFAULT_TIP


class ForceWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.positions = np.array([[1., 2., 3.1], [2.4, -.3, 3.6], [-11., 101., 2.8]])
        self.sim = Simulation.__new__(Simulation)
        self.sim.n = 3
        self.sim.charges = np.array([.1, -.2, .1])
        self.sim.alphas = np.array([.7, 1., .3])
        self.sim.surface = None
        self.sim.settings = dict(DEFAULT_SETTINGS)
        self.sim.tip = dict(DEFAULT_TIP)
        self.sim.drive = {"mode":"gradient", "strength":.75}
        self.sim.drive_enabled = True
        self.sim.pulse_until = -1
        self.sim.callback_error = None
        forcefield = Mock()
        forcefield.CalcEnergy.side_effect = lambda p: .5*float(np.dot(p, p))
        forcefield.CalcGrad.side_effect = lambda p: p
        self.sim.prepared = SimpleNamespace(forcefield=forcefield)
        self.sim.lmp = Mock()

    def check_force(self, *, step=0, positions=None, order=None):
        positions = self.positions if positions is None else positions
        order = np.arange(3) if order is None else np.asarray(order)
        output = np.empty_like(positions)
        self.sim._forces(None, step, 3, order+1, positions[order], output)
        self.sim._check_callback()
        drive = self.sim.effective_drive(step)
        energy, force = fields.energy_and_forces(positions, self.sim.charges,
            self.sim.alphas, self.sim.tip, self.sim.settings["corrugation"], drive)
        np.testing.assert_allclose(output, (force-positions)[order], rtol=2e-14, atol=2e-12)
        actual = self.sim.lmp.fix_external_set_energy_global.call_args.args[1]
        self.assertAlmostEqual(actual, energy+.5*float(np.sum(positions**2)), places=9)

    def test_direct_dictionary_controls_and_update_preserve_force(self):
        self.check_force()
        self.sim.drive.update(strength=.35, growth=.6)
        self.check_force()
        self.sim.settings["corrugation"] = .23
        self.check_force()
        self.sim.update({"drive":{"mode":"gradient", "strength":.25}})
        self.sim.tip.clear()  # A gradient has no dependence on tip parameters.
        self.check_force()

    def test_fixed_drive_and_disabled_preparation_drive(self):
        self.check_force()
        self.check_force(step=99)
        self.check_force(step=100)
        self.sim.drive_enabled = False
        self.check_force(step=99)

    def test_order_changes_and_multiple_calls_preserve_atom_identity(self):
        self.check_force()
        self.check_force(order=[2, 0, 1])
        self.check_force(order=[1, 2, 0], positions=self.positions+[.8, -.2, .1])
        self.check_force()

    def test_duplicate_or_invalid_atom_ids_stop_before_uff(self):
        self.check_force()
        for tags in ([1, 1, 3], [0, 2, 3], [1, 2, 4]):
            self.sim.callback_error = None
            self.sim.prepared.forcefield.reset_mock()
            output = np.ones_like(self.positions)
            self.sim._forces(None, 0, 3, np.asarray(tags), self.positions, output)
            np.testing.assert_array_equal(output, np.zeros_like(output))
            self.sim.prepared.forcefield.CalcEnergy.assert_not_called()
            with self.assertRaisesRegex(RuntimeError, "原子标识"):
                self.sim._check_callback()

    def test_numpy_reference_backend_preserves_the_callback(self):
        self.sim._prepare_force_evaluation("numpy")
        self.assertEqual(self.sim.field_evaluator.metadata()["accelerator"], "numpy")
        self.check_force()
        self.sim.update({"drive":{"mode":"gradient", "strength":.4, "growth":.5}})
        self.check_force(order=[2, 0, 1])


if __name__ == "__main__":
    unittest.main()
