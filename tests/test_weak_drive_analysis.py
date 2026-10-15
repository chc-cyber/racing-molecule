"""Scientific window metrics distinguish backward motion, pausing and censoring."""
import unittest
import numpy as np
from scripts.scan_weak_drive import trajectory_metrics


class WeakDriveAnalysisTests(unittest.TestCase):
    def test_longitudinal_pause_can_include_lateral_motion(self):
        centers = np.array([[0., 0., .3], [.03, .1, .3], [-.01, .2, .3], [0., .3, .3]])
        metrics = trajectory_metrics([0., 1., 2., 3.], centers)
        self.assertAlmostEqual(metrics["forward_fraction"], 1/3)
        self.assertAlmostEqual(metrics["backward_fraction"], 1/3)
        self.assertAlmostEqual(metrics["pause_fraction"], 1/3)
        self.assertAlmostEqual(metrics["sampled_lateral_path_nm"], .3)
        self.assertTrue(metrics["right_censored"])
        self.assertIsNone(metrics["finish_sample_time_ps"])

    def test_atomic_scale_oscillations_do_not_become_backward_windows(self):
        time = np.linspace(0., 4., 41)
        centers = np.column_stack([.04*time + .03*np.sin(2*np.pi*time), time*0, time*0+.3])
        metrics = trajectory_metrics(time, centers)
        self.assertEqual(metrics["coarse_window_count"], 4)
        self.assertEqual(metrics["forward_fraction"], 1.)
        self.assertEqual(metrics["backward_fraction"], 0.)
        self.assertGreater(np.sum(np.diff(centers[:, 0]) < 0), 0)

    def test_first_sampled_finish_and_partial_window(self):
        metrics = trajectory_metrics([0., 1., 1.2], [[19.9, 0., .3], [19.98, 0., .3], [20.01, 0., .3]])
        self.assertTrue(metrics["finish_observed"])
        self.assertFalse(metrics["right_censored"])
        self.assertEqual(metrics["finish_sample_time_ps"], 1.2)
        self.assertEqual(metrics["coarse_window_count"], 1)


if __name__ == "__main__":
    unittest.main()
