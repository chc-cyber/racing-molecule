"""Prevent incomplete, censored or airborne trials from becoming race means."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from analyze_crystal_contrast import CARS, CONDITIONS, panel_shape, summarize


class CrystalAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.jobs=[dict(molecule=car,temperature_K=150.,strength=.6,lattice_scale=2.,growth=1.,seed=seed)
                   for car in CARS for seed in (1,2,3)]
        self.rows=[dict(config=c,status='passed',competitive_valid=True,outcome='finish',
                        finish_time_ps=100.+c['seed'],final_forward_nm=20.,sustained_lift=False)
                   for c in self.jobs]

    def test_missing_seed_cannot_be_reported_as_complete(self):
        group=summarize(self.rows[:-1],self.jobs)[0]
        self.assertFalse(group['all_finished_valid'])
        car=group['cars'][-1]
        self.assertEqual((car['planned'],car['observed'],car['missing']),(3,2,1))
        self.assertIsNone(car['finish_mean_ps'])

    def test_timeout_is_not_averaged_as_a_finish(self):
        self.rows[-1].update(outcome='time-limit',finish_time_ps=None,final_forward_nm=.2)
        group=summarize(self.rows,self.jobs)[0]
        self.assertFalse(group['all_finished_valid'])
        self.assertEqual(group['cars'][-1]['censored'],1)
        self.assertEqual(group['cars'][-1]['finished'],2)
        self.assertIsNone(group['cars'][-1]['finish_mean_ps'])
        self.assertNotIn('slow_fast_ratio',group)

    def test_airborne_finish_cannot_improve_surface_score(self):
        self.rows[-1].update(competitive_valid=False,sustained_lift=True,finish_time_ps=1.)
        group=summarize(self.rows,self.jobs)[0]
        self.assertEqual(group['cars'][-1]['lifted'],1)
        self.assertEqual(group['cars'][-1]['valid'],2)
        self.assertIsNone(group['cars'][-1]['finish_mean_ps'])

    def test_midpoint_conditions_are_ordered_and_not_dropped(self):
        jobs=[];rows=[]
        for scale,growth in reversed(CONDITIONS):
            for row in self.rows:
                config=dict(row['config'],lattice_scale=scale,growth=growth)
                jobs.append(config)
                rows.append(dict(row,config=config))
        groups=summarize(rows,jobs)
        self.assertEqual(tuple((g['lattice_scale'],g['growth']) for g in groups),CONDITIONS)
        self.assertEqual(len(groups),7)
        self.assertTrue(all(g['all_finished_valid'] for g in groups))
        self.assertEqual(sum(c['observed'] for g in groups for c in g['cars']),147)

    def test_midpoint_uses_sample_sd_not_standard_error(self):
        for row in self.rows:row['config']['lattice_scale']=1.5
        group=summarize(self.rows,self.jobs)[0]
        self.assertEqual(group['lattice_scale'],1.5)
        self.assertEqual(group['cars'][0]['finish_mean_ps'],102.)
        self.assertEqual(group['cars'][0]['finish_sd_ps'],1.)
        self.assertAlmostEqual(group['mean_within_cv_percent'],100./102.)

    def test_panel_grid_holds_every_condition(self):
        for count in (1,5,7,10):
            rows,columns=panel_shape(count)
            self.assertGreaterEqual(rows*columns,count)
            self.assertLess(rows*columns-count,columns)
        self.assertEqual(panel_shape(7),(3,3))
        with self.assertRaises(ValueError):panel_shape(0)


if __name__=='__main__':unittest.main()
