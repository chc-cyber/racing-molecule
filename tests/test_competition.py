"""Rules must reach preparation atomically, and never mutate a running race."""
import copy
import unittest
from unittest.mock import patch, Mock
from molecular_car.competition import validate_competition, engine_parameters
from molecular_car.server import GameController
from tests.test_editor_race_handoff import HandoffTests

class CompetitionTests(unittest.TestCase):
    def test_standard_is_fixed_and_independent(self):
        rules=validate_competition()
        self.assertEqual((rules['temperature'],rules['strength'],rules['lattice_scale'],rules['timestep_fs']), (150,.6,1.5,1))
        rules['start_pose']['axis_atoms'].append(2)
        self.assertEqual(validate_competition()['start_pose']['axis_atoms'], [])
        for patch in ({'strength':.7},{'surface':'continuum'},{'timestep_fs':.5}, {'length_nm':50}, {'tip':{}}, {'pulse':True}):
            with self.subTest(patch=patch),self.assertRaises(ValueError):validate_competition(patch)

    def test_invalid_rules_do_not_replace_race(self):
        game=GameController();before=game.snapshot()
        for rules in ({'mode':'custom','temperature':True},{'mode':'custom','strength':float('nan')}, {'mode':'custom','growth':0}):
            with self.assertRaises(ValueError):game.load('molecule','car.mol',competition=rules)
        self.assertEqual(before,game.snapshot())

    def test_custom_rules_are_passed_before_worker_start_and_survive_reset(self):
        game=GameController();context=Mock()
        rules=validate_competition({'mode':'custom','temperature':300,'strength':.35,'growth':.5,'surface':'continuum','start_pose':{'mode':'preserve'}})
        with patch.object(game,'_stop'),patch('molecular_car.server.mp.get_context',return_value=context),patch('molecular_car.server.threading.Thread'):
            game.load('molecule','car.mol',competition=rules)
            args=context.Process.call_args.kwargs['args']
            self.assertEqual(args[2],engine_parameters(rules)[0]);self.assertEqual(args[5],engine_parameters(rules)[1])
            self.assertEqual(args[6]['mode'],'preserve')
            game.action('reset')
            self.assertEqual(context.Process.call_args.kwargs['args'][2],args[2])
        self.assertEqual(game.replay()['competition'],rules)
        self.assertNotIn('tip',game.snapshot()['control'])

    def test_live_controls_and_pulse_are_removed(self):
        game=GameController();game.state['status']='paused';game.commands=Mock()
        for changes in ({'tip':{'x':5}},{'autopilot':True},{'drive':{'strength':.9}},{'temperature':300}):
            with self.assertRaises(ValueError):game.control(changes)
        with self.assertRaises(ValueError):game.action('pulse')
        game.commands.put.assert_not_called()

class CompetitionHandoffTests(HandoffTests):
    def test_candidate_rules_do_not_depend_on_existing_race(self):
        self.payload['competition']={'mode':'custom','temperature':350,'strength':.3,'growth':.7,'surface':'continuum'}
        job=self.ready();job['ready']['settings']=dict(job['settings']);job['ready']['control']={'drive':dict(job['drive'])}
        self.assertEqual(job['settings']['temperature'],350)
        self.assertEqual(self.game.competition['mode'],'standard')
        with patch('molecular_car.editor.race_handoff.stop_worker'),patch('molecular_car.editor.race_handoff.threading.Thread'):
            self.game.handoff.commit(job['id'])
        self.assertEqual(self.game.competition['mode'],'custom')
        self.assertEqual(self.game.replay()['competition']['strength'],.3)
        self.assertEqual(self.game.previous_replay['competition']['mode'],'standard')

    def test_retry_cannot_change_rules(self):
        self.begin()
        with self.assertRaises(ValueError):self.game.handoff.begin({**self.payload,'competition':{'mode':'custom','temperature':300}})
