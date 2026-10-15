import assert from 'node:assert/strict';
import {test} from 'node:test';
import {raceCourse, raceFinished} from '../../frontend/race-course.mjs';

const track = {start_x: 0, finish_x: 200};
test('finish line follows the prepared origin, including truncated replay frames', () => {
  for (const origin of [.7, -.7]) {
    const frame = {race_start_x: origin, race_length_angstrom: 200, distance_angstrom: 70, seq: 800};
    const course = raceCourse(track, frame);
    assert.equal(course.start_x, origin);
    assert.equal(course.finish_x, origin+200);
    assert.equal(course.rule, 'longitudinal-displacement-v1');
  }
});
test('old replays retain their fixed finish line and do not infer a new origin', () => {
  assert.equal(raceCourse(track, {center: [70.7, 0, 3], distance_angstrom: 70}).finish_x, 200);
  assert.equal(raceCourse(track, null).start_x, 0);
});
test('worker completion follows displacement progress, never absolute center', () => {
  assert.equal(raceFinished({center: [200, 0, 3], progress: 199.3/200}), false);
  assert.equal(raceFinished({center: [199.3, 0, 3], progress: 1}), true);
  assert.equal(raceFinished({center: [0, 200, 3], progress: 0}), false);
});
