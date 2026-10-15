import test from 'node:test';
import assert from 'node:assert/strict';
import {raceResult} from '../frontend/race-result.mjs';
const state={status:'finished',finish_reason:'finish',frame:{time_fs:100000,distance_angstrom:200,progress:1},molecule:{name:'car'},competition:{mode:'standard'}};
test('20 nm / 100 ps is 200 m/s, independent of instantaneous speed and playback',()=>{
  const r=raceResult({...state,frame:{...state.frame,speed_angstrom_ps:90},playbackRate:4});
  assert.equal(r.averageSpeed,200);assert.equal(r.distanceNm,20);assert.equal(r.timePs,100);assert(r.finished);
});
test('uses net displacement including final-step overshoot; time limit is not a finish',()=>{
  const r=raceResult({...state,finish_reason:'time_limit',frame:{...state.frame,distance_angstrom:123,progress:.615}});
  assert.equal(r.averageSpeed,123);assert(!r.finished);
  assert.equal(raceResult({...state,frame:{...state.frame,distance_angstrom:-10,progress:0}}).averageSpeed,-10);
  assert(Math.abs(raceResult({...state,frame:{...state.frame,distance_angstrom:200.5}}).averageSpeed-200.5)<1e-10);
});
test('running races and invalid times do not produce results',()=>{
  assert.equal(raceResult({...state,status:'running'}),null);
  for(const time_fs of [0,NaN,Infinity,-1])assert.equal(raceResult({...state,frame:{...state.frame,time_fs}}),null);
});
