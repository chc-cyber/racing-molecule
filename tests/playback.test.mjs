import {test} from 'node:test';
import assert from 'node:assert/strict';
import {PlaybackClock} from '../frontend/playback.mjs';

test('100 ps takes 50 seconds at 1x and 25 seconds at 2x',()=>{
  for (const rate of [1,2]) {
    const clock=new PlaybackClock();clock.setRate(rate);clock.playing=true;
    for(let t=0;t<=50000/rate;t+=10) clock.advance(t,100000);
    assert.equal(clock.timeFs,100000);
  }
});
test('switching rate preserves position, and reset preserves the chosen rate',()=>{
  const clock=new PlaybackClock();clock.playing=true;
  clock.advance(0,100000);clock.advance(100,100000);
  clock.setRate(2);assert.equal(clock.timeFs,200);
  clock.advance(100,100000);clock.advance(200,100000);assert.equal(clock.timeFs,600);
  clock.reset();assert.equal(clock.timeFs,0);assert.equal(clock.rate,2);assert.equal(clock.playing,false);
});
test('pause and hidden frames freeze presentation; underrun never accumulates catch-up',()=>{
  const clock=new PlaybackClock();clock.playing=true;
  clock.advance(0,50);clock.advance(100,50);assert.equal(clock.timeFs,50);assert(clock.buffering);
  clock.advance(10000,50);clock.advance(10010,1000);assert.equal(clock.timeFs,70);
  clock.playing=false;clock.advance(11000,1000);assert.equal(clock.timeFs,70);
  clock.playing=true;clock.advance(12000,1000,false);assert.equal(clock.timeFs,70);
  clock.advance(12010,1000);assert.equal(clock.timeFs,90);
});
