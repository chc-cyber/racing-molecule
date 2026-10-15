// Real rebuilt WASM: race origin survives preparation and absolute-line crossings.
import fs from 'node:fs/promises';
import zlib from 'node:zlib';
import assert from 'node:assert/strict';
import createMolecularCar from '../dist/browser/runtime/molecular-car.mjs';
import {raceCourse, raceFinished} from '../../frontend/race-course.mjs';

const assets = new URL('../dist/browser/assets/', import.meta.url);
const metadata = JSON.parse(await fs.readFile(new URL('model.json', assets), 'utf8'));
async function load(name) {
  const data = zlib.gunzipSync(await fs.readFile(new URL(name, assets)));
  return new Float64Array(data.buffer.slice(data.byteOffset, data.byteOffset+data.byteLength));
}
const module = await createMolecularCar({print:()=>{}, printErr:()=>{}});
const engine = new module.BrowserEngine(metadata, ...await Promise.all(
  ['tile.bin.gz','saip-grid.bin.gz','hetero-grid.bin.gz'].map(load)));
const info = JSON.parse(await fs.readFile(new URL('fourier-atlas.json', assets), 'utf8'));
engine.installFourier(info);
await engine.loadFourierSpecies(0, await load(info.arrays[0].path.replace('assets/','')));
const content = await fs.readFile(new URL('../../examples/benzene.mol', import.meta.url), 'utf8');
const prepared = engine.prepare(content, 'benzene.mol', {...metadata.settings, timestep_fs:1});
engine.update({drive:metadata.drive});
const initial = engine.frame();
assert.equal(initial.distance_angstrom, 0);
assert.equal(initial.progress, 0);
assert.equal(initial.race_start_x, initial.center[0]);
assert.equal(prepared.runtime.race.finish_x, initial.center[0]+200);
const report = {initial, absolute_crossings:[]};
const phase = engine.phaseSpace();
// This existing diagnostic resets step zero and establishes a new race origin.
// Start close to the old fixed line: crossing x=200 must NOT complete 20 nm.
const translated = phase.positions.map((v,i)=>i%3===0?v-initial.center[0]+199.9:v);
const velocity = phase.velocities.map((_,i)=>i%3===0?.1:0);
engine.setPhaseSpace(translated,velocity);
const crossed = engine.advance(5);
assert(crossed.center[0]>200);
assert(crossed.distance_angstrom<2);
assert(!raceFinished(crossed));
assert.equal(raceCourse(metadata.track,crossed).finish_x,crossed.race_start_x+200);
report.absolute_crossings.push(crossed);
// Restore the prepared phase space and run one actual full race at 50 fs samples.
engine.setPhaseSpace(phase.positions,phase.velocities);
let previous=engine.frame(),frame=previous;
const started=performance.now();
while(!raceFinished(frame)&&frame.time_fs<250000) {
  previous=frame;frame=engine.advance(50);
  assert(Math.abs(frame.distance_angstrom-(frame.center[0]-frame.race_start_x))<1e-10);
}
assert(raceFinished(frame));assert(previous.distance_angstrom<200);assert(frame.distance_angstrom>=200);
report.finish={previous,frame,wall_seconds:(performance.now()-started)/1000};
engine.close();
const output=new URL('../../outputs/race-distance-check/',import.meta.url);
await fs.mkdir(output,{recursive:true});
await fs.writeFile(new URL('wasm.json',output),JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify({passed:true,origin_A:frame.race_start_x,previous_distance_nm:previous.distance_angstrom/10,
  finish_distance_nm:frame.distance_angstrom/10,time_ps:frame.time_fs/1000,wall_seconds:report.finish.wall_seconds}));
