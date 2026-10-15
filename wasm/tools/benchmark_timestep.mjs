// Run serially while other CPU-heavy validations are idle.
import fs from 'node:fs/promises';
import path from 'node:path';
import os from 'node:os';
import assert from 'node:assert/strict';
import {root,folder,cars,initialize,sources,checkSources,read,save,stepsFor,hash} from './timestep_support.mjs';
const {engine:e,meta}=await initialize(),gate=await read(path.join(folder,'nve.json'));
assert.equal(gate.status,'complete');await checkSources(gate.sources);
const report={status:'running',sources:await sources(),hardware:{cpu:os.cpus()[0].model,node:process.version},protocol:{duration_ps:2,physical_sample_fs:50,repeats:3,temperature_K:150,lattice_scales:[1.5,2],damping_fs:350,strength:.6,timing:'serial MD wall time, initialization/relaxation excluded; no concurrent benchmark workloads'},rows:[]};
report.sources[path.join(root,'wasm/tools/benchmark_timestep.mjs')]=await hash(path.join(root,'wasm/tools/benchmark_timestep.mjs'));
try{
 for(const car of Object.keys(cars))for(const scale of [1.5,2]){
  const content=await fs.readFile(path.join(root,cars[car]),'utf8');
  e.prepareWithPose(content,car+'.mol',{...meta.settings,surface:'periodic-au',temperature:150,damping_fs:350,lattice_scale:scale,timestep_fs:.5},{mode:'auto'});
  const initial=e.phaseSpace();
  for(const surface of ['periodic-au','fourier-au'])for(let repeat=0;repeat<3;repeat++)for(const dt of repeat%2?[1,.5]:[.5,1]){
   const screened=gate.rows.find(r=>r.car===car&&r.temperature_K===150&&r.lattice_scale===scale&&r.surface===surface&&r.dt_fs===dt);
   if(!screened?.stability_guard_pass)continue;
   e.fourierEnabled(surface==='fourier-au');e.setPhaseSpace(initial.positions,initial.velocities);e.setTimestep(dt);
   e.update({temperature:150,damping_fs:350,drive:{mode:'gradient',strength:.6,growth:1}});e.refresh();
   const started=performance.now();let frame;for(let k=0;k<40;k++)frame=e.advance(stepsFor(50)(dt));
   const seconds=(performance.now()-started)/1000;assert.equal(frame.time_fs,2000);assert(frame.positions.every(Number.isFinite));
   report.rows.push({car,lattice_scale:scale,surface,dt_fs:dt,repeat,seconds,seconds_per_100ps:seconds*50,precision_target_pass:screened.numerical_screen_pass});
  }
  await save(path.join(folder,'benchmark.json'),report);console.log(JSON.stringify({car,scale}));
 }
 await checkSources(report.sources);report.status='complete';
}catch(err){report.status='failed';report.error=String(err.stack||err);process.exitCode=1;}
finally{e.close();e.delete();await save(path.join(folder,'benchmark.json'),report);}
