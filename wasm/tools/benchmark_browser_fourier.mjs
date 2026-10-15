import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import os from 'node:os';
import {root,out,read,initialize,sources,checkSources} from './fourier_support.mjs';
const {engine:e,meta}=await initialize(),gate=await read(path.join(out,'numerics.json'));
assert.equal(gate.status,'passed');await checkSources(gate.sources);
const fixtures=(await read(path.join(out,'native-fixtures.json'))).filter(f=>f.label==='start'&&f.scale===2);
const report={status:'running',sources:await sources(),hardware:{cpu:os.cpus()[0].model,node:process.version},
  protocol:{dt_fs:.5,sample_ps:2,repeats:3,serial:true,temperature:150,strength:.6,initial_state:'same relaxed periodic-surface phase space',timing:'100 ps costs extrapolate 2 ps serial timing; excludes initialization and UI'},md:[],nve:[]};
const save=()=>fs.writeFile(path.join(out,'benchmark.json'),JSON.stringify(report,null,2)+'\n');
try {
  for(const f of fixtures)for(const scale of [1.5,2]) {
    e.prepareWithPose(f.content,f.car+'.mol',{...meta.settings,temperature:150,lattice_scale:scale},{mode:'auto'});
    e.update({drive:{mode:'gradient',strength:.6,growth:1}});const initial=e.phaseSpace();
    for(let repeat=0;repeat<3;repeat++)for(const proxy of repeat%2?[true,false]:[false,true]) {
      e.fourierEnabled(proxy);e.setPhaseSpace(initial.positions,initial.velocities);e.update({temperature:150});
      const counts=[0,0,0],started=performance.now();let frame;
      for(let k=0;k<40;k++){frame=e.advance(100);const c=e.proxyCounts();for(let j=0;j<3;j++)counts[j]+=c[j];}
      const seconds=(performance.now()-started)/1000;assert(frame.positions.every(Number.isFinite));
      report.md.push({car:f.car,scale,repeat,proxy,seconds,seconds_per_100ps:seconds*50,counts});
    }
    // Equal initial positions/velocities; measure each Hamiltonian's own drift.
    for(const proxy of [false,true]) {
      e.fourierEnabled(proxy);e.setPhaseSpace(initial.positions,initial.velocities);e.conservative(true);
      const initialFrame=e.refresh();let drift=0;
      for(let k=0;k<40;k++){const frame=e.advance(100);drift=Math.max(drift,Math.abs(frame.energy-initialFrame.energy));}
      report.nve.push({car:f.car,scale,proxy,max_total_energy_drift:drift});
      e.conservative(false);
    }
    console.log(JSON.stringify({car:f.car,scale}));await save();
  }
  await checkSources(report.sources);report.status='passed';
} catch(error){report.status='failed';report.error=error.stack||String(error);process.exitCode=1;}
finally{e.close();e.delete();await save();}
console.log(JSON.stringify({status:report.status,error:report.error}));
