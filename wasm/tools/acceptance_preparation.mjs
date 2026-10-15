// Check that the paired periodic preparation is representative of actual v2 preparation.
import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import {root,folder,cars,initialize,sources,checkSources,save,hash} from './acceptance_support.mjs';
const report={status:'running',sources:await sources(),rows:[],scope:'Actual v2 preparation versus paired periodic preparation at 150 K / scale 1.5; positions, not a full additional race'};
report.sources[path.join(root,'wasm/tools/acceptance_preparation.mjs')]=await hash(path.join(root,'wasm/tools/acceptance_preparation.mjs'));
const {engine:e,meta}=await initialize();
try {
  for(const car of Object.keys(cars)) {
    const content=await fs.readFile(path.join(root,cars[car]),'utf8'),states=[],preparations=[];
    for(const surface of ['periodic-au','fourier-au']) {
      const info=e.prepareWithPose(content,car+'.mol',{...meta.settings,surface,temperature:150,damping_fs:350,lattice_scale:1.5,timestep_fs:1},{mode:'auto'});
      assert(info.runtime.preparation.accepted);states.push(e.phaseSpace());preparations.push(info.runtime.preparation);
    }
    const dx=states[0].positions.map((x,i)=>x-states[1].positions[i]);
    report.rows.push({car,preparations,maximum_component_difference_A:Math.max(...dx.map(Math.abs)),
      atom_rms_displacement_A:Math.sqrt(dx.reduce((s,x)=>s+x*x,0)/(dx.length/3))});
  }
  await checkSources(report.sources);report.status='complete';
} catch(error) {report.status='failed';report.error=String(error.stack||error);process.exitCode=1;}
finally {e.close();e.delete();await save(path.join(folder,'preparation.json'),report);}
console.log(JSON.stringify(report.rows));
