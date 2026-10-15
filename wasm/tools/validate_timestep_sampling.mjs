// Check whether a 10 fs diagnostic output interval hides energy-error peaks.
import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import {root,folder,cars,initialize,sources,checkSources,save,hash} from './timestep_support.mjs';
const {engine:e,meta}=await initialize(),provenance=await sources(),rows=[];
provenance[path.join(root,'wasm/tools/validate_timestep_sampling.mjs')]=await hash(path.join(root,'wasm/tools/validate_timestep_sampling.mjs'));
try{
 for(const [car,temperature] of [['fullerene-c60',150],['pyridine-11',600]]){
  const content=await fs.readFile(path.join(root,cars[car]),'utf8');
  e.prepareWithPose(content,car+'.mol',{...meta.settings,surface:'periodic-au',temperature,lattice_scale:2,timestep_fs:.5},{mode:'auto'});
  const initial=e.phaseSpace();
  for(const surface of ['periodic-au','fourier-au'])for(const dt of [.5,1,2]){
   e.fourierEnabled(surface==='fourier-au');e.setPhaseSpace(initial.positions,initial.velocities);e.setTimestep(dt);e.conservative(true);
   const first=e.refresh(),kinetic=first.energy-e.energyForces(first.positions,1,0,false).energy;
   let everyStep=0,every10=0;
   for(let step=1;step<=500/dt;step++){
    const f=e.advance(1);assert.equal(f.time_fs,step*dt);const error=Math.abs(f.energy-first.energy);assert(Number.isFinite(error));
    everyStep=Math.max(everyStep,error);if(step%(10/dt)===0)every10=Math.max(every10,error);
   }
   rows.push({car,temperature_K:temperature,lattice_scale:2,surface,dt_fs:dt,duration_ps:.5,initial_kinetic:kinetic,max_energy_error_every_step:everyStep,max_energy_error_every_10fs:every10,relative_every_step:everyStep/kinetic,relative_every_10fs:every10/kinetic,underestimate_pct:100*(1-every10/everyStep),precision_target_pass:everyStep/kinetic<.01});
   console.log(JSON.stringify(rows.at(-1)));
  }
 }
 await checkSources(provenance);await save(path.join(folder,'sampling.json'),{status:'passed',sources:provenance,rows,qualification:'500 fs every-integration-step observation on same trajectory; checks peak sampling, does not establish long-time stability'});
}finally{e.close();e.delete();}
