import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import {root,folder,cars,initialize,sources,checkSources,save} from './timestep_support.mjs';
const {engine:e,meta}=await initialize();const provenance=await sources(),rows=[];let preparedReference;
const content=await fs.readFile(path.join(root,cars['benzene-12']),'utf8');
try{
 for(const dt of [.25,.5,1,2]){
  const info=e.prepareWithPose(content,'benzene.mol',{...meta.settings,temperature:150,timestep_fs:dt},{mode:'auto'});
  assert.equal(info.runtime.timestep_fs,dt);const initial=e.phaseSpace();
  if(preparedReference)assert.deepEqual(initial,preparedReference,'changing race timestep must not alter preparation/velocities');else preparedReference=initial;
  e.pulse();
  assert.equal(e.frame().pulse_until_step,400/dt);assert.equal(e.frame().pulse_active,true);
  let f;for(let i=0;i<8;i++)f=e.advance(50/dt);
  assert.equal(f.time_fs,400);assert.equal(f.pulse_active,false);assert.equal(f.timestep_fs,dt);
  assert.throws(()=>e.setTimestep(.5));e.setPhaseSpace(initial.positions,initial.velocities);e.setTimestep(dt);assert.equal(e.frame().time_fs,0);
  assert.throws(()=>e.setTimestep(0));assert.throws(()=>e.setTimestep(NaN));assert.throws(()=>e.setTimestep(.75));rows.push({dt_fs:dt,passed:true});
 }
 const info=e.prepare(content,'benzene.mol',meta.settings);assert.equal(info.runtime.timestep_fs,.5);await checkSources(provenance);
 await save(path.join(folder,'interface.json'),{status:'passed',sources:provenance,rows,default_timestep_fs:.5,pulse_duration_fs:400,identical_prepared_phase_space:true});
}finally{e.close();e.delete();}
