// Explicit job list prevents accidentally starting a large unbounded scan.
import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import zlib from 'node:zlib';
import {root,folder,cars,initialize,sources,checkSources,read,save,stepsFor,geometry,hash} from './acceptance_support.mjs';
const index=process.argv.indexOf('--jobs');assert(index>=0,'pass --jobs JSON');
const jobs=await read(process.argv[index+1]),gate=await read(path.join(folder,'nve.json'));
assert.equal(gate.status,'complete');await checkSources(gate.sources);
const provenance=await sources();provenance[path.join(root,'wasm/tools/acceptance_races.mjs')]=await hash(path.join(root,'wasm/tools/acceptance_races.mjs'));
provenance[path.join(folder,'nve.json')]=await hash(path.join(folder,'nve.json'));
const label=process.env.TIMESTEP_RUN_LABEL||'races';
const report={status:'running',sources:provenance,protocol:{sample_fs:50,coordinate_fs:500,damping_fs:350,velocity_seed:17329,track_nm:[20,20],same_initial_phase_space:'periodic preparation at 0.5 fs for every backend and race timestep',finish:'absolute center x >=200 A, first sampled crossing; no extrapolation',lift_guard:'all atoms above 8 A for >=5ps',timing:'Node serial MD excluding preparation; do not compare when concurrent CPU jobs are present'},jobs,rows:[]};
const target=path.join(folder,label+'.json');
try{
 for(const c of jobs){
  assert([.5,1,2].includes(c.dt_fs));assert(cars[c.car]);
  const screened=gate.rows.find(r=>r.car===c.car&&r.temperature_K===c.temperature_K&&r.lattice_scale===c.lattice_scale&&r.surface===c.surface&&r.dt_fs===c.dt_fs);
  assert(screened?.numerical_screen_pass||(c.sensitivity_only===true&&screened?.stability_guard_pass&&c.dt_fs===1),'NVE gate rejected '+JSON.stringify(c));
  const key=[c.car,c.surface,'T'+c.temperature_K,'L'+c.lattice_scale,'dt'+c.dt_fs,'seed'+c.seed].join('-');
  const trial=path.join(folder,label,key+'.json');
  try{const old=await read(trial);if(old.status==='complete'&&JSON.stringify(old.config)===JSON.stringify(c)&&JSON.stringify(old.sources)===JSON.stringify(provenance)){report.rows.push(old);continue;}}catch{}
  const row={key,config:c,sources:provenance,status:'running',samples:[]};let e,wasmModule;
  try{
   const initialized=await initialize(c.seed);e=initialized.engine;wasmModule=initialized.module;
   const {meta}=initialized,content=await fs.readFile(path.join(root,cars[c.car]),'utf8');
   const began=performance.now();const info=e.prepareWithPose(content,c.car+'.mol',{...meta.settings,surface:'periodic-au',temperature:c.temperature_K,damping_fs:350,lattice_scale:c.lattice_scale,timestep_fs:.5},{mode:'auto'});
   assert(info.runtime.preparation.accepted);row.preparation_seconds=(performance.now()-began)/1000;
   e.fourierEnabled(c.surface==='fourier-au');e.setTimestep(c.dt_fs);e.update({drive:{mode:'gradient',strength:c.strength??.6,growth:1}});
   let f=e.refresh();const measure=geometry(info,f.positions),initial=f.center.map(x=>x/10),positions=[],steps=[];
   row.initial_center_nm=initial;row.runtime=e.info().runtime;assert(row.runtime.surface.proxy_charts===3);row.initial_phase_space=e.phaseSpace();
   const iterations=Math.round((c.limit_ps??400)*20),batch=stepsFor(50)(c.dt_fs),started=performance.now();
   let low=1,high=1,lift=0,maxLift=0,outcome='time-limit';
   for(let k=0;k<=iterations;k++){
    if(k)f=e.advance(batch);assert.equal(f.time_fs,k*50);const g=measure(f);
    low=Math.min(low,g.bond_min);high=Math.max(high,g.bond_max);assert(low>=.65&&high<=1.45,'bond guard');
    if(k){lift=g.min_z_A>8?lift+1:0;maxLift=Math.max(maxLift,lift);}
    row.samples.push({step:f.step,time_ps:f.time_fs/1000,center_nm:f.center.map(x=>x/10),energy_kcal_mol:f.energy,internal_temperature_K:f.internal_temperature,...g});
    const finish=f.center[0]>=200;if(k%10===0||finish||k===iterations){positions.push(f.positions);steps.push(f.step);}
    if(finish){outcome='finish';break;}
   }
   const last=row.samples.at(-1);Object.assign(row,{status:'complete',outcome,finish_time_ps:outcome==='finish'?last.time_ps:null,finish_crossing_interval_ps:outcome==='finish'?[last.time_ps-.05,last.time_ps]:null,right_censored:outcome!=='finish',observed_time_ps:last.time_ps,bond_ratio_range:[low,high],maximum_lift_duration_ps:maxLift*.05,competitive_valid:maxLift<100,propagation_seconds:(performance.now()-started)/1000,maximum_lateral_excursion_nm:Math.max(...row.samples.map(s=>Math.abs(s.center_nm[1]-initial[1])))});
   const drift=row.samples.filter((_,k)=>k>=20).map((s,k)=>s.center_nm[0]-row.samples[k].center_nm[0]);
   row.low_progress_fraction=drift.filter(v=>v<.02).length/drift.length;
   row.lateral_rms_nm=Math.sqrt(row.samples.reduce((sum,s)=>sum+(s.center_nm[1]-initial[1])**2,0)/row.samples.length);
   await fs.mkdir(path.dirname(trial),{recursive:true});await fs.writeFile(trial+'.coords.gz',zlib.gzipSync(JSON.stringify({positions,steps,masses:info.molecule.atoms.map(a=>a.mass),dt_fs:c.dt_fs})));
  }catch(err){row.status='failed';try{row.error=wasmModule.getExceptionMessage(err);wasmModule.decrementExceptionRefcount(err);}catch{row.error=String(err.stack||err);}}finally{if(e){e.close();e.delete();}}
  await save(trial,row);report.rows.push(row);await save(target,report);console.log(JSON.stringify({key,status:row.status,finish_ps:row.finish_time_ps,error:row.error}));
 }
 await checkSources(provenance);report.status='complete';
}catch(err){report.status='failed';report.error=String(err.stack||err);process.exitCode=1;}
finally{await save(target,report);}
