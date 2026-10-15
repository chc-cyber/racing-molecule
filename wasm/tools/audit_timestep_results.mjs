// Read-only integrity audit: no simulation or assumptions of statistical equivalence.
import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import assert from 'node:assert/strict';
import {root,folder,read,save,hash,checkSources} from './timestep_support.mjs';
const baseJobs=await read(path.join(folder,'race-jobs-all.json'));
const supplementJobs=await read(path.join(folder,'race-jobs-supplement.json'));
const jobs=[...baseJobs,...supplementJobs];
const key=c=>[c.car,c.surface,'T'+c.temperature_K,'L'+c.lattice_scale,'dt'+c.dt_fs,'seed'+c.seed].join('-');
const partial=process.env.TIMESTEP_AUDIT_PARTIAL==='1';
const rows=new Map(),files=new Map();
for(const label of ['races-1','races-2','races-supplement']){
 const report=await read(path.join(folder,label+'.json'));if(!partial)assert.equal(report.status,'complete');await checkSources(report.sources);
 for(const row of report.rows){assert(!rows.has(row.key),'duplicate job');rows.set(row.key,row);files.set(row.key,path.join(folder,label,row.key+'.json'));}
}
if(partial){for(const label of ['races-1','races-2','races-supplement'])for(const name of await fs.readdir(path.join(folder,label))){if(!name.endsWith('.json'))continue;const file=path.join(folder,label,name),row=await read(file);if(!rows.has(row.key)){rows.set(row.key,row);files.set(row.key,file);}}}
if(!partial)assert.equal(rows.size,jobs.length);for(const row of rows.values())assert(jobs.some(c=>key(c)===row.key),'unplanned job');const results=[],initialCoordinates=new Map();let snapshots=0,maxCOMError=0,maxHeightError=0;
for(const c of jobs){
 const r=rows.get(key(c));if(!r&&partial)continue;assert(r,'missing '+key(c));assert.deepEqual(r.config,c);
 const raw=await read(files.get(r.key));assert.deepEqual(raw,r,'owner aggregate differs from trial cache');
 if(r.status!=='complete'){results.push({key:r.key,status:r.status,error:r.error});continue;}
 assert.equal(r.runtime.timestep_fs,c.dt_fs);assert(r.runtime.preparation.accepted);assert(r.runtime.preparation.global_force_norm_kcal_mol_A<=1e-4);
 for(let i=0;i<r.samples.length;i++){
  const s=r.samples[i];assert.equal(s.step,i*50/c.dt_fs);assert(Math.abs(s.time_ps-i*.05)<1e-10);assert(s.center_nm.every(Number.isFinite));assert(Number.isFinite(s.energy_kcal_mol));
  if(i<r.samples.length-1)assert(s.center_nm[0]<20,'ignored early finish');
 }
 const last=r.samples.at(-1);assert.equal(r.observed_time_ps,last.time_ps);
 assert.equal(r.right_censored,r.outcome!=='finish');
 if(r.outcome==='finish'){
  assert(last.center_nm[0]>=20);assert.equal(r.finish_time_ps,last.time_ps);assert.deepEqual(r.finish_crossing_interval_ps,[last.time_ps-.05,last.time_ps]);
 }else{assert.equal(r.outcome,'time-limit');assert.equal(last.time_ps,c.limit_ps);assert.equal(r.finish_time_ps,null);}
 let lift=0,maxLift=0;for(const s of r.samples.slice(1)){lift=s.min_z_A>8?lift+1:0;maxLift=Math.max(maxLift,lift);}
 assert.equal(r.maximum_lift_duration_ps,maxLift*.05);assert.equal(r.competitive_valid,maxLift<100);assert(r.bond_ratio_range[0]>=.65&&r.bond_ratio_range[1]<=1.45);
 const coordinates=JSON.parse(zlib.gunzipSync(await fs.readFile(files.get(r.key)+'.coords.gz')));
 const initialKey=[c.car,c.temperature_K,c.lattice_scale].join('/');
 if(initialCoordinates.has(initialKey))assert.deepEqual(coordinates.positions[0],initialCoordinates.get(initialKey),'different prepared phase coordinates');
 else initialCoordinates.set(initialKey,coordinates.positions[0]);
 assert.equal(coordinates.dt_fs,c.dt_fs);assert.equal(coordinates.steps.length,coordinates.positions.length);
 const mass=coordinates.masses.reduce((a,b)=>a+b,0);
 for(let j=0;j<coordinates.steps.length;j++){
  const s=r.samples[coordinates.steps[j]*c.dt_fs/50],p=coordinates.positions[j];assert(s);assert.equal(p.length,3*coordinates.masses.length);
  for(let axis=0;axis<3;axis++){
   const com=coordinates.masses.reduce((sum,m,i)=>sum+m*p[3*i+axis],0)/mass/10;
   maxCOMError=Math.max(maxCOMError,Math.abs(com-s.center_nm[axis]));
  }
  const z=p.filter((_,i)=>i%3===2);maxHeightError=Math.max(maxHeightError,Math.abs(Math.min(...z)-s.min_z_A),Math.abs(Math.max(...z)-s.max_z_A));snapshots++;
 }
 if(r.assisted_by){assert.equal(await hash(path.join(root,r.assisted_by.runner)),r.assisted_by.runner_sha256);const [guard,buffer]=[r.assisted_by.guard_index,r.assisted_by.unclaimed_buffer_index];assert((guard===3&&buffer===4)||(guard===4&&buffer===7));assert(r.assisted_by.owner_queue_index>buffer);}
 results.push({key:r.key,status:'verified',outcome:r.outcome,valid_finish:r.outcome==='finish'&&r.competitive_valid});
}
assert(maxCOMError<1e-10);assert(maxHeightError<1e-10);
const report={status:partial?'partial-verified':'passed',qualification:'Integrity/protocol verification only, not a claim of step-size equivalence',planned:jobs.length,base_planned:baseJobs.length,supplement_planned:supplementJobs.length,verified:results.filter(r=>r.status==='verified').length,valid_finishes:results.filter(r=>r.valid_finish).length,snapshots,matched_initial_coordinate_groups:initialCoordinates.size,maxCOMError_nm:maxCOMError,maxHeightError_A:maxHeightError,results};
await save(path.join(folder,partial?'integrity-partial.json':'integrity.json'),report);console.log(JSON.stringify({...report,results:undefined}));
