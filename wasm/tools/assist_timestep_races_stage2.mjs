// Resume the helper queue after the sandbox rejected hardlink publication.
// Original MD/source hashes remain untouched; complete helper trials are reused.
import fs from 'node:fs/promises';
import path from 'node:path';
import {spawn} from 'node:child_process';
import assert from 'node:assert/strict';
import {root,folder,read,save,hash} from './timestep_support.mjs';
const runner='wasm/tools/assist_timestep_races_stage2.mjs';
const script=path.join(root,'wasm/tools/scan_timestep_races.mjs'),owner='races-2';
const jobs=await read(path.join(folder,'race-jobs-2.json'));
const key=c=>[c.car,c.surface,'T'+c.temperature_K,'L'+c.lattice_scale,'dt'+c.dt_fs,'seed'+c.seed].join('-');
const trial=(label,c)=>path.join(folder,label,key(c)+'.json');
const exists=async p=>{try{await fs.access(p);return true;}catch{return false;}};
const guard=trial(owner,jobs[4]),buffer=trial(owner,jobs[7]),sources={};
for(const f of [runner,'wasm/tools/assist_timestep_races.mjs','wasm/tools/scan_timestep_races.mjs','outputs/timestep-study/race-jobs-2.json','outputs/timestep-study/nve.json'])sources[f]=await hash(path.join(root,f));
const report={status:'running',sources,recovery:{original_report:'helper-runner.json',error:'Sandbox denies NTFS hardlink creation (WinError 5)',repair:'Coordinates first and JSON last, each via a temporary file and atomic rename; reuses completed helper trials without recomputation'},protocol:{maximum_helpers:2,maximum_total_md_processes:4,owner_label:owner,guard_index:4,unclaimed_buffer_index:7,assisted_indices:'8 and later only',poll_ms:25,stop:'terminate active helper as soon as guard result exists',wall_times_are_not_benchmarks:true},rows:[]};
let cursor=8,writeQueue=Promise.resolve();
const checkpoint=()=>{writeQueue=writeQueue.then(()=>save(path.join(folder,'helper-runner-stage2.json'),report));return writeQueue;};
async function publish(source,c,index,label){
 const row=await read(source);assert.equal(row.status,'complete');assert.deepEqual(row.config,c);
 for(const [file,h] of Object.entries(row.sources))assert.equal(await hash(file),h,'changed helper source '+file);
 const destination=trial(owner,c);
 if(await exists(destination)){
  const existing=await read(destination);assert.deepEqual(existing.config,row.config);assert.deepEqual(existing.sources,row.sources);assert.deepEqual(existing.samples,row.samples);return {status:'existing-verified'};
 }
 assert(!(await exists(buffer)),'owner buffer completed before cache publication');
 await fs.mkdir(path.dirname(destination),{recursive:true});
 const coordSource=source+'.coords.gz',coordTarget=destination+'.coords.gz',coordTemp=coordTarget+'.safe-tmp';
 if(await exists(coordTarget))assert.equal(await hash(coordTarget),await hash(coordSource));
 else{await fs.copyFile(coordSource,coordTemp);await fs.rename(coordTemp,coordTarget);}
 row.assisted_by={runner,runner_sha256:sources[runner],helper_label:label,owner_queue_index:index,guard_index:4,unclaimed_buffer_index:7,publication:'coords then JSON, temporary file + atomic rename',recovered_from_hardlink_error:false};
 const temporary=destination+'.safe-tmp';await save(temporary,row);
 assert(!(await exists(destination)),'unexpected competing publisher');await fs.rename(temporary,destination);
 return {status:'cached',finish_time_ps:row.finish_time_ps,outcome:row.outcome,source_trial:path.relative(root,source),row_sha256:await hash(destination),coords_sha256:await hash(coordTarget)};
}
// Recover both already-running original helper jobs before starting new MD.
for(const [slot,index] of []){
 const c=jobs[index],label=`assist-${slot}-${index}`,source=trial(label,c);
 if(await exists(source)){
  const r={index,slot,key:key(c),recovered_existing_trial:true};report.rows.push(r);
  try{Object.assign(r,await publish(source,c,index,label));}catch(err){r.status='recovery-failed';r.error=String(err.stack||err);}
 }
}
await checkpoint();
async function worker(slot){
 while(true){
  if(await exists(guard))return;
  const index=cursor++;if(index>=jobs.length)return;const c=jobs[index];if(await exists(trial(owner,c)))continue;
  const label=`assist-stage2-${slot}-${index}`,jobFile=path.join(folder,label+'-job.json'),source=trial(label,c);
  await save(jobFile,[c]);const r={index,slot,key:key(c),status:'running',started:new Date().toISOString()};report.rows.push(r);await checkpoint();
  if(await exists(guard)){r.status='cancelled-before-start';await checkpoint();return;}
  const child=spawn(process.execPath,[script,'--jobs',jobFile],{cwd:root,env:{...process.env,TIMESTEP_RUN_LABEL:label},stdio:['ignore','pipe','pipe']});
  let log='',cancelled=false,checking=false;child.stdout.on('data',v=>{log+=v;});child.stderr.on('data',v=>{log+=v;});
  const timer=setInterval(async()=>{if(checking)return;checking=true;try{if(await exists(guard)){cancelled=true;child.kill();}}finally{checking=false;}},25);
  const code=await new Promise((resolve,reject)=>{child.on('error',reject);child.on('close',resolve);});clearInterval(timer);await fs.writeFile(path.join(folder,label+'.log'),log);
  let row;try{row=await read(source);}catch{}
  if(row?.status==='complete')Object.assign(r,await publish(source,c,index,label));
  else if(cancelled)r.status='cancelled-at-guard';
  else{r.status='failed';r.exit_code=code;r.error=row?.error||log;}
  r.ended=new Date().toISOString();await checkpoint();console.log(JSON.stringify(r));if(cancelled||await exists(guard))return;
 }
}
try{await Promise.all([worker(1),worker(2)]);for(const [f,h] of Object.entries(sources))assert.equal(await hash(path.join(root,f)),h);report.status='complete';}
catch(err){report.status='failed';report.error=String(err.stack||err);process.exitCode=1;}
finally{await checkpoint();await writeQueue;}

