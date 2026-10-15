// Assist only distant jobs of the unchanged legacy queue. Two unclaimed heavy
// jobs separate its frontier from assisted jobs; stop helpers at the FIRST guard.
import fs from 'node:fs/promises';
import path from 'node:path';
import {spawn} from 'node:child_process';
import assert from 'node:assert/strict';
import {root,folder,read,save,hash} from './timestep_support.mjs';
const script=path.join(root,'wasm/tools/scan_timestep_races.mjs');
const owner='races-2',jobs=await read(path.join(folder,'race-jobs-2.json'));
const key=c=>[c.car,c.surface,'T'+c.temperature_K,'L'+c.lattice_scale,'dt'+c.dt_fs,'seed'+c.seed].join('-');
const trial=(label,c)=>path.join(folder,label,key(c)+'.json');
const exists=async p=>{try{await fs.access(p);return true;}catch{return false;}};
assert.equal(jobs[3].car,'folded-eicosane');assert.equal(jobs[4].car,'folded-eicosane');
assert.equal(jobs[3].dt_fs,.5);assert.equal(jobs[4].dt_fs,.5);
const guard=trial(owner,jobs[3]),sources={};
for(const f of ['wasm/tools/assist_timestep_races.mjs','wasm/tools/scan_timestep_races.mjs','outputs/timestep-study/race-jobs-2.json','outputs/timestep-study/nve.json'])sources[f]=await hash(path.join(root,f));
const report={status:'running',sources,protocol:{maximum_helpers:2,maximum_total_md_processes:4,owner_label:owner,guard_index:3,unclaimed_buffer_index:4,assisted_indices:'5 and later only',poll_ms:25,stop:'terminate active helper as soon as guard result exists; never assist either guard/buffer job',cache:'identical config/source hashes; complete row and coordinates; atomic no-overwrite hardlink publish; owner report untouched',wall_times_are_not_benchmarks:true},rows:[]};
let cursor=5,writeQueue=Promise.resolve();
const checkpoint=()=>{writeQueue=writeQueue.then(()=>save(path.join(folder,'helper-runner.json'),report));return writeQueue;};
async function worker(slot){
 while(true){
  if(await exists(guard))return;
  const index=cursor++;if(index>=jobs.length)return;
  const c=jobs[index],destination=trial(owner,c);
  if(await exists(destination))continue;
  const label=`assist-${slot}-${index}`,jobFile=path.join(folder,`${label}-job.json`),source=trial(label,c);
  await save(jobFile,[c]);
  const r={index,slot,key:key(c),status:'running',started:new Date().toISOString(),guard_absent_at_start:true};report.rows.push(r);await checkpoint();
  if(await exists(guard)){r.status='cancelled-before-start';await checkpoint();return;}
  const child=spawn(process.execPath,[script,'--jobs',jobFile],{cwd:root,env:{...process.env,TIMESTEP_RUN_LABEL:label},stdio:['ignore','pipe','pipe']});
  let log='',cancelled=false,checking=false;
  child.stdout.on('data',v=>{log+=v;});child.stderr.on('data',v=>{log+=v;});
  const timer=setInterval(async()=>{if(checking)return;checking=true;try{if(await exists(guard)){cancelled=true;child.kill();}}finally{checking=false;}},25);
  const code=await new Promise((resolve,reject)=>{child.on('error',reject);child.on('close',resolve);});clearInterval(timer);
  await fs.writeFile(path.join(folder,label+'.log'),log);
  let row;try{row=await read(source);}catch{}
  if(row?.status==='complete'){
   assert.deepEqual(row.config,c);for(const [file,h] of Object.entries(row.sources))assert.equal(await hash(file),h,'changed helper source '+file);
   assert(!(await exists(destination)),'owner reached an assisted job before guarded publication');
   row.assisted_by={runner:'wasm/tools/assist_timestep_races.mjs',runner_sha256:sources['wasm/tools/assist_timestep_races.mjs'],helper_label:label,owner_queue_index:index,guard_index:3,unclaimed_buffer_index:4};
   await fs.mkdir(path.dirname(destination),{recursive:true});
   await fs.copyFile(source+'.coords.gz',destination+'.coords.gz',1);
   const temporary=destination+'.helper-tmp';await save(temporary,row);
   await fs.link(temporary,destination);await fs.unlink(temporary);
   r.status='cached';r.finish_time_ps=row.finish_time_ps;r.outcome=row.outcome;r.source_trial=path.relative(root,source);r.row_sha256=await hash(destination);r.coords_sha256=await hash(destination+'.coords.gz');
  }else if(cancelled){r.status='cancelled-at-guard';}
  else{r.status='failed';r.exit_code=code;r.error=row?.error||log;}
  r.ended=new Date().toISOString();await checkpoint();console.log(JSON.stringify(r));
  if(cancelled||await exists(guard))return;
 }
}
try{await Promise.all([worker(1),worker(2)]);for(const [f,h] of Object.entries(sources))assert.equal(await hash(path.join(root,f)),h);report.status='complete';}
catch(err){report.status='failed';report.error=String(err.stack||err);process.exitCode=1;}
finally{await checkpoint();await writeQueue;}
