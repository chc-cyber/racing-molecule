// Revalidate archived scale=1 race results against the current WASM module.
// This is a bounded compatibility check, not a new race or parameter scan.
import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {fileURLToPath,pathToFileURL} from 'node:url';

const root=fileURLToPath(new URL('../../',import.meta.url));
const folder=path.join(root,'outputs/crystal-contrast-guarded');
const archived=path.join(folder,'development-v2/source');
const assets=path.join(root,'wasm/dist/browser/assets');
const runtime='wasm/dist/browser/runtime/';
const cars=['benzene-12','coronene-36','ethynyl-rigid','test-car-52','fullerene-c60','pyridine-11','benzonitrile-13'];
const files=['wasm/src/potentials.cpp','wasm/src/engine.cpp','wasm/tools/scan_crystal_races.mjs',runtime+'molecular-car.wasm'];
const hash=async file=>crypto.createHash('sha256').update(await fs.readFile(file)).digest('hex');
const difference=(a,b)=>{
  assert.equal(a.length,b.length);
  let result=0;
  for(let i=0;i<a.length;i++) {
    assert(Number.isFinite(a[i])&&Number.isFinite(b[i]));
    result=Math.max(result,Math.abs(a[i]-b[i]));
  }
  return result;
};
const report={status:'running',scope:'Archived scale=1 baseline compatibility only; no full races rerun',
  protocol:{temperature_K:150,strength:.6,lattice_scale:1,growth:[1,.25],
    langevin_seed:2718281,velocity_seed:17329,damping_fs:350,timestep_fs:.5,
    steps:2000,sample_every_steps:100,auto_start_pose:true,warmup_ps:0},
  tolerances:{position_A:1e-9,velocity_A_fs:1e-11,energy_kcal_mol:1e-7,force_kcal_mol_A:1e-7},
  archived_source_sha256:{},source_sha256:{},baseline_trials_sha256:{},cases:[]};
let oldModule,currentModule;
try {
  const names=(await fs.readdir(path.join(folder,'browser-trials'))).filter(x=>x.endsWith('.json'));
  const baseline=[];
  for(const name of names) {
    const file=path.join(folder,'browser-trials',name),row=JSON.parse(await fs.readFile(file,'utf8'));
    if(row.config.lattice_scale!==1)continue;
    baseline.push(row);
    report.baseline_trials_sha256['browser-trials/'+name]=await hash(file);
  }
  assert.equal(baseline.length,42,'Expected all 42 archived baseline races');
  for(const file of files) {
    const sha=await hash(path.join(archived,file));
    report.archived_source_sha256[file]=sha;
    for(const row of baseline)assert.equal(row.source_sha256[file],sha,`${row.key}: archived ${file}`);
  }
  for(const file of [...files,runtime+'molecular-car.mjs','wasm/dist/browser/assets/model.json',
    ...['tile','saip-grid','hetero-grid'].map(x=>`wasm/dist/browser/assets/${x}.bin.gz`),
    'wasm/tools/verify_crystal_baseline.mjs'])report.source_sha256[file]=await hash(path.join(root,file));
  for(const car of cars) {
    const file=`inputs/${car}.mol`,sha=await hash(path.join(folder,file));
    report.source_sha256[`outputs/crystal-contrast-guarded/${file}`]=sha;
    for(const row of baseline)assert.equal(row.source_sha256[file],sha,`Input changed: ${file}`);
  }
  for(const row of baseline)assert.equal(row.source_sha256['wasm/dist/browser/assets/model.json'],
    report.source_sha256['wasm/dist/browser/assets/model.json']);
  const metadata=JSON.parse(await fs.readFile(path.join(assets,'model.json'),'utf8'));
  const arrays=await Promise.all(['tile','saip-grid','hetero-grid'].map(async name=>{
    const data=zlib.gunzipSync(await fs.readFile(path.join(assets,`${name}.bin.gz`)));
    return new Float64Array(data.buffer.slice(data.byteOffset,data.byteOffset+data.byteLength));
  }));
  const {default:create}=await import(pathToFileURL(path.join(root,runtime+'molecular-car.mjs')));
  oldModule=await create({wasmBinary:await fs.readFile(path.join(archived,runtime+'molecular-car.wasm')),print:()=>{},printErr:()=>{}});
  currentModule=await create({wasmBinary:await fs.readFile(path.join(root,runtime+'molecular-car.wasm')),print:()=>{},printErr:()=>{}});
  for(const car of cars)for(const growth of report.protocol.growth) {
    const old=new oldModule.BrowserEngine({...metadata,seed:2718281},...arrays);
    const current=new currentModule.BrowserEngine({...metadata,seed:2718281},...arrays);
    const result={molecule:car,growth,sampled_frames:0,preparation_position_error_A:0,
      initial_velocity_error_A_fs:0,preparation_force_norm_error:0,static_energy_error:0,
      static_force_error:0,trajectory_position_error_A:0,trajectory_velocity_error_A_fs:0,
      trajectory_energy_error:0,saved_baseline_position_error_A:0,saved_baseline_center_error_nm:0,
      saved_baseline_energy_error:0};
    try {
      const file=car+'.mol',content=await fs.readFile(path.join(folder,'inputs',file),'utf8');
      const settings={...metadata.settings,temperature:150,damping_fs:350,lattice_scale:1};
      const oldInfo=old.prepareWithPose(content,file,settings,{mode:'auto'});
      const newInfo=current.prepareWithPose(content,file,settings,{mode:'auto'});
      assert.deepEqual(oldInfo.molecule,newInfo.molecule);
      assert.deepEqual(oldInfo.normal_counts,newInfo.normal_counts);
      assert.deepEqual(oldInfo.normal_neighbours,newInfo.normal_neighbours);
      assert.deepEqual(oldInfo.runtime.start_pose,newInfo.runtime.start_pose);
      for(const info of [oldInfo,newInfo]) {
        assert(info.runtime.preparation.accepted&&info.runtime.preparation.global_force_norm_kcal_mol_A<=1e-4);
        assert.equal(info.runtime.velocity_seed,17329);assert.equal(info.runtime.langevin_seed,2718281);
      }
      result.preparation_position_error_A=difference(oldInfo.positions,newInfo.positions);
      result.preparation_force_norm_error=Math.abs(oldInfo.runtime.preparation.global_force_norm_kcal_mol_A-newInfo.runtime.preparation.global_force_norm_kcal_mol_A);
      result.initial_velocity_error_A_fs=difference(old.phaseSpace().velocities,current.phaseSpace().velocities);
      for(const engine of [old,current])engine.update({drive:{mode:'gradient',strength:.6,growth}});
      const oe=old.energyForces(oldInfo.positions,1,.6,false),ne=current.energyForces(oldInfo.positions,1,.6,false);
      for(const key of ['energy','uff_energy','surface_energy','field_energy'])result.static_energy_error=Math.max(result.static_energy_error,Math.abs(oe[key]-ne[key]));
      for(const key of ['forces','uff_forces','surface_forces','field_forces'])result.static_force_error=Math.max(result.static_force_error,difference(oe[key],ne[key]));
      const saved=baseline.find(x=>x.config.molecule===car&&x.config.growth===growth&&x.config.seed===2718281);
      assert(saved,'Missing corresponding saved baseline');
      const coords=JSON.parse(zlib.gunzipSync(await fs.readFile(path.join(folder,'browser-trials',saved.coordinates_file))));
      for(let step=0;step<=2000;step+=100) {
        const a=step?old.advance(100):old.refresh(),b=step?current.advance(100):current.refresh();
        assert.equal(a.step,step);assert.equal(b.step,step);result.sampled_frames++;
        result.trajectory_position_error_A=Math.max(result.trajectory_position_error_A,difference(a.positions,b.positions));
        result.trajectory_velocity_error_A_fs=Math.max(result.trajectory_velocity_error_A_fs,difference(old.phaseSpace().velocities,current.phaseSpace().velocities));
        result.trajectory_energy_error=Math.max(result.trajectory_energy_error,Math.abs(a.energy-b.energy));
        const s=saved.samples[step/100];assert.equal(s.step,step);
        result.saved_baseline_center_error_nm=Math.max(result.saved_baseline_center_error_nm,difference(a.center.map(x=>x/10),s.center_nm));
        result.saved_baseline_energy_error=Math.max(result.saved_baseline_energy_error,Math.abs(a.energy-s.energy_kcal_mol));
        const ci=coords.steps.indexOf(step);
        if(ci>=0)result.saved_baseline_position_error_A=Math.max(result.saved_baseline_position_error_A,difference(a.positions,coords.positions[ci]));
      }
      for(const key of ['preparation_position_error_A','trajectory_position_error_A','saved_baseline_position_error_A'])assert(result[key]<1e-9,`${car} ${key}=${result[key]}`);
      for(const key of ['initial_velocity_error_A_fs','trajectory_velocity_error_A_fs'])assert(result[key]<1e-11,`${car} ${key}=${result[key]}`);
      for(const key of ['preparation_force_norm_error','static_energy_error','static_force_error','trajectory_energy_error','saved_baseline_energy_error'])assert(result[key]<1e-7,`${car} ${key}=${result[key]}`);
      assert(result.saved_baseline_center_error_nm<1e-10);
      report.cases.push({...result,status:'passed'});
      console.log(JSON.stringify({molecule:car,growth,status:'passed',position_error_A:result.trajectory_position_error_A,saved_baseline_energy_error:result.saved_baseline_energy_error}));
    } finally {old.close();old.delete();current.close();current.delete();}
  }
  report.status='passed';report.checked_baseline_trials=baseline.length;
  report.interpretation='Both binaries reproduce identical preparation, static response and 1 ps seeded Langevin trajectories for seven molecules and both baseline growth values within recorded tolerances. Saved baseline samples are independently reproduced. This validates reuse of the 42 specifically hashed scale=1 trials; it does not license arbitrary old binaries.';
} catch(error) {
  report.status='failed';
  if(typeof error==='number') {
    for(const module of [oldModule,currentModule])try{report.error=module.getExceptionMessage(error);break;}catch{}
  }
  report.error??=error.stack||String(error);
} finally {
  await fs.writeFile(path.join(folder,'baseline-compatibility.json'),JSON.stringify(report,null,2)+'\n');
}
console.log(JSON.stringify({status:report.status,cases:report.cases.length,error:report.error}));
if(report.status!=='passed')process.exitCode=1;
