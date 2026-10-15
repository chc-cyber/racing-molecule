// Run only in a reserved serial CPU window. Never promotes a timestep or asset.
import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import crypto from 'node:crypto';
import os from 'node:os';
import assert from 'node:assert/strict';
import {fileURLToPath,pathToFileURL} from 'node:url';

const root=fileURLToPath(new URL('../../',import.meta.url));
const browser=path.join(root,'wasm/dist/browser');
const v1Browser=path.join(root,'outputs/parallel-md-validation/baseline/wasm/dist/browser');
const atlasRoot=path.join(root,'outputs/c60-proxy-study');
const out=path.join(root,'outputs/parallel-md-validation');
const read=async f=>JSON.parse(await fs.readFile(f,'utf8'));
const sha=async f=>crypto.createHash('sha256').update(await fs.readFile(f)).digest('hex');
const cars=['benzene-12','pyridine-11','test-car-52','fullerene-c60'];
const load=async(base,d)=>{
  const b=zlib.gunzipSync(await fs.readFile(path.join(base,d.path)));
  assert.equal(b.length,d.bytes);assert.equal(crypto.createHash('sha256').update(b).digest('hex'),d.sha256);
  return new Float64Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength));
};
const meta=await read(path.join(browser,'assets/model.json'));
const old=await read(path.join(v1Browser,'assets/fourier.json'));
const atlas=await read(path.join(atlasRoot,'assets/fourier-atlas.json'));
for(const p of [old,atlas])assert.equal(p.base_model_sha256,await sha(path.join(browser,'assets/model.json')));
const build=await read(path.join(browser,'runtime/build.json'));
for(const [f,h] of Object.entries(build.own_sources))assert.equal(await sha(path.join(root,f)),h,'stale main build '+f);
const sourceFiles=[...Object.keys(build.own_sources).map(f=>path.join(root,f)),
  ...['runtime/molecular-car.wasm','runtime/molecular-car.mjs','runtime/build.json','assets/model.json'].map(f=>path.join(browser,f)),
  ...Object.values(meta.surface.arrays).map(d=>path.join(browser,d.path)),path.join(v1Browser,'assets/fourier.json'),path.join(v1Browser,old.arrays[0].path),
  path.join(atlasRoot,'assets/fourier-atlas.json'),path.join(atlasRoot,atlas.arrays[0].path),
  fileURLToPath(import.meta.url),...cars.map(c=>path.join(root,'outputs/wasm-fourier-validation/races/inputs',c+'.mol'))];
const sources={};for(const f of new Set(sourceFiles))sources[f]=await sha(f);
const checkSources=async()=>{for(const [f,h] of Object.entries(sources))assert.equal(await sha(f),h,'source changed during serial benchmark '+f);};
const arrays=await Promise.all(['tile','saip-grid','hetero-grid'].map(k=>load(browser,meta.surface.arrays[k])));
const oldCoefficients=await load(v1Browser,old.arrays[0]),newCoefficients=await load(atlasRoot,atlas.arrays[0]);
const {default:create}=await import(pathToFileURL(path.join(browser,'runtime/molecular-car.mjs')));
const module=await create({print:()=>{},printErr:()=>{}}),e=new module.BrowserEngine(meta,...arrays);
const report={status:'running',sources,hardware:{cpu:os.cpus()[0].model,node:process.version},protocol:{
  temperature_K:150,strength:.6,growth:1,lattice_scales:[1.5,2],duration_ps:2,repeats:3,
  serial:true,load_and_preparation_excluded:true,sample_fs:50,initial_state:'identical periodic-relaxed positions and seeded velocities',
  cases:[{mode:'periodic',dt_fs:.5},{mode:'v1',dt_fs:.5},{mode:'v2',dt_fs:.5},{mode:'v2',dt_fs:1,exploratory:true}],
  candidate_notice:'1 fs is exploratory even if its MD remains finite; this script does not override timestep screening',
  nve:{mode:'v2',timesteps_fs:[.25,.5,1],duration_ps:2,sample_fs:10,drive:false,thermostat:false},
  extrapolation:'seconds_per_100ps is 50 times measured 2 ps MD, not a timed full race; concurrent CPU work invalidates speedup interpretation'},md:[],nve:[]};
const save=async()=>{await fs.mkdir(out,{recursive:true});await fs.writeFile(path.join(out,'benchmark.json'),JSON.stringify(report,null,2)+'\n');};
let installed;
function mode(name){
  if(name==='periodic'){e.fourierEnabled(false);return;}
  if(installed!==name){e.installFourier(name==='v1'?old:atlas);e.loadFourierSpecies(0,name==='v1'?oldCoefficients:newCoefficients);installed=name;}
  e.fourierEnabled(true);
}
function guards(info,initial){
  const distance=(p,b)=>Math.hypot(...[0,1,2].map(k=>p[3*(b.a-1)+k]-p[3*(b.b-1)+k]));
  const bonds=info.molecule.bonds,reference=bonds.map(b=>distance(initial.positions,b));
  return frame=>{
    assert(Number.isFinite(frame.energy)&&frame.positions.every(Number.isFinite),'nonfinite trajectory');
    const r=bonds.map((b,i)=>distance(frame.positions,b)/reference[i]);
    const bounds=[Math.min(...r),Math.max(...r)];assert(bounds[0]>=.65&&bounds[1]<=1.45,'bond guard');return bounds;
  };
}
try{
  for(const car of cars)for(const scale of [1.5,2]){
    const content=await fs.readFile(path.join(root,'outputs/wasm-fourier-validation/races/inputs',car+'.mol'),'utf8');
    const info=e.prepareWithPose(content,car+'.mol',{...meta.settings,surface:'periodic-au',temperature:150,damping_fs:350,lattice_scale:scale,timestep_fs:.5},{mode:'auto'});
    assert(info.runtime.preparation.accepted);e.update({drive:{mode:'gradient',strength:.6,growth:1}});
    const initial=e.phaseSpace(),guard=guards(info,initial);
    const configs=[['periodic',.5],['v1',.5],['v2',.5],['v2',1]];
    for(let repeat=0;repeat<3;repeat++){
      const order=[...configs.slice(repeat),...configs.slice(0,repeat)];
      for(const [name,dt] of order){
        mode(name);e.setPhaseSpace(initial.positions,initial.velocities);e.setTimestep(dt);e.conservative(false);e.update({temperature:150});e.refresh();
        let frame;const counts=[0,0,0];let lo=1,hi=1;const began=performance.now();
        for(let i=1;i<=40;i++){
          frame=e.advance(50/dt);assert.equal(frame.time_fs,50*i);const b=guard(frame);lo=Math.min(lo,b[0]);hi=Math.max(hi,b[1]);
          const c=frame.fourier_counts;for(let j=0;j<3;j++)counts[j]+=c[j];
        }
        const seconds=(performance.now()-began)/1000;
        report.md.push({car,scale,repeat,mode:name,dt_fs:dt,exploratory:dt===1,seconds,seconds_per_100ps:seconds*50,
          counts,bond_ratio_range:[lo,hi],final_center_A:frame.center,final_time_fs:frame.time_fs});
      }
    }
    for(const dt of [.25,.5,1]){
      mode('v2');e.setPhaseSpace(initial.positions,initial.velocities);e.setTimestep(dt);e.conservative(true);
      const first=e.refresh(),ke=first.energy-e.energyForces(first.positions,1,0,false).energy;
      assert(Number.isFinite(ke)&&ke>0);let drift=0,lo=1,hi=1;const samples=[];
      for(let i=1;i<=200;i++){
        const f=e.advance(10/dt);assert.equal(f.time_fs,10*i);const b=guard(f);lo=Math.min(lo,b[0]);hi=Math.max(hi,b[1]);
        drift=Math.max(drift,Math.abs(f.energy-first.energy));samples.push({time_fs:f.time_fs,energy_delta:f.energy-first.energy});
      }
      report.nve.push({car,scale,mode:'v2',dt_fs:dt,exploratory:dt===1,initial_kinetic_energy:ke,
        max_total_energy_error:drift,relative_to_initial_kinetic:drift/ke,bond_ratio_range:[lo,hi],samples});
      e.conservative(false);
    }
    await save();console.log(JSON.stringify({car,scale,md:report.md.length,nve:report.nve.length}));
  }
  await checkSources();report.status='complete';
}catch(error){report.status='failed';report.error=String(error.stack||error);process.exitCode=1;}
finally{e.close();e.delete();await save();}
console.log(JSON.stringify({status:report.status,md:report.md.length,nve:report.nve.length,error:report.error}));
