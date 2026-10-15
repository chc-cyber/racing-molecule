// Serial, alternating old/new benchmarks; no race workers run concurrently.
import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import os from 'node:os';
import {fileURLToPath,pathToFileURL} from 'node:url';
const root=fileURLToPath(new URL('../../',import.meta.url));
const out=path.join(root,'outputs/surface-kernel-optimization');
const study=path.join(root,'outputs/crystal-contrast-guarded');
const profileOnly=process.argv.includes('--profile-only');
const assets=path.join(root,'wasm/dist/browser/assets');
const read=async f=>JSON.parse(await fs.readFile(f,'utf8'));
const metadata=await read(path.join(assets,'model.json'));
const arrays=await Promise.all(['tile','saip-grid','hetero-grid'].map(async n=>{
  const b=zlib.gunzipSync(await fs.readFile(path.join(assets,n+'.bin.gz')));
  return new Float64Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength));
}));
const {default:create}=await import(pathToFileURL(path.join(root,'wasm/dist/browser/runtime/molecular-car.mjs')));
const module=await create({print:()=>{},printErr:()=>{}});
const engine=new module.BrowserEngine(metadata,...arrays);
const fixtures=(await read(path.join(study,'1p5-numerics.json'))).wasm_fixtures.filter(f=>f.growth===1);
const report={status:'running',hardware:{cpu:os.cpus()[0].model,logical_cpus:os.cpus().length,node:process.version},
  protocol:{serial:true,repeats:3,md_sample_ps:2,dt_fs:.5,temperature_K:150,strength:.6,growth:1,
    description:'One engine; identical initial phase space per paired run; old/new alternate order; MD scores every .05 ps. 100 ps times are extrapolated from 2 ps intervals, not full race timings.',
    profiling:'Kernel-only component clocks; warmed, reusable buffers; not a full per-step profiler'},
  static_edges:[],components:[],md:[],source_sha256:{}};
report.protocol.profile_only=profileOnly;
const save=async()=>fs.writeFile(path.join(out,profileOnly?'profile.json':'benchmark.json'),JSON.stringify(report,null,2)+'\n');
const digest=async f=>crypto.createHash('sha256').update(await fs.readFile(path.join(root,f))).digest('hex');
for(const f of ['wasm/tools/benchmark_surface_optimization.mjs','wasm/src/potentials.cpp','wasm/src/potentials.hpp','wasm/src/engine.cpp',
  'wasm/dist/browser/runtime/molecular-car.wasm','wasm/dist/browser/runtime/molecular-car.mjs','wasm/dist/browser/assets/model.json'])report.source_sha256[f]=await digest(f);
try {
  for(const fixture of fixtures) {
    engine.prepareChemistry(fixture.content,fixture.filename);
    for(const scale of [0,.5,1,1.5,2]) {
      engine.update({lattice_scale:scale,drive:{mode:'gradient',strength:.6,growth:1}});
      engine.surfaceOptimization(0);const reference=engine.energyForces(fixture.positions,1,.6,false);
      engine.surfaceOptimization(3);const current=engine.energyForces(fixture.positions,1,.6,false);
      const force=Math.max(...current.forces.map((v,i)=>Math.abs(v-reference.forces[i]))),energy=Math.abs(current.energy-reference.energy);
      assert(force<1e-7&&energy<1e-7);report.static_edges.push({car:fixture.filename,scale,force_error:force,energy_error:energy});
      if(![1,1.5,2].includes(scale))continue;
      for(let repeat=0;repeat<3;repeat++)for(const level of repeat%2?[3,2,1,0]:[0,1,2,3]) {
        engine.surfaceOptimization(level);engine.forceTiming(fixture.positions,5);
        const row=engine.forceTiming(fixture.positions,100);
        report.components.push({car:fixture.filename,scale,level,repeat,...row});
      }
    }
    console.log(JSON.stringify({profiled:fixture.filename}));await save();
  }
  for(const fixture of profileOnly?[]:fixtures)for(const scale of [1.5,2]) {
    engine.surfaceOptimization(0);
    engine.prepareWithPose(fixture.content,fixture.filename,{...metadata.settings,temperature:150,lattice_scale:scale},{mode:'auto'});
    engine.update({drive:{mode:'gradient',strength:.6,growth:1}});
    const phase=engine.phaseSpace();
    for(let repeat=0;repeat<3;repeat++)for(const level of repeat%2?[3,0]:[0,3]) {
      engine.surfaceOptimization(level);engine.setPhaseSpace(phase.positions,phase.velocities);
      engine.update({temperature:150}); // Reset Langevin state identically if supported by current control path.
      const start=performance.now();let frame;
      for(let i=0;i<40;i++)frame=engine.advance(100);
      const seconds=(performance.now()-start)/1000;
      report.md.push({car:fixture.filename,scale,level,repeat,seconds,seconds_per_100ps:seconds*50,final_center:frame.center,finite:frame.positions.every(Number.isFinite)});
      assert(report.md.at(-1).finite);
    }
    console.log(JSON.stringify({md:fixture.filename,scale}));await save();
  }
  for(const [f,h] of Object.entries(report.source_sha256))assert.equal(await digest(f),h,`Changed source ${f}`);
  report.source_integrity='passed';report.status='passed';
} catch(error){report.status='failed';report.error=error.message||String(error);console.error(report.error);process.exitCode=1;}
finally{engine.close();engine.delete();await save();}
console.log(JSON.stringify({status:report.status,components:report.components.length,md:report.md.length,error:report.error}));
