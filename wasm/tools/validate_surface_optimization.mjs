// Same-energy optimization acceptance against the frozen pre-optimization WASM.
import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {fileURLToPath,pathToFileURL} from 'node:url';
const root=fileURLToPath(new URL('../../',import.meta.url));
const out=path.join(root,'outputs/surface-kernel-optimization');
const study=path.join(root,'outputs/crystal-contrast-guarded');
const assets=path.join(root,'wasm/dist/browser/assets');
const read=async f=>JSON.parse(await fs.readFile(f,'utf8'));
const digest=async f=>crypto.createHash('sha256').update(await fs.readFile(f)).digest('hex');
const save=async report=>fs.writeFile(path.join(out,'numerical-validation.json'),JSON.stringify(report,null,2)+'\n');
const metadata=await read(path.join(assets,'model.json'));
const arrays=await Promise.all(['tile','saip-grid','hetero-grid'].map(async n=>{
  const b=zlib.gunzipSync(await fs.readFile(path.join(assets,n+'.bin.gz')));
  assert.equal(crypto.createHash('sha256').update(b).digest('hex'),metadata.surface.arrays[n].sha256);
  return new Float64Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength));
}));
const report={status:'running',static:[],finite_difference:[],paired_nve:[],sources:{},
  protocol:{precision:'float64',dt_fs:.5,physics:'six-phase-atomic-crystal-contrast-v3 unchanged',
    static_tolerance:'1e-7 absolute plus 2e-12 times component magnitude for high-energy boundary probes',
    nve_tolerance:{position_A:1e-6,energy_kcal_mol:1e-6}}};
const original=await read(path.join(out,'reference/manifest.json'));
for(const [f,h] of Object.entries(original))assert.equal(await digest(path.join(out,'reference',f)),h,f);
const build=await read(path.join(root,'wasm/dist/browser/runtime/build.json'));
for(const [f,h] of Object.entries(build.own_sources))assert.equal(await digest(path.join(root,f)),h,f);
const sourceFiles=[...Object.keys(build.own_sources),'wasm/tools/validate_surface_optimization.mjs',
  'wasm/dist/browser/runtime/build.json','wasm/dist/browser/runtime/molecular-car.mjs','wasm/dist/browser/runtime/molecular-car.wasm',
  'wasm/dist/browser/assets/model.json',...['tile','saip-grid','hetero-grid'].map(n=>`wasm/dist/browser/assets/${n}.bin.gz`),
  'outputs/crystal-contrast-guarded/numerics.json','outputs/crystal-contrast-guarded/1p5-numerics.json'];
for(const f of sourceFiles)report.sources[f]=await digest(path.join(root,f));
const modules=[];const engines=[];
const maxdiff=(a,b)=>Math.max(0,...a.map((v,i)=>Math.abs(v-b[i])));
try {
  for(const directory of [path.join(out,'reference'),root]) {
    const {default:create}=await import(pathToFileURL(path.join(directory,'wasm/dist/browser/runtime/molecular-car.mjs')));
    const module=await create({print:()=>{},printErr:()=>{}});modules.push(module);
    engines.push(new module.BrowserEngine(metadata,...arrays));
  }
  const [old,current]=engines;
  const base=await read(path.join(study,'numerics.json')),extra=await read(path.join(study,'1p5-numerics.json'));
  const fixtures=[...base.wasm_fixtures,...extra.wasm_fixtures];
  const compare=(positions,label,levels=[0,1,2,3])=>{
    const reference=old.energyForces(positions,1,.6,false);
    for(const level of levels) {
      current.surfaceOptimization(level);const actual=current.energyForces(positions,1,.6,false);
      const e=Math.abs(actual.surface_energy-reference.surface_energy),f=maxdiff(actual.surface_forces,reference.surface_forces);
      assert(e<1e-7+Math.abs(reference.surface_energy)*2e-12,`${label} level ${level} energy ${e}`);
      assert(f<1e-7+Math.max(...reference.surface_forces.map(Math.abs))*2e-12,`${label} level ${level} force ${f}`);
      assert.equal(actual.uff_energy,reference.uff_energy);assert.deepEqual(actual.uff_forces,reference.uff_forces);
      report.static.push({label,level,energy_error:e,force_error:f});
    }
  };
  for(const fixture of fixtures) {
    for(const engine of engines){engine.prepareChemistry(fixture.content,fixture.filename);engine.update({lattice_scale:fixture.scale,drive:{growth:fixture.growth}});}
    compare(fixture.positions,`${fixture.filename}/L${fixture.scale}/G${fixture.growth}`);
  }
  const [px,py]=metadata.surface.period;
  for(const car of ['benzene-12','test-car-52','fullerene-c60','pyridine-11','benzonitrile-13']) {
    const fixture=extra.wasm_fixtures.find(f=>f.filename===car+'.mol');
    for(const engine of engines){engine.prepareChemistry(fixture.content,fixture.filename);engine.update({lattice_scale:1.5,drive:{growth:1}});}
    const original=fixture.positions;
    // Tilt the complete molecule: bonded normals span three-dimensional directions.
    const tilted=original.map((v,i)=>i%3===0?v*.8+original[i+2]*.6:i%3===2?v*.8-original[i-2]*.6:v);
    const zmin=Math.min(...tilted.filter((v,i)=>i%3===2));
    for(const height of [1.5,1.999999,2.,2.000001,2.15,2.3,2.6,3.,11.7,11.999999,12.,12.000001,16.1,35.]) {
      const p=tilted.map((v,i)=>v+(i%3===2?height-zmin:0));
      compare(p,`${car}/tilted/minZ${height}`,[1,2,3]);
    }
    for(const offset of [-1e-9,0,1e-9])for(const cell of [0,1,3,4]) {
      const dx=cell*px/4+offset-original[0],dy=cell*py/6+offset-original[1];
      const p=original.map((v,i)=>v+(i%3===0?dx:i%3===1?dy:0));
      compare(p,`${car}/cell${cell}/${offset}`,[3]);
      compare(p.map((v,i)=>v+(i%3===0?-103*px:i%3===1?97*py:0)),`${car}/translated-cell${cell}/${offset}`,[3]);
    }
    current.surfaceOptimization(3);
    // All Cartesian derivatives in the short-range gate and arbitrary tilt.
    const p=tilted.map((v,i)=>v+(i%3===2?2.15-zmin:0)),force=current.energyForces(p,1,.6,false).surface_forces;
    const convergence=[];
    // Catmull-Rom is C1, not C2: an atom exactly on a knot has O(h), not
    // O(h^2), central-difference convergence. Check the old kernel too.
    for(const h of [1e-5,1e-6]) {
      let error=0,reference_error=0;
      const reference_force=old.energyForces(p,1,.6,false).surface_forces;
      for(let i=0;i<p.length;i++) {
        const a=p.slice(),b=p.slice();a[i]+=h;b[i]-=h;
        const fd=-(current.energyForces(a,1,.6,false).surface_energy-current.energyForces(b,1,.6,false).surface_energy)/(2*h);
        const reference_fd=-(old.energyForces(a,1,.6,false).surface_energy-old.energyForces(b,1,.6,false).surface_energy)/(2*h);
        error=Math.max(error,Math.abs(fd-force[i]));reference_error=Math.max(reference_error,Math.abs(reference_fd-reference_force[i]));
      }
      convergence.push({h,error,reference_error});
    }
    const error=convergence.at(-1).error;
    assert(error<2e-5,`${car} finite difference ${error}`);report.finite_difference.push({car,error,coordinates:p.length,convergence});
  }
  await save(report);console.log(JSON.stringify({static:report.static.length,finite_difference:report.finite_difference.length}));
  current.surfaceOptimization(3);
  for(const fixture of [...base.nve,...extra.nve]) {
    const content=await fs.readFile(path.join(study,'inputs',fixture.molecule+'.mol'),'utf8');
    current.prepareWithPose(content,fixture.molecule+'.mol',{...metadata.settings,temperature:150,lattice_scale:1},{mode:'auto'});
    current.update({lattice_scale:fixture.scale,drive:{mode:'gradient',strength:.6,growth:fixture.growth}});
    current.nveWithDrive();current.setPhaseSpace(fixture.phase_space.positions,fixture.phase_space.velocities);
    let x=0,e=0;
    for(const row of fixture.frames) {
      const actual=row.step?current.advance(100):current.refresh();
      assert.equal(actual.step,row.step);x=Math.max(x,maxdiff(actual.positions,row.positions));e=Math.max(e,Math.abs(actual.energy-row.energy));
    }
    assert(x<1e-6&&e<1e-6,`${fixture.molecule}/L${fixture.scale} NVE x=${x} e=${e}`);
    report.paired_nve.push({car:fixture.molecule,scale:fixture.scale,growth:fixture.growth,position_error_A:x,energy_error:e});
    await save(report);
  }
  for(const [f,h] of Object.entries(report.sources))assert.equal(await digest(path.join(root,f)),h,`Changed source ${f}`);
  report.source_integrity='passed';report.status='passed';
} catch(error) {
  report.status='failed';report.error=error.message||String(error);
  console.error(report.error);process.exitCode=1;
} finally {
  for(const engine of engines){engine.close();engine.delete();}
  await save(report);
}
console.log(JSON.stringify({status:report.status,static:report.static.length,fd:report.finite_difference,paired_nve:report.paired_nve.length,error:report.error}));
