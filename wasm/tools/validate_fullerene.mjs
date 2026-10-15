// Execute the real browser binary on the previously rejected workshop molecule.
import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import {fileURLToPath,pathToFileURL} from 'node:url';
import assert from 'node:assert/strict';
const root=fileURLToPath(new URL('../../',import.meta.url));
const assets=path.join(root,'wasm/dist/browser/assets');
const metadata=JSON.parse(await fs.readFile(path.join(assets,'model.json'),'utf8'));
const load=async key=>{const data=zlib.gunzipSync(await fs.readFile(path.join(assets,`${key}.bin.gz`)));
  return new Float64Array(data.buffer.slice(data.byteOffset,data.byteOffset+data.byteLength));};
const [tile,saip,hetero]=await Promise.all(['tile','saip-grid','hetero-grid'].map(load));
const {default:create}=await import(pathToFileURL(path.join(root,'wasm/dist/browser/runtime/molecular-car.mjs')));
const module=await create({print:()=>{},printErr:()=>{}});
const content=await fs.readFile(path.join(root,'tests/data/connected-fullerene.mol'),'utf8');
const output=path.join(root,'outputs/joint-fullerene-diagnosis/wasm-report.json');
const report={passed:false,backend_required:false,cases:[]};
function message(error){try {const value=module.getExceptionMessage(error);module.decrementExceptionRefcount(error);
  return Array.isArray(value)?value.join(': '):String(value);}catch{return error.message||String(error);}}
try {
  for(const surface of ['periodic-au','continuum']) {
    const engine=new module.BrowserEngine(metadata,tile,saip,hetero);
    try {
      const started=performance.now();
      const info=engine.prepare(content,'connected-fullerene.mol',{...metadata.settings,surface});
      const row={surface,atom_count:info.molecule.atoms.length,preparation:info.runtime.preparation,
        prepare_wall_seconds:(performance.now()-started)/1000};
      assert.equal(row.atom_count,72);assert(info.runtime.preparation.accepted);
      assert(info.runtime.preparation.final_coordinate_force_verified);
      assert(info.runtime.preparation.global_force_norm_kcal_mol_A<=1e-4);
      assert(info.runtime.preparation.force_callback_evaluations<=25010);
      assert.equal(info.runtime.preparation.protocol_version,'joint-force-relaxation-v2');
      let frame=engine.frame();assert.equal(frame.step,0);
      const conservative=engine.energyForces(frame.positions,1,0,surface==='continuum');
      row.independent_initial_force_norm=Math.hypot(...conservative.forces);
      assert(row.independent_initial_force_norm<=1e-4);
      if(surface==='periodic-au') {
        const native=JSON.parse(await fs.readFile(path.join(root,'outputs/joint-fullerene-diagnosis/native-static.json'),'utf8'));
        const paired=engine.energyForces(native.positions,1,0,false);
        row.paired_native_force_error=Math.max(...paired.forces.map((f,index)=>Math.abs(f-native.forces[index])));
        row.paired_native_energy_error=Math.abs(paired.energy-native.energy);
        assert(row.paired_native_force_error<1e-7&&row.paired_native_energy_error<1e-7);
      }
      row.initial=frame;
      const run=performance.now();
      for(let index=0;index<100;index++)frame=engine.advance(100);
      row.propagation_wall_seconds=(performance.now()-run)/1000;row.final=frame;
      assert.equal(frame.step,10000);assert.equal(frame.time_fs,5000);
      assert(frame.positions.every(Number.isFinite));assert(frame.center[0]>row.initial.center[0]);
      if(surface==='periodic-au') {
        const native=JSON.parse(await fs.readFile(path.join(root,'outputs/joint-fullerene-diagnosis/native-5ps/report.json'),'utf8'));
        assert.equal(native.smoke.steps,10000);
        row.paired_5ps_center_error_A=Math.max(...frame.center.map((v,index)=>Math.abs(v-native.smoke.final_center_A[index])));
        assert(row.paired_5ps_center_error_A<1e-5);
      }
      report.cases.push(row);
    }finally{engine.delete();}
  }
  // Deliberately exhaust a shared budget in a BFGS line search, then verify that
  // a fresh preparation on the same owner is unaffected by the old soft timer.
  const limited=structuredClone(metadata);limited.preparation.stages=[];
  limited.preparation.evaluation_limit=450;
  const engine=new module.BrowserEngine(limited,tile,saip,hetero);
  try {
    let reason;try{engine.prepare(content,'connected-fullerene.mol',metadata.settings);}catch(error){reason=message(error);}
    assert(reason?.includes('force-evaluation-budget'),reason);
    report.budget_rejection=reason;
    const benzene=await fs.readFile(path.join(assets,'benzene.mol'),'utf8');
    const recovery=engine.prepare(benzene,'benzene.mol',metadata.settings);
    assert(recovery.runtime.preparation.accepted);
    const recovered=engine.advance(100);assert.equal(recovered.step,100);
    assert(recovered.positions.every(Number.isFinite));
    report.recovery_after_budget_stop=true;
  }finally{engine.delete();}
  report.passed=true;
}catch(error){report.error=message(error);}
await fs.mkdir(path.dirname(output),{recursive:true});
await fs.writeFile(output,JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify({passed:report.passed,cases:report.cases.map(row=>({surface:row.surface,
  preparation:row.preparation,final_step:row.final.step,final_center:row.final.center})),
  budget_rejection:report.budget_rejection,error:report.error,output},null,2));
if(!report.passed)process.exitCode=1;
