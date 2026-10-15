// Rigid geometry parity and actual browser startup, using the built WASM module.
import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import {fileURLToPath,pathToFileURL} from 'node:url';
import assert from 'node:assert/strict';
const root=fileURLToPath(new URL('../../',import.meta.url)),assets=path.join(root,'wasm/dist/browser/assets');
const metadata=JSON.parse(await fs.readFile(path.join(assets,'model.json'),'utf8'));
const load=async key=>{const bytes=zlib.gunzipSync(await fs.readFile(path.join(assets,`${key}.bin.gz`)));return new Float64Array(bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength));};
const [tile,saip,hetero]=await Promise.all(['tile','saip-grid','hetero-grid'].map(load));
const {default:create}=await import(pathToFileURL(path.join(root,'wasm/dist/browser/runtime/molecular-car.mjs')));
const module=await create({print:()=>{},printErr:()=>{}}),engine=new module.BrowserEngine(metadata,tile,saip,hetero);
const fixturesDirectory=path.resolve(root,process.argv[2]||'outputs/start-pose-validation/final');
const fixtures=JSON.parse(await fs.readFile(path.join(fixturesDirectory,'pose-fixtures.json'),'utf8'));
const report={passed:false,geometry_cases:[],race_cases:[]};
function message(error){try{const m=module.getExceptionMessage(error);module.decrementExceptionRefcount(error);return Array.isArray(m)?m.join(': '):String(m);}catch{return error.message||String(error);}}
try {
  for(const fixture of fixtures) {
    const result=engine.prepareChemistryWithPose(fixture.content,fixture.filename,fixture.config),pose=result.runtime.start_pose;
    const error=Math.max(...result.positions.map((x,i)=>Math.abs(x-fixture.positions[i])));
    const rotationError=Math.max(...pose.rotation_matrix.flat().map((x,i)=>Math.abs(x-fixture.report.rotation_matrix.flat()[i])));
    assert(error<1e-7&&rotationError<1e-7,`${fixture.filename}: ${error} / ${rotationError}`);
    assert.equal(pose.strategy,fixture.report.strategy);assert.deepEqual(pose.heading_atom_ids,fixture.report.heading_atom_ids);
    report.geometry_cases.push({filename:fixture.filename,mode:fixture.config.mode,position_error_A:error,rotation_error:rotationError,strategy:pose.strategy});
  }
  const bad=[{mode:'invalid'},{mode:'manual',yaw_deg:NaN},{mode:'auto',axis_atoms:[1,1]},{mode:'auto',axis_atoms:[1,256]},{extra:1}];
  report.invalid_requests=[];
  for(const config of bad) {let reason;try{engine.prepareChemistryWithPose(fixtures[0].content,fixtures[0].filename,config);}catch(e){reason=message(e);}assert(reason);report.invalid_requests.push(reason);}
  for(const [name,config] of [['test-car-upright.mol',{mode:'auto'}],['test-car-original.mol',{mode:'manual',yaw_deg:90}],
                            ['connected-fullerene.mol',{mode:'auto'}]]) {
    const fixture=fixtures.find(f=>f.filename===name),start=performance.now();
    const info=engine.prepareWithPose(fixture.content,name,metadata.settings,config);
    assert(info.runtime.preparation.accepted&&info.runtime.preparation.global_force_norm_kcal_mol_A<=1e-4);
    const initial=engine.frame();assert.equal(initial.step,0);let frame;
    for(let i=0;i<20;i++)frame=engine.advance(100);
    assert.equal(frame.step,2000);assert.equal(frame.time_fs,1000);assert(frame.positions.every(Number.isFinite));
    report.race_cases.push({filename:name,config,pose:info.runtime.start_pose,preparation:info.runtime.preparation,
      wall_seconds:(performance.now()-start)/1000,final_center:frame.center});
  }
  report.passed=true;
}catch(error){report.error=message(error);}finally{engine.delete();}
await fs.writeFile(path.join(root,'outputs/start-pose-validation/wasm-report.json'),JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify({passed:report.passed,geometry_cases:report.geometry_cases.length,race_cases:report.race_cases.length,error:report.error}));
if(!report.passed)process.exitCode=1;
