// Execute the new controls in the real WASM module, against native MD data.
import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import {fileURLToPath,pathToFileURL} from 'node:url';
import assert from 'node:assert/strict';
import crypto from 'node:crypto';
import {execFile} from 'node:child_process';
import {promisify} from 'node:util';
const root=fileURLToPath(new URL('../../',import.meta.url));
const folder=path.join(root,'outputs/crystal-contrast-guarded');
const assets=path.join(root,'wasm/dist/browser/assets');
const report={status:'running',control_version:'six-phase-atomic-crystal-contrast-v3',source_sha256:{},static:[],paired_nve:[]};
const digest=data=>crypto.createHash('sha256').update(data).digest('hex');
const digestFile=async file=>digest(await fs.readFile(file));
const lockedFiles=new Map();
let module,engine;
try {
  const referenceFile='outputs/crystal-contrast-guarded/numerics.json';
  const referenceBytes=await fs.readFile(path.join(root,referenceFile));
  const fixtures=JSON.parse(referenceBytes.toString('utf8'));
  report.native_reference={path:referenceFile,sha256:digest(referenceBytes),status:fixtures.status,
    source_sha256:fixtures.source_sha256,source_integrity:fixtures.source_integrity};
  assert.equal(fixtures.status,'passed','Native numerical reference must pass before WASM validation');
  assert.equal(fixtures.source_integrity?.status,'passed','Native reference lacks passing source integrity');
  assert.deepEqual(fixtures.source_integrity.changed_files,[],'Native sources changed during reference generation');
  assert.equal(fixtures.source_integrity.files,Object.keys(fixtures.source_sha256||{}).length);
  assert.equal(fixtures.wasm_fixtures?.length,28,'Expected all 28 fixed-coordinate references');
  assert.equal(fixtures.nve?.length,12,'Expected all 12 native NVE references');
  assert(fixtures.nve.every(row=>Array.isArray(row.frames)&&row.frames.length>0),'NVE references must contain trajectories');
  // The native helper resolves external CHAu, RDKit and LAMMPS dependencies.
  // Check the helper itself before importing it, then verify every resolved file.
  const helper='scripts/validate_crystal_contrast.py';
  assert.equal(await digestFile(path.join(root,helper)),fixtures.source_sha256[helper],'Native resolver changed; regenerate numerics.json');
  const python=path.join(root,process.platform==='win32'?'.venv/Scripts/python.exe':'.venv/bin/python');
  const resolver=[
    'import json,sys',
    'from pathlib import Path',
    'sys.path.insert(0,str(Path.cwd()/"scripts"))',
    'from validate_crystal_contrast import validation_sources',
    'print(json.dumps({name:str(file.resolve()) for name,file in validation_sources(Path.cwd()/"outputs/crystal-contrast-guarded").items()}))'
  ].join(';');
  const {stdout}=await promisify(execFile)(python,['-E','-B','-c',resolver],{cwd:root,windowsHide:true,maxBuffer:1024*1024});
  const nativePaths=JSON.parse(stdout);
  assert.deepEqual(Object.keys(nativePaths).sort(),Object.keys(fixtures.source_sha256).sort(),'Native dependency set changed');
  report.native_reference.resolved_source_paths=nativePaths;
  for(const [name,file] of Object.entries(nativePaths)) {
    const expected=fixtures.source_sha256[name];
    assert.equal(await digestFile(file),expected,`Stale native reference: ${name}`);
    lockedFiles.set(file,expected);
  }
  const ownFiles=[referenceFile,'wasm/tools/validate_crystal_contrast.mjs',
    'wasm/dist/browser/runtime/molecular-car.mjs','wasm/dist/browser/runtime/molecular-car.wasm',
    'wasm/dist/browser/runtime/build.json','wasm/dist/browser/assets/model.json',
    ...['tile','saip-grid','hetero-grid'].map(name=>`wasm/dist/browser/assets/${name}.bin.gz`)];
  const build=JSON.parse(await fs.readFile(path.join(root,'wasm/dist/browser/runtime/build.json'),'utf8'));
  for(const [file,expected] of Object.entries(build.own_sources)) {
    assert.equal(await digestFile(path.join(root,file)),expected,`WASM build source is stale: ${file}`);
    ownFiles.push(file);
  }
  for(const file of ownFiles) {
    const absolute=path.join(root,file),sha=await digestFile(absolute);
    report.source_sha256[file]=sha;lockedFiles.set(absolute,sha);
  }
  assert.equal(report.source_sha256[referenceFile],report.native_reference.sha256,'Native report changed during setup');
  const metadata=JSON.parse(await fs.readFile(path.join(assets,'model.json'),'utf8'));
  const load=async name=>{
    const descriptor=metadata.surface.arrays[name];
    const file=path.join(root,'wasm/dist/browser',descriptor.path);
    assert(lockedFiles.has(file),`Untracked potential asset: ${descriptor.path}`);
    const data=zlib.gunzipSync(await fs.readFile(file));
    assert.equal(data.byteLength,descriptor.bytes,`Potential asset size: ${name}`);
    assert.equal(digest(data),descriptor.sha256,`Potential asset checksum: ${name}`);
    return new Float64Array(data.buffer.slice(data.byteOffset,data.byteOffset+data.byteLength));
  };
  const arrays=await Promise.all(['tile','saip-grid','hetero-grid'].map(load));
  const {default:create}=await import(pathToFileURL(path.join(root,'wasm/dist/browser/runtime/molecular-car.mjs')));
  module=await create({print:()=>{},printErr:()=>{}});
  engine=new module.BrowserEngine(metadata,...arrays);
  for(const fixture of fixtures.wasm_fixtures) {
    engine.prepareChemistry(fixture.content,fixture.filename);
    engine.update({lattice_scale:fixture.scale,drive:{growth:fixture.growth}});
    const current=engine.energyForces(fixture.positions,1,fixture.strength,false);
    const energy=Math.abs(current.energy-fixture.energy),force=Math.max(...current.forces.map((f,i)=>Math.abs(f-fixture.forces[i])));
    assert(energy<1e-7 && force<1e-7,`${fixture.filename} E=${energy} F=${force}`);
    report.static.push({filename:fixture.filename,scale:fixture.scale,growth:fixture.growth,energy_error:energy,force_error:force});
  }
  for(const fixture of fixtures.nve.filter(x=>x.frames)) {
    const content=await fs.readFile(path.join(folder,'inputs',fixture.molecule+'.mol'),'utf8');
    // Native has already accepted this zero-drive geometry. The baseline
    // preparation establishes an independent engine, then installs x/v exactly.
    engine.prepareWithPose(content,fixture.molecule+'.mol',{...metadata.settings,temperature:150,lattice_scale:1},{mode:'auto'});
    engine.update({lattice_scale:fixture.scale,drive:{mode:'gradient',strength:.6,growth:fixture.growth}});
    engine.nveWithDrive();
    engine.setPhaseSpace(fixture.phase_space.positions,fixture.phase_space.velocities);
    let position=0,energy=0;
    for(const row of fixture.frames) {
      const current=row.step?engine.advance(100):engine.refresh();
      position=Math.max(position,...current.positions.map((v,i)=>Math.abs(v-row.positions[i])));
      energy=Math.max(energy,Math.abs(current.energy-row.energy));
      assert(current.step===row.step);
    }
    assert(position<1e-7 && energy<1e-7,`${fixture.molecule} NVE x=${position} E=${energy}`);
    report.paired_nve.push({molecule:fixture.molecule,scale:fixture.scale,growth:fixture.growth,position_error_A:position,energy_error:energy});
  }
  const changed=[];
  for(const [file,expected] of lockedFiles) {
    try{if(await digestFile(file)!==expected)changed.push(file);}catch{changed.push(file);}
  }
  report.source_integrity={status:changed.length?'failed':'passed',files:lockedFiles.size,
    check:'native reference, resolved native dependencies, validator, WASM source/build/runtime and assets hashed before and after validation',changed_files:changed};
  assert.deepEqual(changed,[],'Validation inputs changed during WASM execution');
  report.status='passed';
} catch(error) {
  report.status='failed';
  if(error instanceof Error)report.error=error.message;
  else try{report.error=module.getExceptionMessage(error);module.decrementExceptionRefcount(error);}catch{report.error=String(error);}
} finally {
  if(engine){engine.close();engine.delete();}
  await fs.writeFile(path.join(folder,'wasm-numerics.json'),JSON.stringify(report,null,2)+'\n');
}
console.log(JSON.stringify(report));
if(report.status!=='passed')process.exitCode=1;
