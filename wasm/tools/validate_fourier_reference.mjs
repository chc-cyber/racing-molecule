// Ensure enabling the optional proxy did not change the default periodic model.
import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import assert from 'node:assert/strict';
import {pathToFileURL} from 'node:url';
import {root,out,read,initialize,sources,checkSources,hash} from './fourier_support.mjs';
const current=await initialize(),base=path.join(out,'reference');
const manifest=await read(path.join(base,'manifest.json'));
for(const [f,h] of Object.entries(manifest))assert.equal(await hash(path.join(base,f)),h);
const {default:create}=await import(pathToFileURL(path.join(base,'wasm/dist/browser/runtime/molecular-car.mjs')));
const oldModule=await create({print:()=>{},printErr:()=>{}});
const arrays=await Promise.all(['tile','saip-grid','hetero-grid'].map(async name=>{const b=zlib.gunzipSync(await fs.readFile(path.join(root,'wasm/dist/browser/assets',name+'.bin.gz')));return new Float64Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength));}));
const old=new oldModule.BrowserEngine(current.meta,...arrays),e=current.engine;
const report={status:'running',sources:await sources(),cases:0,max_force_error:0,max_energy_error:0};
try {
  let car;
  for(const f of await read(path.join(out,'native-fixtures.json'))) {
    if(f.car!==car){car=f.car;old.prepareChemistry(f.content,car+'.mol');e.prepareChemistry(f.content,car+'.mol');}
    old.update({lattice_scale:f.scale});e.update({lattice_scale:f.scale});e.fourierEnabled(false);
    const a=old.energyForces(f.positions,1,.6,false),b=e.energyForces(f.positions,1,.6,false);
    report.max_force_error=Math.max(report.max_force_error,...a.surface_forces.map((v,i)=>Math.abs(v-b.surface_forces[i])));
    report.max_energy_error=Math.max(report.max_energy_error,Math.abs(a.surface_energy-b.surface_energy));report.cases++;
  }
  assert(report.max_force_error<1e-10&&report.max_energy_error<1e-10);await checkSources(report.sources);report.status='passed';
}catch(error){report.status='failed';report.error=error.stack||String(error);process.exitCode=1;}
finally{old.close();old.delete();e.close();e.delete();await fs.writeFile(path.join(out,'periodic-regression.json'),JSON.stringify(report,null,2)+'\n');}
console.log(JSON.stringify(report));
