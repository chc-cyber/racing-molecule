// Exercise the optimized hetero-element cache combination against the archived module.
import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {fileURLToPath,pathToFileURL} from 'node:url';
const root=fileURLToPath(new URL('../../',import.meta.url));
const out=path.join(root,'outputs/surface-kernel-optimization');
const assets=path.join(root,'wasm/dist/browser/assets');
const read=async p=>JSON.parse(await fs.readFile(p,'utf8'));
const digest=async p=>crypto.createHash('sha256').update(await fs.readFile(p)).digest('hex');
const gate=await read(path.join(out,'numerical-validation.json'));
assert.equal(gate.status,'passed');
for(const [f,h] of Object.entries(gate.sources))assert.equal(await digest(path.join(root,f)),h);
const meta=await read(path.join(assets,'model.json'));
const arrays=await Promise.all(['tile','saip-grid','hetero-grid'].map(async n=>{
  const b=zlib.gunzipSync(await fs.readFile(path.join(assets,n+'.bin.gz')));
  return new Float64Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength));
}));
const fixtures=(await read(path.join(root,'outputs/wasm-preparation/force-fixtures.json'))).filter(f=>!f.continuum&&f.mode===1);
const report={status:'running',cases:[],sources:{...gate.sources},elements:[]},engines=[];
for(const f of ['wasm/tools/validate_surface_elements.mjs','outputs/wasm-preparation/force-fixtures.json'])report.sources[f]=await digest(path.join(root,f));
try {
  for(const dir of [path.join(out,'reference'),root]) {
    const {default:create}=await import(pathToFileURL(path.join(dir,'wasm/dist/browser/runtime/molecular-car.mjs')));
    const module=await create({print:()=>{},printErr:()=>{}});
    engines.push(new module.BrowserEngine(meta,...arrays));
  }
  const elements=new Set();let filename;
  for(const fixture of fixtures) {
    if(filename!==fixture.filename) {
      filename=fixture.filename;const file=path.join(assets,filename),content=await fs.readFile(file,'utf8');
      report.sources[path.relative(root,file).replaceAll('\\','/')]=await digest(file);
      for(const e of engines){const info=e.prepareChemistry(content,filename);for(const a of info.molecule.atoms)elements.add(a.element);}
    }
    for(const scale of [0,1,1.5,2,12]) {
      for(const e of engines)e.update({lattice_scale:scale});
      const [old,current]=engines.map(e=>e.energyForces(fixture.positions,1,.75,false));
      const energy=Math.abs(old.surface_energy-current.surface_energy);
      const force=Math.max(...old.surface_forces.map((v,i)=>Math.abs(v-current.surface_forces[i])));
      assert(energy<1e-7&&force<1e-7,filename+' '+scale);
      assert.equal(old.uff_energy,current.uff_energy);assert.deepEqual(old.uff_forces,current.uff_forces);
      report.cases.push({filename,scale,energy_error:energy,force_error:force});
    }
  }
  report.elements=[...elements].sort();assert.equal(report.elements.length,12);
  for(const [f,h] of Object.entries(report.sources))assert.equal(await digest(path.join(root,f)),h);
  report.status='passed';report.source_integrity='passed';
} catch(error){report.status='failed';report.error=error.stack||String(error);process.exitCode=1;}
finally{for(const e of engines){e.close();e.delete();}await fs.writeFile(path.join(out,'element-validation.json'),JSON.stringify(report,null,2)+'\n');}
console.log(JSON.stringify({status:report.status,cases:report.cases.length,elements:report.elements,error:report.error}));
