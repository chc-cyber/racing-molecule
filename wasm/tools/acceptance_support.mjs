import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {fileURLToPath,pathToFileURL} from 'node:url';
export const root=fileURLToPath(new URL('../../',import.meta.url));
export const folder=path.join(root,'outputs/default-candidate-acceptance');
export const browser=path.resolve(path.join(root,'wasm/dist/browser'));
export const read=async file=>JSON.parse(await fs.readFile(file,'utf8'));
export const hash=async file=>crypto.createHash('sha256').update(await fs.readFile(file)).digest('hex');
export const save=async(file,value)=>{await fs.mkdir(path.dirname(file),{recursive:true});await fs.writeFile(file,JSON.stringify(value,null,2)+'\n');};
export const cars={
  'test-car-52':'outputs/wasm-fourier-validation/races/inputs/test-car-52.mol',
  'pyridine-11':'outputs/wasm-fourier-validation/races/inputs/pyridine-11.mol',
  'fullerene-c60':'outputs/wasm-fourier-validation/races/inputs/fullerene-c60.mol',
  'folded-eicosane':'outputs/md-release-validation/inputs/folded-eicosane.mol',
};
export async function initialize(seed=2718281){
  const meta=await read(path.join(browser,'assets/model.json'));
  const load=async d=>{const b=zlib.gunzipSync(await fs.readFile(path.join(browser,d.path)));assert.equal(crypto.createHash('sha256').update(b).digest('hex'),d.sha256);return new Float64Array(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength));};
  const arrays=await Promise.all(['tile','saip-grid','hetero-grid'].map(k=>load(meta.surface.arrays[k])));
  const {default:create}=await import(pathToFileURL(path.join(browser,'runtime/molecular-car.mjs')));
  const module=await create({print:()=>{},printErr:()=>{}}),engine=new module.BrowserEngine({...meta,seed},...arrays);
  const proxy=await read(path.join(browser,'assets/fourier-atlas.json'));
  assert.equal(await hash(path.join(browser,'assets/model.json')),proxy.base_model_sha256);
  engine.installFourier(proxy);engine.loadFourierSpecies(0,await load(proxy.arrays[0]));
  return {meta,module,engine};
}
export async function sources(){
  const files=['runtime/molecular-car.mjs','runtime/molecular-car.wasm','runtime/build.json','assets/model.json','assets/fourier-atlas.json','assets/fourier-atlas-0.bin.gz'];
  const build=await read(path.join(browser,'runtime/build.json'));for(const [f,h] of Object.entries(build.own_sources))assert.equal(await hash(path.join(root,f)),h,'stale build '+f);const meta=await read(path.join(browser,'assets/model.json'));files.push(...Object.values(meta.surface.arrays).map(d=>d.path));const result={};for(const file of files)result[path.join(browser,file)]=await hash(path.join(browser,file));
  for(const file of [...Object.values(cars),'wasm/tools/acceptance_support.mjs'])result[path.join(root,file)]=await hash(path.join(root,file));return result;
}
export async function checkSources(s){for(const [file,h] of Object.entries(s))assert.equal(await hash(file),h,'source changed: '+file);}
export const stepsFor=fs=>dt=>{const n=fs/dt;assert.equal(n,Math.round(n));return n;};
export function geometry(info,initial){
  const bonds=info.molecule.bonds.map(b=>[b.a-1,b.b-1]);
  const length=(p,b)=>Math.hypot(...[0,1,2].map(a=>p[3*b[0]+a]-p[3*b[1]+a]));
  const reference=bonds.map(b=>length(initial,b));
  return frame=>{
    assert(Number.isFinite(frame.energy)&&frame.positions.every(Number.isFinite),'nonfinite trajectory');
    const ratios=bonds.map((b,i)=>length(frame.positions,b)/reference[i]);
    const z=frame.positions.filter((_,i)=>i%3===2);
    return {bond_min:Math.min(...ratios),bond_max:Math.max(...ratios),min_z_A:Math.min(...z),max_z_A:Math.max(...z)};
  };
}
