import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {fileURLToPath,pathToFileURL} from 'node:url';
export const root=fileURLToPath(new URL('../../',import.meta.url));
export const out=path.join(root,'outputs/wasm-fourier-validation');
export const read=async p=>JSON.parse(await fs.readFile(p,'utf8'));
export const hash=async p=>crypto.createHash('sha256').update(await fs.readFile(p)).digest('hex');
export async function initialize() {
  const assets=path.join(root,'wasm/dist/browser/assets'),meta=await read(path.join(assets,'model.json'));
  const load=async d=>{
    const bytes=zlib.gunzipSync(await fs.readFile(path.join(root,'wasm/dist/browser',d.path)));
    assert.equal(crypto.createHash('sha256').update(bytes).digest('hex'),d.sha256);
    return new Float64Array(bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength));
  };
  const arrays=await Promise.all(['tile','saip-grid','hetero-grid'].map(k=>load(meta.surface.arrays[k])));
  const {default:create}=await import(pathToFileURL(path.join(root,'wasm/dist/browser/runtime/molecular-car.mjs')));
  const module=await create({print:()=>{},printErr:()=>{}});
  const engine=new module.BrowserEngine(meta,...arrays);
  const proxy=await read(path.join(assets,'fourier.json'));
  assert.equal(await hash(path.join(assets,'model.json')),proxy.base_model_sha256);
  engine.installFourier(proxy);engine.loadFourierSpecies(0,await load(proxy.arrays[0]));
  return {meta,module,engine,proxy};
}
export async function sources() {
  const build=await read(path.join(root,'wasm/dist/browser/runtime/build.json'));
  const result={};
  for(const [f,h] of Object.entries(build.own_sources))assert.equal(await hash(path.join(root,f)),h,'stale build '+f);
  for(const f of [...Object.keys(build.own_sources),'wasm/dist/browser/runtime/molecular-car.wasm',
    'wasm/dist/browser/runtime/molecular-car.mjs','wasm/dist/browser/assets/model.json',
    'wasm/dist/browser/assets/fourier.json','wasm/dist/browser/assets/fourier-0.bin.gz',
    'wasm/tools/fourier_support.mjs'])result[f]=await hash(path.join(root,f));
  return result;
}
export async function checkSources(s) {for(const [f,h] of Object.entries(s))assert.equal(await hash(path.join(root,f)),h,'changed '+f);}
