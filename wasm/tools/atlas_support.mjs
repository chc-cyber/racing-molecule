import fs from 'node:fs/promises';
import path from 'node:path';
import zlib from 'node:zlib';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {fileURLToPath,pathToFileURL} from 'node:url';
export const root=fileURLToPath(new URL('../../',import.meta.url));
export const out=path.join(root,'outputs/c60-proxy-study');
export const read=async p=>JSON.parse(await fs.readFile(p,'utf8'));
export const hash=async p=>crypto.createHash('sha256').update(await fs.readFile(p)).digest('hex');
export async function initialize() {
  const browser=path.join(root,'wasm/dist/browser'),meta=await read(path.join(browser,'assets/model.json'));
  const load=async(base,d)=>{
    const bytes=zlib.gunzipSync(await fs.readFile(path.join(base,d.path)));
    assert.equal(bytes.byteLength,d.bytes);
    assert.equal(crypto.createHash('sha256').update(bytes).digest('hex'),d.sha256);
    return new Float64Array(bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength));
  };
  const arrays=await Promise.all(['tile','saip-grid','hetero-grid'].map(k=>load(browser,meta.surface.arrays[k])));
  const {default:create}=await import(pathToFileURL(path.join(browser,'runtime/molecular-car.mjs')));
  const module=await create({print:()=>{},printErr:()=>{}});
  const engine=new module.BrowserEngine(meta,...arrays);
  const proxy=await read(path.join(out,'assets/fourier-atlas.json'));
  assert.equal(await hash(path.join(browser,'assets/model.json')),proxy.base_model_sha256);
  engine.installFourier(proxy);engine.loadFourierSpecies(0,await load(out,proxy.arrays[0]));
  return {meta,module,engine,proxy};
}
export async function sources() {
  const build=await read(path.join(root,'wasm/dist/browser/runtime/build.json'));
  const result={};
  for(const [f,h] of Object.entries(build.own_sources))assert.equal(await hash(path.join(root,f)),h,'stale build '+f);
  for(const f of [...Object.keys(build.own_sources),'wasm/dist/browser/runtime/molecular-car.wasm',
    'wasm/dist/browser/runtime/molecular-car.mjs','wasm/dist/browser/assets/model.json',
    'outputs/c60-proxy-study/assets/fourier-atlas.json','outputs/c60-proxy-study/assets/fourier-atlas-0.bin.gz',
    'wasm/tools/atlas_support.mjs'])result[f]=await hash(path.join(root,f));
  return result;
}
export async function checkSources(s) {for(const [f,h] of Object.entries(s))assert.equal(await hash(path.join(root,f)),h,'changed '+f);}
