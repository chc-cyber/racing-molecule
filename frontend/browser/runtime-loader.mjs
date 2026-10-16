import createModule from './runtime/molecular-car.mjs';
import {loadJSON,loadTable} from './assets.mjs';
export async function loadRuntime(onProgress=()=>{}){
  const descriptor=await loadJSON(new URL('./runtime/download.json',import.meta.url));
  const wasmBinary=await loadTable({...descriptor,encoding:'raw'},{base:import.meta.url,onProgress:n=>onProgress(Math.min(100,Math.floor(n/descriptor.bytes*100)))});
  const compiled=await WebAssembly.compile(wasmBinary);
  return createModule({instantiateWasm(imports,receive){const instance=new WebAssembly.Instance(compiled,imports);receive(instance,compiled);return instance.exports;},print:()=>{},printErr:()=>{}});
}
