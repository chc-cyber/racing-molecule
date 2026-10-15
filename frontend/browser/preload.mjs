import {loadJSON,loadTable} from './assets.mjs';

// No engine or molecular preparation runs in this worker.
try {
  const cache=await caches.open('racing-molecule-tables-v1');
  const probe=new URL('./cache-probe',import.meta.url).href;
  await cache.put(probe,new Response('ok'));await cache.delete(probe);
  postMessage({status:'loading'});
  const model=await loadJSON(new URL('./assets/model.json',import.meta.url));
  const atlas=await loadJSON(new URL('./assets/fourier-atlas.json',import.meta.url));
  if(atlas.parameter_sha256!==model.surface.parameter_sha256)throw Error('Asset version mismatch');
  const arrays=model.surface.arrays;
  for(const descriptor of [arrays.tile,arrays['saip-grid'],arrays['hetero-grid'],atlas.arrays[0]]) {
    await loadTable(descriptor,{base:import.meta.url,priority:'low',requireCache:true});
    await new Promise(resolve=>setTimeout(resolve,100));
  }
  postMessage({status:'ready'});
} catch {
  postMessage({status:'unavailable'});
}
close();
