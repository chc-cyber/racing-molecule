import {loadJSON,loadTable,tableDownloadBytes} from './assets.mjs';
let requestedHetero=false,active=false,completedHetero=false,commonReady=false;
onmessage=({data})=>{if(data.hetero){requestedHetero=true;run();}};
async function run(){
 if(active)return;active=true;
 try{
  const cache=await caches.open('racing-molecule-tables-v1'),probe=new URL('./cache-probe',import.meta.url).href;
  await cache.put(probe,new Response('ok'));await cache.delete(probe);
  const model=await loadJSON(new URL('./assets/model.json',import.meta.url));
  const atlas=await loadJSON(new URL('./assets/fourier-atlas.json',import.meta.url));
  if(atlas.parameter_sha256!==model.surface.parameter_sha256)throw Error('Asset version mismatch');
  const arrays=model.surface.arrays;
  do{
   const withHetero=requestedHetero&&!completedHetero;
   const list=[...(!commonReady?[arrays.tile,arrays['saip-grid'],atlas.arrays[0]]:[]),...(withHetero?[arrays['hetero-grid']]:[])];
   const total=list.reduce((n,d)=>n+tableDownloadBytes(d),0);let done=0;
   for(const descriptor of list){
    await loadTable(descriptor,{base:import.meta.url,priority:'low',requireCache:true,onProgress:n=>postMessage({status:'loading',percent:Math.min(100,Math.floor((done+n)/total*100))})});
    done+=tableDownloadBytes(descriptor);await new Promise(r=>setTimeout(r,100));
   }
   commonReady=true;if(withHetero)completedHetero=true;
  }while(requestedHetero&&!completedHetero);
  postMessage({status:'ready',percent:100});
 }catch(error){postMessage({status:'unavailable',error:error.message});}
 finally{active=false;}
}
run();
