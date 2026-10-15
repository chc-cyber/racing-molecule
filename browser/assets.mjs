// Cache only checksum-verified tables; a cancelled worker cannot publish partial data.
export async function loadTable(descriptor, {base=import.meta.url,onProgress=()=>{},fetcher=fetch,cacheStorage=globalThis.caches,stallMs=45000}={}) {
  const url=new URL('./'+descriptor.path,base),key=new URL(url);key.searchParams.set('sha256',descriptor.sha256);
  let cache;
  try{cache=await cacheStorage?.open('racing-molecule-tables-v1');}catch{}
  const decode=async bytes=>{
    const raw=await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
    if(raw.byteLength!==descriptor.bytes)throw Error('赛道资源大小不正确，请重试。');
    const hash=await crypto.subtle.digest('SHA-256',raw);
    if(Array.from(new Uint8Array(hash),n=>n.toString(16).padStart(2,'0')).join('')!==descriptor.sha256)throw Error('赛道资源校验失败，请重试。');
    return new Float64Array(raw);
  };
  try{const saved=await cache?.match(key.href);if(saved){const data=await decode(await saved.arrayBuffer());onProgress(descriptor.download_bytes,true);return data;}}catch{try{await cache?.delete(key.href);}catch{}}
  for(let attempt=0;attempt<2;attempt++){
    const controller=new AbortController();let timer;
    const arm=()=>{clearTimeout(timer);timer=setTimeout(()=>controller.abort(),stallMs);};
    try{
      onProgress(0,false);arm();const response=await fetcher(url,{signal:controller.signal});
      if(!response.ok)throw Error('赛道资源下载失败（'+response.status+'），请重试。');
      const reader=response.body.getReader(),chunks=[];let loaded=0;
      for(;;){const {done,value}=await reader.read();if(done)break;chunks.push(value);loaded+=value.byteLength;onProgress(loaded,false);arm();}
      clearTimeout(timer);const bytes=await new Blob(chunks).arrayBuffer(),data=await decode(bytes);
      try{await cache?.put(key.href,new Response(bytes));}catch{}
      return data;
    }catch(error){if(attempt===1)throw Error(controller.signal.aborted?'赛道资源下载中断，请检查网络后重试。':error.message||'赛道资源下载失败，请重试。');}
    finally{clearTimeout(timer);}
  }
}

export async function loadJSON(url){
  const response=await fetch(url,{signal:AbortSignal.timeout(45000)});
  if(!response.ok)throw Error('赛道参数下载失败，请重试。');return response.json();
}
