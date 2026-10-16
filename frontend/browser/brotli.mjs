let initialized;
export async function decodeBrotli(bytes){
  initialized??=(async()=>{
    const {default:init,decompress}=await import('./vendor/brotli/brotli_wasm.js');
    const url=new URL('./vendor/brotli/brotli_wasm_bg.wasm.gz',import.meta.url);
    url.searchParams.set('sha256','904ef2045c64aeb1cc15e2df6ab0f6fdd79dc05fba7512be339fc519e06ba1c1');
    const read=async()=>{
      let cache;try{cache=await caches.open('racing-molecule-tables-v1');}catch{}
      let response=await cache?.match(url.href),packed=response?await response.arrayBuffer():null;
      const valid=async b=>b&&Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',b)),n=>n.toString(16).padStart(2,'0')).join('')==='904ef2045c64aeb1cc15e2df6ab0f6fdd79dc05fba7512be339fc519e06ba1c1';
      if(!await valid(packed)){
        response=await fetch(url,{cache:'force-cache',signal:AbortSignal.timeout(45000)});
        if(!response.ok)throw Error('解码器下载失败，请重试。');packed=await response.arrayBuffer();
        if(!await valid(packed))throw Error('解码器校验失败，请重试。');
        try{await cache?.put(url.href,new Response(packed));}catch{}
      }
      return new Response(new Blob([packed]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
    };
    const binary=await (globalThis.navigator?.locks?navigator.locks.request('race-brotli-decoder',read):read());
    await init(binary);return decompress;
  })().catch(error=>{initialized=null;throw error;});
  const decode=await initialized,result=decode(new Uint8Array(bytes));
  return result.buffer.slice(result.byteOffset,result.byteOffset+result.byteLength);
}
