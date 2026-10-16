import test from 'node:test';
import assert from 'node:assert/strict';
import {gzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {loadTable} from '../frontend/browser/assets.mjs';
const raw=new Float64Array([1,2,3]),bytes=Buffer.from(raw.buffer),packed=gzipSync(bytes);
const descriptor={path:'assets/test.bin.gz',bytes:bytes.length,download_bytes:packed.length,sha256:createHash('sha256').update(bytes).digest('hex')};
const makeCache=()=>{const data=new Map();return{data,open:async()=>({match:async k=>data.get(k)?.clone(),put:async(k,r)=>data.set(k,r.clone()),delete:async k=>data.delete(k)})};};
test('verified cached tables survive new loaders without another network request',async()=>{
 const cacheStorage=makeCache();let calls=0;const fetcher=async()=>{calls++;return new Response(packed);};
 const options={base:'https://game.invalid/browser/worker.mjs',cacheStorage,fetcher};
 assert.deepEqual(await loadTable(descriptor,options),raw);assert.deepEqual(await loadTable(descriptor,options),raw);assert.equal(calls,1);
});
test('corrupt cache and bad download are never accepted',async()=>{
 const cacheStorage=makeCache();cacheStorage.data.set('https://game.invalid/browser/assets/test.bin.gz?sha256='+descriptor.sha256,new Response(gzipSync('bad')));
 let calls=0;const result=await loadTable(descriptor,{base:'https://game.invalid/browser/worker.mjs',cacheStorage,fetcher:async()=>new Response(++calls===1?gzipSync('bad'):packed)});
 assert.deepEqual(result,raw);assert.equal(calls,2);
});
test('blocked cache storage does not prevent preparation',async()=>{
 assert.deepEqual(await loadTable(descriptor,{base:'https://game.invalid/',cacheStorage:{open:async()=>{throw Error('quota');}},fetcher:async()=>new Response(packed)}),raw);
});
test('stalled downloads abort and report a retryable download error',async()=>{
 let calls=0;await assert.rejects(loadTable(descriptor,{base:'https://game.invalid/',cacheStorage:null,stallMs:10,fetcher:async(_,{signal})=>{calls++;return new Promise((_,reject)=>signal.addEventListener('abort',()=>reject(Error('aborted'))));}}),/下载中断/);assert.equal(calls,2);
});
test('active slow streaming is not treated as a fixed total timeout',async()=>{
 const progress=[];const result=await loadTable(descriptor,{base:'https://game.invalid/',cacheStorage:null,stallMs:100,onProgress:n=>progress.push(n),fetcher:async()=>new Response(new ReadableStream({async start(c){for(let i=0;i<packed.length;i+=4){await new Promise(r=>setTimeout(r,20));c.enqueue(packed.subarray(i,i+4));}c.close();}}))});
 assert.deepEqual(result,raw);assert.equal(progress.at(-1),packed.length);assert(progress.length>2);
});

test('overlapping preload and race share one verified download under the cross-worker lock',async()=>{
 const cacheStorage=makeCache();let tail=Promise.resolve(),calls=0;
 const locks={request(_key,fn){const next=tail.then(fn);tail=next.catch(()=>{});return next;}};
 const options={base:'https://game.invalid/browser/worker.mjs',cacheStorage,locks,fetcher:async()=>{calls++;await new Promise(r=>setTimeout(r,20));return new Response(packed);}};
 const results=await Promise.all([loadTable(descriptor,{...options,priority:'low'}),loadTable(descriptor,options)]);
 assert.equal(calls,1);assert.deepEqual(results,[raw,raw]);
});

test('interrupted chunk download resumes completed chunks without changing float64 bytes',async()=>{
 const cut=12,parts=[packed.subarray(0,cut),packed.subarray(cut)];
 const chunks=parts.map((b,i)=>({path:`assets/chunks/${i}.bin`,bytes:b.length,download_bytes:b.length,sha256:createHash('sha256').update(b).digest('hex')}));
 const cacheStorage=makeCache(),calls=[0,0];let fail=true;
 const options={base:'https://game.invalid/browser/worker.mjs',cacheStorage,fetcher:async url=>{const i=Number(url.pathname.match(/(\d)\.bin/)[1]);calls[i]++;if(i===1&&fail)return new Response('failure',{status:503});return new Response(parts[i]);}};
 await assert.rejects(loadTable({...descriptor,chunks},options),/503/);
 fail=false;const result=await loadTable({...descriptor,chunks},options);
 assert.deepEqual(result,raw);assert.deepEqual(calls,[1,3]);
});

test('chunk manifests retain previously cached full tables',async()=>{
 const cacheStorage=makeCache(),options={base:'https://game.invalid/browser/worker.mjs',cacheStorage,fetcher:async()=>new Response(packed)};
 await loadTable(descriptor,options);
 const value=await loadTable({...descriptor,brotli:{path:'must-not-fetch.br',download_bytes:12,chunks:[]},chunks:[{path:'must-not-fetch.bin'}]},{...options,fetcher:async()=>{throw Error('unexpected download');}});
 assert.deepEqual(value,raw);
});
