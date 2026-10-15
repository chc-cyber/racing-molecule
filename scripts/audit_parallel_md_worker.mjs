// Independent scheduling checks for the actual production Worker source.
// Only the imported MD module and network are replaced; worker control flow runs unchanged.
import fs from 'node:fs/promises';
import path from 'node:path';
import vm from 'node:vm';
import assert from 'node:assert/strict';
import crypto from 'node:crypto';
import zlib from 'node:zlib';
import {fileURLToPath} from 'node:url';

const root=fileURLToPath(new URL('../',import.meta.url));
const workerFile=path.join(root,'frontend/browser/worker.mjs');
const original=await fs.readFile(workerFile,'utf8');
assert.match(original,/^import createMolecularCar from /);
const source=original.replace(/^import createMolecularCar from .*;\r?\n/,'')
  .replaceAll('import.meta.url',JSON.stringify('https://test.invalid/browser/worker.mjs'));
const bytes=Buffer.alloc(8),packed=zlib.gzipSync(bytes);
const descriptor={path:'assets/test.bin.gz',bytes:8,sha256:crypto.createHash('sha256').update(bytes).digest('hex')};
const report={status:'running',scope:'production Worker scheduling with mocked MD and assets; not a dynamics validation',
  source_sha256:crypto.createHash('sha256').update(original).digest('hex'),cases:[]};

async function scenario(dt) {
  let step=0,pulseUntil=-1,nextTimer=0;
  const messages=[],queue=new Map(),calls=[];
  const metadata={timestep_fs:.5,frame_steps:100,max_steps:10000000,track:{finish_x:200},
    surface:{arrays:{tile:descriptor,'saip-grid':descriptor,'hetero-grid':descriptor}}};
  const frame=()=>({step,time_fs:step*dt,center:[0,0,3],pulse_active:step<pulseUntil});
  class Engine {
    prepareWithPose(){return {molecule:{},runtime:{timestep_fs:dt,start_pose:{requested:{mode:'auto'}}}};}
    controls(){return {};}
    update(){}
    frame(){return frame();}
    advance(n){calls.push(n);step+=n;return frame();}
    pulse(){pulseUntil=step+400/dt;}
    refresh(){return frame();}
  }
  const context=vm.createContext({URL,Response,DecompressionStream,Uint8Array,Float64Array,Number,Math,
    crypto:crypto.webcrypto,performance,console,onmessage:null,
    createMolecularCar:async()=>({BrowserEngine:Engine}),
    postMessage:m=>messages.push(m),
    setTimeout:fn=>{const id=++nextTimer;queue.set(id,fn);return id;},
    clearTimeout:id=>queue.delete(id),
    fetch:async url=>String(url).endsWith('model.json')?{ok:true,json:async()=>metadata}:new Response(packed)});
  vm.runInContext(source,context,{filename:workerFile});
  const action=async action=>context.onmessage({data:{kind:'action',action}});
  const run=()=>{const entry=queue.entries().next().value;assert(entry,'Expected scheduled tick');queue.delete(entry[0]);entry[1]();};
  await context.onmessage({data:{kind:'load',settings:{surface:'periodic-au'},content:'fixture',filename:'fixture.mol',drive:{},start_pose:{mode:'auto'}}});
  assert.equal(messages.at(-1).kind,'ready');
  await action('start');run();
  const first=messages.find(m=>m.kind==='frame');
  assert.equal(calls[0],50/dt);assert.equal(first.frame.time_fs,50);
  await action('pause');assert.equal(queue.size,0);assert.equal(messages.at(-1).status,'paused');
  await action('start');run();assert.equal(messages.filter(m=>m.kind==='frame').at(-1).frame.time_fs,100);
  await action('pause');
  // Time limit is the metadata physical horizon, independent of MD step count.
  step=(metadata.max_steps*metadata.timestep_fs-50)/dt;
  await action('start');run();
  assert.equal(messages.at(-1).status,'finished');assert.equal(messages.at(-1).reason,'time_limit');
  assert.equal(messages.filter(m=>m.kind==='frame').at(-1).frame.time_fs,5000000);
  assert.equal(queue.size,0);
  // Pulse-only pause follows actual pulse_active, not a fixed number of steps.
  step=0;await action('pulse');for(let i=0;i<8;i++)run();
  assert.equal(messages.at(-1).status,'paused');assert.equal(queue.size,0);
  assert.equal(messages.filter(m=>m.kind==='frame').at(-1).frame.time_fs,400);
  assert(!messages.some(m=>m.kind==='error'));
  return {dt_fs:dt,frame_steps:calls[0],sample_fs:50,max_time_fs:5000000,pulse_only_fs:400,start_pause_resume:true};
}
try {for(const dt of [.25,.5,1,2])report.cases.push(await scenario(dt));report.status='passed';}
catch(error){report.status='failed';report.error=error.stack;process.exitCode=1;}
await fs.mkdir(path.join(root,'outputs/parallel-md-validation'),{recursive:true});
await fs.writeFile(path.join(root,'outputs/parallel-md-validation/worker-scheduling.json'),JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify(report));
