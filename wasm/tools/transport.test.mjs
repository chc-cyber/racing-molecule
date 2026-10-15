import assert from 'node:assert/strict';
import {test} from 'node:test';
import {BrowserTransport} from '../../frontend/browser/transport.mjs';

globalThis.document={addEventListener(){}, hidden:false};
globalThis.addEventListener=()=>{};
globalThis.fetch=async()=>({ok:true, text:async()=> 'test molecule'});
class Worker {
  constructor(){this.sent=[];this.terminated=false;}
  postMessage(data){this.sent.push(data);}
  terminate(){this.terminated=true;}
}
globalThis.Worker=Worker;
const ready = async (competition={mode:"custom"}) => {
  const client=new BrowserTransport();await client.load({example:"benzene.mol",competition});
  client.worker.onmessage({data:{kind:'ready',data:{frame:{seq:0,step:0,positions:[0,0,3]}}}});
  return client;
};

test('Fourier selection replaces the worker and survives reset',async()=> {
  const client=await ready();
  await client.load({reuse:true,competition:{mode:'custom',surface:'periodic-au'}});
  client.worker.onmessage({data:{kind:'ready',data:{frame:{seq:0,step:0,positions:[0,0,3]}}}});
  const old=client.worker;
  await client.load({reuse:true,competition:{mode:'custom',surface:'fourier-au'}});
  assert(old.terminated);assert.equal(client.worker.sent[0].settings.surface,'fourier-au');
  client.worker.onmessage({data:{kind:'ready',data:{frame:{seq:0,step:0,positions:[0,0,3]}}}});
  await client.request('/api/action',{action:'reset'});
  assert.equal(client.worker.sent[0].settings.surface,'fourier-au');
});

test('invalid controls leave requested settings and worker untouched', async()=> {
  const client=await ready();const before=JSON.stringify({settings:client.settings,drive:client.drive});
  for(const changes of [{temperature:NaN},{temperature:true},{drive:{mode:''}}, {tip:{z:0}},
    {tip:{unknown:1}}, {drive:{strength:.5},tip:{strength:.4}}, {surface:'continuum',temperature:600},
    {autopilot:'false'}, {surface:'invalid'}, {corrugation:.1}, {other:2},
    {lattice_scale:NaN},{lattice_scale:true},{lattice_scale:13},{drive:{growth:0}},{drive:{growth:true}}]) {
    await assert.rejects(client.request('/api/control',changes));
    assert.equal(JSON.stringify({settings:client.settings,drive:client.drive}),before);
  }
  assert.equal(client.worker.sent.length,1);
});

test('crystal and field controls survive a reset with a fresh generation',async()=> {
  const client=await ready(),generation=client.generation;
  await client.load({reuse:true,competition:{mode:'custom',lattice_scale:2,growth:.25}});
  await client.request('/api/action',{action:'reset'});
  assert.equal(client.settings.lattice_scale,2);assert.equal(client.drive.growth,.25);
  assert(client.generation>generation);
  const message=client.worker.sent[0];assert.equal(message.settings.lattice_scale,2);assert.equal(message.drive.growth,.25);
});

test('a stale example download cannot replace the latest imported molecule', async()=> {
  const original=globalThis.fetch;
  try {
    let finish;globalThis.fetch=()=>new Promise(resolve=>{finish=resolve;});
    const client=new BrowserTransport();const pending=client.load({example:'benzene.mol'});
    await client.load({filename:'new.mol',content:'new structure'});
    const worker=client.worker,generation=client.generation;
    finish({ok:true,text:async()=> 'old example'});await pending;
    assert.equal(client.worker,worker);assert.equal(client.generation,generation);
    assert.equal(client.input.filename,'new.mol');assert(!worker.terminated);
  } finally {globalThis.fetch=original;}
});

test('replay retention is bounded and truncation is explicit', async()=> {
  const client=await ready();
  for(let seq=1;seq<=5000;seq++)client.appendFrame({seq,step:seq*100,positions:new Array(768).fill(0)});
  assert(client.frames.length<4000);assert(client.droppedFrames>0);
  const replay=await client.replay();assert.equal(replay.replay_complete,false);
  assert.equal(replay.first_retained_step,client.frames[0].step);
  const stream=await client.request('/api/frames?after=-1');assert.equal(stream.dropped,true);assert.equal(stream.frames.length,250);
});

test('invalid imports stop preparation and remain recoverable', async()=> {
  const client=await ready();const previous=client.worker;
  await assert.rejects(client.load({filename:'bad.mol',content:''}));
  assert.equal(client.state.status,'error');assert(previous.terminated);
  await client.load({filename:'valid.mol',content:'another molecule'});
  assert.equal(client.state.status,'preparing');assert.equal(client.input.filename,'valid.mol');
});

test('a finished race permits a surface change and starts fresh preparation', async()=> {
  const client=await ready();client.state.status='finished';const generation=client.generation,worker=client.worker;
  await assert.rejects(client.request('/api/control',{temperature:500}));
  await client.load({reuse:true,competition:{mode:'custom',surface:'continuum'}});
  assert.equal(client.state.status,'preparing');assert.equal(client.generation,generation+1);
  assert.equal(client.settings.surface,'continuum');assert(worker.terminated);
});

test('pose applies by a fresh generation, survives reset, and is recorded in replay',async()=>{
  const client=await ready(),old=client.worker,generation=client.generation;
  await client.request('/api/start-pose',{mode:'manual',yaw_deg:90});
  assert(old.terminated);assert.equal(client.generation,generation+1);
  assert.equal(client.worker.sent[0].start_pose.yaw_deg,90);
  await client.request('/api/action',{action:'reset'});
  assert.equal(client.worker.sent[0].start_pose.yaw_deg,90);
  const replay=await client.replay();assert.equal(replay.start_pose.yaw_deg,90);
  await client.request('/api/load',{content:'new structure',filename:'new.mol'});
  assert.equal(client.worker.sent[0].start_pose.mode,'auto');assert.equal(client.startPose.yaw_deg,0);
});

test('replay retains finish rule, per-frame origin and finish reason', async()=>{
  const client=await ready();
  const race={rule:'longitudinal-displacement-v1',start_x:.7,finish_x:200.7,length_angstrom:200};
  client.state.runtime={race};
  client.frames=[{seq:100,step:5000,race_start_x:.7,race_length_angstrom:200,distance_angstrom:200,progress:1}];
  client.droppedFrames=100;
  client.worker.onmessage({data:{kind:'status',status:'finished',reason:'finish'}});
  const replay=await client.replay();
  assert.equal(replay.finish_reason,'finish');assert.equal(replay.runtime.race.finish_x,200.7);
  assert.equal(replay.frames[0].race_start_x,.7);assert.equal(replay.replay_complete,false);
});

test('invalid or running pose request leaves the worker and current race intact',async()=>{
  const client=await ready(),worker=client.worker,generation=client.generation;
  for(const pose of [{mode:'wrong'},{mode:'manual',yaw_deg:NaN},{axis_atoms:[1,1]},{other:1}])
    await assert.rejects(client.request('/api/start-pose',pose));
  client.state.status='running';await assert.rejects(client.request('/api/start-pose',{mode:'auto'}));
  assert.equal(client.worker,worker);assert.equal(client.generation,generation);assert(!worker.terminated);
});


test('playback stream pages chronologically and only grants current-generation credit', async()=>{
  const client=await ready();
  for(let seq=1;seq<600;seq++)client.appendFrame({seq,step:seq*50,positions:[0,0,3]});
  const first=await client.request(`/api/frames?after=-1&played=50&generation=${client.generation}`);
  assert.deepEqual(first.frames.map(f=>f.seq),Array.from({length:250},(_,i)=>i));
  assert.equal(first.dropped,false);
  assert.deepEqual(client.worker.sent.at(-1),{kind:'playback',until_fs:10050});
  const count=client.worker.sent.length;
  const second=await client.request(`/api/frames?after=249&played=100&generation=${client.generation-1}`);
  assert.equal(second.frames[0].seq,250);assert.equal(client.worker.sent.length,count);
});
