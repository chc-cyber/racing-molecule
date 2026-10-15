import test from 'node:test';
import assert from 'node:assert/strict';
import {RaceCompletion} from '../frontend/editor/race-completion.js';

const deferred=()=>{let resolve;const promise=new Promise(r=>resolve=r);return {promise,resolve};};
const network=()=>Object.assign(new Error('连接已中断'),{recoverable:true});
function setup(overrides={}){
  const calls=[],states=[],project={name:'car',atoms:[{id:'stable',pos3:[1,2,3]}]};
  const deps={
    uuid:()=> 'test-id', settle:async()=>structuredClone(project),
    save:async p=>calls.push(['save',p]), export:async p=>{calls.push(['export',p]);return{filename:'car.mol',content:'MOL'};},
    prepare:async j=>{calls.push(['prepare',j]);return {status:'checking'};},
    status:async id=>{calls.push(['status',id]);return {status:'ready'};},
    commit:async id=>{calls.push(['commit',id]);return {status:'committed',generation:8};},
    active:async()=>({generation:8,status:'paused',molecule:{atoms:[]},frame:{step:0}}),
    cancel:async id=>{calls.push(['cancel',id]);return {status:'cancelled'};},
    remember:j=>calls.push(['remember',structuredClone(j)]),delay:async()=>{},
    notify:s=>states.push(s),navigate:g=>calls.push(['navigate',g]),...overrides};
  return {flow:new RaceCompletion(deps),calls,states,project,deps};
}
test('saves exact design before preparing once and navigates only after real active generation',async()=>{
  const {flow,calls,states,project}=setup();await flow.start();
  assert.deepEqual(calls.map(c=>c[0]),['save','export','remember','prepare','status','commit','remember','navigate']);
  assert.deepEqual(calls[0][1],project);assert.deepEqual(calls[1][1],project);
  assert.equal(flow.state,'success');assert.equal(states.at(-1).message,'赛车已就绪。');
});
test('rapid double click shares one operation before editor settling',async()=>{
  const wait=deferred();const {flow,calls}=setup({settle:()=>wait.promise});
  const first=flow.start();await flow.start();wait.resolve({name:'car'});await first;
  assert.equal(calls.filter(c=>c[0]==='prepare').length,1);
});

test('chosen race rules are captured before asynchronous export and retained for retries',async()=>{
  const wait=deferred(),rules={mode:'custom',temperature:300};
  const {flow,calls}=setup({competition:()=>rules,settle:()=>wait.promise});
  const start=flow.start();rules.temperature=600;wait.resolve({name:'car'});await start;
  assert.equal(calls.find(c=>c[0]==='prepare')[1].competition.temperature,300);
});
test('unapplied preview or save failure never reaches candidate preparation',async()=>{
  for(const stage of ['settle','save']){
    const {flow,calls}=setup({[stage]:async()=>{throw new Error('resolve design first');}});
    await flow.start();assert.equal(flow.state,'failed');assert(!calls.some(c=>c[0]==='prepare'));
    await flow.cancel();assert.equal(flow.state,'idle');
  }
});
test('admission/preparation failure exposes one return action and no promotion',async()=>{
  const {flow,calls,states}=setup({status:async()=>({status:'error',error:'请连接片段'})});
  await flow.start();assert.equal(flow.state,'failed');assert.equal(states.at(-1).retry,false);
  assert.equal(states.at(-1).message,'请连接片段');assert(!calls.some(c=>c[0]==='commit'));
  await flow.cancel();assert.equal(flow.state,'idle');
});
test('lost prepare response is recovered via same job without saving or preparing again',async()=>{
  const {flow,calls}=setup({prepare:async j=>{calls.push(['prepare',j]);throw network();}});
  await flow.start();assert.equal(flow.state,'failed');assert(flow.locked);
  await flow.retry();assert.equal(flow.state,'success');
  assert.equal(calls.filter(c=>c[0]==='prepare').length,1);assert.equal(calls.filter(c=>c[0]==='save').length,1);
});
test('if original prepare never arrived retry resends same id and data',async()=>{
  let n=0,requests=[];const {flow}=setup({prepare:async j=>{requests.push(j);if(!n++)throw network();return{status:'ready'};},status:async()=>{throw Object.assign(new Error('missing'),{missing:true});}});
  await flow.start();await flow.retry();assert.equal(flow.state,'success');assert.deepEqual(requests[0],requests[1]);
});
test('lost commit response recovers committed generation without another promotion',async()=>{
  const {flow,calls,deps}=setup({commit:async id=>{calls.push(['commit',id]);throw network();}});
  await flow.start();deps.status=async()=>({status:'committed',generation:8});await flow.retry();
  assert.equal(flow.state,'success');assert.equal(calls.filter(c=>c[0]==='commit').length,1);
});
test('reload recovery uses persisted job without another save or preparation',async()=>{
  const {flow,calls}=setup();await flow.restore({id:'restored',content:'MOL',filename:'a.mol'});
  assert.equal(flow.state,'success');assert(!calls.some(c=>['save','export','prepare'].includes(c[0])));
});
test('cancel during polling prevents late ready result from promoting',async()=>{
  const wait=deferred();let signal;const started=new Promise(r=>signal=r);
  const {flow,calls}=setup({status:async()=>{signal();return wait.promise;}});
  const run=flow.start();await started;await flow.cancel();wait.resolve({status:'ready'});await run;
  assert.equal(flow.state,'idle');assert(!calls.some(c=>c[0]==='commit'));
});
test('cancel during save prevents any preparation after late save completion',async()=>{
  const wait=deferred();let signal;const started=new Promise(r=>signal=r);
  const {flow,calls}=setup({save:async()=>{signal();return wait.promise;}});
  const run=flow.start();await started;await flow.cancel();wait.resolve();await run;
  assert.equal(flow.state,'idle');assert(!calls.some(c=>c[0]==='prepare'));
});
test('cancel is disabled during commit; ambiguous committed cancellation can be recovered',async()=>{
  const wait=deferred();let signal;const started=new Promise(r=>signal=r);
  const {flow,calls}=setup({commit:async()=>{signal();return wait.promise;}});
  const run=flow.start();await started;await flow.cancel();assert(!calls.some(c=>c[0]==='cancel'));
  wait.resolve({status:'committed',generation:8});await run;assert.equal(flow.state,'success');
});
test('stale active race never navigates; cancelled status never promotes',async()=>{
  for(const override of [{active:async()=>({generation:9,status:'paused',molecule:{},frame:{}})}, {status:async()=>({status:'cancelled'})}]){
    const {flow,calls}=setup(override);await flow.start();assert.equal(flow.state,'failed');assert(!calls.some(c=>c[0]==='navigate'));
  }
});
test('network failure during cancellation returns to editing immediately and retries cleanup before a new job',async()=>{
  let cancellations=0;const {flow,calls,deps}=setup({status:async()=>{throw network();},cancel:async id=>{calls.push(['cancel',id]);if(!cancellations++)throw network();return{status:'cancelled'};}});
  await flow.start();assert.equal(flow.state,'failed');await flow.cancel();assert.equal(flow.state,'idle');assert(!flow.locked);assert(flow.cleanupJob.cancelRequested);
  deps.status=async()=>({status:'ready'});await flow.start();assert.equal(flow.state,'success');assert.equal(cancellations,2);assert.equal(flow.cleanupJob,null);
  const secondPrepare=calls.map(c=>c[0]).lastIndexOf('prepare'),cleanup=calls.map(c=>c[0]).lastIndexOf('cancel');assert(cleanup<secondPrepare);
});
test('reload of a cancelled offline job only attempts cleanup and never promotes or navigates',async()=>{
  const {flow,calls}=setup({cancel:async()=>{throw network();}});await flow.restore({id:'cancelled',cancelRequested:true});
  assert.equal(flow.state,'idle');assert(!flow.locked);assert(!calls.some(c=>['commit','prepare','navigate'].includes(c[0])));
});
