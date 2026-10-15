import test from 'node:test';
import assert from 'node:assert/strict';
import {spawnSync} from 'node:child_process';
import {validateCompetition, engineParameters, isRaceInProgress} from '../frontend/competition.mjs';
import {BrowserTransport} from '../frontend/browser/transport.mjs';

test('Python and JS canonical rules agree across standard and custom cases',()=>{
  const cases=[{}, {mode:'custom',temperature:600,strength:.25,lattice_scale:2},
    {mode:'custom',surface:'continuum',corrugation:.3,growth:.5,start_pose:{mode:'manual',yaw_deg:90}},
    {mode:'custom',surface:'periodic-au',start_pose:{mode:'auto',axis_atoms:[1,5]}}];
  const result=spawnSync('.venv/Scripts/python.exe',['-X','utf8','-c',
    'import json,sys;from molecular_car.competition import validate_competition;print(json.dumps([validate_competition(x) for x in json.load(sys.stdin)]))'],{input:JSON.stringify(cases),encoding:'utf8'});
  assert.equal(result.status,0,result.stderr);
  assert.deepEqual(cases.map(validateCompetition),JSON.parse(result.stdout));
});
test('standard cannot silently accept altered rules or removed drive modes',()=>{
  for(const rules of [{strength:.7},{surface:'continuum'},{timestep_fs:.5},{mode:'custom',tip:{}},{mode:'custom',pulse:true},{length_nm:50},{mode:'custom',temperature:true},{mode:'custom',strength:NaN}])
    assert.throws(()=>validateCompetition(rules));
  const c=validateCompetition();assert.equal(engineParameters(c).drive.mode,'gradient');
});
test('browser begins idle; invalid rules and live controls preserve current race',async()=>{
  globalThis.document={addEventListener(){}};globalThis.addEventListener=()=>{};
  globalThis.Worker=class {constructor(){this.sent=[];}postMessage(m){this.sent.push(m);}terminate(){this.terminated=true;}};
  const client=new BrowserTransport();await client.initial;
  assert.equal(client.worker,null);assert.equal(client.state.status,'idle');
  await client.load({filename:'car.mol',content:'MOL'});const worker=client.worker;
  worker.onmessage({data:{kind:'ready',data:{frame:{seq:0,positions:[0,0,3]}}}});
  const before=structuredClone(client.state);
  for(const rules of [{temperature:600},{mode:'custom',strength:2}])await assert.rejects(client.load({reuse:true,competition:rules}));
  await assert.rejects(client.request('/api/action',{action:'pulse'}));
  await assert.rejects(client.request('/api/control',{drive:{strength:.9}}));
  assert.equal(client.worker,worker);assert(!worker.terminated);assert.deepEqual(client.state,before);
  assert.equal((await client.replay()).competition.mode,'standard');
});

test('race mode stays locked while paused and unlocks after completion or reset',()=>{
  assert.equal(isRaceInProgress({status:'running',frame:{step:0}}),true);
  assert.equal(isRaceInProgress({status:'paused',frame:{step:0}},true),true);
  assert.equal(isRaceInProgress({status:'paused',frame:{step:1}}),true);
  assert.equal(isRaceInProgress({status:'paused',frame:{step:0}},false),false);
  for(const status of ['finished','error','idle'])
    assert.equal(isRaceInProgress({status,frame:{step:100}},true),false);
});
