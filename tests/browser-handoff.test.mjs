import test from 'node:test';
import assert from 'node:assert/strict';
const messages=[],events=new EventTarget();
Object.assign(globalThis,{document:{documentElement:{dataset:{engine:'wasm'}}},window:{},parent:{postMessage:data=>messages.push(data)},location:{origin:'https://game.invalid'},addEventListener:events.addEventListener.bind(events),removeEventListener:events.removeEventListener.bind(events)});
const {prepareBrowserRace,cancelBrowserRace}=await import('../frontend/editor/browser-service.mjs');
function deliver(data,{source=parent,origin=location.origin}={}){const event=new Event('message');Object.defineProperties(event,{data:{value:data},source:{value:source},origin:{value:origin}});events.dispatchEvent(event);}
test('browser handoff accepts only the direct parent, origin, and matching job',async()=>{
  let progress='';const promise=prepareBrowserRace({filename:'car.mol',content:'molecule'},text=>progress=text);
  const {id}=messages.at(-1);
  deliver({channel:'racing-molecule',type:'prepare-result',id},{origin:'https://elsewhere.invalid'});
  deliver({channel:'racing-molecule',type:'prepare-result',id},{source:{}});
  deliver({channel:'racing-molecule',type:'prepare-progress',id,message:'准备赛车…'});
  assert.equal(progress,'准备赛车…');
  deliver({channel:'racing-molecule',type:'prepare-result',id,generation:3});assert.equal((await promise).generation,3);
});
test('cancel rejects immediately and late results cannot complete the next job',async()=>{
  const old=prepareBrowserRace({},()=>{}),oldID=messages.at(-1).id;
  cancelBrowserRace();await assert.rejects(old,/返回编辑/);
  assert.equal(messages.at(-1).type,'cancel-prepare');
  const current=prepareBrowserRace({},()=>{}),id=messages.at(-1).id;
  deliver({channel:'racing-molecule',type:'prepare-result',id:oldID});
  deliver({channel:'racing-molecule',type:'prepare-result',id,error:'当前结构未收敛'});
  await assert.rejects(current,/未收敛/);
});
