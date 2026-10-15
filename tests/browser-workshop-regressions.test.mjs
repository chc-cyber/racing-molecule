import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import createModule from '../wasm/dist/browser/runtime/molecular-car.mjs';
import {validateProject,createProject,reconcileC60Addition} from '../frontend/editor/model.js';
const module=await createModule({print:()=>{},printErr:()=>{}});
function request(operation,payload){try{return module.workshopRequest(operation,payload);}catch(e){const message=module.getExceptionMessage(e);module.decrementExceptionRefcount(e);throw Error(String(message));}}
function prepare(project,operation='prepare'){const t=performance.now();const r=request('geometry',{project,operation,revision:3,fixedIds:[]});assert.equal(r.revision,3);validateProject(r.project);console.log(project.name,Math.round(performance.now()-t)+' ms');return r;}
const data=JSON.parse(await fs.readFile(new URL('../frontend/editor/ketcher-source/c60-source.json',import.meta.url)));
function cage(){return validateProject({...createProject('C60'),atoms:structuredClone(data.atoms),bonds:structuredClone(data.bonds),metadata:{ketcherGroups:[{id:'cage',type:'SUP',kind:'c60',name:'C60',expanded:false,atomIds:data.atoms.map(a=>a.id),attachmentPoints:[{atomId:data.atoms[0].id,attachmentId:'1'}]}]}});}
const coords=p=>p.atoms.map(a=>[a.id,a.pos3]);
function checkBonds(p){const atoms=new Map(p.atoms.map(a=>[a.id,a]));for(const b of p.bonds){const a=atoms.get(b.a),c=atoms.get(b.b),d=Math.hypot(...a.pos3.map((v,i)=>v-c.pos3[i]));assert(d>.7&&d<1.85,`${a.element}-${c.element}: ${d}`);}}
for(const name of ['hexahelicene','spiral'])test(name+' builds a nonplanar conformer without changing topology or 2D',async()=>{
  const content=await fs.readFile(new URL(`./data/workshop/${name}.mol`,import.meta.url),'utf8');
  const p=request('import',{content,filename:name+'.mol'}).project,before=structuredClone(p);
  const {project:q,warnings}=prepare(p);assert(q.has3d);assert.deepEqual(p,before);
  assert.deepEqual(q.atoms.slice(0,p.atoms.length).map(a=>[a.id,a.pos2]),p.atoms.map(a=>[a.id,a.pos2]));
  assert.deepEqual(q.bonds.slice(0,p.bonds.length),p.bonds);checkBonds(q);
  assert(warnings.some(w=>w.includes('弯曲')));
  assert.deepEqual(coords(prepare(q).project),coords(q));
});
test('native C60 uses its indexed spherical seed, repeated switching is exact',()=>{
  const p=cage(),q=prepare(p).project;assert.equal(q.atoms.length,60);assert.equal(q.bonds.length,90);assert.deepEqual(coords(q),coords(p));checkBonds(q);
  assert.deepEqual(coords(prepare(q).project),coords(q));assert.deepEqual(coords(prepare(p,'rebuild').project),coords(p));
});
test('C60 methyl addition constructs only methyl and H; sphere remains exact',()=>{
  const p=cage(),q=structuredClone(p);q.name='C60-CH3';q.atoms.push({...q.atoms[0],id:'methyl',pos2:[8,0],pos3:[8,0,0]});q.bonds.push({id:'link',a:q.atoms[0].id,b:'methyl',order:1});reconcileC60Addition(p,q);
  const r=prepare(q).project;assert.equal(r.atoms.length,65);assert.equal(r.bonds.length,95);assert.deepEqual(coords(r).slice(0,60),coords(p));checkBonds(r);
  assert.deepEqual(coords(prepare(r).project),coords(r));
  // Replace one generated methyl H with a new carbon after manual translation.
  for(const a of r.atoms)a.pos3=a.pos3.map((v,i)=>v+[3,-2,1][i]);
  const hydrogen=r.bonds.find(b=>(b.a==='methyl'&&r.atoms.find(a=>a.id===b.b).element==='H')||(b.b==='methyl'&&r.atoms.find(a=>a.id===b.a).element==='H'));
  const hid=hydrogen.a==='methyl'?hydrogen.b:hydrogen.a;
  const next=structuredClone(r);next.atoms=next.atoms.filter(a=>a.id!==hid);next.bonds=next.bonds.filter(b=>b.id!==hydrogen.id);
  const stable=coords(next);next.atoms.push({...q.atoms.at(-1),id:'extension',pos3:[100,0,0]});next.bonds.push({id:'extend',a:'methyl',b:'extension',order:1});next.metadata.geometryPending={atomIds:['extension']};
  const extended=prepare(next).project;assert.deepEqual(coords(extended).slice(0,stable.length),stable);checkBonds(extended);
});
test('adding a ring preserves every existing atom, including H, and canonical IDs',async()=>{
  const p=prepare(cage()).project;const content=await fs.readFile(new URL('../examples/benzene.mol',import.meta.url),'utf8');
  const ring=request('import',{filename:'ring.mol',content}).project;
  // The fixture includes H: connect to a carbon after removing its attached H.
  const heavy=ring.atoms.filter(a=>a.element!=='H'),members=new Set(heavy.map(a=>a.id));
  const q=structuredClone(p);q.atoms.push(...heavy);q.bonds.push(...ring.bonds.filter(b=>members.has(b.a)&&members.has(b.b)));q.bonds.push({id:'link-ring',a:q.atoms[0].id,b:heavy[0].id,order:1});reconcileC60Addition(p,q);q.metadata.geometryPending={atomIds:heavy.map(a=>a.id)};
  const r=prepare(q).project;assert.deepEqual(coords(r).slice(0,60),coords(p));checkBonds(r);
});
test('invalid neutral imine N is rejected by geometry; charged iminium is accepted',()=>{
  const p=validateProject({...createProject('imine'),atoms:['C','N','C','C'].map((element,i)=>({id:'a'+i,element,pos2:[i,0],pos3:[i,0,0]})),bonds:[{id:'b0',a:'a0',b:'a1',order:2},{id:'b1',a:'a1',b:'a2',order:1},{id:'b2',a:'a1',b:'a3',order:1}]});
  assert.throws(()=>prepare(p),/valence|价态/i);p.atoms[1].charge=1;const q=prepare(p).project;assert(!q.bonds.some(b=>(b.a==='a1'&&q.atoms.find(a=>a.id===b.b).element==='H')||(b.b==='a1'&&q.atoms.find(a=>a.id===b.a).element==='H')));
});
