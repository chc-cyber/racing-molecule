import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import createModule from '../wasm/dist/browser/runtime/molecular-car.mjs';
import {validateProject} from '../frontend/editor/model.js';
const module=await createModule({print:()=>{},printErr:()=>{}});
function request(op,data){try{return module.workshopRequest(op,data);}catch(error){const text=module.getExceptionMessage(error);module.decrementExceptionRefcount(error);throw new Error(Array.isArray(text)?text.join(': '):String(text));}}
const mol=await fs.readFile(new URL('../frontend/editor/examples/tetraphenylbenzene.mol',import.meta.url),'utf8');
const input=request('import',{filename:'car.mol',content:mol}).project;
test('browser MOL/SDF roundtrip preserves atoms, bonds, coordinates and stereo',()=>{
  for(const filename of ['mol','sdf']){
    const file=request('export',{project:input,dimension:3,format:filename});
    const result=request('import',file).project;validateProject(result);
    assert.deepEqual(result.atoms.map(a=>[a.element,a.charge,a.chiralTag]),input.atoms.map(a=>[a.element,a.charge,a.chiralTag]));
    assert.equal(result.bonds.length,input.bonds.length);
    assert(result.atoms.every((a,i)=>a.pos3.every((v,k)=>Math.abs(v-input.atoms[i].pos3[k])<1e-5)));
  }
});
test('2D construction, hydrogen completion, layout and fixed UFF stay independent',()=>{
  const p=structuredClone(input);p.has3d=false;
  const g=request('geometry',{project:p,operation:'prepare',fixedIds:[],revision:7});validateProject(g.project);
  assert.equal(g.revision,7);assert(g.project.has3d);
  assert.deepEqual(g.project.atoms.slice(0,p.atoms.length).map(a=>[a.id,a.pos2]),p.atoms.map(a=>[a.id,a.pos2]));
  const h=request('geometry',{project:g.project,operation:'hydrogenate',fixedIds:[],revision:8}).project;
  assert.deepEqual(h.atoms.map(a=>a.pos3),g.project.atoms.map(a=>a.pos3));
  const layout=request('geometry',{project:h,operation:'layout2d',fixedIds:[],revision:9}).project;
  assert.deepEqual(layout.atoms.map(a=>[a.id,a.pos3]),h.atoms.map(a=>[a.id,a.pos3]));
  const fixed=h.atoms.slice(0,3).map(a=>a.id);
  const relaxed=request('geometry',{project:h,operation:'relax',fixedIds:fixed,revision:10}).project;
  assert.deepEqual(relaxed.atoms.slice(0,3).map(a=>a.pos3),h.atoms.slice(0,3).map(a=>a.pos3));
  const rebuilt=request('geometry',{project:h,operation:'rebuild',fixedIds:[],revision:11}).project;
  assert.equal(rebuilt.atoms.length,h.atoms.length);assert.deepEqual(rebuilt.bonds,h.bonds);
});
test('workshop rejects multiple records and unsupported query chemistry',()=>{
  assert.throws(()=>request('import',{content:mol+'\n$$$$\n'+mol,filename:'two.sdf'}),/单分子/);
  assert.throws(()=>request('geometry',{project:input,operation:'unknown',revision:0}),/未知/);
});
test('chiral example retains handedness and stable IDs through browser export',async()=>{
  const text=await fs.readFile(new URL('../frontend/editor/examples/chiral.mol',import.meta.url),'utf8');
  const p=request('import',{content:text,filename:'chiral.mol'}).project;
  const round=request('import',request('export',{project:p,dimension:3,format:'mol'})).project;
  assert.deepEqual(round.atoms.map(a=>a.chiralTag),p.atoms.map(a=>a.chiralTag));
  const g=request('geometry',{project:p,operation:'hydrogenate',fixedIds:[],revision:1}).project;
  assert.deepEqual(g.atoms.slice(0,p.atoms.length).map(a=>a.id),p.atoms.map(a=>a.id));
  const planar=request('import',request('export',{project:p,dimension:2,format:'mol'})).project;
  assert.deepEqual(planar.atoms.map(a=>a.chiralTag),p.atoms.map(a=>a.chiralTag));
  const built=request('geometry',{project:planar,operation:'prepare',fixedIds:[],revision:2}).project;
  const verified=request('import',request('export',{project:built,dimension:3,format:'mol'})).project;
  assert.deepEqual(verified.atoms.map(a=>a.chiralTag),p.atoms.map(a=>a.chiralTag));
});
test('SUP groups retain real C60 atoms and attachment references across files',async()=>{
  const fragments=JSON.parse(await fs.readFile(new URL('../frontend/editor/fragments.json',import.meta.url),'utf8')).fragments;
  const cage=structuredClone(fragments[2]);cage.format='molecular-car-workshop-v1';
  cage.metadata={ketcherGroups:[{id:'g-c60',type:'SUP',name:'C60',expanded:false,atomIds:cage.atoms.map(a=>a.id),attachmentPoints:[{atomId:cage.atoms[0].id,attachmentId:'1'}]}]};
  validateProject(cage);
  const p=request('import',request('export',{project:cage,dimension:3,format:'mol'})).project;
  assert.equal(p.atoms.length,60);assert.equal(p.bonds.length,90);
  assert.equal(p.metadata.ketcherGroups[0].atomIds.length,60);
  assert.equal(p.metadata.ketcherGroups[0].attachmentPoints[0].atomId,p.atoms[0].id);
});
