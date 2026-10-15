import test from 'node:test';
import assert from 'node:assert/strict';
import {createProject,WorkshopModel,clone,reconcileC60Addition,reconcileGeneratedHydrogens,validateBondEndpoints} from '../frontend/editor/model.js';
import {activeChemistryWarnings,C60_ADDITION_NOTICE,LEGACY_C60_ADDITION_NOTICE} from '../frontend/editor/editor-notices.js';
import {C60_POS3,C60_BONDS} from '../frontend/editor/native-c60-data.js';
import {KetcherAdapter} from '../frontend/editor/ketcher-adapter.js';

test('A completed C60 operation does not masquerade as an unresolved warning',()=>{
 for(const text of [C60_ADDITION_NOTICE,LEGACY_C60_ADDITION_NOTICE]){
  const p=createProject();p.metadata={stereoNotice:text};assert.deepEqual(activeChemistryWarnings(p),[]);
  p.metadata.stereoNotice+=' 手性需要确认。';assert.deepEqual(activeChemistryWarnings(p),['手性需要确认。']);
 }
});
test('Actual valence warnings remain visible',()=>{
 const p=createProject();p.atoms=[{id:'a-c',element:'C',charge:0,pos2:[0,0],pos3:[0,0,0]},...Array.from({length:5},(_,i)=>({id:'a-h'+i,element:'H',charge:0,pos2:[i,1],pos3:[i,1,0]}))];p.bonds=p.atoms.slice(1).map(a=>({id:'b-'+a.id,a:'a-c',b:a.id,order:1}));assert.equal(activeChemistryWarnings(p).length,1);
});
test('Interleaved disconnected atoms use a declared component transport order without changing the project',()=>{
 const p=createProject();p.atoms=['c','o','hc','ho'].map((id,i)=>({id,element:['C','O','H','H'][i],charge:0,pos2:[i,0],pos3:[i,0,i/10]}));p.bonds=[{id:'oh',a:'o',b:'ho',order:1},{id:'ch',a:'c',b:'hc',order:1}];const before=clone(p);
 const {ket,reference}=new KetcherAdapter().toKet(p);assert.deepEqual(reference.atomIds,['c','hc','o','ho']);assert.deepEqual(reference.bondIds,['ch','oh']);assert.deepEqual(ket.mol0.bonds.map(b=>b.atoms),[[0,1],[2,3]]);assert.deepEqual(p,before);
});
test('Declared native C60 addition lowers one double bond and shares exact undo',()=>{
 const p=createProject();p.has3d=true;p.atoms=C60_POS3.map((v,i)=>({id:'c-'+i,element:'C',charge:0,pos2:v.slice(0,2),pos3:[...v]}));p.bonds=C60_BONDS.map((b,i)=>({id:'c-b-'+i,a:'c-'+b.a,b:'c-'+b.b,order:b.order}));
 p.metadata={ketcherGroups:[{id:'g',kind:'c60',type:'SUP',name:'C60',atomIds:p.atoms.map(a=>a.id),expanded:false}]};
 const m=new WorkshopModel(p),before=clone(m.project),draft=clone(m.project);draft.atoms.push({id:'branch',element:'C',charge:0,pos2:[5,0],pos3:[5,0,0]});draft.bonds.push({id:'connect',a:'c-0',b:'branch',order:1});
 reconcileGeneratedHydrogens(before,draft);assert(reconcileC60Addition(before,draft));validateBondEndpoints(draft,'c-0','branch');
 m.commit('连接',p=>Object.assign(p,draft),{geometryConfirmedIds:['branch']});assert.equal(m._undo.length,1);assert.deepEqual(m.project.atoms.slice(0,60).map(a=>a.pos3),before.atoms.map(a=>a.pos3));assert(m.undo());assert.deepEqual(m.project,before);
});
