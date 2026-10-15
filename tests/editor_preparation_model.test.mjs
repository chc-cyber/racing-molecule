import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
const source = readFileSync(new URL('../frontend/editor/model.js', import.meta.url), 'utf8');
const {WorkshopModel,createProject,clone,validateProject} = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
const atom=(id,element='C',extra={})=>({id,element,charge:0,pos2:[0,0],pos3:[0,0,0],...extra});
const bond=(id,a,b,order=1)=>({id,a,b,order});
function existing(){return {...createProject('existing'),has3d:true,atoms:[atom('c'),atom('d','C',{pos2:[1,0],pos3:[1.5,0,0]})],bonds:[bond('cd','c','d')]};}

test('native 3D construction confirms new stable IDs in the same reversible transaction',()=>{
  const m=new WorkshopModel(existing()),before=clone(m.project),history=m._undo.length;
  m.commit('native 3D carbon',p=>{p.atoms.push(atom('new','C',{pos3:[3.02,1,0]}));p.bonds.push(bond('dn','d','new'));},{geometryConfirmedIds:['new']});
  assert.deepEqual(m.project.metadata.geometryPending,{atomIds:[],changedAtomIds:['d']});
  assert.equal(m._undo.length,history+1);const constructed=clone(m.project);
  m.undo();assert.deepEqual(m.project,before);m.redo();assert.deepEqual(m.project,constructed);
});

test('only explicitly placed 3D vertices are confirmed; a 2D drag leaves new geometry pending',()=>{
  const m=new WorkshopModel(existing());
  m.commit('2D new chain',p=>{p.atoms.push(atom('n1','C',{pos3:[2.4,1,0]}),atom('n2','C',{pos3:[3.4,1,1]}));p.bonds.push(bond('dn','d','n1'),bond('nn','n1','n2'));});
  m.moveSelection([2,3],{dimension:2,rigid:false,ids:['n1']});
  assert.deepEqual(m.project.metadata.geometryPending.atomIds,['n1','n2']);const before=clone(m.project);
  m.moveSelection([0,1,1],{dimension:3,rigid:false,ids:['n1']});
  assert.deepEqual(m.project.metadata.geometryPending.atomIds,['n2']);assert.deepEqual(m.project.atoms.find(a=>a.id==='n1').pos3,[2.4,2,1]);
  const placed=clone(m.project);m.undo();assert.deepEqual(m.project,before);m.redo();assert.deepEqual(m.project,placed);
});

test('exact distance, angle, torsion and fragment rotation confirm the moved new branch',()=>{
  for(const action of [m=>m.setDistance('d','n1',1.8),m=>m.setAngle('c','d','n1',120),m=>m.setDihedral('c','d','n1','n2',70),m=>m.rotateSelection('x',35,{rigid:false,ids:['n1','n2']})]){
    const m=new WorkshopModel(existing());
    m.commit('new branch',p=>{p.atoms.push(atom('n1','C',{pos3:[2.4,1,0]}),atom('n2','C',{pos3:[3.4,1,1]}));p.bonds.push(bond('dn','d','n1'),bond('nn','n1','n2'));});
    const before=clone(m.project);action(m);assert.deepEqual(m.project.metadata.geometryPending.atomIds,[]);
    assert.deepEqual(m.project.atoms.slice(0,2),before.atoms.slice(0,2));const positioned=clone(m.project);
    m.undo();assert.deepEqual(m.project,before);m.redo();assert.deepEqual(m.project,positioned);
  }
});

test('unknown confirmed geometry IDs reject atomically without a new history entry',()=>{
  const m=new WorkshopModel(existing()),before=clone(m.project),history=m._undo.length;
  assert.throws(()=>m.commit('invalid placement',p=>{p.name='bad';},{geometryConfirmedIds:['missing']}),/不存在/);
  assert.deepEqual(m.project,before);assert.equal(m._undo.length,history);
});

test('layout, true 3D motion, and selection never mark topology for automatic optimization',()=>{
  const m=new WorkshopModel(existing());m.selection=new Set(['c']);
  m.moveSelection([3,2],{dimension:2,rigid:false});assert.equal(m.project.metadata?.geometryPending,undefined);
  m.moveSelection([0,1,2],{dimension:3,rigid:false});assert.equal(m.project.metadata?.geometryPending,undefined);
  const original=clone(m.project);m.undo();m.redo();assert.deepEqual(m.project,original);
});

test('new chemical region and its old attachment are tracked by stable ID and restored with geometry undo',()=>{
  const m=new WorkshopModel(existing()),fixed=clone(m.project.atoms);
  m.commit('add ring',p=>{p.atoms.push(atom('new1'),atom('new2'));p.bonds.push(bond('cn','c','new1'),bond('nn','new1','new2',2));});
  assert.deepEqual(m.project.metadata.geometryPending,{atomIds:['new1','new2'],changedAtomIds:['c']});
  const before=clone(m.project),prepared=clone(before);delete prepared.metadata.geometryPending;
  prepared.atoms.find(a=>a.id==='new1').pos3=[0,1.4,0];
  prepared.atoms.find(a=>a.id==='new2').pos3=[0,2.8,0];
  assert.equal(m.applyAsync(m.revision,prepared,'prepare'),true);
  assert.equal(m.project.metadata.geometryPending,undefined);
  assert.deepEqual(m.project.atoms.slice(0,2),fixed);
  m.undo();assert.deepEqual(m.project,before);m.redo();assert.deepEqual(m.project,prepared);
});

test('explicit substitution consumes only editor-generated H in the same undo transaction',()=>{
  const p=existing();p.atoms.push(atom('h1','H',{generatedHydrogen:true}),atom('h2','H',{generatedHydrogen:true}),atom('h3','H',{generatedHydrogen:true}));
  p.bonds.push(bond('ch1','c','h1'),bond('ch2','c','h2'),bond('ch3','c','h3'));
  const m=new WorkshopModel(p),original=clone(m.project);
  m.commit('add substituent',q=>{q.atoms.push(atom('n','N'));q.bonds.push(bond('cn','c','n'));});
  assert.equal(m.project.atoms.filter(a=>a.element==='H').length,2);
  assert.equal(m.project.atoms.some(a=>a.id==='h1'),false);
  assert.deepEqual(m.project.atoms.find(a=>a.id==='d'),original.atoms.find(a=>a.id==='d'));
  assert.ok(m.project.metadata.hydrogenNotice.includes('1'));
  const connected=clone(m.project);m.undo();assert.deepEqual(m.project,original);m.redo();assert.deepEqual(m.project,connected);
});

test('imported or manually drawn H is kept even when the new structure is overvalent',()=>{
  const p=existing();for(let i=0;i<3;i++){p.atoms.push(atom('h'+i,'H'));p.bonds.push(bond('ch'+i,'c','h'+i));}
  const m=new WorkshopModel(p);m.commit('draft bond',q=>{q.atoms.push(atom('n','N'));q.bonds.push(bond('cn','c','n'));});
  assert.equal(m.project.atoms.filter(a=>a.element==='H').length,3);
  assert.equal(validateProject(m.project).atoms.length,6);
});

test('deleting a SUP member prunes stable group references and invalidates template provenance',()=>{
  const p=existing();p.metadata={ketcherGroups:[{id:'g',type:'SUP',name:'C60',kind:'c60',expanded:false,atomIds:['c','d'],attachmentPoints:[{atomId:'c',attachmentId:'1'}],template:{indexToAtomId:['c','d']}}]};
  const m=new WorkshopModel(p);m.deleteSelection({ids:['c']});
  assert.deepEqual(m.project.metadata.ketcherGroups[0].atomIds,['d']);
  assert.deepEqual(m.project.metadata.ketcherGroups[0].attachmentPoints,[]);
  assert.equal(m.project.metadata.ketcherGroups[0].template,undefined);
  m.undo();assert.deepEqual(m.project,p);
});

test('stale prepare results cannot erase a newer topology edit or its pending region',()=>{
  const m=new WorkshopModel(existing()),revision=m.revision,old=clone(m.project);
  const n=m.addAtom({element:'N'});assert.equal(m.applyAsync(revision,old),false);
  assert.ok(m.project.atoms.some(a=>a.id===n));assert.ok(m.project.metadata.geometryPending.atomIds.includes(n));
});

test('old chemistry edits may move their explicit targets but an added substituent keeps its old anchor fixed',()=>{
  const m=new WorkshopModel(existing());m.updateAtom('c',{element:'N'});
  assert.deepEqual(m.project.metadata.geometryPending.editedAtomIds,['c']);
  m.undo();m.addAtom({id:'n',element:'N'});m.addBond('c','n');
  assert.equal(m.project.metadata.geometryPending.editedAtomIds,undefined);
  assert.deepEqual(m.project.metadata.geometryPending.changedAtomIds,['c']);
});

test('a native C60 crossing bond checks valence in every tool and rejects transactionally',()=>{
  const p=existing();p.atoms.push(atom('x'),atom('y'),atom('external','N'));
  p.bonds.push(bond('cx','c','x',2),bond('cy','c','y'));
  p.metadata={ketcherGroups:[{id:'g',type:'SUP',kind:'c60',name:'C60',expanded:false,atomIds:['c','d','x','y'],attachmentPoints:[]}]};
  const m=new WorkshopModel(p),before=clone(m.project);
  assert.throws(()=>m.addBond('c','external'),/C60 连接未应用/);
  assert.deepEqual(m.project,before);assert.equal(m.canUndo,false);
  m.updateBond('cx',{order:1});m.addBond('c','external');
  assert.equal(m.project.bonds.filter(b=>b.a==='c'||b.b==='c').reduce((n,b)=>n+b.order,0),4);
  m.undo();assert.equal(m.project.bonds.some(b=>b.a==='external'||b.b==='external'),false);
});

test('automatic stereo invalidation during H substitution never unfreezes the old chiral anchor',()=>{
  const p=existing();p.atoms[0].chiralTag='CHI_TETRAHEDRAL_CW';
  p.atoms.push(atom('h','H',{generatedHydrogen:true,pos3:[-1,-1,-1]}),atom('f','F',{pos3:[0,1.4,0]}),atom('cl','Cl',{pos3:[0,0,1.7]}));
  p.bonds.push(bond('ch','c','h'),bond('cf','c','f'),bond('ccl','c','cl'));
  const m=new WorkshopModel(p),anchor=clone(p.atoms[0].pos3);
  m.addAtom({id:'br',element:'Br'});m.addBond('c','br');
  assert.equal(m.project.atoms.find(a=>a.id==='c').chiralTag,'CHI_UNSPECIFIED');
  assert.equal(m.project.metadata.geometryPending.editedAtomIds,undefined);
  assert.ok(m.project.metadata.geometryPending.changedAtomIds.includes('c'));
  assert.deepEqual(m.project.atoms.find(a=>a.id==='c').pos3,anchor);
});

test('generated H replacement uses the same charged-carbon valence as canonical chemistry checks',()=>{
  const p={...createProject(),has3d:true,atoms:[atom('c','C',{charge:-1})],bonds:[]};
  for(let i=0;i<3;i++){p.atoms.push(atom('h'+i,'H',{generatedHydrogen:true}));p.bonds.push(bond('ch'+i,'c','h'+i));}
  const m=new WorkshopModel(p);m.addAtom({id:'f',element:'F'});m.addBond('c','f');
  assert.equal(m.project.atoms.filter(a=>a.element==='H').length,2);
  assert.equal(m.project.bonds.filter(b=>b.a==='c'||b.b==='c').reduce((n,b)=>n+b.order,0),3);
});

test('declared bracket-H counts are included when an explicit connection substitutes generated H',()=>{
  const p={...createProject(),has3d:true,atoms:[atom('c','C',{explicitHs:2})],bonds:[]};
  for(let i=0;i<2;i++){p.atoms.push(atom('h'+i,'H',{generatedHydrogen:true}));p.bonds.push(bond('ch'+i,'c','h'+i));}
  const m=new WorkshopModel(p);m.addAtom({id:'f',element:'F'});m.addBond('c','f');
  assert.equal(m.project.atoms.filter(a=>a.element==='H').length,1);
  assert.equal(m.project.atoms.find(a=>a.id==='c').explicitHs,2);
});
