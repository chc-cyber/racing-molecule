import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {KetcherAdapter, KETCHER_SIDECAR_FORMAT} from '../frontend/editor/ketcher-adapter.js';
import {WorkshopModel, createProject, clone, serializeProject, parseProject} from '../frontend/editor/model.js';
import {C60_POS3,C60_BONDS} from '../frontend/editor/native-c60-data.js';
// Execute the pinned real Ketcher Pool implementation (not a replacement Map).
import {Pool} from '../frontend/editor/research/upstream/ketcher/packages__ketcher-core__src__domain__entities__pool.ts';

function importer(ket, offset = 0) {
  const atoms = new Pool(), bonds = new Pool();
  for (let i = 0; i < offset; i++) { atoms.add({}); bonds.add({}); }
  atoms.clear(); bonds.clear();
  const molecule = ket.mol0 ?? {atoms: [], bonds: []}, poolIds = [];
  for (const atom of molecule.atoms) poolIds.push(atoms.add({label: atom.label, charge: atom.charge ?? null, isotope: atom.isotope ?? null,
    aam: atom.mapping ?? 0, alias: atom.alias ?? null, radical: atom.radical ?? 0, implicitHCount: atom.implicitHCount ?? null,
    implicitH: 2, stereoParity: atom.stereoParity ?? 0, stereoLabel: atom.stereoLabel ?? null,
    pp: {x: atom.location[0], y: -atom.location[1], z: 0}}));
  for (const bond of molecule.bonds) bonds.add({type: bond.type, begin: poolIds[bond.atoms[0]], end: poolIds[bond.atoms[1]], stereo: bond.stereo ?? 0});
  const sgroups=new Pool();
  for(const group of molecule.sgroups??[])sgroups.add({type:group.type,isMonomer:false,data:{name:group.name,expanded:group.expanded},atoms:group.atoms.map(i=>poolIds[i]),getAttachmentPoints(){return(group.attachmentPoints??[]).map(p=>({atomId:poolIds[p.attachmentAtom],leaveAtomId:p.leavingAtom===undefined?undefined:poolIds[p.leavingAtom],attachmentId:p.attachmentId}));}});
  return {atoms, bonds, frags: new Pool(), sgroups, rxnArrows: new Pool(), rxnPluses: new Pool()};
}
function fixture() {
  return {...createProject('adapter test'), has3d: true,
    atoms: [
      {id: 'stable-a', element: 'C', charge: 0, isotope: 13, explicitHs: 1, noImplicit: true, chiralTag: 'CHI_TETRAHEDRAL_CW', pos2: [0, 0], pos3: [9, 8, 7]},
      {id: 'stable-b', element: 'N', charge: 1, pos2: [1.5, 0], pos3: [10.5, 8, 7]},
      {id: 'stable-c', element: 'O', charge: 0, pos2: [3, 0], pos3: [12, 8, 7]},
    ], bonds: [{id: 'stable-ab', a: 'stable-a', b: 'stable-b', order: 1}, {id: 'stable-bc', a: 'stable-b', b: 'stable-c', order: 1}],
    fragments: [{id: 'fragment', name: 'custom', atomIds: ['stable-a', 'stable-b'], ports: [{atomId: 'stable-b', maxOrder: 1}], rigid: true}]};
}
function setup(project = fixture(), offset = 0) {
  let sequence = 0;
  const model = new WorkshopModel(project), adapter = new KetcherAdapter({makeId: prefix => `${prefix}-new-${++sequence}`});
  const payload = adapter.toKet(project), struct = importer(payload.ket, offset);
  adapter.bindRebuild(struct, project, payload.reference);
  return {model, adapter, struct};
}
function addAtom(struct, x, y = 0, extra = {}) {
  return struct.atoms.add({label: 'C', charge: 0, isotope: 0, radical: 0, implicitHCount: null, aam: 0, stereoParity: 0, stereoLabel: null, pp: {x, y, z: 999}, ...extra});
}

test('explicit KET reference binds nonzero pool IDs without using atom maps or graph matching', () => {
  const {model, adapter, struct} = setup(fixture(), 17);
  assert.equal(adapter.atoms.get(17).id, 'stable-a'); assert.equal(adapter.bonds.get(18).id, 'stable-bc');
  assert.deepEqual([...adapter.selectionIds(struct, {atoms: [17], bonds: [18]})], ['stable-a', 'stable-b', 'stable-c']);
  assert.deepEqual(adapter.toolSelection(struct, new Set(['stable-a', 'stable-c']), {bondId: 'stable-ab'}), {atoms: [17, 19], bonds: [17]});
  assert.equal(struct.atoms.get(17).aam, 0);
  assert.equal(adapter.commit(struct, model).changed, false); assert.equal(model.revision, 0);
});

test('Ketcher-only 2D layout leaves 3D, isotope, explicit H, stereo, fragment, and stable IDs intact', () => {
  const {model, adapter, struct} = setup(), original = clone(model.project);
  struct.atoms.get(0).pp = {x: 22, y: -13, z: 500};
  const result = adapter.commit(struct, model);
  assert.equal(result.changed, true); assert.equal(result.requiresPush, false);
  assert.deepEqual(model.project.atoms[0].pos2, [22, 13]);
  assert.deepEqual(model.project.atoms.map(a => a.pos3), original.atoms.map(a => a.pos3));
  assert.equal(model.project.atoms[0].explicitHs, 1); assert.equal(model.project.atoms[0].noImplicit, true);
  assert.equal(model.project.atoms[0].chiralTag, 'CHI_TETRAHEDRAL_CW'); assert.equal(model.project.atoms[0].isotope, 13);
  assert.deepEqual(model.project.fragments, original.fragments);
  assert.equal(model.project.metadata.ketcherAdapter.format, KETCHER_SIDECAR_FORMAT);
  assert.equal(model.project.metadata.ketcherAdapter.atomIdsByPool['0'], 'stable-a');
  assert.deepEqual(parseProject(serializeProject(model.project)), model.project);
  assert.equal(adapter.commit(struct, model).changed, false); assert.equal(model.revision, 1);
});

test('chemical atom map aam is preserved as chemistry data, never overwritten with editor UUIDs', () => {
  const {model, adapter, struct} = setup();
  struct.atoms.get(1).aam = 317;
  // Atom-map-only updates are native data and must survive a project save.
  adapter.commit(struct, model);
  assert.equal(model.project.metadata.ketcherAdapter.atomData['stable-b'].aam, 317);
  const rebuilt = adapter.toKet(parseProject(serializeProject(model.project)));
  assert.equal(rebuilt.ket.mol0.atoms[1].mapping, 317);
  assert.equal(struct.atoms.get(1).aam, 317);
  assert.equal(model.project.atoms[1].id, 'stable-b');
});

test('computed implicit H does not overwrite canonical explicit hydrogen data', () => {
  const {model, adapter, struct} = setup();
  struct.atoms.get(0).implicitH = 7;
  assert.equal(adapter.commit(struct, model).changed, false);
  struct.atoms.get(0).implicitHCount = 2;
  adapter.commit(struct, model);
  assert.equal(model.project.atoms[0].explicitHs, 2); assert.equal(model.project.atoms[0].noImplicit, true);
  struct.atoms.get(0).implicitHCount = null;
  adapter.commit(struct, model);
  assert.equal(model.project.atoms[0].explicitHs, 0); assert.equal(model.project.atoms[0].noImplicit, false);
});

test('natural abundance omits KET isotope zero while declared isotopes and canonical coordinates round-trip', () => {
  const project = fixture();
  project.atoms[1].isotope = 0;
  // The third atom intentionally has no isotope field; both forms are natural.
  const {model, adapter, struct} = setup(project), before = clone(model.project);
  const payload = adapter.toKet(model.project);
  assert.equal(payload.ket.mol0.atoms[0].isotope, 13);
  assert.equal(Object.hasOwn(payload.ket.mol0.atoms[1], 'isotope'), false);
  assert.equal(Object.hasOwn(payload.ket.mol0.atoms[2], 'isotope'), false);
  assert.equal(struct.atoms.get(0).isotope, 13);
  assert.equal(struct.atoms.get(1).isotope, null);
  assert.equal(struct.atoms.get(2).isotope, null);
  assert.equal(adapter.commit(struct, model).changed, false);
  assert.deepEqual(model.project, before);
  struct.atoms.get(0).pp.x += 0.2;
  adapter.commit(struct, model);
  assert.deepEqual(model.project.atoms.map(a => a.pos3), before.atoms.map(a => a.pos3));
  assert.equal(model.project.atoms[0].isotope, 13);
  assert.equal(model.project.atoms[1].isotope, 0);
  assert.equal(Object.hasOwn(model.project.atoms[2], 'isotope'), false);
  assert.deepEqual(model.project.bonds, before.bonds);
});

test('local addition initializes new 3D positions via explicit neighbors while all surviving coordinates stay exact', () => {
  const {model, adapter, struct} = setup(), previous = clone(model.project);
  const poolId = addAtom(struct, 4.5);
  struct.bonds.add({begin: 2, end: poolId, type: 1, stereo: 0});
  adapter.commit(struct, model);
  assert.equal(model.project.atoms.length, 4);
  assert.deepEqual(model.project.atoms.slice(0, 3).map(a => a.pos3), previous.atoms.map(a => a.pos3));
  assert.deepEqual(model.project.atoms[3].pos3, [13.5, 8, 7]);
  assert.notEqual(model.project.atoms[3].id, 'stable-c');
  assert.match(model.project.metadata.geometryNotice, /其余玩家构象保持不变/);
  struct.bonds.delete(1); struct.bonds.delete(2); struct.atoms.delete(2);
  adapter.commit(struct, model);
  assert.deepEqual(model.project.atoms.map(a => a.id), ['stable-a', 'stable-b', 'a-new-1']);
});

test('actual Pool delete/set AtomAdd-style undo preserves a tombstoned UUID and its prior 3D coordinate', () => {
  const {model, adapter, struct} = setup();
  const oldAtom = clone(struct.atoms.get(2)), oldBond = clone(struct.bonds.get(1));
  struct.atoms.delete(2); struct.bonds.delete(1); adapter.commit(struct, model);
  struct.atoms.set(2, clone(oldAtom)); struct.bonds.set(1, clone(oldBond));
  adapter.commit(struct, model);
  assert.equal(model.project.atoms[2].id, 'stable-c'); assert.deepEqual(model.project.atoms[2].pos3, [12, 8, 7]);
  assert.equal(model.project.bonds[1].id, 'stable-bc');
});

test('copy/paste of identical atoms ignores copied markers and allocates distinct UUIDs', () => {
  const {model, adapter, struct} = setup();
  const copy = clone(struct.atoms.get(2)); copy.__molecularCarEditorUUID = 'stable-c'; copy.pp.x += 2;
  const newPoolId = struct.atoms.add(copy);
  adapter.commit(struct, model);
  assert.notEqual(adapter.atoms.get(newPoolId).id, 'stable-c');
  assert.equal(new Set(model.project.atoms.map(a => a.id)).size, 4);
  assert.deepEqual(model.project.atoms[2].pos3, [12, 8, 7]);
});

test('external canonical undo/redo rebuilds the tool pool with explicit references and restores fragment state', () => {
  const {model, adapter, struct} = setup(), original = clone(model.project);
  struct.atoms.delete(1); struct.bonds.clear(); adapter.commit(struct, model);
  assert.deepEqual(model.project.fragments[0].atomIds, ['stable-a']);
  const edited = clone(model.project);
  model.undo(); assert.deepEqual(model.project, original);
  let payload = adapter.toKet(model.project), rebuilt = importer(payload.ket, 30);
  adapter.bindRebuild(rebuilt, model.project, payload.reference);
  assert.equal(adapter.atoms.get(31).id, 'stable-b');
  assert.equal(adapter.commit(rebuilt, model).changed, false);
  model.redo(); assert.deepEqual(model.project, edited);
  payload = adapter.toKet(model.project); rebuilt = importer(payload.ket, 70);
  adapter.bindRebuild(rebuilt, model.project, payload.reference);
  assert.equal(adapter.atoms.get(71).id, 'stable-c');
  assert.equal(adapter.commit(rebuilt, model).changed, false);
});

test('unexpected new pools are rejected even if chemistry and every atom coordinate are identical', () => {
  const {model, adapter, struct} = setup(), other = importer(adapter.toKet(model.project).ket);
  assert.throws(() => adapter.read(other, model.project), /必须从共享项目重新推送/);
  const payload = adapter.toKet(model.project);
  const duplicate = clone(payload.reference); duplicate.atomIds[2] = duplicate.atomIds[0];
  assert.throws(() => adapter.bindRebuild(other, model.project, duplicate), /原子数量或导入引用/);
  other.bonds.get(0).end = 2;
  assert.throws(() => adapter.bindRebuild(other, model.project, payload.reference), /键端点/);
  assert.equal(adapter.commit(struct, model).changed, false);
});

test('unsupported molecular features cause explicit rejection without changing canonical state', () => {
  for (const mutate of [s => s.sgroups.add({type: 'SUP'}), s => { s.atoms.get(0).hCount = 2; },
    s => { s.atoms.get(0).queryProperties = {aromaticity: 'aromatic'}; }, s => { s.bonds.get(0).type = 9; },
    s => { s.atoms.get(0).radical = 1; }]) {
    const {model, adapter, struct} = setup(), before = clone(model.project);
    mutate(struct); assert.throws(() => adapter.commit(struct, model), /不支持|不能无损|未应用/);
    assert.deepEqual(model.project, before); assert.equal(model.revision, 0);
  }
});

test('AND/OR enhanced stereo edits are rejected atomically instead of being kept only in a lossy sidecar', () => {
  // Actual Ketcher 3.18 labels are '&<number>' and 'or<number>', not 'AND'/'OR'.
  for (const label of ['&1', '&12', 'or1', 'or42', '&', 'or', 'unknown']) {
    const {model, adapter, struct} = setup(), before = clone(model.project);
    const sidecarBefore = adapter.snapshotSidecar(struct);
    struct.atoms.get(0).stereoLabel = label;
    // Even simultaneous coordinate/new-atom edits must remain uncommitted.
    struct.atoms.get(0).pp.x = 98;
    addAtom(struct, 5, 2);
    assert.throws(() => adapter.commit(struct, model), /不能无损保存 AND\/OR 增强立体组/);
    assert.deepEqual(model.project, before);
    assert.equal(model.revision, 0);
    assert.equal(adapter.atoms.size, before.atoms.length);
    struct.atoms.delete(3);
    assert.deepEqual(adapter.snapshotSidecar(struct), sidecarBefore);
    struct.atoms.get(0).stereoLabel = null;
    struct.atoms.get(0).pp.x = 0;
    assert.equal(adapter.commit(struct, model).changed, false);
  }
});

test('ABS and empty stereo labels remain accepted through native serialization and explicit pool rebuild', () => {
  for (const label of ['abs', '', null, undefined]) {
    const {model, adapter, struct} = setup(), positions = clone(model.project.atoms.map(a => a.pos3));
    struct.atoms.get(0).stereoLabel = label;
    struct.atoms.get(0).pp.x += 0.25;
    assert.equal(adapter.commit(struct, model).changed, true);
    assert.deepEqual(model.project.atoms.map(a => a.pos3), positions);
    assert.equal(model.project.metadata.ketcherAdapter.atomData['stable-a'].stereoLabel, label ?? null);
    const saved = parseProject(serializeProject(model.project)), payload = adapter.toKet(saved);
    assert.equal(payload.ket.mol0.atoms[0].stereoLabel, label === 'abs' ? 'abs' : undefined);
    const rebuilt = importer(payload.ket, 21);
    adapter.bindRebuild(rebuilt, saved, payload.reference);
    assert.equal(adapter.commit(rebuilt, model).changed, false);
  }
});

test('enhanced stereo in an old sidecar or an unexpected imported pool cannot bypass live-read rejection', () => {
  for (const label of ['&2', 'or3']) {
    const {model, adapter, struct} = setup(), saved = clone(model.project);
    saved.metadata = {ketcherAdapter: {format: KETCHER_SIDECAR_FORMAT, atomData: {
      'stable-a': {stereoLabel: label, canonicalChiralTag: 'CHI_TETRAHEDRAL_CW'},
    }}};
    assert.throws(() => adapter.toKet(saved), /不能无损保存 AND\/OR/);
    // A tag mismatch must not silently strip an unsupported saved declaration.
    saved.atoms[0].chiralTag = 'CHI_UNSPECIFIED';
    assert.throws(() => adapter.toKet(saved), /不能无损保存 AND\/OR/);
    const payload = adapter.toKet(model.project), imported = importer(payload.ket, 30);
    imported.atoms.get(30).stereoLabel = label;
    assert.throws(() => adapter.bindRebuild(imported, model.project, payload.reference), /不能无损保存 AND\/OR/);
    assert.equal(adapter.atomPool, struct.atoms);
    assert.equal(adapter.commit(struct, model).changed, false);
    assert.equal(model.revision, 0);
  }
});

test('explicit wedge edits invalidate stale atom tags and retain the new drawn declaration', () => {
  const {model, adapter, struct} = setup();
  struct.bonds.get(0).stereo = 6;
  adapter.commit(struct, model);
  assert.equal(model.project.bonds[0].direction, 'BEGINDASH');
  assert.equal(model.project.atoms[0].chiralTag, 'CHI_UNSPECIFIED');
  assert.match(model.project.metadata.stereoNotice, /楔形键已修改/);
  assert.equal(adapter.commit(struct, model).changed, false);
});

test('successive toolkit stereo declarations remain durable sidecar data even after the canonical tag was cleared', () => {
  const {model, adapter, struct} = setup();
  struct.atoms.get(0).stereoParity = 1;
  adapter.commit(struct, model);
  assert.equal(model.project.atoms[0].chiralTag, 'CHI_UNSPECIFIED');
  struct.atoms.get(0).stereoParity = 2;
  assert.equal(adapter.commit(struct, model).changed, true);
  assert.equal(model.project.metadata.ketcherAdapter.atomData['stable-a'].stereoParity, 2);
  assert.equal(adapter.toKet(model.project).ket.mol0.atoms[0].stereoParity, 2);
  assert.equal(adapter.commit(struct, model).changed, false);
  model.undo(); assert.equal(model.project.metadata.ketcherAdapter.atomData['stable-a'].stereoParity, 1);
});

test('C60 copying gives 60 fresh IDs, while changing its original 2D layout preserves every 3D vertex', () => {
  const library = JSON.parse(readFileSync(new URL('../frontend/editor/fragments.json', import.meta.url), 'utf8'));
  const c60 = library.fragments.find(f => f.id === 'c60');
  const canonical = new WorkshopModel(); canonical.addFragment(c60);
  const {model, adapter, struct} = setup(canonical.project), old = clone(model.project);
  const remap = new Map();
  for (const [poolId, atom] of [...struct.atoms]) { const copied = clone(atom); copied.pp.x += 20; remap.set(poolId, struct.atoms.add(copied)); }
  for (const [, bond] of [...struct.bonds]) struct.bonds.add({...bond, begin: remap.get(bond.begin), end: remap.get(bond.end)});
  struct.atoms.get(0).pp.x += .2;
  adapter.commit(struct, model);
  assert.equal(model.project.atoms.length, 120); assert.equal(model.project.bonds.length, 180);
  assert.equal(new Set(model.project.atoms.map(a => a.id)).size, 120);
  assert.deepEqual(model.project.atoms.slice(0, 60).map(a => a.pos3), old.atoms.map(a => a.pos3));
  assert.deepEqual(model.project.fragments, old.fragments);
});

test('native drawing cannot bypass C60 port valence; explicit reduction opens the attachment', () => {
  const library = JSON.parse(readFileSync(new URL('../frontend/editor/fragments.json', import.meta.url), 'utf8'));
  const canonical = new WorkshopModel(); canonical.addFragment(library.fragments.find(f => f.id === 'c60'));
  const {model, adapter, struct} = setup(canonical.project), original = clone(model.project);
  const newPoolId = addAtom(struct, 12, 12); struct.bonds.add({begin: 0, end: newPoolId, type: 1, stereo: 0});
  assert.throws(() => adapter.commit(struct, model), /C60 等饱和位点/); assert.deepEqual(model.project, original);
  const [, double] = [...struct.bonds].find(([, b]) => b.type === 2 && [b.begin, b.end].includes(0)); double.type = 1;
  adapter.commit(struct, model);
  assert.equal(model.project.atoms.length, 61); assert.equal(model.project.bonds.length, 91);
  assert.deepEqual(model.project.fragments, original.fragments);
});

test('declared no-port fragment blocks external connection while free drawing accepts incomplete valence', () => {
  const {model, adapter, struct} = setup();
  const next = addAtom(struct, -1.5); struct.bonds.add({begin: 0, end: next, type: 1, stereo: 0});
  assert.throws(() => adapter.commit(struct, model), /不是声明的连接位点/);
  const free = fixture(); free.fragments = [];
  const f = setup(free); f.struct.bonds.get(0).type = 3; f.struct.bonds.get(1).type = 3;
  assert.equal(f.adapter.commit(f.struct, f.model).changed, true);
});

test('a bound empty tool pool allocates new UUIDs and never pretends it already has a generated conformer', () => {
  const {model, adapter, struct} = setup(createProject());
  addAtom(struct, 5, -2); adapter.commit(struct, model);
  assert.deepEqual(model.project.atoms[0].pos2, [5, 2]);
  assert.deepEqual(model.project.atoms[0].pos3, [5, 2, 0]); assert.equal(model.project.has3d, false);
});

test('programmatic asynchronous pushes are serialized, change events suppressed, and final pool explicitly rebound', async () => {
  const {model, adapter} = setup(), entered = [], releases = [];
  let live;
  const ketcher = {editor: {struct: () => live}, async setMolecule(text, options) {
    assert.deepEqual(options, {preserveCanvasPosition: true, needZoom: false});
    const ket = JSON.parse(text); entered.push(ket.mol0.atoms.length);
    await new Promise(resolve => releases.push(resolve)); live = importer(ket, entered.length * 10);
    assert.equal(adapter.commit(live, model).suppressed, true);
  }};
  const first = adapter.push(ketcher, model.project), changed = clone(model.project); changed.atoms[0].pos2 = [33, 44];
  const second = adapter.push(ketcher, changed);
  await new Promise(resolve => setImmediate(resolve)); assert.deepEqual(entered, [3]); assert.equal(adapter.suppressing, true);
  releases.shift()(); await first; await new Promise(resolve => setImmediate(resolve)); assert.deepEqual(entered, [3, 3]);
  releases.shift()(); await second;
  assert.equal(adapter.suppressing, false); assert.equal(adapter.atoms.get(20).id, 'stable-a');
  assert.deepEqual(adapter.read(live, changed).atoms[0].pos2, [33, 44]); assert.equal(model.revision, 0);
});

function c60Project(){
  const atomIds=C60_POS3.map((_,i)=>'sphere-'+i);
  return {...createProject('native SUP'),has3d:true,atoms:C60_POS3.map((point,i)=>({id:atomIds[i],element:'C',charge:0,pos2:[i*0.1,i%5],pos3:[...point]})),bonds:C60_BONDS.map((b,i)=>({id:'sphere-bond-'+i,a:atomIds[b.a],b:atomIds[b.b],order:b.order})),metadata:{ketcherGroups:[{id:'native-group',type:'SUP',name:'C60',kind:'c60',atomIds,expanded:false,attachmentPoints:[{atomId:atomIds[0],attachmentId:'1'}],template:{id:'c60',version:2,indexToAtomId:atomIds}}]}};
}
function insertNativeC60(struct,{marker=true}={}){
  const pools=C60_POS3.map((_,i)=>addAtom(struct,20+i*0.1,i%5));
  for(const bond of C60_BONDS)struct.bonds.add({begin:pools[bond.a],end:pools[bond.b],type:bond.order,stereo:0});
  const data={name:'C60',expanded:false};if(marker)data.molecularCarTemplate={id:'c60',version:2};
  struct.sgroups.add({type:'SUP',isMonomer:false,data,atoms:pools,getAttachmentPoints(){return[{atomId:pools[0],attachmentId:'1'}];}});
  return pools;
}
test('native SUP round-trip binds ordered member UUIDs, selection, collapse state and 3D seed without graph matching',()=>{
  const project=c60Project(),{adapter,model,struct}=setup(project,24);
  const payload=adapter.toKet(project);
  assert.deepEqual(payload.ket.mol0.sgroups[0].atoms,Array.from({length:60},(_,i)=>i));
  assert.equal(struct.sgroups.get(0).data.expanded,false);
  assert.equal(adapter.commit(struct,model).changed,false);
  assert.deepEqual([...adapter.selectionIds(struct,{sgroups:[0]})],project.metadata.ketcherGroups[0].atomIds);
  struct.sgroups.get(0).data.expanded=true;adapter.commit(struct,model);
  assert.equal(model.project.metadata.ketcherGroups[0].expanded,true);
  assert.deepEqual(model.project.atoms.map(a=>a.pos3),project.atoms.map(a=>a.pos3));
  assert.equal(model.undo(),true);assert.equal(model.project.metadata.ketcherGroups[0].expanded,false);
  assert.equal(model.redo(),true);assert.equal(model.project.metadata.ketcherGroups[0].expanded,true);
});
test('new native C60 source marker plus explicit ordered Pool membership initializes a real spherical cage',()=>{
  const {adapter,model,struct}=setup(createProject()),pools=insertNativeC60(struct);
  adapter.commit(struct,model);const group=model.project.metadata.ketcherGroups[0];
  assert.equal(group.kind,'c60');assert.deepEqual(group.template.indexToAtomId,pools.map(id=>adapter.atoms.get(id).id));
  const lookup=new Map(model.project.atoms.map(a=>[a.id,a]));
  assert.ok(Math.max(...model.project.atoms.map(a=>a.pos3[2]))-Math.min(...model.project.atoms.map(a=>a.pos3[2]))>6);
  for(const bond of model.project.bonds)assert.ok(Math.abs(Math.hypot(...lookup.get(bond.a).pos3.map((v,i)=>v-lookup.get(bond.b).pos3[i]))-1.42)<1e-12);
  assert.equal(model.project.has3d,false,'a physical seed is not a completed generated conformer');
});
test('a user SUP named C60 never gains template identity or a guessed spherical mapping',()=>{
  const {adapter,model,struct}=setup(createProject());insertNativeC60(struct,{marker:false});adapter.commit(struct,model);
  assert.equal(model.project.metadata.ketcherGroups[0].kind,undefined);
  assert.equal(model.project.metadata.ketcherGroups[0].template,undefined);
  assert.ok(model.project.atoms.every(a=>a.pos3[2]===0));
});
test('copying native C60 creates a distinct group and sixty UUIDs while retaining every old conformer coordinate',()=>{
  const {adapter,model,struct}=setup(c60Project()),before=clone(model.project),pools=insertNativeC60(struct);adapter.commit(struct,model);
  assert.equal(model.project.metadata.ketcherGroups.length,2);
  assert.notEqual(model.project.metadata.ketcherGroups[1].id,'native-group');
  assert.ok(pools.every(i=>!before.atoms.some(a=>a.id===adapter.atoms.get(i).id)));
  assert.deepEqual(model.project.atoms.slice(0,60),before.atoms);
});
test('native C60 single-bond attachment is one announced undoable addition and never leaves a five-valent cage carbon',()=>{
  const {adapter,model,struct}=setup(c60Project()),before=clone(model.project),outside=addAtom(struct,15);
  struct.bonds.add({begin:0,end:outside,type:1,stereo:0});
  const result=adapter.commit(struct,model),anchor=model.project.atoms[0].id;
  assert.equal(result.requiresPush,true);assert.match(result.warnings.join(' '),/加成连接/);
  const sum=model.project.bonds.filter(b=>b.a===anchor||b.b===anchor).reduce((n,b)=>n+b.order,0);
  assert.equal(sum,4);assert.equal(model.project.bonds.filter(b=>b.order===2).length,29);
  assert.deepEqual(model.project.atoms.slice(0,60).map(a=>a.pos3),before.atoms.map(a=>a.pos3));
  model.undo();assert.deepEqual(model.project,before);model.redo();assert.equal(model.project.atoms.length,61);
});
test('a C60 connection into another saturated atom is rejected atomically with a useful reason',()=>{
  const project=c60Project();project.atoms.push({id:'outside',element:'C',charge:0,pos2:[20,0],pos3:[20,0,0]});
  for(let i=0;i<4;i++){project.atoms.push({id:'h'+i,element:'H',charge:0,pos2:[20,i+1],pos3:[20,i+1,0]});project.bonds.push({id:'oh'+i,a:'outside',b:'h'+i,order:1});}
  const {adapter,model,struct}=setup(project),before=clone(model.project);struct.bonds.add({begin:0,end:60,type:1,stereo:0});
  assert.throws(()=>adapter.commit(struct,model),/位点已满价.*显式氢/);assert.deepEqual(model.project,before);assert.equal(model.revision,0);
});
test('C60 connection may replace only an automatically generated H and undo restores it together with cage bond order',()=>{
  const project=c60Project();project.atoms.push({id:'outside',element:'C',charge:0,pos2:[20,0],pos3:[20,0,0]});
  for(let i=0;i<4;i++){project.atoms.push({id:'h'+i,element:'H',charge:0,generatedHydrogen:true,pos2:[20,i+1],pos3:[20,i+1,0]});project.bonds.push({id:'oh'+i,a:'outside',b:'h'+i,order:1});}
  const {adapter,model,struct}=setup(project),before=clone(model.project);struct.bonds.add({begin:0,end:60,type:1,stereo:0});
  const result=adapter.commit(struct,model);assert.equal(result.requiresPush,true);assert.equal(model.project.atoms.length,64);assert.match(model.project.metadata.hydrogenNotice,/替代了 1 个/);
  assert.equal(model.project.bonds.filter(b=>b.a==='outside'||b.b==='outside').length,4);
  model.undo();assert.deepEqual(model.project,before);
});
