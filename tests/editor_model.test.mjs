import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync, existsSync} from 'node:fs';
import {spawnSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import {join} from 'node:path';

// The browser owns model.js as an ES module; no root package.json or app-wide
// Node configuration is required just to test this standalone editor module.
const source = readFileSync(new URL('../frontend/editor/model.js', import.meta.url), 'utf8') + '\n//# sourceURL=workshop-model.js';
const {WorkshopModel, createProject, clone, validateProject, parseProject, serializeProject,
  chemistryWarnings, selectedIds, distance, angle, dihedral} = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);

const close = (actual, expected, tolerance = 1e-8) => assert.ok(Math.abs(actual - expected) < tolerance, `${actual} != ${expected}`);
const pointsClose = (a, b) => a.forEach((v, i) => close(v, b[i]));
function graph(points, edges, elements = []) {
  return {...createProject('test molecule'), has3d: true,
    atoms: points.map((pos3, i) => ({id: `a${i}`, element: elements[i] ?? 'C', charge: 0, pos2: [i * 1.5, i % 2], pos3: [...pos3]})),
    bonds: edges.map(([a, b, order = 1], i) => ({id: `b${i}`, a: `a${a}`, b: `a${b}`, order}))};
}
function chain() { return graph([[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 2, 0]], [[0, 1], [1, 2], [2, 3], [3, 4]]); }
function hexagon() {
  return graph(Array.from({length: 6}, (_, i) => [1.4 * Math.cos(i * Math.PI / 3), 1.4 * Math.sin(i * Math.PI / 3), 0]),
    Array.from({length: 6}, (_, i) => [i, (i + 1) % 6, 1.5]));
}

test('native project preserves stable IDs, both layouts, stereo, isotope, explicit H and editor metadata', () => {
  const p = chain();
  Object.assign(p.atoms[0], {charge: 1, isotope: 13, chiralTag: 'CHI_TETRAHEDRAL_CW', explicitHs: 1, noImplicit: true, radicalElectrons: 0});
  Object.assign(p.bonds[1], {order: 2, stereo: 'STEREOZ', direction: 'EITHERDOUBLE', stereoAtoms: ['a0', 'a3']});
  p.fragments = [{id: 'cage', name: 'custom', atomIds: ['a2', 'a3', 'a4'], ports: [{atomId: 'a2', maxOrder: 1}], rigid: true}];
  p.metadata = {toolMappings: {rdkit: ['a4', 'a2', 'a3', 'a1', 'a0']}, notes: '保存设计'};
  const saved = parseProject(serializeProject(p));
  assert.deepEqual(saved, p);
  assert.notEqual(saved.atoms[0], p.atoms[0]);
});

test('validation rejects broken identity mappings while accepting incomplete chemical structures', () => {
  const p = chain(); p.atoms[1].id = p.atoms[0].id;
  assert.throws(() => validateProject(p), /ID 重复/);
  const badBond = chain(); badBond.bonds[0].a = 'missing';
  assert.throws(() => validateProject(badBond), /不存在/);
  const invalidCoordinates = chain(); invalidCoordinates.atoms[0].pos3[0] = NaN;
  assert.throws(() => validateProject(invalidCoordinates), /有限坐标/);
  const m = new WorkshopModel(graph([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [[0, 1, 3], [0, 2, 3]]));
  assert.ok(chemistryWarnings(m.project).some(message => message.includes('超过')));
  assert.equal(parseProject(serializeProject(m.project)).atoms.length, 3);
  assert.equal(validateProject(createProject()).atoms.length, 0);
});

test('RDKit empty stereo reference lists and explicitly disabled fragment ports remain representable', () => {
  const p = chain(); p.bonds.forEach(b => { b.stereoAtoms = []; });
  p.fragments = [{id: 'f', name: 'disabled port', atomIds: ['a0'], ports: [{atomId: 'a0', maxOrder: 0}], rigid: true}];
  assert.deepEqual(validateProject(p), p);
  const m = new WorkshopModel(p);
  assert.throws(() => m.connect('a0', 'a4'), /最多允许 0/);
});

test('2D rearrangement and 3D atom motion never overwrite the other layout or perceive bonds from distance', () => {
  const m = new WorkshopModel(chain());
  const original = clone(m.project);
  m.selection.add('a0');
  m.moveSelection([50, -10], {dimension: 2});
  assert.deepEqual(m.project.atoms.map(a => a.pos3), original.atoms.map(a => a.pos3));
  const arranged2d = m.project.atoms.map(a => a.pos2);
  m.moveSelection([100, 0, 17]);
  assert.deepEqual(m.project.atoms.map(a => a.pos2), arranged2d);
  assert.deepEqual(m.project.bonds, original.bonds);
  assert.deepEqual(m.project.atoms.slice(1).map(a => a.pos3), original.atoms.slice(1).map(a => a.pos3));
});

test('whole rigid fragment motion expands selection; explicit override permits individual atom editing', () => {
  const p = chain();
  p.fragments = [{id: 'rigid', name: 'cage', atomIds: ['a2', 'a3', 'a4'], ports: [{atomId: 'a2'}], rigid: true}];
  const m = new WorkshopModel(p); m.selection.add('a3');
  assert.deepEqual(m.expandedSelection(), new Set(['a3', 'a2', 'a4']));
  m.moveSelection([1, 2, 3]);
  for (let i = 0; i < 5; i++) pointsClose(m.project.atoms[i].pos3, i >= 2 ? p.atoms[i].pos3.map((v, k) => v + [1, 2, 3][k]) : p.atoms[i].pos3);
  const before = clone(m.project);
  m.moveSelection([0, 0, 1], {rigid: false});
  assert.deepEqual(m.project.atoms[2].pos3, before.atoms[2].pos3);
  close(m.project.atoms[3].pos3[2], before.atoms[3].pos3[2] + 1);
});

test('rotation preserves rigid fragment distances and only changes 3D coordinates', () => {
  const m = new WorkshopModel(chain()), before = clone(m.project);
  m.selection = new Set(['a2', 'a3', 'a4']);
  m.rotateSelection([1, 0, 0], 90, {origin: [1, 0, 0]});
  pointsClose(m.project.atoms[3].pos3, [1, 0, 1]);
  pointsClose(m.project.atoms[4].pos3, [1, 0, 2]);
  close(distance(m.project, 'a3', 'a4'), distance(before, 'a3', 'a4'));
  assert.deepEqual(m.project.atoms.map(a => a.pos2), before.atoms.map(a => a.pos2));
  assert.deepEqual(m.project.bonds, before.bonds);
});

test('distance edit moves exactly the connected side without changing adjacent distances', () => {
  const m = new WorkshopModel(chain()), before = clone(m.project);
  m.setDistance('a1', 'a2', 1.8);
  close(distance(m.project, 'a1', 'a2'), 1.8);
  pointsClose(m.project.atoms[0].pos3, before.atoms[0].pos3);
  pointsClose(m.project.atoms[1].pos3, before.atoms[1].pos3);
  close(distance(m.project, 'a2', 'a4'), distance(before, 'a2', 'a4'));
  assert.deepEqual(m.project.bonds, before.bonds);
});

test('angle edit reaches target and leaves fixed side plus moving fragment internals intact', () => {
  const m = new WorkshopModel(chain()), before = clone(m.project);
  m.setAngle('a0', 'a1', 'a2', 125);
  close(angle(m.project, 'a0', 'a1', 'a2'), 125);
  close(distance(m.project, 'a2', 'a4'), distance(before, 'a2', 'a4'));
  close(distance(m.project, 'a1', 'a2'), distance(before, 'a1', 'a2'));
  pointsClose(m.project.atoms[0].pos3, before.atoms[0].pos3);
  pointsClose(m.project.atoms[1].pos3, before.atoms[1].pos3);
});

test('angle adjustment handles a collinear initial geometry without NaNs', () => {
  const m = new WorkshopModel(graph([[-1, 0, 0], [0, 0, 0], [1, 0, 0]], [[0, 1], [1, 2]]));
  m.setAngle('a0', 'a1', 'a2', 109.5);
  close(angle(m.project, 'a0', 'a1', 'a2'), 109.5);
});

test('dihedral edit is real 3D motion, preserves topology, survives 2D edit, undo and redo', () => {
  const p = chain(); Object.assign(p.atoms[2], {chiralTag: 'CHI_TETRAHEDRAL_CCW'});
  const m = new WorkshopModel(p);
  for (const target of [65, -130, 179, -5]) {
    m.setDihedral('a0', 'a1', 'a2', 'a3', target);
    close(dihedral(m.project, 'a0', 'a1', 'a2', 'a3'), target);
    close(distance(m.project, 'a2', 'a4'), distance(p, 'a2', 'a4'));
  }
  assert.deepEqual(m.project.bonds, p.bonds);
  const torsion = clone(m.project);
  m.selection.add('a4'); m.moveSelection([10, 20], {dimension: 2});
  assert.deepEqual(m.project.atoms.map(a => a.pos3), torsion.atoms.map(a => a.pos3));
  m.undo(); assert.deepEqual(m.project, torsion);
  m.redo(); close(dihedral(m.project, 'a0', 'a1', 'a2', 'a3'), -5);
  assert.equal(m.project.atoms[2].chiralTag, 'CHI_TETRAHEDRAL_CCW');
});

test('ring bond edits and degenerate torsions fail transactionally with useful messages', () => {
  const m = new WorkshopModel(hexagon()), before = serializeProject(m.project);
  assert.throws(() => m.setDistance('a0', 'a1', 2), /环/);
  assert.throws(() => m.setAngle('a0', 'a1', 'a2', 90), /环/);
  assert.throws(() => m.setDihedral('a0', 'a1', 'a2', 'a3', 30), /环/);
  assert.equal(serializeProject(m.project), before); assert.equal(m.revision, 0);
  const linear = new WorkshopModel(graph([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], [[0, 1], [1, 2], [2, 3]]));
  assert.throws(() => linear.setDihedral('a0', 'a1', 'a2', 'a3', 45), /几何方向/);
});

test('fragment insertion remaps stable IDs and stereo references without graph matching', () => {
  const fragment = chain(); fragment.name = 'custom'; fragment.ports = [{atomId: 'a0', maxOrder: 1}];
  fragment.bonds[1].order = 2; fragment.bonds[1].stereo = 'STEREOE'; fragment.bonds[1].stereoAtoms = ['a0', 'a3'];
  const m = new WorkshopModel();
  const first = m.addFragment(fragment, {pos3: [5, 0, 0]});
  const second = m.addFragment(fragment, {pos3: [-5, 0, 0]});
  assert.equal(new Set(m.project.atoms.map(a => a.id)).size, 10);
  assert.equal(first.ports[0].atomId, first.idMap.a0);
  assert.deepEqual(m.project.bonds[1].stereoAtoms, [first.idMap.a0, first.idMap.a3]);
  assert.deepEqual(m.project.bonds[5].stereoAtoms, [second.idMap.a0, second.idMap.a3]);
  pointsClose(m.project.atoms[0].pos3, [5, 1, 0]);
  assert.equal(m.project.fragments.length, 2);
});

test('explicit fragment connection enforces declared ports and local valence; free drawing remains available', () => {
  const benzene = hexagon(); benzene.name = '苯环'; benzene.ports = [{atomId: 'a0', maxOrder: 1}];
  const m = new WorkshopModel(); const ring = m.addFragment(benzene);
  const carbon = m.addAtom({pos3: [4, 0, 0]});
  assert.throws(() => m.connect(ring.idMap.a1, carbon), /连接位点/);
  m.connect(ring.idMap.a0, carbon);
  const another = m.addAtom({pos3: [5, 0, 0]});
  assert.throws(() => m.connect(ring.idMap.a0, another), /足够价态/);
  m.addBond(ring.idMap.a0, another);
  assert.ok(chemistryWarnings(m.project).length > 0);
});

test('full history restores deletion, stereo, fragment ports and coordinates together', () => {
  const p = chain();
  p.bonds[1].stereo = 'STEREOZ'; p.bonds[1].stereoAtoms = ['a0', 'a3'];
  p.fragments = [{id: 'f', name: 'whole', atomIds: p.atoms.map(a => a.id), ports: [{atomId: 'a3', maxOrder: 1}], rigid: true}];
  const m = new WorkshopModel(p); m.selection = new Set(['a3']);
  m.deleteSelection(); const deleted = clone(m.project);
  assert.equal(m.project.atoms.length, 4); assert.deepEqual(m.project.fragments[0].ports, []);
  assert.equal(m.project.bonds[1].stereo, 'STEREONONE'); assert.match(m.project.metadata.stereoNotice, /重新确认/);
  m.undo(); assert.deepEqual(m.project, p); assert.deepEqual(m.selection, new Set(['a3']));
  m.redo(); assert.deepEqual(m.project, deleted);
  m.undo(); m.addAtom(); assert.equal(m.canRedo, false);
});

test('revision gates reject obsolete async results even after undo reaches old coordinates', () => {
  const m = new WorkshopModel(chain()), original = clone(m.project), start = m.revision;
  const optimized = clone(original); optimized.atoms[0].pos3 = [7, 8, 9];
  m.selection.add('a4'); m.moveSelection([1, 0, 0]); m.undo();
  assert.deepEqual(m.project, original); assert.ok(m.revision > start);
  assert.equal(m.applyAsync(start, optimized), false); assert.deepEqual(m.project, original);
  assert.equal(m.applyAsync(m.revision, optimized), true); assert.deepEqual(m.project, optimized);
  m.undo(); assert.deepEqual(m.project, original);
  m.redo(); assert.deepEqual(m.project, optimized);
});

test('transaction failures never mutate state/history; change events produce no feedback loops', () => {
  const m = new WorkshopModel(chain()), before = clone(m.project); let count = 0;
  m.addEventListener('change', event => { count++; assert.equal(event.detail.revision, m.revision); });
  assert.throws(() => m.commit('broken', p => { p.atoms[0].pos3[0] = NaN; }), /有限坐标/);
  assert.deepEqual(m.project, before); assert.equal(m.revision, 0); assert.equal(m.canUndo, false);
  m.commit('no change', () => {}); assert.equal(count, 0);
  m.updateAtom('a0', {charge: 1}); m.undo(); m.redo(); assert.equal(count, 3);
  assert.throws(() => m.updateAtom('a0', {id: 'replacement'}), /稳定原子 ID/);
});

test('local stereo inversion is announced without silently changing the saved tetrahedral tag', () => {
  const p = graph([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 1], [0, 2], [0, 3]]);
  p.atoms[0].chiralTag = 'CHI_TETRAHEDRAL_CW';
  const m = new WorkshopModel(p); m.selection.add('a3');
  m.moveSelection([0, 0, -2]);
  assert.equal(m.project.atoms[0].chiralTag, 'CHI_TETRAHEDRAL_CW');
  assert.ok(chemistryWarnings(m.project).some(w => w.includes('手性方向已翻转')));
  m.undo(); assert.deepEqual(m.project, p);
  m.selection = new Set(p.atoms.map(a => a.id)); m.rotateSelection('y', 180);
  assert.equal(chemistryWarnings(m.project).length, 0);
});

test('deleting a stereo reference bond clears stale E/Z with a visible notice; invalid new references are rejected', () => {
  const p = chain(); Object.assign(p.bonds[1], {order: 2, stereo: 'STEREOZ', stereoAtoms: ['a0', 'a3']});
  const m = new WorkshopModel(p);
  m.commit('delete reference bond', draft => { draft.bonds = draft.bonds.filter(b => b.id !== 'b0'); });
  assert.equal(m.project.bonds[0].stereo, 'STEREONONE');
  assert.deepEqual(m.project.bonds[0].stereoAtoms, []);
  assert.ok(chemistryWarnings(m.project).some(w => w.includes('双键立体参考失效')));
  m.undo(); assert.deepEqual(m.project, p);
  assert.throws(() => m.updateBond('b1', {stereoAtoms: ['a4', 'a3']}), /连续原子/);
  m.updateBond('b1', {order: 1});
  assert.equal(m.project.bonds[1].stereo, 'STEREONONE');
});

test('changing a tetrahedral neighborhood clears the order-dependent tag, including delete then re-add', () => {
  const p = graph([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, -1, -1]], [[0, 1], [0, 2], [0, 3], [0, 4]], ['C', 'F', 'Cl', 'Br', 'I']);
  p.atoms[0].chiralTag = 'CHI_TETRAHEDRAL_CW'; p.bonds[0].direction = 'BEGINWEDGE';
  const m = new WorkshopModel(p);
  m.deleteBond('b0'); m.addBond('a0', 'a1');
  assert.equal(m.project.atoms[0].chiralTag, 'CHI_UNSPECIFIED');
  assert.ok(chemistryWarnings(m.project).some(w => w.includes('请重新指定手性')));
  m.undo(); m.undo(); assert.deepEqual(m.project, p);
  m.updateAtom('a1', {element: 'H'});
  assert.equal(m.project.atoms[0].chiralTag, 'CHI_UNSPECIFIED');
  assert.equal(m.project.bonds[0].direction, 'NONE');
  m.undo(); m.updateBond('b1', {order: 2});
  assert.equal(m.project.atoms[0].chiralTag, 'CHI_UNSPECIFIED');
  m.undo(); m.commit('reorder incident bonds', draft => { draft.bonds.reverse(); });
  assert.equal(m.project.atoms[0].chiralTag, 'CHI_UNSPECIFIED');
});

test('rigid selection expands transitively but filters stale/deleted IDs', () => {
  const p = chain();
  p.fragments = [
    {id: 'f1', name: 'first', atomIds: ['a0', 'a1'], ports: [], rigid: true},
    {id: 'f2', name: 'second', atomIds: ['a1', 'a2'], ports: [], rigid: true},
    {id: 'f3', name: 'not rigid', atomIds: ['a2', 'a3'], ports: [], rigid: false},
  ];
  assert.deepEqual(selectedIds(p, new Set(['a0', 'missing'])), new Set(['a0', 'a1', 'a2']));
  assert.deepEqual(selectedIds(p, new Set(['a0', 'missing']), false), new Set(['a0']));
});

test('all shipped JSON fragments load with independent atom IDs and real graph counts', () => {
  const library = JSON.parse(readFileSync(new URL('../frontend/editor/fragments.json', import.meta.url), 'utf8'));
  const m = new WorkshopModel();
  for (const template of library.fragments) {
    const first = m.addFragment(template), second = m.addFragment(template);
    assert.equal(first.atomIds.length, template.atoms.length);
    assert.equal(second.atomIds.length, template.atoms.length);
    assert.equal(new Set([...first.atomIds, ...second.atomIds]).size, template.atoms.length * 2);
    if (template.id === 'c60') {
      assert.equal(template.atoms.length, 60); assert.equal(template.bonds.length, 90);
      assert.ok(first.ports.every(p => p.maxOrder === 1));
      assert.equal(selectedIds(m.project, [first.atomIds[0]]).size, 60);
    }
  }
  assert.equal(m.project.fragments.length, library.fragments.length * 2);
  assert.deepEqual(parseProject(serializeProject(m.project)), m.project);
});

const root = fileURLToPath(new URL('../', import.meta.url));
const python = [join(root, '.venv', 'Scripts', 'python.exe'), join(root, '.venv', 'bin', 'python')].find(existsSync);
function callGeometryService(code, input) {
  const result = spawnSync(python, ['-E', '-B', '-c', code], {cwd: root, encoding: 'utf8', input: input === undefined ? undefined : JSON.stringify(input), maxBuffer: 4 * 1024 * 1024});
  assert.equal(result.status, 0, result.stderr);
  return JSON.parse(result.stdout);
}

test('integration: imported tetraphenylbenzene torsion survives model edit, 2D movement and actual MOL reimport', {skip: !python}, () => {
  const fixture = callGeometryService(`
import json
from pathlib import Path
from molecular_car.editor import dispatch
from molecular_car.editor.service import to_mol
p = dispatch('/api/editor/import', {'filename': 'tetraphenylbenzene.mol', 'content': Path('examples/tetraphenylbenzene.mol').read_text()})['project']
m = to_mol(p, sanitize=True)
torsion = None
for bond in m.GetBonds():
    if bond.GetBondTypeAsDouble() == 1 and not bond.IsInRing():
        b, c = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        aa = [a.GetIdx() for a in m.GetAtomWithIdx(b).GetNeighbors() if a.GetIdx() != c and a.GetAtomicNum() != 1]
        dd = [a.GetIdx() for a in m.GetAtomWithIdx(c).GetNeighbors() if a.GetIdx() != b and a.GetAtomicNum() != 1]
        if aa and dd:
            torsion = [p['atoms'][i]['id'] for i in (aa[0], b, c, dd[0])]
            break
print(json.dumps({'project': p, 'torsion': torsion}))
`);
  assert.ok(fixture.torsion);
  const m = new WorkshopModel(fixture.project);
  m.setDihedral(...fixture.torsion, 73);
  close(dihedral(m.project, ...fixture.torsion), 73);
  m.selection.add(fixture.torsion[3]); m.moveSelection([7, -6], {dimension: 2});
  close(dihedral(m.project, ...fixture.torsion), 73);
  const exported = callGeometryService(`
import json, sys
from molecular_car.editor import dispatch
from molecular_car.editor.service import to_mol
from rdkit.Chem import rdMolTransforms
request = json.load(sys.stdin)
p = request['project']
indices = {a['id']: i for i, a in enumerate(p['atoms'])}
torsion = [indices[i] for i in request['torsion']]
rdkit_angle = rdMolTransforms.GetDihedralDeg(to_mol(p).GetConformer(), *torsion)
out = dispatch('/api/editor/export', {'project': p, 'dimension': 3})
q = dispatch('/api/editor/import', out)['project']
print(json.dumps({'project': q, 'header': out['content'].splitlines()[1], 'rdkitAngle': rdkit_angle}))
`, {project: m.project, torsion: fixture.torsion});
  assert.match(exported.header, /3D/);
  close(exported.rdkitAngle, 73);
  const q = validateProject(exported.project);
  assert.equal(q.atoms.length, m.project.atoms.length);
  q.atoms.forEach((atom, i) => atom.pos3.forEach((v, k) => close(v, m.project.atoms[i].pos3[k], 0.000051)));
  assert.deepEqual(q.bonds.map(b => b.order), m.project.bonds.map(b => b.order));
});

test('integration: actual C60 fragment graph has 60 atoms, 90 bonds, rigid movement and dynamic valence checks', {skip: !python}, () => {
  const fullerene = callGeometryService(`
import json
from molecular_car.editor.fragments import fullerene
print(json.dumps(fullerene()))
`);
  const m = new WorkshopModel(), added = m.addFragment(fullerene);
  assert.equal(added.atomIds.length, 60); assert.equal(m.project.bonds.length, 90);
  m.selection.add(added.atomIds[0]); assert.equal(m.expandedSelection().size, 60);
  const before = clone(m.project); m.moveSelection([3, 2, 1]);
  for (let i = 0; i < 60; ++i) pointsClose(m.project.atoms[i].pos3, before.atoms[i].pos3.map((v, k) => v + [3, 2, 1][k]));
  const outer = m.addAtom({pos3: [10, 10, 10]});
  assert.throws(() => m.connect(added.atomIds[0], outer), /足够价态/);
  const double = m.project.bonds.find(b => b.order === 2 && [b.a, b.b].includes(added.atomIds[0]));
  m.updateBond(double.id, {order: 1});
  m.connect(added.atomIds[0], outer);
  assert.equal(chemistryWarnings(m.project).length, 0);
  m.undo(); assert.equal(m.project.bonds.length, 90);
  m.undo(); assert.equal(m.project.bonds.find(b => b.id === double.id).order, 2);
});

test('integration: deleting/re-adding or reordering F[C@](Cl)(Br)I bonds cannot silently flip S to R', {skip: !python}, () => {
  const fixture = callGeometryService(`
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from molecular_car.editor import dispatch
m = Chem.MolFromSmiles('F[C@](Cl)(Br)I')
AllChem.Compute2DCoords(m)
p = dispatch('/api/editor/import', {'filename': 'stereo.mol', 'content': Chem.MolToMolBlock(m)})['project']
print(json.dumps(p))
`);
  const center = fixture.atoms.find(a => a.chiralTag && a.chiralTag !== 'CHI_UNSPECIFIED');
  assert.ok(center);
  for (const reorderOnly of [false, true]) {
    const m = new WorkshopModel(fixture);
    if (reorderOnly) m.commit('reorder', draft => { draft.bonds.push(draft.bonds.shift()); });
    else { const b = m.project.bonds[0]; m.deleteBond(b.id); m.addBond(b.a, b.b, b.order); }
    assert.equal(m.project.atoms.find(a => a.id === center.id).chiralTag, 'CHI_UNSPECIFIED');
    assert.ok(chemistryWarnings(m.project).some(w => w.includes('请重新指定手性')));
    const result = callGeometryService(`
import json, sys
from rdkit import Chem
from molecular_car.editor import dispatch
p = json.load(sys.stdin)
out = dispatch('/api/editor/export', {'project': p, 'dimension': 2})
observed = Chem.MolFromMolBlock(out['content'], removeHs=False)
Chem.AssignStereochemistry(observed, cleanIt=True, force=True)
print(json.dumps({'CIP': [a.GetProp('_CIPCode') for a in observed.GetAtoms() if a.HasProp('_CIPCode')]}))
`, m.project);
    assert.deepEqual(result.CIP, []);
  }
});
