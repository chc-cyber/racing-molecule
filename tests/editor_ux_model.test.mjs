import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';

const source = readFileSync(new URL('../frontend/editor/model.js', import.meta.url), 'utf8');
const {WorkshopModel, createProject, clone, geometryPreview, fragmentPlacementPreview,
  movableSide, distance, angle, dihedral, selectedIds} = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
const library = JSON.parse(readFileSync(new URL('../frontend/editor/fragments.json', import.meta.url), 'utf8')).fragments;
const ethynyl = library.find(f => f.id === 'ethynyl');
const cage = library.find(f => f.id === 'c60');
const close = (a, b, tolerance = 1e-8) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
const pointClose = (a, b) => a.forEach((v, i) => close(v, b[i]));
const atom = (p, id) => p.atoms.find(a => a.id === id);
function graph(points, edges = [], elements = []) {
  return {...createProject('interaction fixture'), has3d: true,
    atoms: points.map((pos3, i) => ({id: `a${i}`, element: elements[i] ?? 'C', charge: 0,
      pos2: [i * 3 + 20, i % 2 + 30], pos3: [...pos3]})),
    bonds: edges.map(([a, b, order = 1], i) => ({id: `b${i}`, a: `a${a}`, b: `a${b}`, order}))};
}
const chain = () => graph([[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 2, 0]], [[0, 1], [1, 2], [2, 3], [3, 4]]);
function freezeDeep(value) { Object.freeze(value); for (const child of Object.values(value)) if (child && typeof child === 'object') freezeDeep(child); return value; }

test('a blank-to-blank drawing gesture creates the actual graph as one reversible transaction', () => {
  const m = new WorkshopModel(); let events = 0;
  m.addEventListener('change', () => ++events);
  const result = m.drawBondGesture({startPoint: [2, 3], endPoint: [3.3, 3.5], element: 'N', order: 2});
  assert.equal(m.project.atoms.length, 2); assert.equal(m.project.bonds.length, 1);
  assert.equal(m.project.bonds[0].order, 2); assert.equal(m.project.atoms[0].element, 'N');
  assert.deepEqual(m.project.bonds[0].a, result.startId); assert.equal(result.createdAtomIds.length, 2);
  assert.equal(m.revision, 1); assert.equal(events, 1);
  const finished = clone(m.project);
  m.undo(); assert.equal(m.project.atoms.length, 0); assert.equal(m.canUndo, false);
  m.redo(); assert.deepEqual(m.project, finished);
});

test('extending an existing atom preserves its identity/properties and uses local offsets in the other layout', () => {
  const p = graph([[5, 6, 7]]); Object.assign(p.atoms[0], {element: 'O', isotope: 18, charge: -1});
  const m = new WorkshopModel(p);
  const made = m.drawBondGesture({startId: 'a0', endPoint: [21.5, 31], element: 'C', dimension: 2});
  assert.deepEqual(m.project.atoms[0], p.atoms[0]);
  pointClose(atom(m.project, made.endId).pos3, [6.5, 7, 7]);
  assert.deepEqual(atom(m.project, made.endId).pos2, [21.5, 31]);
  assert.equal(m.revision, 1);
  m.undo(); assert.deepEqual(m.project, p);
  const reverse = m.drawBondGesture({startPoint: [4, 7, 9], endId: 'a0', dimension: 3});
  pointClose(atom(m.project, reverse.startId).pos2, [19, 31]);
  assert.deepEqual(m.project.atoms[0], p.atoms[0]);
});

test('drawing between existing atoms updates the same bond and honors wedge drag orientation', () => {
  const p = graph([[0, 0, 0], [1, 0, 0]], [[0, 1]]), m = new WorkshopModel(p);
  const result = m.drawBondGesture({startId: 'a1', endId: 'a0', order: 1, direction: 'BEGINDASH'});
  assert.equal(result.bondId, 'b0'); assert.deepEqual(result.createdAtomIds, []);
  assert.equal(m.project.bonds.length, 1); assert.equal(m.project.atoms.length, 2);
  assert.equal(m.project.bonds[0].a, 'a1'); assert.equal(m.project.bonds[0].b, 'a0');
  assert.equal(m.project.bonds[0].direction, 'BEGINDASH');
  assert.deepEqual(m.project.atoms, p.atoms);
  m.undo(); assert.deepEqual(m.project, p);
});

test('invalid drawing and failed new-atom transactions leave no orphan atoms or history', () => {
  const p = chain(), m = new WorkshopModel(p);
  for (const request of [
    {startId: 'a0', endId: 'a0'},
    {startPoint: [0, 0], endPoint: [NaN, 0]},
    {startId: 'a0', endPoint: [3, 4], order: 9},
    {startId: 'missing', endPoint: [3, 4]},
  ]) assert.throws(() => m.drawBondGesture(request));
  assert.deepEqual(m.project, p); assert.equal(m.revision, 0); assert.equal(m.canUndo, false);
});

test('replacing an element uses the same atom and preserves both layouts through undo/redo', () => {
  const p = chain(), m = new WorkshopModel(p);
  m.updateAtom('a2', {element: 'N'});
  assert.equal(m.project.atoms[2].id, 'a2');
  assert.equal(m.project.atoms[2].element, 'N');
  assert.deepEqual(m.project.atoms.map(a => a.pos2), p.atoms.map(a => a.pos2));
  assert.deepEqual(m.project.atoms.map(a => a.pos3), p.atoms.map(a => a.pos3));
  assert.deepEqual(m.project.bonds, p.bonds);
  m.undo(); assert.deepEqual(m.project, p); m.redo(); assert.equal(m.project.atoms[2].element, 'N');
});

test('geometry preview is pure, identifies whole movable components, and does not accumulate drag samples', () => {
  const p = freezeDeep(chain()), m = new WorkshopModel(p);
  for (const value of [1.2, 2.6, 1.8]) {
    const preview = geometryPreview(p, 'distance', ['a1', 'a2'], value);
    close(distance(preview.project, 'a1', 'a2'), value);
    assert.deepEqual(new Set(preview.movingIds), new Set(['a2', 'a3', 'a4']));
    assert.deepEqual(preview.fixedIds, ['a0', 'a1']);
    assert.deepEqual(preview.project.atoms.slice(0, 2), p.atoms.slice(0, 2));
    assert.deepEqual(preview.project.atoms.map(a => a.pos2), p.atoms.map(a => a.pos2));
    close(distance(preview.project, 'a2', 'a4'), distance(p, 'a2', 'a4'));
  }
  assert.equal(m.revision, 0); assert.equal(m.canUndo, false); assert.deepEqual(m.project, p);
  m.setDistance('a1', 'a2', 1.8);
  assert.equal(m.revision, 1); m.undo(); assert.deepEqual(m.project, p);
});

test('reverse geometry moves only the opposite side and reaches independently measured targets', () => {
  const p = chain();
  const cases = [
    ['distance', ['a1', 'a2'], 2.3, distance, ['a2', 'a3', 'a4']],
    ['angle', ['a0', 'a1', 'a2'], 116, angle, ['a1', 'a2', 'a3', 'a4']],
    ['dihedral', ['a0', 'a1', 'a2', 'a3'], -127, dihedral, ['a2', 'a3', 'a4']],
  ];
  for (const [type, ids, target, measure, fixed] of cases) {
    const preview = geometryPreview(p, type, ids, target, {reverse: true});
    close(measure(preview.project, ...ids), target);
    for (const atomId of fixed) pointClose(atom(preview.project, atomId).pos3, atom(p, atomId).pos3);
    assert.deepEqual(preview.project.bonds, p.bonds);
    assert.deepEqual(preview.project.atoms.map(a => a.pos2), p.atoms.map(a => a.pos2));
    const m = new WorkshopModel(p), method = `set${type[0].toUpperCase()}${type.slice(1)}`;
    m[method](...ids, target, {reverse: true});
    assert.deepEqual(m.project, preview.project); assert.equal(m.revision, 1);
    m.undo(); assert.deepEqual(m.project, p);
  }
});

test('ring-side preflight rejects ambiguity before preview or apply and retains the exact source', () => {
  const p = graph([[0, 0, 0], [1, 0, 0], [1, 1, 0]], [[0, 1], [1, 2], [2, 0]]), before = clone(p);
  assert.throws(() => movableSide(p, 'a0', 'a1'), /环/);
  assert.throws(() => geometryPreview(p, 'angle', ['a0', 'a1', 'a2'], 100), /环/);
  const m = new WorkshopModel(p); assert.throws(() => m.setDistance('a0', 'a1', 1.5, {reverse: true}), /环/);
  assert.deepEqual(p, before); assert.deepEqual(m.project, before); assert.equal(m.revision, 0);
});

test('fragment drag previews retain identical IDs; cancel has no history and apply matches one snapshot', () => {
  const p = freezeDeep(graph([[0, 0, 0]])), m = new WorkshopModel(p), template = freezeDeep(clone(ethynyl));
  const first = fragmentPlacementPreview(p, template, {pos3: [3, 4, 5]});
  const second = fragmentPlacementPreview(p, template, {pos3: [9, 4, 1]});
  assert.deepEqual(first.atomIds, second.atomIds); assert.equal(first.fragmentId, second.fragmentId);
  assert.deepEqual(first.idMap, second.idMap); assert.equal(m.revision, 0); assert.equal(m.canUndo, false);
  assert.deepEqual(m.project, p);
  const applied = m.placeFragment(template, {pos3: [9, 4, 1], placementId: second.placementId});
  assert.deepEqual(applied.project, second.project); assert.deepEqual(m.project, second.project);
  assert.equal(m.revision, 1); assert.equal(m.project.fragments.length, 1);
  m.undo(); assert.deepEqual(m.project, p); m.redo(); assert.deepEqual(m.project, second.project);
  const current = clone(m.project);
  assert.throws(() => m.placeFragment(template, {placementId: second.placementId}), /已经应用|已被使用/);
  assert.deepEqual(m.project, current); assert.equal(m.revision, 3);
  const subsequent = fragmentPlacementPreview(m.project, template);
  assert.equal(new Set([...subsequent.atomIds, ...second.atomIds]).size, 4);
});

test('2D fragment attachment has a real connecting bond, reasonable bond length and outward orientation', () => {
  const p = graph([[4, 8, 2], [3, 8, 2]], [[0, 1]]);
  p.atoms[0].pos2 = [4, 7]; p.atoms[1].pos2 = [3, 7];
  const m = new WorkshopModel(p);
  const result = m.placeFragment(ethynyl, {targetId: 'a0', portIndex: 0, dimension: 2, pos2: [12, 7]});
  const port = atom(m.project, result.ports[0].atomId), inside = atom(m.project, result.atomIds[1]);
  pointClose(port.pos2, [5.52, 7]); pointClose(inside.pos2, [6.72, 7]);
  assert.equal(m.project.bonds.length, p.bonds.length + ethynyl.bonds.length + 1);
  assert.equal(m.project.bonds.find(b => b.id === result.bondId).b, port.id);
  assert.deepEqual(m.project.atoms.slice(0, 2), p.atoms);
  assert.deepEqual(result.atomIds.map(id => atom(m.project, id).pos3), ethynyl.atoms.map(a => a.pos3));
  assert.equal(m.revision, 1); m.undo(); assert.deepEqual(m.project, p);
});

test('3D attachment rotates the entire new fragment while preserving existing coordinates and inactive layout', () => {
  const p = graph([[1, 2, 3], [0, 2, 3]], [[0, 1]]), before = clone(p);
  const result = fragmentPlacementPreview(p, ethynyl, {targetId: 'a0', dimension: 3, pos3: [1, 12, 3]});
  const port = atom(result.project, result.ports[0].atomId), inside = atom(result.project, result.atomIds[1]);
  pointClose(port.pos3, [1, 3.52, 3]); pointClose(inside.pos3, [1, 4.72, 3]);
  close(distance(result.project, port.id, inside.id), 1.2);
  assert.deepEqual(result.project.atoms.slice(0, 2), p.atoms);
  assert.deepEqual(result.atomIds.map(id => atom(result.project, id).pos2), ethynyl.atoms.map(a => a.pos2));
  assert.deepEqual(p, before);
});

test('connection failure checks both port restrictions and explicit valence without deleting hydrogen', () => {
  const saturated = graph([[0, 0, 0]]); saturated.atoms[0].explicitHs = 4;
  const m = new WorkshopModel(saturated);
  assert.throws(() => m.placeFragment(ethynyl, {targetId: 'a0'}), /足够价态/);
  assert.deepEqual(m.project, saturated); assert.equal(m.revision, 0); assert.equal(m.canUndo, false);
  const restricted = graph([[0, 0, 0], [1, 0, 0]], [[0, 1]]);
  restricted.fragments = [{id: 'f', name: 'restricted', atomIds: ['a0', 'a1'], ports: [{atomId: 'a1', maxOrder: 1}], rigid: true}];
  assert.throws(() => fragmentPlacementPreview(restricted, ethynyl, {targetId: 'a0'}), /连接位点/);
  assert.throws(() => fragmentPlacementPreview(restricted, ethynyl, {targetId: 'a1', order: 2}), /足够价态|最多允许/);
  assert.throws(() => fragmentPlacementPreview(restricted, ethynyl, {portIndex: 99}), /存在的片段连接位点/);
});

test('real C60 placement keeps 60 atoms, 90 bonds and rigid structure; saturated attachment is rejected transactionally', () => {
  const p = graph([[0, 0, 0]]), m = new WorkshopModel(p);
  assert.throws(() => m.placeFragment(cage, {targetId: 'a0', pos3: [5, 0, 0]}), /足够价态/);
  assert.deepEqual(m.project, p); assert.equal(m.revision, 0);
  const placed = m.placeFragment(cage, {pos3: [12, 6, 2]});
  assert.equal(placed.atomIds.length, 60); assert.equal(m.project.bonds.length, 90);
  assert.equal(selectedIds(m.project, [placed.atomIds[0]]).size, 60);
  for (const bond of cage.bonds) close(distance(m.project, placed.idMap[bond.a], placed.idMap[bond.b]), distance(cage, bond.a, bond.b));
  m.undo(); assert.deepEqual(m.project, p); assert.equal(m.canUndo, false);
  const prepared = clone(cage), portId = prepared.ports[0].atomId;
  prepared.bonds.find(b => b.order === 2 && [b.a, b.b].includes(portId)).order = 1;
  const attached = m.placeFragment(prepared, {targetId: 'a0', pos3: [5, 2, 7]});
  assert.equal(m.project.bonds.length, 91); close(distance(m.project, 'a0', attached.ports[0].atomId), 1.52);
  for (const bond of prepared.bonds) close(distance(m.project, attached.idMap[bond.a], attached.idMap[bond.b]), distance(prepared, bond.a, bond.b));
});

test('portless custom fragment can be placed at its centroid but cannot be silently attached', () => {
  const template = {...clone(ethynyl), ports: []}, p = graph([[0, 0, 0]]);
  const preview = fragmentPlacementPreview(p, template, {pos3: [8, 9, 10]});
  const points = preview.atomIds.map(id => atom(preview.project, id).pos3);
  pointClose(points[0].map((v, i) => (v + points[1][i]) / 2), [8, 9, 10]);
  assert.throws(() => fragmentPlacementPreview(p, template, {targetId: 'a0'}), /没有声明连接位点/);
});

test('preview/apply remains compatible with revision-gated async geometry and one-step cancellation', () => {
  const m = new WorkshopModel(chain()), base = clone(m.project), revision = m.revision;
  geometryPreview(base, 'dihedral', ['a0', 'a1', 'a2', 'a3'], 42);
  fragmentPlacementPreview(base, ethynyl, {pos3: [8, 0, 0]});
  assert.equal(m.revision, revision); assert.equal(m.canUndo, false);
  m.drawBondGesture({startId: 'a4', endPoint: [35, 33]});
  assert.equal(m.applyAsync(revision, base), false);
  m.undo(); assert.deepEqual(m.project, base);
  assert.equal(m.applyAsync(revision, base), false);
});
