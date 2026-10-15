/**
 * Reproducible, offline feasibility probe of pinned upstream source.
 * Requires Node >= 24 (built-in TypeScript stripping); no npm installation.
 * This is a source/function probe, not a claim that the full browser UIs ran.
 */
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { stripTypeScriptTypes } from 'node:module';
import { runInNewContext } from 'node:vm';
import { fileURLToPath } from 'node:url';
import { Pool } from './upstream/ketcher/packages__ketcher-core__src__domain__entities__pool.ts';
import { ifDef } from './upstream/ketcher/packages__ketcher-core__src__utilities__ifDef.ts';

const root = new URL('./', import.meta.url);
const read = name => readFileSync(new URL(name, root), 'utf8');
const manifest = JSON.parse(read('upstream-manifest.json'));
const observations = [];
function check(name, fn) {
  const details = fn();
  observations.push({ name, passed: true, details });
}
const ks = path => read('upstream/ketcher/packages__ketcher-core__src__' + path.replaceAll('/', '__'));
const js = path => read('upstream/jsmol/src__org__jmol__' + path.replaceAll('/', '__'));

check('Pinned source integrity', () => {
  for (const entry of manifest.source_files) {
    const bytes = readFileSync(new URL(entry.path, root));
    assert.equal(createHash('sha256').update(bytes).digest('hex'), entry.sha256);
  }
  return { files: manifest.source_files.length, allSHA256Matched: true };
});

// Execute the original Pool and ifDef. For the serializer, execute the exact
// upstream functions after stripping TypeScript. Import wiring is replaced by
// an isolated context; Vec2 and coordinate helpers are also upstream code.
const vecSource = ks('domain/entities/vec2.ts').replace(/^import .*;\r?\n/gm, '');
const Vec2 = runInNewContext(stripTypeScriptTypes(vecSource).replaceAll('export ', '') + '\nVec2;', { assert });
const helpers = ks('domain/serializers/ket/helpers.ts');
const coordSource = helpers.slice(helpers.indexOf('const rotateCoordAxisBy180Degrees'), helpers.indexOf('export const modifyTransformation'));
const ketSource = ks('domain/serializers/ket/toKet/moleculeToKet.ts');
const serializerSource = ketSource.slice(ketSource.indexOf('export function moleculeToKet'), ketSource.indexOf('function rglabelToKet'));
const serialize = runInNewContext(stripTypeScriptTypes(coordSource + serializerSource).replaceAll('export ', '') + '\nmoleculeToKet;', {
  Vec2, ifDef, Axis: { y: 'y' },
});

let pool, serialized;
function atom(aam, x) {
  return { label: 'C', aam, pp: {x, y: 2, z: 7}, charge: 1, stereoParity: 1,
    stereoLabel: 'abs', queryProperties: {}, getInitiallySelected: () => false };
}
check('Ketcher actual Pool preserves IDs during local delete/add', () => {
  pool = new Pool();
  pool.add(atom(101, 1)); pool.add(atom(102, 2)); pool.add(atom(103, 3));
  pool.delete(1);
  const added = pool.add(atom(104, 4));
  assert.equal(added, 3);
  assert.deepEqual([...pool.keys()], [0, 2, 3]);
  return { livePoolIDs: [...pool.keys()] };
});
check('Ketcher actual KET serializer omits internal Pool IDs', () => {
  serialized = serialize({ atoms: pool, bonds: new Map(), sgroups: new Map(), frags: new Map() });
  assert.equal(serialized.atoms.length, 3);
  assert(serialized.atoms.every(a => !Object.hasOwn(a, 'id')));
  const reloaded = new Pool();
  serialized.atoms.forEach(a => reloaded.add(a));
  assert.deepEqual([...reloaded.keys()], [0, 1, 2]);
  return { beforePoolIDs: [...pool.keys()], newPoolIDs: [...reloaded.keys()], atomRecord: serialized.atoms[1] };
});
check('Ketcher mapping, formal charge, parity and z survive its atom serializer', () => {
  assert.deepEqual(Array.from(serialized.atoms, a => a.mapping), [101, 103, 104]);
  assert.equal(serialized.atoms[0].charge, 1);
  assert.equal(serialized.atoms[0].stereoParity, 1);
  assert.equal(serialized.atoms[0].location[2], 7);
  return { mappingPossible: true, qualification: 'mapping is the chemical atom-atom mapping field; it is not a permanent editor ID contract' };
});
check('Ketcher import explicitly reconstructs atoms and retains mapping', () => {
  const source = ks('domain/serializers/ket/fromKet/atomBondToStruct.ts');
  assert(source.includes("ifDef(params, 'aam', source.mapping)"));
  assert(source.includes('const newAtom = new Atom(params as AtomAttributes)'));
  return 'Source inspection: importer constructs a new Atom; adapter-owned mapping must be explicit.';
});

const modelkit = js('modelkit/ModelKit.java');
const writer = js('adapter/writers/MOLWriter.java');
check('JSmol ModelKit has real coordinate and topology editing implementations', () => {
  for (const text of ['private int assignAtom(', 'private boolean assignBond(', 'vwr.rotateAboutPointsInternal(', 'vwr.ms.setAtomCoord(']) assert(modelkit.includes(text));
  return 'Source inspection confirms atom/bond mutation and actual coordinate rotation/movement; not merely camera orbit.';
});
check('JSmol normal ModelKit editing has automatic hydrogen/geometry side effects', () => {
  for (const text of ['vwr.ms.removeUnnecessaryBonds(atom, isDelete)', 'vwr.ms.setAtomCoordRelative(atomIndex', 'vwr.addHydrogens(bsA, Viewer.MIN_SILENT)']) assert(modelkit.includes(text));
  return 'Must explicitly configure/adapter-wrap these operations before claiming topology/coordinate preservation.';
});
check('JSmol plain MOL export renumbers records and does not emit external atom IDs', () => {
  assert(writer.includes('atomMap[i] = ++n'));
  assert(writer.includes('PT.sprintf(" %12.5p %12.5p %12.5p 0", "p", o)'));
  return 'MOL atom numbering is export-local; V3000 atom mapping column is written as zero. Jmol SDF has extra atom-name data, but this is not the MOL synchronization path.';
});
check('JSmol stereo export is derived/conditional rather than verbatim metadata', () => {
  for (const text of ['private String getAtomParity(Atom a)', 'Arrays.sort(connections)', 'return (vTemp.dot(vNorm) > 0 ? "1" : "2")', 'if (is2d)']) assert(writer.includes(text));
  return 'V2000 parity comes from 3D coordinates and sorted atom indices; near/far wedge emission is inside is2d. This can be chemically useful but is not a verbatim stereo roundtrip contract.';
});
check('Source license notices are retained', () => {
  assert(read('upstream/ketcher/LICENSE').includes('Apache License'));
  assert(read('upstream/jsmol/LICENSE.txt').includes('GNU LESSER GENERAL PUBLIC LICENSE'));
  assert(modelkit.includes('version 2.1 of the License, or (at your option) any later version.'));
  return { ketcher: 'Apache-2.0', jsmol: 'LGPL-2.1-or-later' };
});

const output = { verifiedDate: '2026-10-11', nodeVersion: process.version,
  scope: 'Ketcher isolated upstream function execution + JSmol pinned-source inspection; not a full browser integration test',
  decision: 'Use a project-owned shared graph with Three.js + RDKit for this prototype; defer mature-editor adapter integration until ID/stereo/history contracts have dedicated tests.',
  observations,
};
writeFileSync(new URL('probe-results.json', root), JSON.stringify(output, null, 2) + '\n');
console.log(JSON.stringify(output, null, 2));
