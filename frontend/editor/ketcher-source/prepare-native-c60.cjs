/** Source-level extension of EPAM Ketcher 3.18.0 (Apache-2.0).
 * Never modifies a compiled JavaScript bundle. See README.md for provenance.
 */
'use strict';
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');
const {execFileSync} = require('node:child_process');
const base = __dirname;
const sourceRoot = path.join(base, 'build/upstream');
const inputPath = path.join(base, 'c60-source.json');
if (!fs.existsSync(inputPath)) {
  const library = JSON.parse(fs.readFileSync(path.join(base, '../fragments.json'), 'utf8'));
  const c60 = library.fragments.find(f => f.id === 'c60');
  if (!c60) throw new Error('C60 source graph missing');
  fs.writeFileSync(inputPath, JSON.stringify({
    source: 'molecular_car/editor/fragments.py: procedural truncated icosahedron; copied at build time',
    atoms: c60.atoms, bonds: c60.bonds,
  }, null, 2) + '\n');
}
const data = JSON.parse(fs.readFileSync(inputPath, 'utf8'));
const ids = new Map(data.atoms.map((a, i) => [a.id, i]));
const bonds = data.bonds.map(b => ({a: ids.get(b.a), b: ids.get(b.b), order: b.order}));
const adjacency = data.atoms.map(() => []);
const valence = data.atoms.map(() => 0);
const bondKeys = new Set();
for (const b of bonds) {
  if (![1, 2].includes(b.order) || b.a === undefined || b.b === undefined || b.a === b.b) throw new Error('Invalid C60 bond');
  const key = [b.a, b.b].sort((a, b) => a-b).join(':');
  if (bondKeys.has(key)) throw new Error('Duplicate C60 bond');
  bondKeys.add(key);
  adjacency[b.a].push(b.b); adjacency[b.b].push(b.a);
  valence[b.a] += b.order; valence[b.b] += b.order;
}
if (data.atoms.length !== 60 || bonds.length !== 90 || adjacency.some(a => a.length !== 3) || valence.some(v => v !== 4) || bonds.filter(b => b.order === 2).length !== 30) throw new Error('C60 must be 60 carbons, 90 bonds, degree 3, valence 4, thirty double bonds');
if (data.atoms.some(a => a.element !== 'C')) throw new Error('C60 must contain only carbon');

// Schlegel projection through a pentagonal face produces a planar chemical
// drawing of the fullerene graph; a sphere projection would hide crossed bonds.
const cycles = [];
function walk(nodes) {
  const end = nodes.at(-1);
  if (nodes.length === 5) {
    if (adjacency[end].includes(nodes[0])) cycles.push(nodes);
    return;
  }
  for (const next of adjacency[end]) if (!nodes.includes(next)) walk([...nodes, next]);
}
walk([0]);
if (!cycles.length) throw new Error('Pentagon face not found');
const outerFace = cycles[0];
const add = (a, b) => a.map((v, i) => v + b[i]);
const sub = (a, b) => a.map((v, i) => v - b[i]);
const scale = (a, k) => a.map(v => v*k);
const dot = (a, b) => a.reduce((v, x, i) => v+x*b[i], 0);
const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
const norm = a => Math.sqrt(dot(a, a));
const unit = a => scale(a, 1/norm(a));
const center = scale(outerFace.reduce((sum, i) => add(sum, data.atoms[i].pos3), [0, 0, 0]), 1/5);
const normal = unit(center);
const height = dot(normal, center);
const eye = scale(normal, height + 0.28);
const u = unit(sub(data.atoms[outerFace[0]].pos3, center));
const v = unit(cross(normal, u));
let points = data.atoms.map(atom => {
  const ray = sub(atom.pos3, eye);
  const t = (height - dot(normal, eye))/dot(normal, ray);
  const p = sub(add(eye, scale(ray, t)), center);
  return [dot(p, u), dot(p, v)];
});
const averageBond = bonds.reduce((sum, b) => sum+norm(sub(points[b.a], points[b.b])), 0)/90;
points = points.map(p => scale(p, 1.5/averageBond));
const atomLine = p => `${p[0].toFixed(4).padStart(10)}${p[1].toFixed(4).padStart(10)}    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0`;
const n = x => String(x).padStart(3);
const groupLines = ['M  STY  1   1 SUP'];
for(let start=0;start<60;start+=15)groupLines.push(`M  SAL   1${n(Math.min(15,60-start))}`+Array.from({length:Math.min(15,60-start)},(_,i)=>' '+n(start+i+1)).join(''));
groupLines.push('M  SMT   1 C₆₀', 'M  SAP   1  1   1   0  1');
const mol = ['C60 Fullerene', '  Ketcher  10132600002D', '', ' 60 90  0  0  0  0            999 V2000', ...points.map(atomLine), ...bonds.map(b => `${n(b.a+1)}${n(b.b+1)}${n(b.order)}  0  0  0  0`), ...groupLines, 'M  END'].join('\n');
fs.writeFileSync(path.join(base, 'c60-template.mol'), mol+'\n');
// The same explicit template index orders all native SUP members and 3D seeds.
// No graph matching or SMILES round-trip participates in editor identity.
fs.writeFileSync(path.join(base,'../native-c60-data.js'), '/** Molecular Car C60 v2; procedural truncated icosahedron, original code. */\nexport const C60_TEMPLATE_VERSION=2;\nexport const C60_POS3='+JSON.stringify(data.atoms.map(a=>a.pos3))+';\nexport const C60_BONDS='+JSON.stringify(bonds)+';\n');

const patchedSources = new Map();
function patch(relative, oldValue, newValue) {
  const target = path.join(sourceRoot, relative);
  const source = patchedSources.get(relative) ?? execFileSync('git',['show','HEAD:'+relative],{cwd:sourceRoot,encoding:'utf8'});
  if (source.split(oldValue).length !== 2) throw new Error(`Pinned source precondition failed: ${relative}`);
  const modified=source.replace(oldValue,newValue);patchedSources.set(relative,modified);fs.writeFileSync(target,modified);
}
// Nitrogen cannot use phosphorus' neutral pentavalent fallback. Returning a
// negative H count activates Ketcher's native red valence-error decoration.
patch('packages/ketcher-core/src/domain/entities/atom.ts',
  "    if (label === 'N' || label === 'P') {",
  `    if (label === 'N' && charge === 0) {
      return { valence: 3, hydrogenCount: 3 - radicalCount - connectionCount };
    }
    if (label === 'N' || label === 'P') {`);
patch('packages/ketcher-core/src/domain/entities/atom.ts',
  '    const hydrogenCount = this.overrideHydrogenCountIfNeeded(',
  `    if (label === 'N' && charge === 0 &&
        connectionCount + radicalCount + (this.implicitHCount ?? 0) > 3) {
      return this.applyValenceResult(3, -1, connectionCount);
    }
    const hydrogenCount = this.overrideHydrogenCountIfNeeded(`);
const react = 'packages/ketcher-react/src';
patch(`${react}/script/ui/data/templates.js`, '\n].map((structStr)', `\n  // Molecular Car extension: real C60 graph; see ketcher-source/README.md.\n  ${JSON.stringify(mol)},\n].map((structStr)`);
patch(`${react}/script/ui/data/templates.js`, '].map((structStr) => molSerializer.deserialize(structStr));', `].map((structStr) => {
  const struct = molSerializer.deserialize(structStr);
  if (struct.name === 'C60 Fullerene') {
    const group = struct.sgroups.values().next().value;
    // This marker originates only in this native template source. fromPaste
    // explicitly maps ordered SUP members to newly allocated live Pool IDs.
    group.data.molecularCarTemplate = {id: 'c60', version: 2};
    group.data.expanded = false;
  }
  return struct;
});`);
patch(`${react}/components/Icon/utils/iconNameToIcon.ts`, "import Template7Icon from '../../../assets/icons/files/template-7.svg';", "import Template7Icon from '../../../assets/icons/files/template-7.svg';\nimport Template8Icon from '../../../assets/icons/files/template-8.svg';");
patch(`${react}/components/Icon/utils/iconNameToIcon.ts`, "  'template-7': Template7Icon,", "  'template-7': Template7Icon,\n  'template-8': Template8Icon,");
patch(`${react}/script/ui/action/action.types.ts`, "  | 'template-7'", "  | 'template-7'\n  | 'template-8'");
// Saturated fullerene must never merge with an existing H/atom/bond. Native
// drawing still owns placement and SUP selection; connections use its bond
// tool and the canonical, undoable addition/valence transaction.
patch(`${react}/script/editor/tool/template.ts`, '    this.event = event;\n\n    this.templatePreview?.hidePreview();', `    this.event = event;

    this.templatePreview?.hidePreview();
    if (this.template.molecule.name === 'C60 Fullerene' && this.editor.findItem(event, this.findItems)) {
      this.editor.event.message.dispatch({info: '请先将 C₆₀ 放在空白处，再用单键连接。'});
      return;
    }`);
const icon = '<svg width="32" height="32" viewBox="0 0 32 32" xmlns="http://www.w3.org/2000/svg"><text x="3" y="22" fill="currentColor" font-family="Arial,sans-serif" font-size="21" font-weight="600">C<tspan dy="4" font-size="13">60</tspan></text></svg>\n';
fs.writeFileSync(path.join(sourceRoot, `${react}/assets/icons/files/template-8.svg`), icon);
fs.writeFileSync(path.join(base, 'template-8.svg'), icon);

const sha256 = value => crypto.createHash('sha256').update(value).digest('hex');
fs.writeFileSync(path.join(base, 'extension-data.json'), JSON.stringify({
  version: '3.18.0+molecular-car.c60.4', upstreamCommit: 'f4f013e1855b71b0cbd37fc79d55b0cbe51c2371',
  nativeButton: 'template-8', title: 'C60 Fullerene', atoms: 60, bonds: 90, doubleBonds: 30,
  sourceGraphSha256: sha256(fs.readFileSync(inputPath)), molSha256: sha256(mol+'\n'),
  outerFace, drawing: 'Native collapsed SUP C60; Schlegel drawing on expand; explicit ordered spherical 3D seed',
  templateProvenance:'SUP data.molecularCarTemplate.id=c60,version=2; ordered native insertion Pool map',
  modifiedUpstreamSources: ['packages/ketcher-core/src/domain/entities/atom.ts', `${react}/script/ui/data/templates.js`, `${react}/components/Icon/utils/iconNameToIcon.ts`, `${react}/script/ui/action/action.types.ts`, `${react}/assets/icons/files/template-8.svg`, `${react}/script/editor/tool/template.ts`],
}, null, 2)+'\n');

// The official Vite source resolver supports building all four packages directly.
// Its dev config imports the entire Webpack/Rollup build configurations merely
// for constant values. Supply those same fixed constants to avoid bundling the
// legacy build tools into Vite's configuration loader (lightningcss exports).
let vite = fs.readFileSync(path.join(sourceRoot, 'example/vite.config.js'), 'utf8');
const values = {
  'process.env.NODE_ENV': JSON.stringify('production'),
  'process.env.VERSION': JSON.stringify('3.18.0'),
  'process.env.BUILD_DATE': JSON.stringify('2026-10-12T00:00:00'),
  'process.env.BUILD_NUMBER': 'undefined',
  'process.env.HELP_LINK': JSON.stringify('v3.18.0'),
  'process.env.INDIGO_VERSION': JSON.stringify('1.46.0'),
  'process.env.INDIGO_MACHINE': JSON.stringify(''),
};
vite = vite.replace("import { valuesToReplace as polymerEditorValues } from '../packages/ketcher-macromolecules/rollup.config.mjs';", `const polymerEditorValues = ${JSON.stringify(values)};`)
  .replace("import { valuesToReplace as ketcherReactValues } from '../packages/ketcher-react/rollup.config.mjs';", `const ketcherReactValues = ${JSON.stringify(values)};`)
  .replace("import { envVariables as exampleEnv } from './config/webpack.config';", "const exampleEnv = { MODE: 'standalone', KETCHER_ENABLE_REDUX_LOGGER: 'false', NODE_ENV: 'production', PUBLIC_URL: './' };")
  .replace("import { INDIGO_WORKER_IMPORTS } from '../packages/ketcher-standalone/rollup.config.mjs';", "const INDIGO_WORKER_IMPORTS = { WASM_LOADER: './indigoWorkerImports/useWasmLoader' };")
  .replace('export default defineConfig({', 'export default defineConfig({\n  root: resolve(__dirname),\n  build: { cssMinify: false },');
const noticePlugin = fs.readFileSync(path.join(base, 'dependency-notices-plugin.cjs'), 'utf8');
vite = "import fs from 'node:fs';\nimport path from 'node:path';\n"+vite;
vite = vite.replace('  plugins: [\n    PackageScopedAliasesPlugin(),', `  plugins: [\n    ${noticePlugin.trim()},\n    PackageScopedAliasesPlugin(),`);
fs.writeFileSync(path.join(sourceRoot, 'example/vite-native.config.js'), vite);
console.log('Native template-8 prepared from source: 60 C / 90 bonds / 30 double bonds.');
