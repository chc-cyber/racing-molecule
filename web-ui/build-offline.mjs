import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import vm from 'node:vm';

// This deliberately small bundler supports the static imports/exports in the
// bundled local Three.js distribution. It does not install or download tools.
// Every transformation is asserted so future source changes fail visibly.
const project = path.dirname(fileURLToPath(import.meta.url));
const dist = path.join(project, 'dist');
const read = name => fs.readFile(path.join(dist, name), 'utf8');
const [html, baseCss, minimalCss, appSource, sceneSource, coreSource, threeSource, moleculeSource, fallbackSource, iconSource, threeLicense, fontLicense, font] = await Promise.all([
  read('index.html'), read('style.css'), read('minimal.css'), read('app.js'), read('scene.js'),
  read('assets/three.core.js'), read('assets/three.module.js'),
  read('assets/molecule.json'), read('assets/molecule-fallback.svg'),
  read('assets/mark.svg'), read('assets/THREE-LICENSE.txt'),
  read('assets/BARLOW-OFL.txt'), fs.readFile(path.join(dist, 'assets/BarlowCondensed-Bold.ttf')),
]);
const cssSource = `${baseCss}\n${minimalCss}`;

function once(source, pattern, replacement, label) {
  if (!pattern.test(source)) throw new Error(`Missing expected ${label}; inspect the source before rebuilding.`);
  return source.replace(pattern, replacement);
}

function parseNames(list) {
  return list.split(',').map(name => name.trim()).filter(Boolean).map(name => {
    const match = /^([A-Za-z_$][\w$]*)(?:\s+as\s+([A-Za-z_$][\w$]*))?$/.exec(name);
    if (!match) throw new Error(`Unsupported module identifier: ${name}`);
    return { local: match[1], exported: match[2] || match[1] };
  });
}

function exportsObject(names) {
  return names.map(({ local, exported }) => `${exported}:${local}`).join(',');
}

function assertClassic(source, label) {
  if (/^\s*(?:import|export)\s/m.test(source) || /\bimport\.meta\b/.test(source)) {
    throw new Error(`${label} contains an unsupported remaining module statement.`);
  }
  new vm.Script(source, { filename: label });
}

let coreExports;
const core = once(coreSource, /^export\s*\{([^}]+)\};?\s*$/m, (_, names) => {
  coreExports = parseNames(names);
  return '';
}, 'Three.js core export');
assertClassic(core, 'three.core.js');
const coreBundle = `const __threeCore = (() => {\n${core}\nreturn {${exportsObject(coreExports)}};\n})();`;

let three = threeSource;
three = once(three, /^import\s*\{([^}]+)\}\s*from\s*['"]\.\/three\.core\.js['"];?\s*$/m, (_, names) => {
  const destructuring = parseNames(names).map(({ local, exported }) => exported === local ? local : `${local}:${exported}`).join(',');
  return `const {${destructuring}} = __threeCore;`;
}, 'Three.js renderer import');
let reexports;
three = once(three, /^export\s*\{([^}]+)\}\s*from\s*['"]\.\/three\.core\.js['"];?\s*$/m, (_, names) => {
  reexports = parseNames(names);
  return '';
}, 'Three.js renderer re-export');
let threeExports;
three = once(three, /^export\s*\{([^}]+)\};?\s*$/m, (_, names) => {
  threeExports = parseNames(names);
  return '';
}, 'Three.js renderer export');
assertClassic(three, 'three.module.js');
const externalProperties = reexports.map(({ local, exported }) => `${exported}:__threeCore.${local}`).join(',');
const threeBundle = `const __three = (() => {\n${three}\nreturn {${externalProperties},${exportsObject(threeExports)}};\n})();`;

const json = JSON.parse(moleculeSource);
if (json.atoms.length !== 52 || json.bonds.length !== 56) throw new Error('Unexpected visual molecule structure.');
const svgData = `data:image/svg+xml;base64,${Buffer.from(fallbackSource).toString('base64')}`;
let scene = once(sceneSource, /^import\s+\*\s+as\s+THREE\s+from\s+['"]\.\/assets\/three\.module\.js['"];?\s*$/m, 'const THREE = __three;', 'scene import');
scene = once(scene, /^const MODEL_URL = new URL\([^\n]+\);\s*$/m, `const MODEL_DATA = ${JSON.stringify(json)};`, 'model URL');
scene = once(scene, /^const FALLBACK_URL = new URL\([^\n]+\);\s*$/m, `const FALLBACK_URL = { href: ${JSON.stringify(svgData)} };`, 'fallback URL');
scene = once(scene, /await fetch\(MODEL_URL(?:\.href)?,\s*\{ signal: abort\.signal \}\)/, 'await Promise.resolve({ ok: true, json: async () => MODEL_DATA })', 'model fetch');
scene = once(scene, /^export function mountMoleculeScene/m, 'function mountMoleculeScene', 'scene export');
assertClassic(scene, 'scene.js');
const sceneBundle = `const __scene = (() => {\n${scene}\nreturn {mountMoleculeScene};\n})();`;

// Execute the library wrappers as well as checking their syntax. DOM-dependent
// rendering is verified separately in the browser; these classes need no DOM.
const smokeContext = { console, AbortController };
vm.runInNewContext(`(() => {${coreBundle}\n${threeBundle}\n${sceneBundle}\nglobalThis.smoke = { renderer: typeof __three.WebGLRenderer, scene: typeof __scene.mountMoleculeScene, vertices: new __three.SphereGeometry(1, 12, 8).attributes.position.count };})();`, smokeContext);
if (smokeContext.smoke.renderer !== 'function' || smokeContext.smoke.scene !== 'function' || smokeContext.smoke.vertices !== 117) {
  throw new Error('Offline module export wiring failed.');
}

const app = once(appSource, /^import\s*\{\s*mountMoleculeScene\s*\}\s*from\s*['"]\.\/scene\.js['"];?\s*$/m, 'const { mountMoleculeScene } = __scene;', 'application import');
assertClassic(app, 'app.js');
const bundle = `(() => {\n'use strict';\n${coreBundle}\n${threeBundle}\n${sceneBundle}\n${app}\n})();`;
assertClassic(bundle, 'offline-bundle.js');

const fontData = `data:font/ttf;base64,${font.toString('base64')}`;
const css = once(cssSource, /url\(['"]?\.\/assets\/BarlowCondensed-Bold\.ttf['"]?\)/, `url('${fontData}')`, 'local font');
if (/url\(\s*['"]?(?!data:|#)[^)'"\s]+/.test(css)) throw new Error('The offline stylesheet still contains an external asset.');
const iconData = `data:image/svg+xml;base64,${Buffer.from(iconSource).toString('base64')}`;
let offline = once(html, /<link rel="icon"[^>]+>/, `<link rel="icon" type="image/svg+xml" href="${iconData}">`, 'favicon link');
offline = once(offline, /\s*<link rel="preload"[^>]+>/, '', 'font preload');
offline = once(offline, /<link rel="stylesheet"[^>]+>/, `<style>\n${css}\n</style>`, 'stylesheet link');
offline = once(offline, /\s*<link rel="stylesheet" href="\.\/minimal\.css">/, '', 'minimal stylesheet link');
// Literal replacement callbacks avoid treating JavaScript `$&` and `$'` text
// as String.replace substitutions. Escaping keeps inline HTML script parsing safe.
const inlineScript = bundle.replace(/<\/script/gi, '<\\/script');
offline = once(offline, /<script type="module" src="\.\/app\.js"><\/script>/, () => `<script>\n${inlineScript}\n</script>`, 'entry module');
const licenses = JSON.stringify({
  format: 'Racing Molecule standalone interface prototype',
  threejs: { license: 'MIT', text: threeLicense },
  barlow: { license: 'SIL Open Font License 1.1', text: fontLicense },
}).replace(/</g, '\\u003c');
offline = once(offline, /<\/body>/, () => `<script type="application/json" id="bundled-open-source-licenses">${licenses}</script>\n</body>`, 'body terminator');
if (/<(?:script|link)\b[^>]+(?:src|href)=["']\.\//i.test(offline)) throw new Error('The offline HTML still references an external source file.');

const output = path.join(dist, 'offline.html');
await fs.writeFile(output, offline, 'utf8');
const bytes = Buffer.byteLength(offline);
console.log(JSON.stringify({ output, bytes, mebibytes: Number((bytes / 1024 / 1024).toFixed(3)), atoms: json.atoms.length, bonds: json.bonds.length, bundledLibraries: ['Three.js'], embeddedLicenses: ['MIT', 'SIL OFL 1.1'], networkRequired: false }, null, 2));
