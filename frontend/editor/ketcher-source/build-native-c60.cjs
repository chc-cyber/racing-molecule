'use strict';
const fs = require('node:fs');
const path = require('node:path');
const {execFileSync} = require('node:child_process');
const crypto = require('node:crypto');
const base = __dirname;
const source = path.join(base, 'build/upstream');
const dist = path.join(base, '../vendor/ketcher-c60');
const original = path.join(base, '../vendor/ketcher');
const expectedCommit = 'f4f013e1855b71b0cbd37fc79d55b0cbe51c2371';
const commit = execFileSync('git', ['rev-parse', 'HEAD'], {cwd: source, encoding: 'utf8'}).trim();
if (commit !== expectedCommit) throw new Error('Only pinned official Ketcher v3.18.0 source is supported');
const run = (script, args, cwd) => execFileSync(process.execPath, [script, ...args], {cwd, stdio: 'inherit', env: {...process.env, MODE: 'standalone', NODE_ENV: 'production', PUBLIC_URL: './'}});
run(path.join(base, 'prepare-native-c60.cjs'), [], base);
run(path.join(source, 'packages/ketcher-core/scripts/compile-ket-schema.cjs'), [], path.join(source, 'packages/ketcher-core'));
run(path.join(source, 'node_modules/vite/bin/vite.js'), ['build', '--base=/editor/vendor/ketcher-c60/', '--outDir='+dist, '--emptyOutDir', '--config=vite-native.config.js'], path.join(source, 'example'));
fs.copyFileSync(path.join(source, 'LICENSE'), path.join(dist, 'LICENSE'));
fs.copyFileSync(path.join(source, 'NOTICE'), path.join(dist, 'NOTICE'));
fs.mkdirSync(path.join(dist, 'licenses'), {recursive: true});
for (const name of fs.readdirSync(path.join(original, 'static/js')).filter(n=>n.endsWith('.LICENSE.txt'))) fs.copyFileSync(path.join(original, 'static/js', name), path.join(dist, 'licenses/official-release-'+name));
fs.copyFileSync(path.join(base, 'README.md'), path.join(dist, 'SOURCE-EXTENSION.md'));
const hash = bytes => crypto.createHash('sha256').update(bytes).digest('hex');
const inventory = directory => fs.readdirSync(directory,{withFileTypes:true}).flatMap(entry=>{
  const full=path.join(directory,entry.name);
  return entry.isDirectory()?inventory(full):[{path:path.relative(dist,full).replaceAll('\\','/'),bytes:fs.statSync(full).size,sha256:hash(fs.readFileSync(full))}];
});
const data = JSON.parse(fs.readFileSync(path.join(base, 'extension-data.json'), 'utf8'));
const dependencies = JSON.parse(fs.readFileSync(path.join(dist, 'BUILD-DEPENDENCIES.json'), 'utf8'));
const patch = execFileSync('git',['diff','--','packages/ketcher-core/src/domain/entities/atom.ts','packages/ketcher-react/src/script/ui/data/templates.js','packages/ketcher-react/src/components/Icon/utils/iconNameToIcon.ts','packages/ketcher-react/src/script/ui/action/action.types.ts','packages/ketcher-react/src/script/editor/tool/template.ts'],{cwd:source,encoding:'utf8'});
fs.writeFileSync(path.join(base,'SOURCE-MODIFICATIONS.patch'),patch);
const manifest={
  component:'Ketcher with Molecular Car native C60 template',version:data.version,upstreamVersion:'3.18.0',
  source:{repository:'https://github.com/epam/ketcher',release:'https://github.com/epam/ketcher/releases/tag/v3.18.0',commit},
  entry:'/editor/vendor/ketcher-c60/index.html',license:'Apache-2.0; included third-party copyright notices retained',
  deployment:{runtimeNpmRequired:false,externalCDNRequired:false,oldFragmentLibraryRequired:false,originalStandaloneDistributionModified:false,nativeReactTemplateButton:true,domInjectedButton:false,minifiedBundlePatched:false},
  extension:{...data,sourceDirectory:'frontend/editor/ketcher-source',buildScript:'frontend/editor/ketcher-source/build-native-c60.cjs',readme:'frontend/editor/ketcher-source/README.md',sourcePatch:'frontend/editor/ketcher-source/SOURCE-MODIFICATIONS.patch'},
  build:{node:process.version,vite:'8.0.16',lessOverride:'4.2.0',upstreamPackageLockSha256:hash(fs.readFileSync(path.join(source,'package-lock.json'))),compatibleLessArchive:'https://registry.npmjs.org/less/-/less-4.2.0.tgz',compatibleLessSha256:hash(fs.readFileSync(path.join(base,'build/tools/less-4.2.0.tgz'))),cssMinified:false,reason:'Keep upstream legacy CSS without requiring lightningcss to parse its legacy container-query syntax'},
  dependencyInventory:'BUILD-DEPENDENCIES.json',includedDependencyCount:dependencies.dependencies.length,
  verification:{probe:'frontend/editor/ketcher-source/probe-native-c60-v2.cjs',report:'outputs/ketcher-native-c60-v2/report.json',rawNativeProbe:'frontend/editor/ketcher-source/probe-native-c60.cjs',minimalSUPProbe:'frontend/editor/ketcher-source/probe-superatom.cjs'},
  files:inventory(dist).sort((a,b)=>a.path.localeCompare(b.path)),
};
fs.writeFileSync(path.join(base,'../vendor/ketcher-c60-manifest.json'),JSON.stringify(manifest,null,2)+'\n');
console.log(`Native C60 distribution ready: ${manifest.files.length} files; ${manifest.includedDependencyCount} included dependency records.`);
