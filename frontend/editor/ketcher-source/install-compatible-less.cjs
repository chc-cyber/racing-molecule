'use strict';
const fs = require('node:fs');
const path = require('node:path');
const {execFileSync} = require('node:child_process');
const root=path.resolve(__dirname,'build');
const archive=path.resolve(root,'tools/less-4.2.0.tgz');
const target=path.resolve(root,'upstream/node_modules/less');
if(!target.startsWith(root+path.sep)||!archive.startsWith(root+path.sep))throw new Error('Build dependency path escaped expected directory');
if(!fs.existsSync(archive))throw new Error('Download the official less-4.2.0.tgz archive first; see README.md');
const pkg=JSON.parse(fs.readFileSync(path.join(target,'package.json'),'utf8'));
if(pkg.version==='4.2.0'){console.log('Compatible Less 4.2.0 already installed');process.exit(0);}
const backup=path.resolve(root,`tools/less-upstream-${pkg.version}-${Date.now()}`);
if(!backup.startsWith(root+path.sep))throw new Error('Backup path escaped expected directory');
fs.renameSync(target,backup);
fs.mkdirSync(target);
execFileSync('tar',['-xzf',archive,'-C',target,'--strip-components=1'],{stdio:'inherit'});
if(JSON.parse(fs.readFileSync(path.join(target,'package.json'),'utf8')).version!=='4.2.0')throw new Error('Official compatible dependency failed to extract');
console.log('Compatible Less 4.2.0 installed inside the isolated build directory; original retained in '+backup);
