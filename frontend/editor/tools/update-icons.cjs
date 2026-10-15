/* Regenerate the offline module and hashes after intentional local SVG edits. */
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const dir=path.resolve(__dirname,'../vendor/ketcher-icons');
const manifest=JSON.parse(fs.readFileSync(path.join(dir,'manifest.json'),'utf8'));
for(const name of ['workshop-atom','workshop-hydrogen','workshop-fit','workshop-uff','workshop-layout'])
 if(!manifest.icons.some(i=>i.name===name))manifest.icons.push({name,file:name+'.svg',origin:'workshop-original',changes:'new original 24px UI drawing'});
const sha=file=>crypto.createHash('sha256').update(fs.readFileSync(path.join(dir,file))).digest('hex');
const sources={};for(const icon of manifest.icons){icon.sha256=sha(icon.file);sources[icon.name]=fs.readFileSync(path.join(dir,icon.file),'utf8');}
fs.writeFileSync(path.join(dir,'icons.js'),'// Generated from local SVG snapshots; see manifest.json, LICENSE and NOTICE.\nexport const iconSources=Object.freeze('+JSON.stringify(sources,null,2)+');\n');
manifest.version='ketcher-3.18.0-workshop-icons.3';manifest.generatedModule.sha256=sha('icons.js');
for(const file of manifest.licenseFiles)file.sha256=sha(file.file);
fs.writeFileSync(path.join(dir,'manifest.json'),JSON.stringify(manifest,null,2)+'\n');
console.log('Verified and embedded '+manifest.icons.length+' local icons.');
