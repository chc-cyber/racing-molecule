/* Own headless browser, static interception only. Never contacts race APIs. */
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),{spawnSync}=require('node:child_process');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root=path.resolve(__dirname,'../../..'),output=path.join(root,'outputs/editor-jsmol-check');
async function launch(){
 const browser=await chromium.launch({executablePath:process.env.CHROME_PATH||'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true});
 const page=await browser.newPage({viewport:{width:1100,height:760}}),errors=[],missing=[],external=[];
 page.on('pageerror',e=>errors.push(e.stack));
 await page.route('**/*',route=>{const url=new URL(route.request().url());if(url.hostname!=='127.0.0.1'||url.pathname.startsWith('/api/')){external.push(url.href);return route.abort();}const target=path.resolve(root,'frontend','.'+decodeURIComponent(url.pathname));assert(target.startsWith(path.join(root,'frontend')+path.sep));if(!fs.existsSync(target)||!fs.statSync(target).isFile()){missing.push(url.pathname);return route.fulfill({status:404,body:'not found'});}const ext=path.extname(target),mime={'.js':'text/javascript','.html':'text/html','.css':'text/css','.json':'application/json','.png':'image/png','.gif':'image/gif','.jpg':'image/jpeg','.mol':'chemical/x-mdl-molfile'}[ext]||'application/octet-stream';return route.fulfill({contentType:mime,body:fs.readFileSync(target)});});
 await page.goto('http://127.0.0.1:8765/editor/research/jsmol-live.html');await page.waitForFunction(()=>window.probeReady,{},{timeout:60000});
 return {browser,page,errors,missing,external};
}
module.exports={launch,root,output};
function python(code,input){const r=spawnSync(process.env.EDITOR_PYTHON||path.join(root,'.venv/Scripts/python.exe'),['-X','utf8','-E','-B','-c',code],{cwd:root,input:input===undefined?undefined:JSON.stringify(input),encoding:'utf8',maxBuffer:4e6});if(r.status!==0)throw Error(r.stderr||r.stdout);return JSON.parse(r.stdout);}
function fixtureData(){return python(`import json,pathlib
from molecular_car.editor.service import import_file,export_file,_stereo_signature
from molecular_car.editor.fragments import fullerene
from rdkit import Chem
out={}
for key,name in [('tpb','tetraphenylbenzene.mol'),('chiral','chiral.mol')]:
 text=(pathlib.Path('frontend/editor/examples')/name).read_text();p=import_file({'filename':name,'content':text})['project'];out[key]={'project':p,'mol':text}
p=fullerene();out['c60']={'project':p,'mol':export_file({'project':p,'dimension':3})['content']}
print(json.dumps(out))`);}
const topology=s=>s.bonds.map(b=>JSON.stringify([Math.min(b.atom1.atomIndex,b.atom2.atomIndex),Math.max(b.atom1.atomIndex,b.atom2.atomIndex),b.order])).sort();
const displacement=(a,b)=>Math.hypot(...a.map((v,k)=>v-b[k]));
function torsion(points){const sub=(a,b)=>a.map((v,k)=>v-b[k]),dot=(a,b)=>a.reduce((n,v,k)=>n+v*b[k],0),scale=(a,n)=>a.map(v=>v*n),norm=a=>Math.hypot(...a),unit=a=>scale(a,1/norm(a)),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];const axis=unit(sub(points[2],points[1])),ab=sub(points[0],points[1]),cd=sub(points[3],points[2]),v=unit(sub(ab,scale(axis,dot(ab,axis)))),w=unit(sub(cd,scale(axis,dot(cd,axis))));return Math.atan2(dot(cross(axis,v),w),dot(v,w))*180/Math.PI;}
async function verify(){
 const fixtures=fixtureData(),{browser,page,errors,missing,external}=await launch(),checks=[],limitations=[];fs.mkdirSync(output,{recursive:true});
 const check=(name,details={})=>{checks.push({name,passed:true,details});console.log('PASS',name);};
 const run=script=>page.evaluate(s=>jsmolProbe.run(s),script),snapshot=()=>page.evaluate(()=>jsmolProbe.inspect());
 const load=async key=>{await page.evaluate(async f=>{const {JSmolAdapter}=await import('/editor/jsmol-adapter.js');const {WorkshopModel}=await import('/editor/model.js');window.hostModel=new WorkshopModel(f.project);window.bridge=new JSmolAdapter(probeApplet);bridge.loadProject(f.project,f.mol,{indexToAtomId:f.project.atoms.map(a=>a.id),revision:hostModel.revision});jsmolProbe.run('spacefill 22%; wireframe 0.12; zoom 85;');},fixtures[key]);return snapshot();};
 let failure=null;
 try{
  let s=await load('tpb');assert.equal(s.version,1603031);assert.equal(s.atoms.length,52);assert.equal(s.bonds.length,56);assert.equal(s.modelkit.addHydrogens,false);assert.equal(s.modelkit.autobond,false);check('Offline runtime 16.3.31 loads TPB and disables automatic H/bond additions',{atoms:52,bonds:56,settings:s.modelkit});
  const loaded=await page.evaluate(()=>bridge.read());assert.deepEqual(loaded.atoms.map(a=>a.id),fixtures.tpb.project.atoms.map(a=>a.id));check('Explicit MOL row mapping assigns project IDs through numeric property_editorid',{count:loaded.atoms.length,maxToken:Math.max(...loaded.atoms.map(a=>a.token))});

  await run('set antialiasDisplay false; set picking dragAtom;');s=await snapshot();
  // Internal screen coordinates are used ONLY to aim an actual test mouse drag,
  // never as a data synchronization/identity interface.
  const target=await page.evaluate(()=>{const v=probeApplet._applet.viewer,r=document.querySelector('canvas').getBoundingClientRect();const a=v.ms.at.find(a=>a.getElementNumber()===1&&v.findNearestAtomIndex(a.sX,a.sY)===a.i);if(!a)throw Error('No visible pickable H');return {index:a.i,x:r.x+a.sX,y:r.y+a.sY,orientation:Jmol.getPropertyAsArray(probeApplet,'orientationInfo')};});
  await page.mouse.move(target.x,target.y);await page.mouse.down();await page.mouse.move(target.x+25,target.y-20,{steps:7});await page.mouse.up();
  let after=await snapshot();assert(displacement(s.atoms[target.index].coord,after.atoms[target.index].coord)>.05);assert.deepEqual(topology(after),topology(s));assert.deepEqual(await page.evaluate(()=>Jmol.getPropertyAsArray(probeApplet,'orientationInfo')),target.orientation);
  const moved=await page.evaluate(()=>bridge.captureGeometry());assert.deepEqual(moved.project.atoms.map(a=>a.pos2),fixtures.tpb.project.atoms.map(a=>a.pos2));assert.deepEqual(moved.project.bonds,fixtures.tpb.project.bonds);assert.deepEqual(moved.project.fragments,fixtures.tpb.project.fragments);assert.equal(moved.baseRevision,0);
  check('Real ModelKit mouse drag changes atom coordinates without camera rotation/topology changes',{atomIndex:target.index,angstroms:displacement(s.atoms[target.index].coord,after.atoms[target.index].coord)});

  s=await load('tpb');let bridgeBond;
  for(const b of s.bonds){if(b.atom1.sym!=='C'||b.atom2.sym!=='C'||b.order!==1)continue;const a=b.atom1.atomIndex,c=b.atom2.atomIndex,branch=new Set([c]),todo=[c];for(let i=0;i<todo.length;i++)for(const d of s.bonds){const x=d.atom1.atomIndex,y=d.atom2.atomIndex;if(x===a&&y===c||x===c&&y===a)continue;const other=x===todo[i]?y:y===todo[i]?x:null;if(other!==null&&!branch.has(other)){branch.add(other);todo.push(other);}}if(!branch.has(a)){const neighbor=(at,except)=>s.bonds.flatMap(x=>[x.atom1.atomIndex===at?x.atom2.atomIndex:null,x.atom2.atomIndex===at?x.atom1.atomIndex:null]).find(i=>i!==null&&i!==except);bridgeBond={a,c,branch:[...branch],torsion:[neighbor(a,c),a,c,neighbor(c,a)]};break;}}
  assert(bridgeBond);const beforeAngle=torsion(bridgeBond.torsion.map(i=>s.atoms[i].coord));
  await run(`modelkit rotate branch {atomindex=${bridgeBond.a}} {atomindex=${bridgeBond.c}} 35;`);after=await snapshot();const afterAngle=torsion(bridgeBond.torsion.map(i=>after.atoms[i].coord)),delta=((afterAngle-beforeAngle+540)%360)-180;
  assert(Math.abs(Math.abs(delta)-35)<1e-5);assert.deepEqual(topology(after),topology(s));for(let i=0;i<s.atoms.length;i++)if(!bridgeBond.branch.includes(i))assert(displacement(s.atoms[i].coord,after.atoms[i].coord)<1e-9);
  check('ModelKit rotate branch performs a real 35 degree TPB torsion edit',{indices:bridgeBond.torsion,before:beforeAngle,after:afterAngle,delta,branchAtoms:bridgeBond.branch.length});
  await page.locator('#inspect').click();await page.screenshot({path:path.join(output,'tpb-torsion-live.png'),fullPage:true});
  const history=await page.evaluate(()=>{const model=hostModel,result=bridge.captureGeometry();const original=JSON.stringify(model.project);if(!model.applyAsync(result.baseRevision,result.project,'JSmol torsion'))throw Error('Apply failed');const edited=JSON.stringify(model.project);model.undo();const undo=JSON.stringify(model.project)===original;model.redo();return{undo,redo:JSON.stringify(model.project)===edited,staleRejected:!model.applyAsync(result.baseRevision,result.project,'obsolete')};});assert(history.undo&&history.redo&&history.staleRejected);check('Coordinate adapter preserves 2D/stereo/fragments and supports host undo/redo plus stale rejection',history);

  s=await load('tpb');const p=s.atoms[0].coord;await run(`modelkit add "C" {${p[0]+.1} ${p[1]} ${p[2]}};`);after=await snapshot();assert.equal(after.atoms.length,53);assert.equal(after.bonds.length,56);assert.equal(after.atoms.filter(a=>a.sym==='H').length,22);let proposal=await page.evaluate(()=>bridge.captureEditProposal());assert.equal(proposal.addedAtoms.length,1);assert.equal(proposal.addedBonds.length,0);const fresh=proposal.addedAtoms[0];assert(fresh.token>0);check('ModelKit atom addition near an existing atom adds no automatic bond or hydrogens',{newId:fresh.id,newJmolIndex:fresh.jmolIndex});
  await run(`{atomindex=${fresh.jmolIndex}}.atomname="custom_uuid"; modelkit assign atom {atomindex=${fresh.jmolIndex}} "O"; modelkit assign atom {atomindex=${fresh.jmolIndex}} "Pl";`);let edits=await page.evaluate(()=>bridge.read()),edited=edits.atoms.find(a=>a.id===fresh.id);assert.equal(edited.element,'O');assert.equal(edited.charge,1);assert.equal(edited.token,fresh.token);const atomname=await page.evaluate(i=>jsmolProbe.value(`{atomindex=${i}}.atomname`),fresh.jmolIndex);assert.notEqual(atomname,'custom_uuid');check('ModelKit element/charge editing preserves custom token; atomName is NOT stable',{element:edited.element,charge:edited.charge,token:edited.token,renamed:atomname});
  await run(`modelkit connect {atomindex=0} {atomindex=${fresh.jmolIndex}} "1";`);assert.equal((await snapshot()).bonds.length,57);
  for(const order of [2,3,1]){await run(`modelkit assign bond {atomindex=0 or atomindex=${fresh.jmolIndex}} "${order}";`);edits=await page.evaluate(()=>bridge.read());assert.equal(edits.bonds.find(b=>b.a===fresh.id||b.b===fresh.id).order,order);}
  await run(`modelkit assign bond {atomindex=0 or atomindex=${fresh.jmolIndex}} "0";`);assert.equal((await snapshot()).bonds.length,56);await run(`modelkit assign atom {atomindex=${fresh.jmolIndex}} "X";`);edits=await page.evaluate(()=>bridge.read());assert.equal(edits.atoms.length,52);assert.deepEqual(edits.atoms.map(a=>a.id),fixtures.tpb.project.atoms.map(a=>a.id));check('ModelKit creates, changes and deletes actual chemical bonds and deletes an atom',{testedBondOrders:[1,2,3],remainingAtoms:52,remainingBonds:56});

  await run('modelkit assign atom {atomindex=33} "X"; modelkit add "C" {12 12 5};');edits=await page.evaluate(()=>bridge.read());const survivorIds=fixtures.tpb.project.atoms.filter((_,i)=>i!==33).map(a=>a.id);assert(survivorIds.every(id=>edits.atoms.some(a=>a.id===id)));assert(!edits.atoms.some(a=>a.id===fixtures.tpb.project.atoms[33].id));assert.equal(edits.atoms.length,52);check('Middle deletion and later insertion keep every surviving ID without index/graph matching',{deletedId:fixtures.tpb.project.atoms[33].id,lastIndex:edits.atoms.at(-1).jmolIndex});
  const beforeTokens=edits.atoms.map(a=>a.token);await run('save state "identity_test"; {*}.property_editorid=0; restore state "identity_test";');
  const nativeRestore=await page.evaluate(()=>{const s=jsmolProbe.inspect();let rejected=false,commandBlocked=false;try{bridge.read({assignNew:false});}catch{rejected=true;}try{bridge.script('restore state "identity_test";');}catch{commandBlocked=true;}return{rejected,commandBlocked,tokens:s.ids,atoms:s.atoms.length};});
  assert(nativeRestore.rejected||JSON.stringify(nativeRestore.tokens)!==JSON.stringify(beforeTokens));
  assert(nativeRestore.commandBlocked);
  limitations.push({name:'Native JSmol state restoration after middle deletion/addition does not reliably restore custom identity',beforeTokens,afterTokens:nativeRestore.tokens,adapterRejected:nativeRestore.rejected,requiredPolicy:'Do not use native ModelKit undo/redo or state restoration for project history. Reload canonical host snapshots through an explicit MOL row map.'});
  check('Native history limitation is reproduced; project history must own canonical snapshots',{adapterRejected:nativeRestore.rejected,beforeAtoms:beforeTokens.length,afterAtoms:nativeRestore.atoms});

  await load('tpb');let raw=await snapshot();const hb=raw.bonds.find(b=>b.atom1.sym==='H'||b.atom2.sym==='H'),carbonIndex=hb.atom1.sym==='H'?hb.atom2.atomIndex:hb.atom1.atomIndex;
  await run(`modelkit assign atom {atomindex=${carbonIndex}} "N";`);proposal=await page.evaluate(()=>bridge.captureEditProposal());assert(proposal.removedAtoms.some(a=>a.element==='H'));assert(proposal.changedAtoms.some(a=>a.jmolIndex===carbonIndex));let rejected=await page.evaluate(()=>{try{bridge.captureGeometry();return false;}catch{return true;}});assert(rejected);
  limitations.push({name:'ModelKit element changes remove explicit attached H even with addHydrogens=false',removedIds:proposal.removedAtoms.map(a=>a.id),requiredPolicy:'Use explicit complete topology diffs; geometry-only bridge rejects this change.'});check('Adapter detects incidental H deletion and rejects topology mutation as a coordinate edit',{removedHydrogens:proposal.removedAtoms.filter(a=>a.element==='H').length});
  await load('tpb');await run('modelkit add "C" {15 15 0}; modelkit connect {atomindex=0} {atomindex=52} "1";');raw=await snapshot();await run('modelkit assign atom {atomindex=52} "O";');after=await snapshot();const changedDistance=displacement(raw.atoms.at(-1).coord,after.atoms.at(-1).coord),bondLength=after.bonds.find(b=>b.atom1.atomIndex===52||b.atom2.atomIndex===52).length_Ang;assert(changedDistance>1);assert(Math.abs(bondLength-1.5)<1e-6);
  const lengthProposal=await page.evaluate(()=>bridge.captureEditProposal());assert(lengthProposal.addedAtoms.some(a=>a.element==='O'));assert(lengthProposal.snapshot.atoms.some(a=>a.jmolIndex===52));
  limitations.push({name:'Changing a singly bonded element resets its distance to 1.5 angstrom',displacement:changedDistance,requiredPolicy:'Include coordinate changes in the proposal or use host-owned exact topology actions.'});check('Documented ModelKit length adjustment is observed and recorded',{angstroms:changedDistance,newBondLength:bondLength});

  s=await load('c60');assert.equal(s.atoms.length,60);assert.equal(s.bonds.length,90);assert(s.atoms.every(a=>a.bondCount===3));assert.deepEqual((await page.evaluate(()=>bridge.read())).atoms.map(a=>a.id),fixtures.c60.project.atoms.map(a=>a.id));
  const reversed=python(`import json,sys,copy
from molecular_car.editor.service import export_file
d=json.load(sys.stdin);p=copy.deepcopy(d);p['atoms'].reverse();print(json.dumps({'project':p,'mol':export_file({'project':p,'dimension':3})['content']}))`,fixtures.c60.project);
  await page.evaluate(f=>bridge.loadProject(f.project,f.mol,{indexToAtomId:f.project.atoms.map(a=>a.id),revision:19}),reversed);edits=await page.evaluate(()=>bridge.read());assert.deepEqual(edits.atoms.map(a=>a.id),reversed.project.atoms.map(a=>a.id));for(const a of edits.atoms){const expected=fixtures.c60.project.atoms.find(b=>b.id===a.id);assert(displacement(a.pos3,expected.pos3)<.0001);}check('C60 is a 60 atom/90 bond graph; reversed MOL records retain correct UUID mapping',{atoms:60,bonds:90,graphMatchingUsed:false});
  await page.locator('#inspect').click();await page.screenshot({path:path.join(output,'c60-live.png'),fullPage:true});

  await load('chiral');const canonical=await page.evaluate(()=>bridge.captureGeometry());await run('select {*}; translateSelected {0.3 -0.2 0.4};');const chiralAfter=await page.evaluate(()=>bridge.captureGeometry());assert.deepEqual(chiralAfter.project.atoms.map(a=>a.chiralTag),canonical.project.atoms.map(a=>a.chiralTag));const jmolMol=await page.evaluate(()=>bridge.exportMol());
  const stereo=python(`import json,sys,numpy as np
from rdkit import Chem
from molecular_car.editor.service import _stereo_signature
d=json.load(sys.stdin);a=Chem.MolFromMolBlock(d['before'],removeHs=False);b=Chem.MolFromMolBlock(d['after'],removeHs=False);assert b is not None;sa=_stereo_signature(a);sb=_stereo_signature(b);assert sa==sb,(sa,sb);assert '3D' in d['after'].splitlines()[1];error=float(np.max(np.abs((b.GetConformer().GetPositions()-a.GetConformer().GetPositions())-[.3,-.2,.4])));assert error<.00011,error;print(json.dumps({'before':sa,'after':sb,'threeDimensionalHeader':True,'maxMolRoundingErrorAngstrom':error}))`,{before:fixtures.chiral.mol,after:jmolMol});check('Real chiral MOL survives JSmol 3D translation and diagnostic MOL export with RDKit verification',stereo);
  assert(!jmolMol.includes('property_editorid'));assert(fixtures.chiral.project.atoms.every(a=>!jmolMol.includes(a.id)));
  await page.locator('#inspect').click();await page.screenshot({path:path.join(output,'chiral-live.png'),fullPage:true});
  assert.equal(errors.length,0);assert.equal(external.length,0);assert(missing.every(x=>x==='/editor/research/'));check('All browser code/resources local; zero external/API requests and zero page errors',{pageErrors:errors,blockedExternal:external,benignDirectoryProbe:missing});
 }catch(e){failure=e;checks.push({name:'Suite completion',passed:false,error:e.stack});throw e;}
 finally{
  const manifest=JSON.parse(fs.readFileSync(path.join(root,'frontend/editor/vendor/jsmol/manifest.json'),'utf8'));
  fs.writeFileSync(path.join(output,'report.json'),JSON.stringify({date:'2026-10-12',runtimeVersion:'16.3.31',source:manifest.source,archiveSha256:manifest.archiveSha256,scope:'Real local Chromium/JSmol runtime, direct mouse plus ModelKit scripts, explicit ID token bridge; NOT a finished two-way Ketcher/JSmol workshop integration.',checks,limitations,pageErrors:errors,blockedExternal:external,missingRequests:missing,passed:failure===null},null,2));
  await browser.close();
 }
}
if(require.main===module)verify().catch(e=>{console.error(e);process.exitCode=1;});
