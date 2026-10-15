/* Real main-page seed-only conversion regression (no automatic UFF). Never starts/stops a server and
   unconditionally blocks /api/load. Native canvas actions are mouse gestures;
   the test does not alter production modules or substitute geometry results. */
'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const base=process.argv[2]||'http://127.0.0.1:8786';
assert.equal(base,'http://127.0.0.1:8786','Use only the independent editor verification service');
const root=path.resolve(__dirname,'..'),out=path.join(root,'outputs/editor-seed-check');
fs.mkdirSync(out,{recursive:true});
const sources=['frontend/editor/index.html','frontend/editor/workshop.css','frontend/editor/workshop.js','frontend/editor/model.js','frontend/editor/ketcher-host.js','frontend/editor/ketcher-adapter.js','frontend/editor/jsmol-host.js','frontend/editor/jsmol-adapter.js','molecular_car/editor/service.py','molecular_car/editor/construction.py','molecular_car/editor/superatoms.py','frontend/editor/vendor/ketcher-c60-manifest.json'];
const hashes=()=>Object.fromEntries(sources.map(file=>[file,crypto.createHash('sha256').update(fs.readFileSync(path.join(root,file))).digest('hex')]));
const report={started:new Date().toISOString(),base,sourceHashes:hashes(),checks:[],geometryRequests:[],geometryResponses:[],pageErrors:[],externalRequests:[],loadCalls:0};
const distance=(a,b)=>Math.hypot(...a.map((value,index)=>value-b[index]));
const lengths=(p,bonds)=>{const atoms=new Map(p.atoms.map(a=>[a.id,a]));return bonds.map(b=>distance(atoms.get(b.a).pos3,atoms.get(b.b).pos3));};
(async()=>{
 const browser=await chromium.launch({executablePath:process.env.CHROME_PATH||'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-unsafe-swiftshader']});
 const page=await browser.newPage({viewport:{width:1480,height:960}});page.setDefaultTimeout(30000);report.browser=browser.version();
 page.on('pageerror',error=>report.pageErrors.push(error.message));
 const responses=[];
 page.on('response',response=>{
  if(response.url()!==base+'/api/editor/geometry')return;
  responses.push((async()=>{const body=await response.json(),request=response.request().postDataJSON();await response.finished();report.geometryResponses.push({operation:request.operation,status:response.status(),construction:body.construction,hasEnergy:Object.hasOwn(body,'energy'),hydrogensAdded:body.hydrogensAdded,responseMs:response.request().timing().responseEnd});})());
 });
 await page.route('**/*',route=>{
  const url=route.request().url();
  if(url===base+'/api/load'){report.loadCalls++;return route.abort();}
  if(!url.startsWith(base+'/')&&!url.startsWith('blob:')&&!url.startsWith('data:')){report.externalRequests.push(url);return route.abort();}
  if(url===base+'/api/editor/geometry')report.geometryRequests.push(route.request().postDataJSON().operation);
  return route.continue();
 });
 const project=()=>page.evaluate(()=>structuredClone(workshop.model.project));
 const settled=async()=>{
  await page.waitForFunction(()=>workshop.mature?.ready&&workshop.views?.ready&&!workshop.mature.suspended&&!workshop.views.suspended&&!workshop.mature.gesture&&!workshop.views.gesture,null,{timeout:90000});
  await page.evaluate(async()=>{await workshop.mature.settled();await workshop.views.settled();});
  await page.waitForFunction(()=>!workshop.mature.suspended&&!workshop.views.suspended&&!workshop.mature.dirty&&!workshop.views.pending);
 };
 const record=async(name,details={})=>{
  await settled();const file=`${String(report.checks.length+1).padStart(2,'0')}.json`;
  const state=await page.evaluate(()=>({project:structuredClone(workshop.model.project),revision:workshop.model.revision,dimension:workshop.dimension,jsmol:workshop.views.adapter.read({assignNew:false})}));
  fs.writeFileSync(path.join(out,file),JSON.stringify({name,details,...state},null,2)+'\n');
  report.checks.push({name,passed:true,evidence:file,details});console.log('PASS',name);
 };
 const assertTools=async()=>{
  await settled();const p=await project(),native=await page.evaluate(()=>({count:workshop.mature.ketcher.editor.struct().atoms.size,bonds:workshop.mature.ketcher.editor.struct().bonds.size,jsmol:workshop.views.adapter.read({assignNew:false})}));
  assert.equal(native.count,p.atoms.length);assert.equal(native.bonds,p.bonds.length);
  assert.deepEqual(native.jsmol.atoms.map(a=>a.id),p.atoms.map(a=>a.id));assert.equal(native.jsmol.bonds.length,p.bonds.length);
  for(const atom of native.jsmol.atoms)assert(distance(atom.pos3,p.atoms.find(a=>a.id===atom.id).pos3)<0.00002,'JSmol coordinate mapping differs');
  return p;
 };
 let frame;
 const fresh=async()=>{await page.locator('#file-menu>summary').click();await page.locator('#new-project').click();await settled();};
 const blank=async(dx=0,dy=0)=>{const b=await frame.locator('[data-testid=canvas]').boundingBox();return{x:b.x+b.width*.45+dx,y:b.y+b.height*.44+dy};};
 const point2=id=>page.evaluate(id=>{
  const h=workshop.mature,r=h.ketcher.editor.render,s=h.ketcher.editor.struct(),entry=[...h.adapter.atoms.entries()].find(([pool,a])=>a.id===id&&s.atoms.has(pool));
  if(!entry)throw new globalThis.Clazz._Error('Missing Ketcher UUID');
  const pp=s.atoms.get(entry[0]).pp,q=r.page2obj({clientX:0,clientY:0}),u=r.page2obj({clientX:100,clientY:0}),v=r.page2obj({clientX:0,clientY:100}),a=(u.x-q.x)/100,b=(v.x-q.x)/100,c=(u.y-q.y)/100,d=(v.y-q.y)/100,det=a*d-b*c,dx=pp.x-q.x,dy=pp.y-q.y,box=h.frame.getBoundingClientRect();
  return{x:box.left+(d*dx-b*dy)/det,y:box.top+(-c*dx+a*dy)/det};
 },id);
 const click=async p=>{await page.mouse.click(p.x,p.y);await page.mouse.move(1460,0);await settled();};
 const drag=async(a,b)=>{await page.mouse.move(a.x,a.y);await page.mouse.down();await page.mouse.move(b.x,b.y,{steps:12});await page.mouse.up();await page.mouse.move(1460,0);await settled();};
 const benzene=async(dx=0,dy=0)=>{await frame.locator('[data-testid=template-0]').filter({visible:true}).click();await click(await blank(dx,dy));};
 const selectNative=async()=>{await frame.locator('[data-testid=select-rectangle]').filter({visible:true}).click();await page.mouse.move(1460,0);await settled();};
 const switchTo=async dim=>{
  await page.locator(`#view-${dim}d`).click();await page.waitForFunction(d=>workshop.dimension===d,dim,{timeout:90000});
  await page.waitForFunction(()=>workshop.pendingGeometry===0,null,{timeout:90000});await settled();
  await page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve))));
 };
 const hydrogenLabels=()=>page.evaluate(()=>{
  const h=workshop.mature,r=h.ketcher.editor.render,s=h.ketcher.editor.struct(),q=r.page2obj({clientX:0,clientY:0}),u=r.page2obj({clientX:100,clientY:0}),v=r.page2obj({clientX:0,clientY:100}),a=(u.x-q.x)/100,b=(v.x-q.x)/100,c=(u.y-q.y)/100,d=(v.y-q.y)/100,det=a*d-b*c;
  return [...h.adapter.atoms.entries()].filter(([pool,atom])=>s.atoms.has(pool)&&workshop.model.project.atoms.find(x=>x.id===atom.id).element==='H').map(([pool,atom])=>{
   const pp=s.atoms.get(pool).pp,node=r.ctab.atoms.get(pool)?.label?.path?.node,box=node?.getBoundingClientRect(),dx=pp.x-q.x,dy=pp.y-q.y,x=(d*dx-b*dy)/det,y=(-c*dx+a*dy)/det;
   return{id:atom.id,text:node?.textContent,visible:!!box&&box.width>0,offset:box?Math.hypot(box.x+box.width/2-x,box.y+box.height/2-y):null};
  });
 });
 try{
  await page.goto(base+'/editor/index.html');await page.waitForFunction(()=>window.workshop?.mature?.ready&&workshop.views?.ready,null,{timeout:90000});
  frame=page.frames().find(f=>/\/vendor\/ketcher-c60\/index\.html(?:\?|$)/.test(f.url()));assert(frame);await settled();
  report.initialGeneration=(await(await page.request.get(base+'/api/state')).json()).generation;
  await fresh();await benzene();let p=await assertTools();assert.equal(p.atoms.length,6);assert.equal(p.bonds.length,6);const drawing=structuredClone(p);
  await switchTo(3);p=await assertTools();assert.equal(p.atoms.length,12);assert.equal(p.bonds.length,12);assert(p.atoms.filter(a=>a.element==='H').every(a=>a.generatedHydrogen===true));
  assert.deepEqual(p.atoms.slice(0,6).map(a=>[a.id,a.pos2]),drawing.atoms.map(a=>[a.id,a.pos2]));assert(!p.metadata?.geometryPending);
  await record('Native benzene switches through real prepare and creates six stable real H atoms',{requests:report.geometryRequests.slice()});
  const first3D=structuredClone(p);await switchTo(2);let labels=await hydrogenLabels();assert.equal(labels.length,6);assert(labels.every(label=>label.visible&&label.text==='H'&&label.offset<7),JSON.stringify(labels));
  assert.deepEqual((await project()).atoms,first3D.atoms);await page.screenshot({path:path.join(out,'hydrogen-2d-labels.png'),fullPage:true});
  await record('Returning from hidden 3D places every H label at its true native atom position',{labels});
  await benzene(210);await selectNative();p=await project();let newIds=p.atoms.filter(a=>!first3D.atoms.some(old=>old.id===a.id)).map(a=>a.id);assert.equal(newIds.length,6);
  const oldCarbons=first3D.atoms.filter(a=>a.element==='C'),newCarbons=p.atoms.filter(a=>newIds.includes(a.id));
  const oldPoints=await Promise.all(oldCarbons.map(async atom=>({id:atom.id,...await point2(atom.id)}))),newPoints=await Promise.all(newCarbons.map(async atom=>({id:atom.id,...await point2(atom.id)})));
  const oldAnchor=oldPoints.sort((a,b)=>b.x-a.x)[0],newAnchor=newPoints.sort((a,b)=>a.x-b.x)[0];
  await frame.locator('[data-testid=bonds]').filter({visible:true}).click();await drag(oldAnchor,newAnchor);await selectNative();p=await assertTools();
  assert(p.bonds.some(b=>[b.a,b.b].includes(oldAnchor.id)&&[b.a,b.b].includes(newAnchor.id)),'Native bond gesture did not connect the rings');
  const removed=first3D.atoms.filter(a=>!p.atoms.some(b=>b.id===a.id));assert.equal(removed.length,1);assert(removed[0].element==='H'&&removed[0].generatedHydrogen,'Only the explicitly substituted automatic H may be removed');
  const beforePrepare=structuredClone(p);await switchTo(3);p=await assertTools();
  for(const atom of first3D.atoms.filter(a=>a.id!==removed[0].id))assert.deepEqual(p.atoms.find(a=>a.id===atom.id),atom);
  for(const atom of beforePrepare.atoms)assert.deepEqual(p.atoms.find(a=>a.id===atom.id).pos2,atom.pos2);
  const ringBonds=p.bonds.filter(b=>newIds.includes(b.a)&&newIds.includes(b.b)),ringLengths=lengths(p,ringBonds);assert.equal(ringBonds.length,6);assert(ringLengths.every(length=>length>1.3&&length<1.55),ringLengths.join(','));
  assert.equal(p.atoms.length,22);assert(!p.metadata?.geometryPending);await page.screenshot({path:path.join(out,'connected-rings-3d.png'),fullPage:true});
  await record('Native added and connected phenyl prepares at chemical scale and preserves every surviving old coordinate',{ringLengths,removedAutomaticHydrogen:removed[0].id});
  const biphenyl=structuredClone(p);await switchTo(2);await switchTo(3);assert.deepEqual(await project(),biphenyl);await record('Repeated view switches without topology changes are idempotent and do not optimize the existing conformer');
  await page.locator('#file-menu>summary').click();await page.locator('#example-tpb').click();await page.waitForFunction(()=>workshop.model.project.atoms.length===52);await settled();
  const ordered=await page.evaluate(async()=>{const {geometryPreview,dihedral}=await import('./model.js'),p=workshop.model.project;for(const bond of p.bonds){const near=id=>p.bonds.filter(b=>b.a===id||b.b===id).map(b=>b.a===id?b.b:b.a).filter(i=>i!==bond.a&&i!==bond.b&&p.atoms.find(a=>a.id===i).element!=='H');for(const a of near(bond.a))for(const d of near(bond.b))try{const ids=[a,bond.a,bond.b,d];geometryPreview(p,'dihedral',ids,dihedral(p,...ids));return ids;}catch{}}throw new globalThis.Clazz._Error('No TPB torsion');});
  await page.evaluate(ids=>workshop.beginGeometry('dihedral',ids),ordered);await page.locator('#geometry-value').fill('71.23456789');await page.locator('#geometry-apply').click();await settled();
  const manualTpb=await project();await switchTo(2);await benzene(350,-120);await selectNative();p=await project();assert.equal(p.atoms.length,58);await switchTo(3);p=await assertTools();
  for(const atom of manualTpb.atoms)assert.deepEqual(p.atoms.find(a=>a.id===atom.id),atom);
  const torsion=await page.evaluate(async ids=>(await import('./model.js')).dihedral(workshop.model.project,...ids),ordered);assert(Math.abs(torsion-71.23456789)<1e-7);
  await record('TPB hand-edited torsion and every old atom survive a new native ring and automatic preparation',{torsion});
  await fresh();await frame.locator('[data-testid=template-8]').filter({visible:true}).click();await click(await blank());await selectNative();p=await assertTools();assert.equal(p.atoms.length,60);assert.equal(p.bonds.length,90);
  const group=p.metadata?.ketcherGroups?.find(g=>g.kind==='c60');assert(group&&group.type==='SUP'&&group.expanded===false&&group.atomIds.length===60);
  const nativeGroups=await frame.evaluate(()=>[...window.ketcher.editor.struct().sgroups].map(([id,g])=>({id,type:g.type,name:g.data?.name,expanded:g.data?.expanded,atoms:g.atoms.length})));
  assert(nativeGroups.some(g=>g.type==='SUP'&&g.atoms===60),'C60 is not a native SUP display group');
  await page.screenshot({path:path.join(out,'c60-native-collapsed-2d.png'),fullPage:true});await record('Native C60 template inserts a collapsed SUP formula with sixty real mapped atoms and ninety bonds',{nativeGroups});
  const c60Drawing=structuredClone(p);await switchTo(3);p=await assertTools();assert.equal(p.atoms.length,60);assert.equal(p.bonds.length,90);assert.deepEqual(p.atoms.map(a=>[a.id,a.pos2]),c60Drawing.atoms.map(a=>[a.id,a.pos2]));
  const center=[0,1,2].map(i=>p.atoms.reduce((sum,a)=>sum+a.pos3[i],0)/60),radii=p.atoms.map(a=>distance(a.pos3,center));assert(radii.every(radius=>radius>3.4&&radius<3.7));
  for(let i=0;i<60;i++)for(let j=i+1;j<60;j++)assert(Math.abs(distance(p.atoms[i].pos3,p.atoms[j].pos3)-distance(c60Drawing.atoms[i].pos3,c60Drawing.atoms[j].pos3))<1e-8,'Sphere seed lost its rigid internal geometry');
  await page.screenshot({path:path.join(out,'c60-sphere-3d.png'),fullPage:true});await record('C60 switches to a genuine spherical conformer while preserving all template IDs and all pair distances',{radiusMin:Math.min(...radii),radiusMax:Math.max(...radii)});
  const sphere=structuredClone(p);await switchTo(2);await switchTo(3);assert.deepEqual(await project(),sphere);await record('Collapsed C60 and prepared 3D sphere survive repeated view switches without regeneration');
  // Draw a real methyl branch using Ketcher's native bond tool. This covers the
  // previously failing entire-cage embedding path, before and after first 3D.
  for(const existing3D of [false,true]){
   await fresh();await frame.locator('[data-testid=template-8]').filter({visible:true}).click();await click(await blank());await selectNative();
   if(existing3D){await switchTo(3);await switchTo(2);}
   const cage=await project(),group=cage.metadata.ketcherGroups.find(g=>g.kind==='c60'),cageIds=group.atomIds;
   const anchor=group.attachmentPoints[0].atomId,a=await point2(anchor);
   await frame.locator('[data-testid=bonds]').filter({visible:true}).click();await drag(a,{x:a.x+125,y:a.y});await selectNative();
   let branch=await assertTools();assert.equal(branch.atoms.length,61,'Native bond tool must add exactly one methyl carbon');assert.equal(branch.bonds.length,91);
   const methyl=branch.atoms.find(atom=>!cageIds.includes(atom.id));assert.equal(methyl.element,'C');assert(branch.bonds.some(bond=>[bond.a,bond.b].includes(anchor)&&[bond.a,bond.b].includes(methyl.id)));
   const drawn=structuredClone(branch),start=Date.now();await switchTo(3);branch=await assertTools();const elapsedMs=Date.now()-start;
   assert.equal(branch.atoms.length,65);assert.equal(branch.bonds.length,95);assert.equal(branch.atoms.filter(atom=>atom.element==='H').length,4);
   for(const atom of cage.atoms)assert.deepEqual(branch.atoms.find(next=>next.id===atom.id).pos3,atom.pos3,'All sixty cage coordinates must be reused exactly');
   for(const atom of drawn.atoms)assert.deepEqual(branch.atoms.find(next=>next.id===atom.id).pos2,atom.pos2);
   const byId=new Map(branch.atoms.map(atom=>[atom.id,atom])),connectionLength=distance(byId.get(anchor).pos3,byId.get(methyl.id).pos3);
   assert(connectionLength>1.4&&connectionLength<1.7,connectionLength);const hydrogenLengths=lengths(branch,branch.bonds.filter(bond=>byId.get(bond.a).element==='H'||byId.get(bond.b).element==='H'));
   assert(hydrogenLengths.every(value=>value>1&&value<1.2),hydrogenLengths.join(','));
   const center=[0,1,2].map(axis=>cageIds.reduce((sum,id)=>sum+byId.get(id).pos3[axis],0)/60);
   assert(distance(byId.get(methyl.id).pos3,center)>distance(byId.get(anchor).pos3,center)+1,'Methyl must point outside the cage');
   await page.screenshot({path:path.join(out,existing3D?'c60-methyl-appended.png':'c60-methyl-first.png'),fullPage:true});
   await record(existing3D?'An existing C60 cage gains a real native methyl branch without re-embedding or moving old atoms':'A native C60-methyl drawing switches directly to a real spherical hydrogenated conformer',{elapsedMs,connectionLength,hydrogenLengths});
   const prepared=structuredClone(branch);await page.locator('#undo').click();await settled();assert.deepEqual(await project(),drawn);await page.locator('#redo').click();await settled();assert.deepEqual(await project(),prepared);
   await switchTo(2);assert((await hydrogenLabels()).every(label=>label.visible&&label.offset<7));await switchTo(3);assert.deepEqual(await project(),prepared);
   await record('Methyl conversion '+(existing3D?'after existing 3D':'from initial drawing')+' shares exact undo/redo and idempotent view history');
  }
  await fresh();await benzene();let release,startedResolve;const held=new Promise(resolve=>release=resolve),started=new Promise(resolve=>startedResolve=resolve);
  await page.route('**/api/editor/geometry',async route=>{const request=route.request().postDataJSON();report.geometryRequests.push(request.operation);const response=await route.fetch();startedResolve();await held;await route.fulfill({response});});
  await page.locator('#view-3d').click();await started;await benzene(240);await selectNative();const whilePending=await project();assert.equal(whilePending.atoms.length,12);assert(!whilePending.atoms.some(atom=>atom.element==='H'));release();
  await page.waitForFunction(()=>workshop.pendingGeometry===0,null,{timeout:90000});await settled();await page.unroute('**/api/editor/geometry');assert.deepEqual(await project(),whilePending);assert.equal(await page.evaluate(()=>workshop.dimension),2);assert(/旧计算结果已丢弃/.test(await page.locator('#feedback').textContent()));
  await record('A delayed real prepare result cannot replace a newer native canvas edit',{currentAtoms:whilePending.atoms.length});
  await Promise.all(responses);assert(report.geometryResponses.length);assert(report.geometryResponses.every(response=>response.operation==='prepare'&&response.status===200&&response.construction.optimized===false&&!response.hasEnergy),'Automatic conversion must only construct initial positions and real H');
  report.finalGeneration=(await(await page.request.get(base+'/api/state')).json()).generation;assert.equal(report.finalGeneration,report.initialGeneration);assert.equal(report.loadCalls,0);assert.deepEqual(report.pageErrors,[]);assert.deepEqual(report.externalRequests,[]);
  await record('All preparation checks leave the race untouched and make no external requests',{loadCalls:report.loadCalls});report.passed=report.checks.length;
 }catch(error){report.failure=error.stack;report.feedback=await page.locator('#feedback').textContent().catch(()=>null);throw error;}
 finally{report.completed=new Date().toISOString();const final=hashes();report.sourcesChangedDuringRun=sources.filter(file=>report.sourceHashes[file]!==final[file]);fs.writeFileSync(path.join(out,'browser-report.json'),JSON.stringify(report,null,2)+'\n');await page.screenshot({path:path.join(out,'last-state.png'),fullPage:true}).catch(()=>{});await browser.close();}
 console.log(JSON.stringify({passed:report.passed,pageErrors:report.pageErrors,loadCalls:report.loadCalls,sourcesChangedDuringRun:report.sourcesChangedDuringRun}));
})().catch(error=>{console.error(error);process.exitCode=1;});

