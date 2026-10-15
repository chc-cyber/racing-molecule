/* Isolated browser visual regression: camera/display only, no MD steps. */
const fs = require('node:fs'), path = require('node:path'), assert = require('node:assert/strict'), crypto = require('node:crypto');
const {chromium} = require('C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root = path.resolve(__dirname, '..'), out = path.join(root, 'outputs/gold-side-check');
const stage = process.argv[2] || 'after', base = 'http://127.0.0.1:8794';
const report = {stage, base, errors: [], warnings: [], shots: []};
const hash = value => crypto.createHash('sha256').update(JSON.stringify(value)).digest('hex');
const hook = `
window.__goldSideQA = {
 pose(p,t) {cameraMode='free';needsCameraFocus=false;orbit.enableDamping=false;orbit.target.set(...t);camera.position.set(...p);orbit.update();updateSceneScale();},
 inspect() {const casters=[];scene.traverse(o=>{if(o.castShadow)casters.push(o.name||o.type);});return {shadows:renderer.shadowMap.enabled,casters,edgeCount:goldSurface.edgeMeshes?.reduce((n,m)=>n+m.count,0)||0,edgeTriangles:goldSurface.edgeMeshes?.reduce((n,m)=>n+m.count*m.geometry.index.count/3,0)||0};}
};
`;
(async()=>{
 fs.mkdirSync(out,{recursive:true});
 const browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-unsafe-swiftshader']});
 const page=await browser.newPage({viewport:{width:1480,height:1080},deviceScaleFactor:1});
 const stable=()=>page.evaluate(()=>new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(()=>requestAnimationFrame(r)))));
 const snap=()=>page.evaluate(()=>window.__molecularCarDebug.getSnapshot());
 await page.route('**/favicon.ico',r=>r.fulfill({status:204}));
 // Freeze only the old surface in the before stage; both stages use the same
 // current interface/transport, so unrelated UI edits cannot change the test.
 // No files on any running user service are replaced.
 await page.route('**/app.js',async r=>{
  const response=await r.fetch();
  const source=await response.text();
  assert(source.includes('\nawait pollState();'));
  await r.fulfill({response,body:source.replace('\nawait pollState();','\n'+hook+'\nawait pollState();')});
 });
 if(stage==='before')await page.route('**/gold-surface.js',r=>r.fulfill({contentType:'text/javascript',body:fs.readFileSync(path.join(out,'before-gold-surface.js'),'utf8')}));
 page.on('pageerror',e=>report.errors.push(e.message));
 page.on('console',m=>{if(m.type()==='error')report.errors.push(m.text());if(m.type()==='warning')report.warnings.push(m.text());});
 const shot=async name=>{
  await stable();const s=await snap(), inspection=await page.evaluate(()=>window.__goldSideQA.inspect());
  assert.equal(s.status,'paused');assert.equal(s.renderedStep,0);assert.equal(hash(s.positions),report.initial.position_hash);
  assert.equal(inspection.shadows,false);assert.deepEqual(inspection.casters,[]);
  if(stage==='after') {assert(inspection.edgeCount>0);assert.equal(inspection.edgeCount,report.edgeCount);assert(s.goldSurface.fallbackHeight>0);assert(s.goldSurface.bottomHeight<0);}
  const filename=stage+'-'+name+'.png';await page.screenshot({path:path.join(out,filename),fullPage:true});
  report.shots.push({name,filename,camera:s.camera,target:s.cameraTarget,goldSurface:s.goldSurface,inspection,renderer:s.renderer,position_hash:hash(s.positions)});
  console.log(name,JSON.stringify({height:s.camera[1],localTriangles:s.goldSurface.detailTriangles ?? s.goldSurface.triangles,edgeTriangles:inspection.edgeTriangles}));
 };
 try {
  await page.goto(base+'/',{waitUntil:'domcontentloaded'});
  await page.waitForFunction(()=>window.__molecularCarDebug&&window.__goldSideQA,undefined,{timeout:15000});
  // The current interface starts idle: prepare the built-in car explicitly,
  // then keep it paused. Starting the race is deliberately never requested.
  if((await snap()).status==='idle')await page.locator('#prepare-race').evaluate(button=>button.click());
  await page.waitForFunction(()=>{const s=window.__molecularCarDebug?.getSnapshot();return s?.status==='paused'&&s.positions?.length>0;},undefined,{timeout:120000});
  const initial=await snap();report.initial={step:initial.renderedStep,atomCount:initial.atomCount,position_hash:hash(initial.positions)};
  report.edgeCount=(await page.evaluate(()=>window.__goldSideQA.inspect())).edgeCount;
  const poses=[['front',[100,35,180],[100,2,95]],['back',[100,35,-180],[100,2,-95]],
   ['right',[285,35,0],[208,2,0]],['left',[-85,35,0],[-8,2,0]],
   ['corner',[240,34,130],[180,2,80]],['grazing',[100,10,180],[100,2,100]],
   ['height-99',[100,99,150],[100,2,0]],['height-101',[100,101,150],[100,2,0]]];
  for(const [name,p,t]of poses){await page.evaluate(({p,t})=>window.__goldSideQA.pose(p,t),{p,t});await shot(name);}
  await page.locator('#home-camera').click();await shot('home');
  await page.locator('#space-filling-model').click();await shot('space-filling');
  await page.evaluate(()=>window.__goldSideQA.pose([100,38,145],[100,2,35]));
  const box=await page.locator('#scene canvas').boundingBox();await page.mouse.move(box.x+box.width*.6,box.y+box.height*.6);
  for(let i=0;i<12;i++){await page.mouse.wheel(0,250);await stable();}
  await shot('wheel-far');
  assert.deepEqual(report.errors,[]);assert.deepEqual(report.warnings,[]);report.passed=true;
 }catch(e){report.passed=false;report.failure=e.stack;throw e;}
 finally{fs.writeFileSync(path.join(out,stage+'-report.json'),JSON.stringify(report,null,2)+'\n');await browser.close();}
 console.log(JSON.stringify({passed:report.passed,shots:report.shots.length,errors:report.errors}));
})().catch(e=>{console.error(e);process.exitCode=1;});
