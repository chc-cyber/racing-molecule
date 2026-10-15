/* Static WASM visual check in an isolated browser context: preparation only,
 * no start/pulse/reset requests and no native service is addressed. */
const fs = require('node:fs'), path = require('node:path'), assert = require('node:assert/strict'), crypto = require('node:crypto');
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const base = 'http://127.0.0.1:8791', out = path.resolve(__dirname, '../outputs/gold-zoom-check');
const report = {base_url: base, created_utc: new Date().toISOString(), purpose: 'static WASM display; no MD steps', errors: [], warnings: [], shots: []};
const hash = p => crypto.createHash('sha256').update(JSON.stringify(p)).digest('hex');
const hook = `
window.__goldZoomQA = {
 setHeight(h) {cameraMode='free'; needsCameraFocus=false; orbit.enableDamping=false; const d=camera.position.clone().sub(orbit.target); d.multiplyScalar((h-orbit.target.y)/d.y); camera.position.copy(orbit.target).add(d); orbit.update(); updateSceneScale();},
 lowAngle() {cameraMode='free'; needsCameraFocus=false; orbit.enableDamping=false; orbit.target.set(15,2.2,0); camera.position.set(25,15,120); orbit.update(); updateSceneScale();},
 inspect() {const casts=[],receives=[]; scene.traverse(o=>{if(o.castShadow)casts.push(o.name||o.type);if(o.receiveShadow)receives.push(o.name||o.type);}); return {enabled:renderer.shadowMap.enabled,casts,receives,atomCastShadow:atomsMesh?.castShadow,bondCastShadow:bondsMesh?.castShadow,farHeight:goldSurface.farPlane.position.y,fallbackHeight:goldSurface.fallbackHeight ?? goldSurface.backing.position.y};}
};
`;
(async () => {
 const browser=await chromium.launch({executablePath:process.env.CHROME_PATH||'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-unsafe-swiftshader']});
 const page=await browser.newPage({viewport:{width:1480,height:1080},deviceScaleFactor:1});
 await page.route('**/favicon.ico',r=>r.fulfill({status:204}));
 await page.route('**/app.js',async r=>{const response=await r.fetch(); await r.fulfill({response,body:(await response.text()).replace(/\nawait pollState\(\);/,'\n'+hook+'\nawait pollState();')});});
 page.on('pageerror',e=>report.errors.push(e.message));
 page.on('console',m=>{if(m.type()==='error')report.errors.push(m.text());if(m.type()==='warning')report.warnings.push(m.text());});
 const stable=()=>page.evaluate(()=>new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(()=>requestAnimationFrame(r)))));
 const snap=()=>page.evaluate(()=>window.__molecularCarDebug.getSnapshot());
 const shot=async name=>{
  await stable(); const s=await snap(), q=await page.evaluate(()=>window.__goldZoomQA.inspect());
  assert.equal(s.status,'paused');assert.equal(s.renderedStep,0);assert.equal(s.atomCount,52);assert.equal(hash(s.positions),report.initial.position_hash);
  assert.equal(q.enabled,false);assert.equal(q.atomCastShadow,false);assert.equal(q.bondCastShadow,false);assert.deepEqual(q.casts,[]);assert.deepEqual(q.receives,[]);
  assert(q.farHeight>=s.goldSurface.radius);assert(q.fallbackHeight>0);
  const filename='wasm-after-'+name+'.png'; await page.screenshot({path:path.join(out,filename),fullPage:true});
  report.shots.push({name,filename,camera:s.camera,cameraTarget:s.cameraTarget,step:s.renderedStep,position_hash:hash(s.positions),goldSurface:s.goldSurface,inspection:q}); console.log(name,s.camera[1],s.goldSurface.triangles);
 };
 try {
  await page.goto(base+'/',{waitUntil:'domcontentloaded'});
  await page.waitForFunction(()=>{const s=window.__molecularCarDebug?.getSnapshot();return s?.status==='paused'&&s.atomCount===52&&s.positions?.length===156;},undefined,{timeout:240000});
  const s=await snap(); report.initial={status:s.status,step:s.renderedStep,generation:s.generation,surface:s.surface,position_hash:hash(s.positions),clipping:s.clipping}; assert.equal(s.renderedStep,0);
  await page.locator('#home-camera').click();await stable();
  await page.evaluate(()=>window.__goldZoomQA.setHeight(101));await shot('height-101');
  await page.evaluate(()=>window.__goldZoomQA.lowAngle());await shot('low-angle');
  await page.locator('#home-camera').click();await stable();
  const box=await page.locator('#scene canvas').boundingBox();await page.mouse.move(box.x+box.width*.6,box.y+box.height*.55);report.wheel_sequence=[];
  for(let i=0;i<18;i++){await page.mouse.wheel(0,250);await stable();const v=await snap();report.wheel_sequence.push({height:v.camera[1],radius:v.goldSurface.patch?.radius,triangles:v.goldSurface.triangles});if(v.camera[1]>=200)break;}
  await shot('wheel-far');
  await page.locator('#home-camera').click();await stable();assert.deepEqual(report.errors,[]);report.passed=true;
 }catch(error){report.passed=false;report.failure=error.stack;throw error;}
 finally{fs.writeFileSync(path.join(out,'wasm-after-report.json'),JSON.stringify(report,null,2)+'\n');await browser.close();}
 console.log(JSON.stringify({passed:report.passed,shots:report.shots.length,errors:report.errors,warnings:report.warnings}));
})().catch(e=>{console.error(e);process.exitCode=1;});
