// Headless acceptance on the isolated, newly built static browser preview.
const fs=require('node:fs');
const path=require('node:path');
const assert=require('node:assert/strict');
const {spawn}=require('node:child_process');
const {chromium}=require('C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const out=path.resolve(__dirname,'../outputs/crystal-contrast-guarded');
(async()=>{
  const root=path.resolve(__dirname,'..');
  const helper=spawn(path.join(root,'.venv/Scripts/python.exe'),['-u',path.join(root,'wasm/tools/serve.py'),'--port','0'],{cwd:root,windowsHide:true,stdio:['ignore','pipe','pipe']});
  const base=await new Promise((resolve,reject)=>{
    let output='';helper.stdout.on('data',b=>{output+=b;const match=output.match(/Molecular Car browser preview: (http:\/\/127\.0\.0\.1:\d+\/)/);if(match)resolve(match[1]);});
    helper.on('error',reject);helper.on('exit',code=>reject(new Error('Preview startup exit '+code)));
  });
  const browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-unsafe-swiftshader']});
  const page=await browser.newPage({viewport:{width:1480,height:1100}});const report={errors:[],checks:[]};
  page.on('pageerror',e=>report.errors.push(e.message));
  await page.route('**/favicon.ico',route=>route.fulfill({status:204}));
  const snapshot=()=>page.evaluate(()=>window.__molecularCarDebug.getSnapshot());
  const ready=()=>page.waitForFunction(()=>window.__molecularCarDebug?.getSnapshot().status==='paused',null,{timeout:180000});
  try {
    await page.goto(base,{waitUntil:'domcontentloaded'});await ready();
    await page.locator('.advanced-section > summary').click();
    assert.equal(await page.locator('#lattice_scale').inputValue(),'1');
    assert.equal(await page.locator('#growth').inputValue(),'1');
    for(const [id,value] of [['lattice_scale','2'],['growth','.25']]) {
      const before=await snapshot();
      await page.locator('#'+id).evaluate((el,v)=>{el.value=v;el.dispatchEvent(new Event('input',{bubbles:true}));el.dispatchEvent(new Event('change',{bubbles:true}));},value);
      await page.waitForFunction(g=>{const s=window.__molecularCarDebug.getSnapshot();return s.generation>g&&s.status==='paused';},before.generation,{timeout:180000});
      const after=await snapshot();assert.equal(after.renderedStep,0);
      if(id==='lattice_scale')assert.equal(after.surfaceMetadata.lattice_scale,2);else assert.equal(after.drive.growth,.25);
      report.checks.push({id,new_generation:after.generation,step:after.renderedStep});
    }
    await page.locator('#start-button').click();
    await page.waitForFunction(()=>window.__molecularCarDebug.getSnapshot().renderedStep>=500,null,{timeout:60000});
    await page.locator('#start-button').click();await ready();
    const final=await snapshot();assert.equal(final.surfaceMetadata.lattice_scale,2);assert.equal(final.drive.growth,.25);
    await page.locator('#growth').scrollIntoViewIfNeeded();
    await page.screenshot({path:path.join(out,'controls-preview.png'),fullPage:true});
    await page.locator('#tip-drive').click();await page.waitForFunction(()=>document.querySelector('#growth').disabled);
    await page.locator('#surface-mode').selectOption('continuum');await ready();assert(await page.locator('#lattice_scale').isDisabled());
    assert.deepEqual(report.errors,[]);report.status='passed';report.final=final;
  }catch(error){report.status='failed';report.error=error.message;}
  finally{await browser.close();helper.kill();fs.writeFileSync(path.join(out,'browser-controls.json'),JSON.stringify(report,null,2));}
  console.log(JSON.stringify({status:report.status,checks:report.checks,error:report.error}));if(report.status!=='passed')process.exitCode=1;
})();
