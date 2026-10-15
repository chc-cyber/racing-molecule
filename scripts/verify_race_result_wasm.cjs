const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),{spawn}=require('node:child_process');
const {chromium}=require('C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root=path.resolve(__dirname,'..'),out=path.join(root,'outputs/race-result-check');
(async()=>{
  const server=spawn(path.join(root,'.venv/Scripts/python.exe'),['wasm/tools/serve.py','--port','0'],{cwd:root,windowsHide:true});
  server.stderr.on('data',()=>{});let browser;const report={passed:false,errors:[]};
  try{
    const base=await new Promise((resolve,reject)=>{const timer=setTimeout(()=>reject(Error('server timeout')),15000);server.once('error',reject);server.stdout.on('data',d=>{const m=String(d).match(/http:\/\/127\.0\.0\.1:\d+\//);if(m){clearTimeout(timer);resolve(m[0]);}});});
    browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-unsafe-swiftshader']});
    const page=await browser.newPage({viewport:{width:1440,height:1000},acceptDownloads:true});page.setDefaultTimeout(120000);
    page.on('pageerror',e=>report.errors.push(e.message));
    await page.goto(base+'race.html');
    await page.waitForFunction(()=>window.__molecularCarDebug?.getSnapshot().status==='idle');
    await page.locator('#example-select').selectOption('benzene.mol');await page.locator('#load-example-button').click();await page.locator('#prepare-race').click();
    await page.waitForFunction(()=>window.__molecularCarDebug.getSnapshot().status==='paused');
    const started=performance.now();await page.locator('#start-button').click();
    await page.locator('#race-result-dialog').waitFor({state:'visible'});report.wall_seconds=(performance.now()-started)/1000;
    const state=await page.evaluate(async()=>{const {request}=await import('./transport.js');return request('/api/state');});
    assert.equal(state.status,'finished');assert(state.frame.distance_angstrom>=200);
    const expected=(state.frame.distance_angstrom/state.frame.time_fs*100000).toFixed(1);
    assert.equal(await page.locator('#race-result-speed').textContent(),expected);
    assert.equal(await page.locator('#race-result-title').textContent(),'冲线成功！');
    report.time_ps=state.frame.time_fs/1000;report.distance_nm=state.frame.distance_angstrom/10;report.speed_ms=Number(expected);
    fs.mkdirSync(out,{recursive:true});await page.waitForTimeout(350);await page.screenshot({path:path.join(out,'wasm-finish.png')});
    const download=page.waitForEvent('download');await page.locator('#result-save').click();const file=await download;
    const replay=JSON.parse(fs.readFileSync(await file.path(),'utf8'));assert.equal(replay.frames.at(-1).time_fs,state.frame.time_fs);
    await page.locator('#result-again').click();await page.waitForFunction(()=>window.__molecularCarDebug.getSnapshot().status==='paused'&&window.__molecularCarDebug.getSnapshot().renderedStep===0);
    assert.equal(await page.locator('#race-result-dialog').isVisible(),false);
    assert.deepEqual(report.errors,[]);report.passed=true;
  }catch(error){report.error=error.stack;process.exitCode=1;}
  finally{if(browser)await browser.close();server.kill();fs.mkdirSync(out,{recursive:true});fs.writeFileSync(path.join(out,'wasm.json'),JSON.stringify(report,null,2));console.log(JSON.stringify(report));}
})();
