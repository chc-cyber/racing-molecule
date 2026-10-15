/* Render fixtures through the production HTTP handler; no user's race is touched. */
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const {spawn}=require('node:child_process');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root=path.resolve(__dirname,'..'),out=path.join(root,'outputs/race-result-check');
(async()=>{
  const server=spawn(path.join(root,'.venv/Scripts/python.exe'),['-c',
    "from molecular_car.server import LocalHTTPServer, Handler, GameController; from functools import partial; s=LocalHTTPServer(('127.0.0.1',0),partial(Handler,game=GameController())); print(f'http://127.0.0.1:{s.server_port}/',flush=True); s.serve_forever()"],{cwd:root,windowsHide:true});
  server.stderr.on('data',()=>{});
  let browser,page;const report={passed:false,errors:[],fixtures:[]};
  try {
    const base=await new Promise((resolve,reject)=>{const t=setTimeout(()=>reject(Error('server timeout')),15000);server.once('error',reject);server.stdout.on('data',d=>{const m=String(d).match(/http:\/\/127\.0\.0\.1:\d+\//);if(m){clearTimeout(t);resolve(m[0]);}});});
    browser=await chromium.launch({executablePath:process.env.CHROME_PATH||'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-unsafe-swiftshader']});
    page=await browser.newPage({viewport:{width:1440,height:1000}});
    page.on('pageerror',e=>report.errors.push(e.message));
    const state={generation:1,status:'paused',settings:{surface:'fourier-au',temperature:150,damping_fs:350,corrugation:.08,lattice_scale:1.5},
      track:{start_x:0,finish_x:200,half_width:100,surface_z:0},molecule:{name:'finish fixture',atoms:[{id:1,element:'C',mass:12}],bonds:[],notes:[]},
      control:{tip:{x:6,y:0,z:7,strength:.6},drive:{mode:'gradient',strength:.6,growth:1},autopilot:true},runtime:{},
      frame:{seq:1,step:100,time_fs:100,center:[200,0,3.4],positions:[200,0,3.4],speed_angstrom_ps:0,energy:0,
        internal_temperature:150,distance_angstrom:199.3,progress:199.3/200,race_start_x:.7,race_length_angstrom:200}};
    await page.route('**/api/**',route=>route.fulfill({contentType:'application/json',body:JSON.stringify(route.request().url().includes('/frames')
      ?{frames:state.frame?[state.frame]:[],latest_seq:state.frame?.seq??-1,generation:state.generation,status:state.status}:state)}));
    await page.goto(base);
    await page.waitForFunction(()=>window.__molecularCarDebug?.getSnapshot().finishX===200.7);
    await page.waitForFunction(()=>document.querySelector('#distance-value').textContent==='19.93');
    report.fixtures.push({origin:.7,finishX:(await page.evaluate(()=>window.__molecularCarDebug.getSnapshot())).finishX});
    state.generation++;state.frame={...state.frame,seq:2,center:[199.31,0,3.4],positions:[199.31,0,3.4],race_start_x:-.7,distance_angstrom:200.01,progress:1};
    state.status='finished';state.finish_reason='finish';state.frame.time_fs=100000;state.frame.step=100000;
    await page.waitForFunction(()=>window.__molecularCarDebug?.getSnapshot().finishX===199.3);
    await page.waitForFunction(()=>Number(document.querySelector('#distance-value').textContent)>=20);
    report.fixtures.push({origin:-.7,finishX:(await page.evaluate(()=>window.__molecularCarDebug.getSnapshot())).finishX,distance_nm:await page.locator('#distance-value').textContent()});
    await page.locator('#race-result-dialog').waitFor({state:'visible'});
    assert.equal(await page.locator('#race-result-speed').textContent(),'200.0');
    assert.equal(await page.locator('#race-result-time').textContent(),'100.00');
    assert.equal(await page.locator('#race-result-distance').textContent(),'20.00');
    fs.mkdirSync(out,{recursive:true});
    await page.screenshot({path:path.join(out,'desktop.png')});
    await page.keyboard.press('Escape');await page.waitForTimeout(1200);
    assert.equal(await page.locator('#race-result-dialog').isVisible(),false);
    await page.locator('#show-race-result').click();
    await page.setViewportSize({width:390,height:844});
    await page.screenshot({path:path.join(out,'mobile.png'),fullPage:true});
    const bounds=await page.locator('#race-result-dialog').boundingBox();assert(bounds.x>=0&&bounds.x+bounds.width<=390);
    state.generation++;state.status='preparing';state.frame=null;state.molecule=null;
    await page.locator('#race-result-dialog').waitFor({state:'hidden'});
    assert.equal(await page.locator('#show-race-result').isVisible(),false);
    state.molecule={name:'Time limit car',atoms:[{id:1,element:'C',mass:12}],bonds:[],notes:[]};
    state.frame={seq:1,step:5000000,time_fs:5000000,center:[10,0,3.4],positions:[10,0,3.4],distance_angstrom:10,progress:.05,race_start_x:0,race_length_angstrom:200};
    state.status='finished';state.finish_reason='time_limit';
    await page.locator('#race-result-dialog').waitFor({state:'visible'});
    assert.equal(await page.locator('#race-result-title').textContent(),'本轮结束');
    assert((await page.locator('#race-result-note').textContent()).includes('尚未完赛'));
    assert.equal(report.errors.length,0);report.passed=true;
  }catch(e){report.errors.push(e.stack||String(e));if(page)report.diagnostics=await page.evaluate(()=>({snapshot:window.__molecularCarDebug?.getSnapshot(),text:document.body.innerText}));process.exitCode=1;}
  finally{if(browser)await browser.close();server.kill();fs.mkdirSync(out,{recursive:true});fs.writeFileSync(path.join(out,'report.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report));}
})();
