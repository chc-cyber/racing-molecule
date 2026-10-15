/* Independent production HTTP + static subpath site; real MD, no API mocks. */
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),{spawn}=require('node:child_process');
const {chromium}=require('C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root=path.resolve(__dirname,'..'),out=path.join(root,'outputs/competition-acceptance');
const sleep=ms=>new Promise(r=>setTimeout(r,ms));
async function serve(wasm){
  const code=wasm?[
    'from wasm.tools.serve import StaticHandler, ROOT','from http.server import ThreadingHTTPServer','from functools import partial',
    "server=ThreadingHTTPServer(('127.0.0.1',0),partial(StaticHandler,directory=str(ROOT/'wasm')))",
    "print(f'http://127.0.0.1:{server.server_port}/dist/',flush=True)",'server.serve_forever()'
  ]:[
    'from molecular_car.server import LocalHTTPServer,Handler,GameController','from functools import partial','game=GameController()',
    "server=LocalHTTPServer(('127.0.0.1',0),partial(Handler,game=game))","print(f'http://127.0.0.1:{server.server_port}/',flush=True)",
    'try: server.serve_forever()','finally: server.server_close();game.close()'];
  const process=spawn(path.join(root,'.venv/Scripts/python.exe'),['-X','utf8','-c',code.join('\n')],{cwd:root,windowsHide:true});
  let stderr='';process.stderr.on('data',d=>stderr+=d);
  const url=await new Promise((resolve,reject)=>{const timeout=setTimeout(()=>reject(Error(stderr||'Server timeout')),15000);process.once('error',reject);process.stdout.on('data',d=>{const m=String(d).match(/http:\/\/127\.0\.0\.1:\d+\/(?:dist\/)?/);if(m){clearTimeout(timeout);resolve(m[0]);}});});
  return {process,url};
}
(async()=>{
  fs.mkdirSync(out,{recursive:true});let browser;const report={passed:false,checks:[],pageErrors:[]};
  try{
    browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-unsafe-swiftshader']});
    for(const wasm of [false,true]){
      const server=await serve(wasm),name=wasm?'wasm':'native';let page;
      try{
        page=await browser.newPage({viewport:{width:1440,height:1000},acceptDownloads:true});page.setDefaultTimeout(90000);
        page.on('pageerror',e=>report.pageErrors.push({engine:name,error:e.message}));
        const requests=[];page.on('request',r=>{if(new URL(r.url()).pathname.startsWith('/api/'))requests.push(r.url());});
        const snapshot=()=>page.evaluate(()=>window.__molecularCarDebug.getSnapshot());
        const ready=()=>page.waitForFunction(()=>{const s=window.__molecularCarDebug?.getSnapshot();return s?.status==='paused'&&s.positions?.length>0;});
        const running=async()=>{await page.locator('#start-button').click();await page.waitForFunction(()=>window.__molecularCarDebug.getSnapshot().renderedStep>=100);assert(await page.locator('#race-setup [data-key="temperature"]').isDisabled());await page.locator('#start-button').click();await ready();};
        await page.goto(server.url);await page.waitForFunction(()=>window.__molecularCarDebug?.getSnapshot().status==='idle');
        assert.equal(await page.locator('#tip-drive,#pulse-button,#move-tip-button').count(),0);
        await page.locator('#example-select').selectOption('benzene.mol');await page.locator('#load-example-button').click();
        const started=performance.now();await page.locator('#prepare-race').click();await ready();
        let s=await snapshot();assert.equal(s.competition.mode,'standard');assert.equal(s.competition.temperature,150);assert.equal(s.atomCount,12);assert.equal(s.competition.timestep_fs,1);
        const standardSeconds=(performance.now()-started)/1000;
        await running();
        // Changes are drafts until explicit preparation; molecule stays the same.
        await page.locator('#race-setup-panel > summary').click();
        await page.locator('#race-setup [data-mode="custom"]').click();
        await page.locator('#race-setup [data-key="temperature"]').fill('300');
        await page.locator('#race-setup [data-key="strength"]').fill('0.35');
        await page.locator('#race-setup [data-key="surface"]').selectOption('continuum');
        await page.locator('#race-setup [data-key="corrugation"]').fill('0.12');
        await page.locator('#race-setup [data-key="growth"]').fill('0.5');
        await page.screenshot({path:path.join(out,name+'-settings.png')});
        assert.equal((await snapshot()).competition.mode,'standard');
        await page.locator('#prepare-race').click();await ready();
        s=await snapshot();assert.equal(s.competition.mode,'custom');assert.equal(s.competition.temperature,300);assert.equal(s.drive.strength,.35);assert.equal(s.drive.growth,.5);assert.equal(s.surface,'continuum');assert.equal(s.atomCount,12);
        await running();
        const download=page.waitForEvent('download');await page.locator('#replay-menu summary').click();await page.locator('#download-button').click();
        const file=await download;await file.saveAs(path.join(out,name+'-replay.json'));const replay=JSON.parse(fs.readFileSync(await file.path(),'utf8'));
        assert.equal(replay.competition.mode,'custom');assert.equal(replay.competition.random_seed,98245165);
        assert(replay.frames.length>1);assert(!('tip' in replay.frames[0]));assert(!('pulse_active' in replay.frames[0]));
        assert.equal(replay.runtime.langevin_seed,98245165);
        await page.locator('#reset-button').click();await ready();assert.equal((await snapshot()).competition.temperature,300);
        await page.screenshot({path:path.join(out,name+'-custom.png')});
        await page.setViewportSize({width:390,height:844});await page.screenshot({path:path.join(out,name+'-mobile.png'),fullPage:true});
        assert(await page.evaluate(()=>document.documentElement.scrollWidth<=391));
        if(wasm)assert.equal(requests.length,0,'Static browser game must make no API calls');
        report.checks.push({engine:name,standard_preparation_seconds:standardSeconds,standard:true,custom:true,reset:true,replay:true,mobile:true,api_calls:requests.length});
        if(!wasm){
          await page.setViewportSize({width:1440,height:1000});await page.goto(server.url+'editor/index.html');
          await page.waitForFunction(()=>window.workshop?.mature?.ready&&workshop.views?.ready);
          await page.evaluate(async content=>{await workshop.importContent(content,'benzene.mol');await workshop.mature.settled();await workshop.views.settled();},fs.readFileSync(path.join(root,'examples/benzene.mol'),'utf8'));
          const prepBefore=requests.filter(u=>u.endsWith('/api/editor/race/prepare')).length;
          await page.locator('#send-track').click();await page.locator('#competition-dialog').waitFor({state:'visible'});
          await page.locator('#competition-setup [data-mode="custom"]').click();
          await page.locator('#competition-setup [data-key="temperature"]').fill('250');
          await page.locator('#competition-setup [data-key="surface"]').selectOption('continuum');
          await page.screenshot({path:path.join(out,'workshop-rules.png')});
          await page.locator('#competition-prepare').click();await page.waitForURL(server.url+'?workshop=*');await ready();
          assert.equal((await snapshot()).competition.temperature,250);
          assert.equal(requests.filter(u=>u.endsWith('/api/editor/race/prepare')).length-prepBefore,1);
          await running();report.checks.push({workshop:true,chosen_temperature:250,preparations:1});
        }
      }finally{
        if(page)await page.close();
        if(!wasm)try{await fetch(server.url+'api/shutdown',{method:'POST',headers:{'Content-Type':'application/json'},body:'{}'});}catch{}
        if(server.process.exitCode===null){if(wasm)server.process.kill();else{await sleep(1500);if(server.process.exitCode===null)server.process.kill();}}
      }
    }
    assert.deepEqual(report.pageErrors,[]);report.passed=true;
  }catch(error){report.error=error.stack;process.exitCode=1;}
  finally{if(browser)await browser.close();fs.writeFileSync(path.join(out,'report.json'),JSON.stringify(report,null,2));console.log(JSON.stringify(report));}
})();
