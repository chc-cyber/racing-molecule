/* Real production HTTP handler + real workshop preparation/commit + Chrome. */
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const {randomUUID}=require('node:crypto'),{spawn}=require('node:child_process');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root=path.resolve(__dirname,'..'),out=path.join(root,'outputs/native-module-handoff');
(async()=>{
  const code=[
    'from molecular_car.server import LocalHTTPServer, Handler, GameController',
    'from functools import partial',
    'game=GameController()',
    "server=LocalHTTPServer(('127.0.0.1',0),partial(Handler,game=game))",
    "print(f'http://127.0.0.1:{server.server_port}/',flush=True)",
    'try: server.serve_forever()',
    'finally: server.server_close(); game.close()'
  ].join('\n');
  const server=spawn(path.join(root,'.venv/Scripts/python.exe'),['-c',code],{cwd:root,windowsHide:true});
  server.stderr.on('data',()=>{});
  let base,browser;const report={passed:false,pageErrors:[],consoleErrors:[]};
  try{
    base=await new Promise((resolve,reject)=>{const timer=setTimeout(()=>reject(Error('server timeout')),15000);server.once('error',reject);server.stdout.on('data',data=>{const m=String(data).match(/http:\/\/127\.0\.0\.1:\d+\//);if(m){clearTimeout(timer);resolve(m[0]);}});});
    browser=await chromium.launch({executablePath:process.env.CHROME_PATH||'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-unsafe-swiftshader']});
    const page=await browser.newPage({viewport:{width:1440,height:1000}});
    page.on('pageerror',e=>report.pageErrors.push(e.message));
    page.on('console',m=>{if(m.type()==='error')report.consoleErrors.push(m.text());});
    const api=async(url,data)=>{const response=data===undefined?await page.request.get(base+url):await page.request.post(base+url,{data});assert(response.ok(),await response.text());return response.json();};
    const module=await page.request.get(base+'race-course.mjs');
    report.module_content_type=module.headers()['content-type'];
    assert.match(report.module_content_type,/^(text|application)\/javascript/);
    const id=randomUUID(),started=performance.now();
    await api('api/editor/race/prepare',{id,filename:'tetraphenylbenzene.mol',content:fs.readFileSync(path.join(root,'examples/tetraphenylbenzene.mol'),'utf8')});
    let job;
    for(let attempt=0;attempt<300;attempt++){
      job=await api('api/editor/race/status?id='+id);
      if(['ready','error'].includes(job.status))break;
      await new Promise(resolve=>setTimeout(resolve,200));
    }
    assert.equal(job.status,'ready',JSON.stringify(job));
    report.preparation_seconds=(performance.now()-started)/1000;
    report.commit=await api('api/editor/race/commit',{id});assert.equal(report.commit.status,'committed');
    await page.goto(base);
    await page.waitForFunction(()=>{const s=window.__molecularCarDebug?.getSnapshot();return s?.connected&&s.status==='paused'&&s.atomCount===52&&s.positions?.length===156;},null,{timeout:30000});
    const snap=await page.evaluate(()=>window.__molecularCarDebug.getSnapshot());
    const state=await api('api/state');
    assert(snap.renderer.triangles>0);assert.equal(snap.finishX,state.runtime.race.finish_x);
    assert.equal(await page.locator('#distance-value').textContent(),'0.00');
    await page.locator('#start-button').click();
    await page.waitForFunction(()=>window.__molecularCarDebug.getSnapshot().renderedStep>=50);
    await api('api/action',{action:'pause'});
    report.render={atomCount:snap.atomCount,triangles:snap.renderer.triangles,finishX:snap.finishX};
    report.running_step=(await api('api/state')).frame.step;
    assert.equal(report.pageErrors.length,0);assert.equal(report.consoleErrors.length,0);report.passed=true;
    fs.mkdirSync(out,{recursive:true});await page.screenshot({path:path.join(out,'race.png')});
  }catch(e){report.error=e.stack||String(e);process.exitCode=1;}
  finally{
    if(base)try{await fetch(base+'api/shutdown',{method:'POST',headers:{'Content-Type':'application/json'},body:'{}'});}catch{}
    if(browser)await browser.close();
    if(server.exitCode===null)await new Promise(resolve=>{const timer=setTimeout(()=>{server.kill();resolve();},4000);server.once('exit',()=>{clearTimeout(timer);resolve();});});
    fs.mkdirSync(out,{recursive:true});fs.writeFileSync(path.join(out,'report.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report));
  }
})();
