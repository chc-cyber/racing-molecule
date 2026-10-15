const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),{spawn}=require('node:child_process'),{createHash}=require('node:crypto');
const {chromium}=require('C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root=path.resolve(__dirname,'..'),mode=process.argv[2]||'wasm',out=path.join(root,'outputs/playback-check');
const serverArgs=mode==='wasm'?['wasm/tools/serve.py','--port','0']:['-c',
 "from molecular_car.server import LocalHTTPServer, Handler, GameController; from functools import partial; g=GameController(); s=LocalHTTPServer(('127.0.0.1',0),partial(Handler,game=g)); print(f'http://127.0.0.1:{s.server_port}/',flush=True); s.serve_forever(); g.close(); s.server_close()"];
(async()=>{
 const server=spawn(path.join(root,'.venv/Scripts/python.exe'),serverArgs,{cwd:root,windowsHide:true});
 let browser,page,base;const report={mode,passed:false,errors:[],runs:[]};let stderr='';server.stderr.on('data',d=>stderr+=String(d));
 try{
  base=await new Promise((resolve,reject)=>{const t=setTimeout(()=>reject(Error('server timeout '+stderr)),15000);server.once('error',reject);server.stdout.on('data',d=>{const m=String(d).match(/http:\/\/127\.0\.0\.1:\d+\//);if(m){clearTimeout(t);resolve(m[0]);}});});
  browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-unsafe-swiftshader']});
  page=await browser.newPage({viewport:{width:1440,height:1000},acceptDownloads:true});page.setDefaultTimeout(160000);
  page.on('pageerror',e=>report.errors.push(e.message));
  await page.goto(base+(mode==='wasm'?'race.html':''));
  const snapshot=()=>page.evaluate(()=>window.__molecularCarDebug.getSnapshot());
  await page.waitForFunction(()=>window.__molecularCarDebug?.getSnapshot().status==='idle');
  await page.locator('#example-select').selectOption('benzene.mol');await page.locator('#load-example-button').click();await page.locator('#prepare-race').click();
  await page.waitForFunction(()=>window.__molecularCarDebug.getSnapshot().status==='paused'&&window.__molecularCarDebug.getSnapshot().renderedStep===0);
  const representation=await page.locator('#space-filling-model').boundingBox(),speed=await page.locator('#playback-1x').boundingBox();
  assert(speed.y>=representation.y+representation.height);
  for(const rate of [1,2]){
   if(rate===2){await page.locator('#result-again').click();await page.waitForFunction(()=>window.__molecularCarDebug.getSnapshot().status==='paused'&&window.__molecularCarDebug.getSnapshot().renderedStep===0);}
   await page.locator(`#playback-${rate}x`).click();
   const started=performance.now();await page.locator('#start-button').click();
   await page.waitForFunction(()=>window.__molecularCarDebug.getSnapshot().playbackTimeFs>1000);
   const a=await snapshot(),wall=performance.now();await page.waitForTimeout(2000);const b=await snapshot();
   const slope=(b.playbackTimeFs-a.playbackTimeFs)/((performance.now()-wall)/1000);
   assert(Math.abs(slope-2000*rate)<600,`rate ${rate}: ${slope}`);
   assert(b.computedTimeFs>b.playbackTimeFs+500,'MD should compute ahead');assert(b.bufferedFrames<=205,'bounded unread buffer');
   assert.equal(await page.locator('#race-result-dialog').isVisible(),false);
   let pausedSeconds=0;
   if(rate===1){
    const pauseStart=performance.now();await page.locator('#start-button').click();const paused=await snapshot();
    await page.waitForTimeout(600);assert.equal((await snapshot()).playbackTimeFs,paused.playbackTimeFs);
    await page.locator('#playback-2x').click();assert.equal((await snapshot()).playbackTimeFs,paused.playbackTimeFs);
    await page.locator('#playback-1x').click();await page.locator('#start-button').click();pausedSeconds+=(performance.now()-pauseStart)/1000;
   }
   await page.waitForFunction(()=>window.__molecularCarDebug.getSnapshot().status==='finished');
   const computed=await snapshot();assert.equal(await page.locator('#race-result-dialog').isVisible(),false,'must wait for displayed finish');
   assert(computed.playbackTimeFs<computed.computedTimeFs);
   if(rate===1){
    const pauseStart=performance.now();await page.locator('#start-button').click();const paused=await snapshot();
    await page.waitForTimeout(600);assert.equal((await snapshot()).playbackTimeFs,paused.playbackTimeFs);
    assert.equal(await page.locator('#race-result-dialog').isVisible(),false);
    await page.locator('#start-button').click();pausedSeconds+=(performance.now()-pauseStart)/1000;
   }
   await page.locator('#race-result-dialog').waitFor({state:'visible'});
   const elapsed=(performance.now()-started)/1000,final=await snapshot();
   const download=page.waitForEvent('download');await page.locator('#result-save').click();const file=await download;
   const replay=JSON.parse(fs.readFileSync(await file.path(),'utf8'));
   for(let i=1;i<replay.frames.length;i++)assert.equal(replay.frames[i].seq,replay.frames[i-1].seq+1);
   const last=replay.frames.at(-1),expected=last.time_fs/(2000*rate);
   assert(Math.abs(elapsed-pausedSeconds-expected)<3,`playback timing: ${elapsed-pausedSeconds} vs ${expected}`);
   assert.equal(last.step,last.time_fs,'MD timestep remains 1 fs');assert(last.distance_angstrom>=200);
   assert.equal(final.renderedSeq,last.seq);assert.equal(final.playbackStatus,'finished');
   const avg=(last.distance_angstrom/last.time_fs*100000).toFixed(1);assert.equal(await page.locator('#race-result-speed').textContent(),avg);
   const hash=createHash('sha256').update(JSON.stringify(replay.frames)).digest('hex');
   report.runs.push({rate,wall_seconds:elapsed,paused_seconds:pausedSeconds,expected_playback_seconds:expected,slope_fs_per_second:slope,time_ps:last.time_fs/1000,distance_nm:last.distance_angstrom/10,average_speed_ms:Number(avg),frame_count:replay.frames.length,trajectory_sha256:hash,buffered:b.bufferedFrames});
   fs.mkdirSync(out,{recursive:true});await page.screenshot({path:path.join(out,`${mode}-${rate}x.png`)});
  }
  assert.equal(report.runs[0].trajectory_sha256,report.runs[1].trajectory_sha256,'rates must produce identical complete MD samples');
  assert.deepEqual(report.errors,[]);report.passed=true;
 }catch(error){report.error=error.stack;if(page)report.diagnostics=await page.evaluate(()=>({snapshot:window.__molecularCarDebug?.getSnapshot(),text:document.body.innerText})).catch(()=>null);report.stderr=stderr;process.exitCode=1;}
 finally{
  if(browser)await browser.close();
  if(mode==='native'&&base)await fetch(base+'api/shutdown',{method:'POST',headers:{'Content-Type':'application/json'},body:'{}'}).catch(()=>{});
  if(mode==='wasm')server.kill();else setTimeout(()=>server.kill(),2000).unref();
  fs.mkdirSync(out,{recursive:true});fs.writeFileSync(path.join(out,`${mode}.json`),JSON.stringify(report,null,2));console.log(JSON.stringify(report));
 }
})();
