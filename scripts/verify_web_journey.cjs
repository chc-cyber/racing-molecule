const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),{spawn}=require('node:child_process');
const {chromium,firefox}=require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const root=path.resolve(__dirname,'..'),browserName=process.argv[2]||'chrome',out=path.join(root,'outputs/web-release-acceptance',browserName);
fs.mkdirSync(out,{recursive:true});
const serverCode=`from http.server import ThreadingHTTPServer
from functools import partial
from wasm.tools.serve import StaticHandler, ROOT
class PrefixHandler(StaticHandler):
 def do_GET(self):
  if not self.path.startswith('/molecular-car/'):
   self.send_error(404); return
  self.path=self.path[len('/molecular-car'):]
  super().do_GET()
s=ThreadingHTTPServer(('127.0.0.1',0),partial(PrefixHandler,directory=str(ROOT/'wasm'/'dist')))
print(f'http://127.0.0.1:{s.server_port}/molecular-car/',flush=True)
s.serve_forever()`;
(async()=>{
 const server=spawn(path.join(root,'.venv/Scripts/python.exe'),['-c',serverCode],{cwd:root,windowsHide:true});server.stderr.on('data',()=>{});
 let browser,page,garage,race,base;const report={browser:browserName,started:new Date().toISOString(),checks:[],errors:[],httpErrors:[],external:[],apiRequests:[],downloads:[],assets:[]};
 const persist=()=>fs.writeFileSync(path.join(out,'report.json'),JSON.stringify(report,null,2));
 const pass=(name,data={})=>{report.checks.push({name,passed:true,...data});persist();console.log('PASS',name,JSON.stringify(data));};
 const check=async(name,fn)=>{try{pass(name,await fn()||{});}catch(e){report.checks.push({name,passed:false,error:e.message});persist();console.log('FAIL',name,e.message);}};
 const shot=name=>page.screenshot({path:path.join(out,name+'.png'),fullPage:true});
 try{
  base=await new Promise((resolve,reject)=>{const t=setTimeout(()=>reject(Error('server startup timeout')),10000);server.stdout.on('data',d=>{const m=String(d).match(/http:\/\/127\.0\.0\.1:\d+\/molecular-car\//);if(m){clearTimeout(t);resolve(m[0]);}});server.once('error',reject);});
  browser=browserName==='firefox'?await firefox.launch({headless:true}):await chromium.launch({executablePath:browserName==='edge'?(process.env.EDGE_PATH||'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe'):(process.env.CHROME_PATH||'C:/Program Files/Google/Chrome/Application/chrome.exe'),headless:true,args:['--enable-unsafe-swiftshader']});
  report.browserVersion=browser.version();const context=await browser.newContext({viewport:{width:1440,height:1000},acceptDownloads:true});
  page=await context.newPage();page.setDefaultTimeout(30000);
  page.on('pageerror',e=>{report.errors.push(e.message);persist();});
  page.on('request',r=>{if(new URL(r.url()).pathname.startsWith('/molecular-car/api/'))report.apiRequests.push(r.url());if(!r.url().startsWith(base)&&!r.url().startsWith('blob:')&&!r.url().startsWith('data:'))report.external.push(r.url());});
  page.on('response',async r=>{if(r.status()>=400)report.httpErrors.push({url:r.url(),status:r.status()});if(/wasm$|bin.gz$/.test(r.url()))report.assets.push({url:r.url().replace(base,''),bytes:Number(r.headers()['content-length']||0)});});
  page.on('download',d=>report.downloads.push(d.suggestedFilename()));
  await page.goto(base);await page.getByRole('link',{name:'加入赛事'}).waitFor();await shot('01-home');pass('首页与子目录静态入口');
  const opened=performance.now();await page.getByRole('link',{name:'加入赛事'}).click();
  await page.locator('#garage-frame').waitFor();garage=await page.locator('#garage-frame').elementHandle().then(e=>e.contentFrame());
  const settled=async()=>{await garage.waitForFunction(()=>window.workshop?.mature?.ready&&workshop.views?.ready&&!workshop.mature.suspended&&!workshop.views.suspended&&!workshop.mature.gesture&&!workshop.views.gesture,null,{timeout:90000});await garage.evaluate(async()=>{await workshop.mature.settled();await workshop.views.settled();});};
  await settled();pass('工坊二维与三维编辑器就绪',{seconds:(performance.now()-opened)/1000});
  const ketcher=garage.frameLocator('#ketcher-frame');
  await ketcher.locator('[data-testid="template-0"]').filter({visible:true}).click();
  const canvas=await ketcher.locator('[data-testid=canvas]').boundingBox();await page.mouse.click(canvas.x+canvas.width*.46,canvas.y+canvas.height*.45);
  await ketcher.locator('[data-testid=select-rectangle]').filter({visible:true}).last().click();await page.mouse.move(1435,985);await settled();
  let project=await garage.evaluate(()=>structuredClone(workshop.model.project));assert.equal(project.atoms.length,6);assert.equal(project.bonds.length,6);
  await garage.locator('#project-name').fill('发布验收 · 苯环赛车');await garage.locator('#project-name').press('Tab');await shot('02-drawn-2d');pass('通过二维模板工具真实绘制苯环',{atoms:6,bonds:6});
  await garage.locator('#undo').click();await settled(); // Undo the name edit, not the ring.
  await garage.locator('#redo').click();await settled();assert.equal((await garage.evaluate(()=>workshop.model.project.name)),'发布验收 · 苯环赛车');pass('工坊撤销与重做');
  const threeStart=performance.now();await garage.locator('#view-3d').click();
  await garage.waitForFunction(()=>workshop.dimension===3&&workshop.pendingGeometry===0,null,{timeout:100000});await settled();
  project=await garage.evaluate(()=>structuredClone(workshop.model.project));assert.equal(project.atoms.length,12);assert(project.has3d);assert(project.atoms.every(a=>a.pos3.every(Number.isFinite)));
  await shot('03-generated-3d');pass('实际 WASM 补氢与三维构象生成',{seconds:(performance.now()-threeStart)/1000,atoms:12});
  let dl=page.waitForEvent('download');await garage.locator('#save-project').click();let file=await dl;await file.saveAs(path.join(out,'drawn-car.mcar'));assert.equal(JSON.parse(fs.readFileSync(path.join(out,'drawn-car.mcar'),'utf8')).atoms.length,12);pass('工程文件下载');
  const prepareStart=performance.now();await garage.locator('#send-track').click();await page.waitForURL('**/#track',{timeout:170000});
  race=await page.locator('#track-frame').elementHandle().then(e=>e.contentFrame());
  await race.waitForFunction(()=>window.__molecularCarDebug?.getSnapshot().status==='paused'&&window.__molecularCarDebug.getSnapshot().renderedStep===0,null,{timeout:90000});
  let snap=await race.evaluate(()=>__molecularCarDebug.getSnapshot());assert.equal(snap.atomCount,12);assert.equal(snap.competition.mode,'standard');assert.equal(snap.playbackRate,1);
  await shot('04-race-ready');pass('工坊准备完成并交接同一分子到比赛',{seconds:(performance.now()-prepareStart)/1000,generation:snap.generation,atoms:snap.atomCount,competition:snap.competition});
  const started=performance.now();await race.locator('#start-button').click();await race.waitForFunction(()=>__molecularCarDebug.getSnapshot().playbackTimeFs>2000);
  const timeA=await race.evaluate(()=>__molecularCarDebug.getSnapshot().playbackTimeFs),wallA=performance.now();await page.waitForTimeout(1000);const timeB=await race.evaluate(()=>__molecularCarDebug.getSnapshot().playbackTimeFs);
  const slope=(timeB-timeA)/((performance.now()-wallA)/1000);assert(Math.abs(slope-2000)<500);pass('默认 1× 按模拟时间播放',{fs_per_second:slope});
  await race.locator('#playback-2x').click();assert.equal(await race.locator('#playback-2x').getAttribute('aria-pressed'),'true');
  const timeC=await race.evaluate(()=>__molecularCarDebug.getSnapshot().playbackTimeFs),wallC=performance.now();await page.waitForTimeout(1000);const timeD=await race.evaluate(()=>__molecularCarDebug.getSnapshot().playbackTimeFs);
  const slope2=(timeD-timeC)/((performance.now()-wallC)/1000);assert(Math.abs(slope2-4000)<700);pass('比赛中切换 2×',{fs_per_second:slope2});
  await race.locator('a[data-app-route="garage"]').click();await page.waitForURL('**/#garage');await page.waitForTimeout(500);
  snap=await race.evaluate(()=>__molecularCarDebug.getSnapshot());const frozen=snap.playbackTimeFs;await page.waitForTimeout(600);assert.equal((await race.evaluate(()=>__molecularCarDebug.getSnapshot())).playbackTimeFs,frozen);
  assert.equal(await garage.evaluate(()=>workshop.model.project.atoms.length),12);
  await garage.locator('a[data-app-route="track"]').click();await page.waitForURL('**/#track');assert.equal((await race.evaluate(()=>__molecularCarDebug.getSnapshot())).playbackTimeFs,frozen);await race.locator('#start-button').click();pass('返回工坊暂停、保留设计、返回比赛继续');
  await race.locator('#race-result-dialog').waitFor({state:'visible',timeout:180000});
  const final=await race.evaluate(async()=>{const {request}=await import('./transport.js');return request('/api/state');});
  assert(final.frame.distance_angstrom>=200);const expected=(final.frame.distance_angstrom/final.frame.time_fs*100000).toFixed(1);
  assert.equal(await race.locator('#race-result-speed').textContent(),expected);assert.equal(await race.locator('#race-result-name').textContent(),'发布验收 · 苯环赛车');
  await shot('05-result-desktop');pass('从工坊绘制到真实完赛与平均速度结果',{wall_seconds:(performance.now()-started)/1000,time_ps:final.frame.time_fs/1000,distance_nm:final.frame.distance_angstrom/10,average_speed_ms:Number(expected)});
  dl=page.waitForEvent('download');await race.locator('#result-save').click();file=await dl;await file.saveAs(path.join(out,'replay.json'));
  const replay=JSON.parse(fs.readFileSync(path.join(out,'replay.json'),'utf8'));assert.equal(replay.frames.at(-1).time_fs,final.frame.time_fs);assert.equal(replay.settings.timestep_fs,1);assert(replay.replay_complete);pass('成绩页导出完整原始轨迹',{frames:replay.frames.length,timestep_fs:replay.settings.timestep_fs});
  await page.setViewportSize({width:390,height:844});await page.waitForTimeout(350);await shot('06-result-mobile');
  await check('手机尺寸成绩弹窗完整可见',async()=>{const b=await race.locator('#race-result-dialog').boundingBox();assert(b.x>=0&&b.x+b.width<=391&&b.y>=0&&b.y+b.height<=845);return{bounds:b};});
  await race.locator('.result-close').click();await page.waitForTimeout(200);await shot('07-track-mobile');
  await check('手机尺寸可操作模型与倍速按钮',async()=>{await race.locator('#playback-1x').scrollIntoViewIfNeeded();await race.locator('#playback-1x').click();assert.equal(await race.locator('#playback-1x').getAttribute('aria-pressed'),'true');const b=await race.locator('#playback-1x').boundingBox();assert(b.x>=0&&b.x+b.width<=391);return{bounds:b};});
  await race.locator('a[data-app-route="garage"]').click();await page.waitForURL('**/#garage');await garage.locator('#view-2d').click();await settled();await page.waitForTimeout(400);await shot('08-workshop-mobile');
  await check('手机尺寸工坊主操作无横向溢出',async()=>{const size=await garage.evaluate(()=>({width:innerWidth,scroll:document.documentElement.scrollWidth}));assert(size.scroll<=size.width+1,JSON.stringify(size));const b=await garage.locator('#send-track').boundingBox();assert(b.x>=0&&b.x+b.width<=391,JSON.stringify(b));return{size,sendBounds:b};});
  await check('缩窄视口并返回二维后分子仍在可见画布内',async()=>{
    const points=await garage.evaluate(()=>{
      const h=workshop.mature,e=h.ketcher.editor,r=e.render,s=e.struct(),q=r.page2obj({clientX:0,clientY:0}),x=r.page2obj({clientX:100,clientY:0}),y=r.page2obj({clientX:0,clientY:100}),a=(x.x-q.x)/100,b=(y.x-q.x)/100,c=(x.y-q.y)/100,d=(y.y-q.y)/100,det=a*d-b*c,box=h.frame.getBoundingClientRect();
      return [...s.atoms.values()].map(atom=>{const dx=atom.pp.x-q.x,dy=atom.pp.y-q.y;return{x:box.left+(d*dx-b*dy)/det,y:box.top+(-c*dx+a*dy)/det};});
    });
    assert(points.some(p=>p.x>=0&&p.x<390&&p.y>=190&&p.y<780),JSON.stringify(points));return{points};
  });
  await garage.locator('#fit-view').click();await page.waitForTimeout(300);await shot('09-workshop-mobile-fit');
  await page.setViewportSize({width:1440,height:1000});await garage.locator('a[data-app-route="track"]').click();await page.waitForURL('**/#track');
  await race.locator('#show-race-result').click();await race.locator('#result-again').click();await race.waitForFunction(()=>__molecularCarDebug.getSnapshot().status==='paused'&&__molecularCarDebug.getSnapshot().renderedStep===0,null,{timeout:90000});pass('成绩页再跑一次正确归零');
  await race.locator('a[data-app-route="garage"]').click();await page.waitForURL('**/#garage');
  await garage.locator('#send-track').click();await garage.locator('#race-dialog').waitFor({state:'visible'});await garage.locator('#race-cancel').click();await page.waitForTimeout(1500);assert(page.url().endsWith('#garage'));assert.equal(await garage.locator('#race-dialog').isVisible(),false);assert.equal(await garage.locator('#send-track').isEnabled(),true);pass('取消准备后留在可编辑工坊');
  await garage.locator('#send-track').click();await page.waitForURL('**/#track',{timeout:170000});await race.waitForFunction(()=>__molecularCarDebug.getSnapshot().status==='paused',null,{timeout:90000});pass('取消后可以重新准备并进入赛道');
  await race.locator('#race-setup-panel > summary').click();
  await race.locator('#race-setup [data-mode="custom"]').click();await race.locator('#race-setup [data-key="temperature"]').fill('200');
  await race.locator('a[data-app-route="garage"]').click();await page.waitForURL('**/#garage');
  await garage.locator('#send-track').click();await page.waitForURL('**/#track',{timeout:170000});
  await race.waitForFunction(()=>__molecularCarDebug.getSnapshot().status==='paused',null,{timeout:90000});
  await check('自定义赛事参数沿工坊交接保留',async()=>{const c=await race.evaluate(()=>__molecularCarDebug.getSnapshot().competition);assert.equal(c.mode,'custom');assert.equal(c.temperature,200);return{competition:c};});
  await race.locator('a[data-app-route="garage"]').click();await page.waitForURL('**/#garage');
  const prior=await garage.evaluate(()=>JSON.stringify(workshop.model.project));await page.reload();
  garage=await page.locator('#garage-frame').elementHandle().then(e=>e.contentFrame());await settled();
  assert.equal(await garage.evaluate(()=>JSON.stringify(workshop.model.project)),prior);pass('刷新后恢复工坊草稿');
  await check('无页面异常、资源丢失或远端计算接口请求',async()=>{assert.deepEqual(report.errors,[]);assert.deepEqual(report.httpErrors,[]);assert.deepEqual(report.external,[]);assert.deepEqual(report.apiRequests,[]);return{asset_requests:report.assets};});
  report.completed=true;
 }catch(error){report.fatal=error.stack;console.log('FATAL',error.stack);if(page){await shot('failure').catch(()=>{});report.diagnostics={url:page.url(),shell:await page.locator('body').innerText().catch(()=>''),garage:await garage?.locator('body').innerText().catch(()=>''),race:await race?.locator('body').innerText().catch(()=>''),frames:page.frames().map(f=>f.url())};}process.exitCode=1;}
 finally{report.passed=Boolean(report.completed)&&report.checks.every(c=>c.passed)&&!report.errors.length;if(!report.passed)process.exitCode=1;report.finished=new Date().toISOString();persist();if(browser)await browser.close();server.kill();console.log(JSON.stringify({completed:report.completed,checks:report.checks.length,failed:report.checks.filter(c=>!c.passed),fatal:report.fatal,errors:report.errors,httpErrors:report.httpErrors}));}
})();

