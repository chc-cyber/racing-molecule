const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),{spawn}=require('node:child_process');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const root=path.resolve(__dirname,'..'),out=path.join(root,'outputs/web-release-acceptance/recovery');fs.mkdirSync(out,{recursive:true});
(async()=>{
 const server=spawn(path.join(root,'.venv/Scripts/python.exe'),['wasm/tools/serve.py','--port','0'],{cwd:root,windowsHide:true});server.stderr.on('data',()=>{});
 let browser,page,garage,race,base;const report={checks:[],errors:[]};
 const record=(name,data)=>{report.checks.push({name,...data});fs.writeFileSync(path.join(out,'report.json'),JSON.stringify(report,null,2));console.log(name,JSON.stringify(data));};
 const shot=name=>page.screenshot({path:path.join(out,name+'.png'),fullPage:true});
 try{
  base=await new Promise((resolve,reject)=>{const t=setTimeout(()=>reject(Error('server timeout')),10000);server.stdout.on('data',d=>{const m=String(d).match(/http:\/\/127\.0\.0\.1:\d+\//);if(m){clearTimeout(t);resolve(m[0]);}});});
  browser=await chromium.launch({executablePath:process.env.CHROME_PATH||'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-unsafe-swiftshader']});
  const context=await browser.newContext({viewport:{width:1440,height:1000}});page=await context.newPage();page.setDefaultTimeout(40000);page.on('pageerror',e=>report.errors.push(e.message));
  await page.goto(base+'#garage');garage=await page.locator('#garage-frame').elementHandle().then(e=>e.contentFrame());
  await garage.waitForFunction(()=>window.workshop?.mature?.ready&&workshop.views?.ready);
  await garage.locator('#open-file').setInputFiles(path.join(root,'outputs/web-release-acceptance/chrome/drawn-car.mcar'));
  await garage.waitForFunction(()=>workshop.model.project.atoms.length===12);
  let failedAsset=false;await page.route('**/fourier-atlas-0.bin.gz',route=>{if(!failedAsset){failedAsset=true;return route.abort('internetdisconnected');}return route.continue();});
  await garage.locator('#send-track').click();await garage.locator('#race-dialog').waitFor({state:'visible'});await garage.locator('#race-dialog').waitFor({state:'hidden',timeout:100000});
  const error=await garage.locator('#feedback').innerText();assert(page.url().endsWith('#garage'));assert(failedAsset);assert(/fetch|失败|错误/i.test(error));
  await shot('01-download-failure');record('势表下载失败显示可读提示并可返回编辑',{passed:true,feedback:error});
  await garage.locator('#send-track').click();await page.waitForURL('**/#track',{timeout:100000});race=await page.locator('#track-frame').elementHandle().then(e=>e.contentFrame());
  await race.waitForFunction(()=>__molecularCarDebug.getSnapshot().status==='paused');record('下载恢复后重新准备成功',{passed:true});
  await race.locator('#start-button').click();await race.waitForFunction(()=>__molecularCarDebug.getSnapshot().playbackTimeFs>2000);
  const before=await race.evaluate(()=>__molecularCarDebug.getSnapshot());let prompts=0;page.on('dialog',async d=>{prompts++;await d.accept();});
  await page.reload();race=await page.locator('#track-frame').elementHandle().then(e=>e.contentFrame());await race.waitForFunction(()=>__molecularCarDebug.getSnapshot().status==='idle');
  const after=await race.evaluate(()=>__molecularCarDebug.getSnapshot());assert.equal(prompts,1);record('比赛中刷新先出现离页提醒，确认离开才清空比赛',{passed:true,beforeStep:before.renderedStep,afterStatus:after.status,prompts});await shot('02-refresh-loses-race');
  await race.locator('a[data-app-route="garage"]').click();await page.waitForURL('**/#garage');garage=await page.locator('#garage-frame').elementHandle().then(e=>e.contentFrame());await garage.waitForFunction(()=>window.workshop?.mature?.ready&&workshop.model.project.atoms.length===12);record('比赛刷新丢失后工坊草稿仍保留',{passed:true});
  // Fault injection only: model browser storage quota failure, then restore it.
  await garage.evaluate(()=>{window.savedSetItem=Storage.prototype.setItem;Storage.prototype.setItem=function(){throw new DOMException('Storage quota exceeded','QuotaExceededError');};});
  await garage.locator('#send-track').click();await page.waitForURL('**/#track',{timeout:100000});
  const storageFeedback=await garage.locator('#feedback').innerText();
  assert(storageFeedback.includes('暂存失败'));record('暂存失败仍能进入比赛，并提供工程下载提示',{passed:true,faultInjection:'Storage.setItem throws QuotaExceededError',feedback:storageFeedback});await shot('03-storage-failure');
  race=await page.locator('#track-frame').elementHandle().then(e=>e.contentFrame());
  await race.locator('a[data-app-route="garage"]').click();await page.waitForURL('**/#garage');
  await garage.evaluate(()=>{Storage.prototype.setItem=window.savedSetItem;delete window.savedSetItem;});
  await garage.locator('#file-menu > summary').click();await garage.locator('#new-project').click();await garage.locator('#send-track').click();await garage.locator('#race-dialog').waitFor({state:'hidden'});const empty=await garage.locator('#feedback').innerText();assert(page.url().endsWith('#garage'));record('空结构送入赛道返回明确错误',{passed:true,feedback:empty});
  await garage.locator('#open-file').setInputFiles({name:'broken.mol',mimeType:'text/plain',buffer:Buffer.from('not a mol file')});await page.waitForTimeout(1000);const bad=await garage.locator('#feedback').innerText();assert.equal(await garage.evaluate(()=>workshop.model.project.atoms.length),0);record('损坏 MOL 导入被拒绝且保持工程',{passed:true,feedback:bad});
  assert.deepEqual(report.errors,[]);report.completed=true;
 }catch(error){report.fatal=error.stack;console.log(error.stack);await shot('failure').catch(()=>{});process.exitCode=1;}
 finally{report.passed=Boolean(report.completed)&&report.checks.every(c=>c.passed)&&!report.errors.length;if(!report.passed)process.exitCode=1;fs.writeFileSync(path.join(out,'report.json'),JSON.stringify(report,null,2));if(browser)await browser.close();server.kill();}
})();

