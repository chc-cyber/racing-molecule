const fs=require('fs'),path=require('path'),assert=require('assert/strict'),{spawn}=require('child_process');const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root=path.resolve(__dirname,'..'),out=path.join(root,'outputs/preload'+(process.env.WEB_TEST_URL?'-public':''));fs.mkdirSync(out,{recursive:true});const report={checks:[],errors:[]};let browser,server;
(async()=>{let base=process.env.WEB_TEST_URL;if(!base){server=spawn(path.join(root,'.venv/Scripts/python.exe'),['-u','wasm/tools/serve.py','--port','0'],{cwd:root,windowsHide:true});server.stderr.on('data',()=>{});base=await new Promise(r=>server.stdout.on('data',d=>{const m=String(d).match(/http:\/\/127.0.0.1:\d+/);if(m)r(m[0]);}));}base=base.replace(/\/$/,'');report.url=base;browser=await chromium.launch({headless:true,executablePath:process.env.CHROME_PATH||'C:/Program Files/Google/Chrome/Application/chrome.exe'});let context=await browser.newContext();
const pass=(name,data={})=>{report.checks.push({name,...data});console.log('PASS',name,JSON.stringify(data));};
async function open(context,content){const page=await context.newPage();page.on('pageerror',e=>report.errors.push(e.message));page.on('dialog',d=>d.accept());await page.goto(base+'/#garage');await page.locator('#garage-frame').waitFor();const g=await page.locator('#garage-frame').elementHandle().then(e=>e.contentFrame());await g.waitForFunction(()=>window.workshop?.mature?.ready&&(workshop.dimension!==3||workshop.views?.ready),null,{timeout:120000});await g.evaluate(t=>workshop.importContent(t,'test.mol'),content);await g.evaluate(()=>{window.prepareMessages=[];new MutationObserver(()=>{const t=document.getElementById('race-message').textContent;if(prepareMessages.at(-1)!==t)prepareMessages.push(t);}).observe(document.getElementById('race-message'),{childList:true,subtree:true,characterData:true});});return{page,g};}
const benzene=fs.readFileSync(path.join(root,'examples/benzene.mol'),'utf8');
async function ready(item){const t=performance.now();await item.g.locator('#send-track').click();await Promise.race([item.page.waitForURL('**/#track',{timeout:610000}),item.g.waitForFunction(()=>!document.getElementById('race-dialog').open&&document.getElementById('feedback').classList.contains('error'),null,{timeout:610000}).then(async()=>{throw Error(await item.g.locator('#feedback').textContent());})]);const race=await item.page.locator('#track-frame').elementHandle().then(e=>e.contentFrame());await race.waitForFunction(()=>window.__molecularCarDebug?.getSnapshot().status==='paused',null,{timeout:30000});return{seconds:(performance.now()-t)/1000,messages:await item.g.evaluate(()=>prepareMessages)};}
let requests={};context.on('request',r=>{if((r.url().endsWith('.bin.gz')||r.url().includes('/assets/chunks/')))requests[r.url()]=(requests[r.url()]||0)+1;});
let item=await open(context,benzene);
await item.page.waitForFunction(()=>document.documentElement.dataset.racePreload==='loading',null,{timeout:120000});
assert.equal(await item.page.locator('#track-frame').count(),0);
await item.g.locator('#view-2d').click();const k=item.g.frameLocator('#ketcher-frame');
const t=performance.now();await k.locator('[data-testid="template-0"]').filter({visible:true}).click();
const b=await k.locator('[data-testid=canvas]').boundingBox();await item.page.mouse.click(b.x+b.width*.8,b.y+b.height*.3);
await k.locator('[data-testid=select-rectangle]').filter({visible:true}).last().click();
await item.g.waitForFunction(()=>workshop.model.project.atoms.length===18,null,{timeout:10000});
pass('draw benzene during background download without starting race',{seconds:(performance.now()-t)/1000});
await item.g.locator('#undo').click();await item.g.waitForFunction(()=>workshop.model.project.atoms.length===12);
await item.g.locator('#redo').click();await item.g.waitForFunction(()=>workshop.model.project.atoms.length===18);
await item.g.locator('#undo').click();await item.g.waitForFunction(()=>workshop.model.project.atoms.length===12);
pass('undo and redo remain responsive during preload');
const active=await ready(item);assert(Object.values(requests).every(n=>n===1));assert.equal(Object.keys(requests).length,39);
pass('send during preload shares downloads',{seconds:active.seconds,requests});await context.close();
context=await browser.newContext();requests={};context.on('request',r=>{if((r.url().endsWith('.bin.gz')||r.url().includes('/assets/chunks/')))requests[r.url()]=(requests[r.url()]||0)+1;});
item=await open(context,benzene);await item.page.waitForFunction(()=>document.documentElement.dataset.racePreload==='ready',null,{timeout:610000});
const complete=await ready(item);assert.equal(Object.keys(requests).length,39);assert(Object.values(requests).every(n=>n===1));
pass('completed preload removes table download from send',{seconds:complete.seconds,requests});
assert.deepEqual(report.errors,[]);report.passed=true;})().catch(e=>{report.failure=e.stack;console.error(e);process.exitCode=1;}).finally(async()=>{fs.writeFileSync(path.join(out,'report.json'),JSON.stringify(report,null,2));await browser?.close();server?.kill();});

