'use strict';
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const front = path.resolve(__dirname, '..', '..');
const base = 'http://127.0.0.1:18781';
const output = path.resolve(front, '..', 'outputs/ketcher-native-c60');
fs.mkdirSync(output, {recursive: true});
const report = {date: new Date().toISOString(), checks: [], pageErrors: [], externalRequests: [], apiRequests: []};
const check = (name, value) => {assert(value, name); report.checks.push({name, passed: true});};
const snapshot = frame => frame.evaluate(() => {
  const e = window.ketcher.editor, s = e.struct();
  return {atoms: [...s.atoms].filter(([,a])=>!a.isPreview).map(([id,a])=>({id,label:a.label,pp:[a.pp.x,a.pp.y],charge:a.charge,isotope:a.isotope})), bonds: [...s.bonds].filter(([,b])=>!b.isPreview).map(([id,b])=>({id,begin:b.begin,end:b.end,type:b.type})), tool:e.tool()?.constructor?.name};
});
const position = (frame, atomId) => frame.evaluate(id => {
  const e=window.ketcher.editor,r=e.render,p=e.struct().atoms.get(id).pp;
  const q=r.page2obj({clientX:0,clientY:0}), x=r.page2obj({clientX:100,clientY:0}),y=r.page2obj({clientX:0,clientY:100});
  const a=(x.x-q.x)/100,b=(y.x-q.x)/100,c=(x.y-q.y)/100,d=(y.y-q.y)/100,det=a*d-b*c,dx=p.x-q.x,dy=p.y-q.y;
  return {x:(d*dx-b*dy)/det,y:(-c*dx+a*dy)/det};
},atomId);
(async()=>{
  const browser = await chromium.launch({executablePath:process.env.CHROME_PATH || 'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true});
  const page = await browser.newPage({viewport:{width:1200,height:860}});
  page.on('pageerror',e=>report.pageErrors.push(e.message));
  await page.route('**/*',async route=>{
    const u=new URL(route.request().url());
    if(u.origin!==base){report.externalRequests.push(u.href);return route.abort();}
    if(u.pathname.startsWith('/api/')){report.apiRequests.push(u.pathname);return route.abort();}
    if(u.pathname==='/probe.html')return route.fulfill({contentType:'text/html',body:'<!doctype html><style>body{margin:0}iframe{border:0;width:1200px;height:860px}</style><iframe src="/editor/vendor/ketcher-c60/index.html"></iframe>'});
    const file=path.resolve(front,'.'+decodeURIComponent(u.pathname));
    assert(file.startsWith(front+path.sep));
    if(!fs.existsSync(file))throw new Error('Missing local asset: '+file);
    const mime={'.html':'text/html','.js':'text/javascript','.css':'text/css','.json':'application/json','.png':'image/png','.svg':'image/svg+xml','.woff2':'font/woff2'};
    return route.fulfill({contentType:mime[path.extname(file)]||'application/octet-stream',body:fs.readFileSync(file)});
  });
  try {
    await page.goto(base+'/probe.html');
    const frame=await page.locator('iframe').contentFrame();
    await frame.locator('[data-testid="template-8"]').waitFor({timeout:60000});
    const f=page.frames().find(f=>f.url().endsWith('/editor/vendor/ketcher-c60/index.html'));
    await f.waitForFunction(()=>Boolean(window.ketcher?.editor),{},{timeout:60000});
    check('same-origin iframe exports initialized Ketcher',await page.evaluate(()=>Boolean(document.querySelector('iframe').contentWindow.ketcher?.editor)));
    const toolbar = f.locator('[data-testid="bottom-toolbar"]');
    report.toolbar = await toolbar.locator('button').evaluateAll(buttons=>buttons.map(b=>({id:b.dataset.testid,title:b.title,text:b.textContent,top:b.getBoundingClientRect().top})));
    check('C60 is a native button in the same row as original templates',report.toolbar.some(b=>b.id==='template-8') && Math.abs(report.toolbar.find(b=>b.id==='template-8').top-report.toolbar.find(b=>b.id==='template-0').top)<1);
    check('all eight original native ring templates remain',Array.from({length:8},(_,i)=>'template-'+i).every(id=>report.toolbar.some(b=>b.id===id)));
    await f.locator('[data-testid="template-8"]').click();
    check('native C60 button enters native template mode with the real C60 structure',await f.evaluate(()=>window.ketcher.editor.tool()?.template?.molecule?.name==='C60 Fullerene'));
    await page.mouse.move(600,390);
    await page.mouse.click(600,390);
    await f.locator('[data-testid="select-rectangle"]').last().click();
    await page.mouse.move(1198,0);
    report.inserted=await snapshot(f);
    check('native canvas placement creates exactly 60 carbon atoms and 90 real bonds',report.inserted.atoms.length===60&&report.inserted.bonds.length===90&&report.inserted.atoms.every(a=>a.label==='C'));
    const degree = new Map(report.inserted.atoms.map(a=>[a.id,0])), valence=new Map(degree);
    for(const b of report.inserted.bonds){degree.set(b.begin,degree.get(b.begin)+1);degree.set(b.end,degree.get(b.end)+1);valence.set(b.begin,valence.get(b.begin)+b.type);valence.set(b.end,valence.get(b.end)+b.type);}
    check('all C60 atoms have degree 3 and valence 4, including thirty double bonds',[...degree.values()].every(n=>n===3)&&[...valence.values()].every(n=>n===4)&&report.inserted.bonds.filter(b=>b.type===2).length===30);
    const mol = await f.evaluate(()=>window.ketcher.getMolfile());
    fs.writeFileSync(path.join(output,'native-c60.mol'),mol);
    check('native MOL export contains 60 atoms and 90 bonds',mol.split(/\r?\n/)[3].slice(0,6)===' 60 90');
    await page.screenshot({path:path.join(output,'native-c60-inserted.png')});
    await f.locator('[data-testid="undo"]').last().click();
    check('native undo removes the entire template',(await snapshot(f)).atoms.length===0);
    await f.locator('[data-testid="redo"]').last().click();
    check('native redo restores the complete template',(await snapshot(f)).bonds.length===90);
    // Actual native bond drawing attaches a new atom; existing C60 graph is
    // preserved. The saturated connecting carbon intentionally shows valence 5.
    const epos=await position(f,report.inserted.atoms[0].id);
    await f.locator('[data-testid="bonds"]').last().click();
    await page.mouse.move(epos.x,epos.y);
    await page.mouse.down();
    await page.mouse.move(epos.x+70,epos.y-38,{steps:12});
    await page.mouse.up();
    await f.locator('[data-testid="select-rectangle"]').last().click();
    await page.mouse.move(1198,0);
    report.connected=await snapshot(f);
    check('native bond tool adds one explicit atom and connecting bond',report.connected.atoms.length===61&&report.connected.bonds.length===91);
    check('C60 original atoms retain their live IDs and ninety bonds',report.inserted.atoms.every(a=>report.connected.atoms.some(c=>c.id===a.id))&&report.inserted.bonds.every(b=>report.connected.bonds.some(c=>c.id===b.id&&c.type===b.type&&c.begin===b.begin&&c.end===b.end)));
    await page.screenshot({path:path.join(output,'native-c60-connected.png')});
    await f.evaluate(async mol=>window.ketcher.setMolecule(mol),mol);
    const reloaded=await snapshot(f);
    check('native exported MOL reimports to complete C60 topology',reloaded.atoms.length===60&&reloaded.bonds.length===90&&reloaded.bonds.filter(b=>b.type===2).length===30);
    check('no external network, API calls, or page errors',report.pageErrors.length===0&&report.externalRequests.length===0&&report.apiRequests.length===0);
    report.passed=true;
    console.log(`${report.checks.length} native C60 browser checks passed`);
  } finally {
    fs.writeFileSync(path.join(output,'browser-report.json'),JSON.stringify(report,null,2)+'\n');
    await browser.close();
  }
})().catch(e=>{console.error(e);process.exitCode=1;});
