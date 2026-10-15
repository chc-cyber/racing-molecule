/* Static-only probe: serves local files through Playwright route fulfilment.
   Never starts a service, calls /api/load, or reads a live editor design.
   node frontend/editor/research/probe-ketcher-browser.cjs */
const fs=require('node:fs');
const path=require('node:path');
const assert=require('node:assert/strict');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root=path.resolve(__dirname,'../../..');
const frontend=path.join(root,'frontend');
const output=path.join(root,'outputs/ketcher-deploy');
const base='http://127.0.0.1:18777';
const mime={'.html':'text/html','.js':'text/javascript','.css':'text/css','.json':'application/json','.svg':'image/svg+xml','.png':'image/png','.ico':'image/x-icon','.wasm':'application/wasm'};
const library=JSON.parse(fs.readFileSync(path.join(frontend,'editor/fragments.json'),'utf8')).fragments;
function ket(fragment){
  const ids=new Map(fragment.atoms.map((a,i)=>[a.id,i]));
  return {root:{nodes:[{$ref:'mol0'}]},mol0:{type:'molecule',atoms:fragment.atoms.map(a=>({label:a.element,location:[a.pos2[0],a.pos2[1],0],charge:a.charge||0})),bonds:fragment.bonds.map(b=>({type:b.order,atoms:[ids.get(b.a),ids.get(b.b)]}))}};
}
(async()=>{
  fs.mkdirSync(output,{recursive:true});
  const browser=await chromium.launch({executablePath:process.env.CHROME_PATH||'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true});
  const page=await browser.newPage({viewport:{width:1100,height:750}});
  const report={date:new Date().toISOString(),version:'3.18.0',checks:[],pageErrors:[],externalRequests:[],requests:[]};
  page.on('pageerror',e=>report.pageErrors.push(e.message));
  await page.route('**/*',async route=>{
    const url=new URL(route.request().url());
    report.requests.push(url.href);
    if(url.origin!==base){report.externalRequests.push(url.href);return route.abort();}
    assert(!url.pathname.startsWith('/api/'),'Probe must never access race APIs');
    const target=path.resolve(frontend,'.'+decodeURIComponent(url.pathname));
    assert(target.startsWith(frontend+path.sep),'Static paths must remain inside frontend');
    if(!fs.existsSync(target))return route.fulfill({status:404,body:'not found'});
    return route.fulfill({contentType:mime[path.extname(target)]||'application/octet-stream',body:fs.readFileSync(target)});
  });
  const pass=(name,evidence={})=>{report.checks.push({name,passed:true,evidence});console.log('PASS',name);};
  try{
    await page.goto(base+'/editor/research/ketcher-probe.html');
    await page.waitForFunction(()=>window.probe?.ready,{},{timeout:60000});
    const api=await page.evaluate(()=>{
      const k=probe.ketcher,e=k.editor,s=e.struct();
      return {api:{setMolecule:typeof k.setMolecule,getKet:typeof k.getKet,getMolfile:typeof k.getMolfile},editorMethods:Object.getOwnPropertyNames(Object.getPrototypeOf(e)),poolMethods:Object.getOwnPropertyNames(Object.getPrototypeOf(s.atoms)),subscribe:typeof e.subscribe,subscribeSource:e.subscribe.toString(),unsubscribeSource:e.unsubscribe.toString(),selectionSource:e.selection.toString(),eventNames:Object.keys(e.event||{}),changeEventMethods:Object.getOwnPropertyNames(Object.getPrototypeOf(k.changeEvent)),editorOwnKeys:Object.keys(e)};
    });
    assert.equal(api.api.setMolecule,'function');assert.equal(api.api.getKet,'function');assert.equal(api.api.getMolfile,'function');
    pass('Same-origin local iframe initializes window.ketcher and real editor.struct()/pool API',api);
    for(const id of ['benzene','c60']){
      const f=library.find(f=>f.id===id),source=ket(f);
      assert(f,'Missing fixture '+id);
      const result=await page.evaluate(async source=>{
        const k=probe.ketcher;
        await k.setMolecule(JSON.stringify(source),{preserveCanvasPosition:true,needZoom:false});
        const s=k.editor.struct();
        const atoms=[],bonds=[];
        s.atoms.forEach((a,id)=>atoms.push({id,label:a.label,charge:a.charge,pos:[a.pp.x,a.pp.y]}));
        s.bonds.forEach((b,id)=>bonds.push({id,begin:b.begin,end:b.end,type:b.type,stereo:b.stereo}));
        return {atoms,bonds,ket:JSON.parse(await k.getKet()),eventCount:probe.events.length};
      },source);
      assert.equal(result.atoms.length,f.atoms.length);assert.equal(result.bonds.length,f.bonds.length);
      const ordinal=new Map(result.atoms.map((a,i)=>[a.id,i]));
      assert.deepEqual(result.atoms.map(a=>a.label),source.mol0.atoms.map(a=>a.label));
      // Match explicit serialized ordinals and coordinates, never symmetry matching or SMILES.
      for(let i=0;i<result.atoms.length;i++){
        const actual=result.atoms[i].pos,expected=source.mol0.atoms[i].location;
        assert(Math.abs(actual[0]-expected[0])<1e-10);assert(Math.abs(actual[1]+expected[1])<1e-10);
      }
      assert.deepEqual(result.bonds.map(b=>[ordinal.get(b.begin),ordinal.get(b.end),b.type]),source.mol0.bonds.map(b=>[...b.atoms,b.type]));
      pass(id+' direct single-mol0 KET import preserves atom ordinals, coordinates and bond endpoint references',{atomCount:result.atoms.length,bondCount:result.bonds.length,poolIds:result.atoms.map(a=>a.id),eventCount:result.eventCount,coordinateConvention:'pp.x=KET.x; pp.y=-KET.y; preserveCanvasPosition:true'});
    }
    const frame=page.frames().find(f=>f.url().endsWith('/vendor/ketcher/index.html'));
    const controls=await frame.locator('button').evaluateAll(bs=>bs.map(b=>({title:b.title,aria:b.getAttribute('aria-label'),'data-testid':b.getAttribute('data-testid'),text:b.textContent.slice(0,40)})));
    report.controls=controls;
    await page.evaluate(async()=>{await probe.ketcher.setMolecule(JSON.stringify({root:{nodes:[{$ref:'mol0'}]},mol0:{type:'molecule',atoms:[],bonds:[]}}),{preserveCanvasPosition:true});probe.events=[];});
    await frame.locator('[data-testid="bonds"]').click();
    const svgs=await frame.locator('svg').evaluateAll(es=>es.filter(e=>e.getBoundingClientRect().width>200&&e.getBoundingClientRect().height>200).map(e=>({outer:e.outerHTML.slice(0,500),box:{x:e.getBoundingClientRect().x,y:e.getBoundingClientRect().y,width:e.getBoundingClientRect().width,height:e.getBoundingClientRect().height}})));
    report.canvas=svgs;
    await page.mouse.move(450,350);await page.mouse.down();await page.mouse.move(530,400,{steps:6});await page.mouse.up();
    await page.waitForFunction(()=>probe.ketcher.editor.struct().atoms.size===2);
    const drawn=await page.evaluate(()=>{const e=probe.ketcher.editor;return {atoms:e.struct().atoms.size,bonds:e.struct().bonds.size,events:probe.events,selection:e.selection()};});
    assert.equal(drawn.bonds,1);assert(drawn.events.length>0,'Real pointer drawing must dispatch editor change');
    pass('Real pointer bond drawing mutates atom/bond pools and dispatches editor change',drawn);
    const selected=await page.evaluate(async()=>{const e=probe.ketcher.editor,ids=[...e.struct().atoms.keys()];probe.selectionEvents=[];e.selection({atoms:ids});await new Promise(resolve=>setTimeout(resolve,50));return {ids,selection:e.selection(),events:probe.selectionEvents,subscription:{dispatch:e.event.selectionChange.dispatch.toString(),add:e.event.selectionChange.add.toString(),has:e.event.selectionChange.hasHandler(probe.selectionToken.handler),tokenKeys:Object.keys(probe.selectionToken)}};});
    assert.deepEqual(selected.selection.atoms,selected.ids);
    assert.deepEqual(selected.events.at(-1).selection.atoms,selected.ids);
    pass('editor.selection({atoms: poolIds}) reads back explicit atom pool identities and dispatches selectionChange',selected);
    const detached=await page.evaluate(()=>{const e=probe.ketcher.editor;e.unsubscribe('selectionChange',probe.selectionToken);const count=probe.selectionEvents.length;e.selection(null);return {count,after:probe.selectionEvents.length};});
    assert.equal(detached.count,detached.after);
    pass('editor.unsubscribe accepts the subscription token and stops selection notifications',detached);
    await page.screenshot({path:path.join(output,'ketcher-initialized.png')});
    assert.deepEqual(report.externalRequests,[]);assert.deepEqual(report.pageErrors,[]);
    pass('Initialization and structure import use only local files; no browser page errors or external requests');
  }finally{
    report.passed=report.checks.length;
    fs.writeFileSync(path.join(output,'browser-probe.json'),JSON.stringify(report,null,2)+'\n');
    await browser.close();
  }
})().catch(e=>{console.error(e);process.exitCode=1;});
