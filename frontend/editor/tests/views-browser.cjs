/* Isolated viewport acceptance: static requests are fulfilled from local files.
   No server is started, no race endpoint is accessed, no live design is changed.
   Run: node frontend/editor/tests/views-browser.cjs */
const fs=require('node:fs');
const path=require('node:path');
const assert=require('node:assert/strict');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root=path.resolve(__dirname,'../../..');
const html=`<!doctype html><style>
body{margin:0}#vp{position:relative;width:980px;height:650px;background:#f8faf5}
#two,#three{position:absolute;inset:0;width:100%;height:100%;touch-action:none}
#three[hidden],#two[hidden]{display:none}.labels-3d{position:absolute;inset:0;pointer-events:none}
.labels-3d span{position:absolute;font-size:12px}.atom-label{font-size:14px}.atom-id{font-size:9px}
</style><script type="importmap">{"imports":{"three":"/vendor/three.module.js"}}</script>
<div id="vp"><svg id="two"></svg><div id="three" hidden></div></div><script type="module">
import {WorkshopModel,selectedIds} from '/editor/model.js'; import {WorkshopViews} from '/editor/views.js';
window.m=new WorkshopModel();window.mode='bond';window.config=null;window.options={numbers:false,lasso:false,anchors:[],movingIds:[],fixedIds:[]};window.drags=[];window.hover=null;window.errors=[];
function selection(ids){m.selection=new Set(ids);v.render(m.project,m.selection);}
function pick(hit,e){const ids=hit.atomId?[hit.atomId]:m.project.bonds.filter(b=>b.id===hit.bondId).flatMap(b=>[b.a,b.b]);if(e.shiftKey){const s=new Set(m.selection);for(const id of ids)s.has(id)?s.delete(id):s.add(id);selection(s);}else selection(ids);}
window.v=new WorkshopViews(document.querySelector('#two'),document.querySelector('#three'),{
mode:()=>mode,axis:()=>[0,0,1],revision:()=>m.revision,moveIds:()=>selectedIds(m.project,m.selection,true),interaction:()=>{},display:()=>options,previewFragment:()=>config,hover:x=>window.hover=x,pick,blank:()=>selection([]),
stroke:s=>{const r=m.drawBondGesture(s);selection([r.endId]);},
region:(ids,e)=>{if(e.shiftKey){const s=new Set(m.selection);ids.forEach(id=>s.has(id)?s.delete(id):s.add(id));selection(s);}else selection(ids);},
doublePick:hit=>{let s=new Set(hit.atomId?[hit.atomId]:m.project.bonds.filter(b=>b.id===hit.bondId).flatMap(b=>[b.a,b.b]));for(let last=-1;last!==s.size;){last=s.size;m.project.bonds.forEach(b=>{if(s.has(b.a)||s.has(b.b)){s.add(b.a);s.add(b.b);}});}selection(s);},
drag:d=>{drags.push({mode:d.mode,axis:d.axis,origin:d.origin,delta:d.delta,degrees:d.degrees});if(d.mode==='rotate')m.rotateSelection(d.axis||[0,0,1],d.degrees,{ids:d.ids,rigid:false,origin:d.origin});else m.moveSelection(d.delta,{ids:d.ids,rigid:false,dimension:d.dimension});},
fragmentPlace:(p,h,e,d)=>{try{m.placeFragment(config.template,{targetId:h.atomId||null,portIndex:config.portIndex||0,pos2:p.slice(0,2),pos3:d===3?p:[...p,0],dimension:d,order:config.order||1});config=null;v.setFragmentPreview(null);}catch(error){errors.push(error.message);}}
});m.addEventListener('change',()=>v.render(m.project,m.selection));window.choose=selection;v.render(m.project,m.selection);
</script>`;

(async()=>{
  const browser=await chromium.launch({executablePath:process.env.CHROME_PATH||'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-unsafe-swiftshader']});
  const page=await browser.newPage({viewport:{width:1000,height:700}}),pageErrors=[],checks=[];
  page.on('pageerror',e=>pageErrors.push(e.stack));
  await page.route('**/*',async route=>{
    const url=new URL(route.request().url());
    if(url.pathname==='/test')return route.fulfill({contentType:'text/html',body:html});
    assert(!url.pathname.startsWith('/api/'),'Viewport suite must never call live APIs');
    const target=path.resolve(root,'frontend','.'+decodeURIComponent(url.pathname));
    assert(target.startsWith(path.join(root,'frontend')+path.sep),'Static test path must remain in frontend');
    if(fs.existsSync(target))return route.fulfill({contentType:target.endsWith('.js')?'text/javascript':'application/json',body:fs.readFileSync(target)});
    return route.fulfill({status:404,body:'not found'});
  });
  const pass=name=>{checks.push(name);console.log('PASS',name);};
  const coord=async i=>page.evaluate(i=>v.screen2(m.project.atoms[i].pos2),i);
  const project=()=>page.evaluate(()=>structuredClone(m.project));
  const drag=async(a,b,button='left')=>{await page.mouse.move(...a);await page.mouse.down({button});await page.mouse.move(...b,{steps:6});await page.mouse.up({button});};
  const setMode=mode=>page.evaluate(mode=>{window.mode=mode;v.updateMode();},mode);
  try{
    await page.goto('http://127.0.0.1:8765/test');await page.waitForFunction(()=>window.v);
    await drag([400,325],[550,340]);let p=await project();assert.equal(p.atoms.length,2);assert.equal(p.bonds.length,1);
    assert(Math.abs(Math.hypot(...p.atoms[1].pos2.map((x,i)=>x-p.atoms[0].pos2[i]))-1.5)<1e-8);assert.equal(await page.evaluate(()=>m._undo.length),1);
    pass('Blank bond drag creates two carbons + one bond in one history step, fixed 1.5 layout');
    await drag(await coord(1),[555,235]);assert.equal((await project()).atoms.length,3);await drag(await coord(2),await coord(0));assert.equal((await project()).bonds.length,3);
    pass('Atom extension and closure snap use existing endpoints');
    assert.equal(await page.locator('.atom-label').count(),0);assert.equal(await page.locator('.atom-id').count(),0);
    await page.evaluate(()=>{options.numbers=true;v.render(m.project,m.selection);});assert.equal(await page.locator('.atom-id').count(),3);
    pass('Skeletal carbon rendering and optional numbering');
    await page.keyboard.down('Alt');const freeStart=await coord(1),freeEnd=[freeStart[0]+180,freeStart[1]+35];await drag(freeStart,freeEnd);await page.keyboard.up('Alt');
    p=await project();assert(Math.abs(Math.hypot(...p.atoms.at(-1).pos2.map((x,i)=>x-p.atoms[1].pos2[i]))-1.5)>.5);
    await page.evaluate(()=>m.undo());assert.equal((await project()).atoms.length,3);pass('Alt disables length/angle snap, undo restores the entire stroke');
    await setMode('select');await page.evaluate(()=>choose([]));const before=await project(),a=await coord(0);await drag(a,[a[0]-25,a[1]+30]);p=await project();
    assert.notDeepEqual(p.atoms[0].pos2,before.atoms[0].pos2);assert.deepEqual(p.atoms.map(a=>a.pos3),before.atoms.map(a=>a.pos3));
    pass('Select tool drags real 2D layout without changing 3D');
    await drag([200,150],[700,550]);assert.equal(await page.evaluate(()=>m.selection.size),3);pass('Marquee selects contained atoms');
    await page.keyboard.down('Shift');await page.mouse.click(...await coord(0));await page.keyboard.up('Shift');assert.equal(await page.evaluate(()=>m.selection.size),2);pass('Shift click toggles atom selection');
    await page.mouse.dblclick(...await coord(0));assert.equal(await page.evaluate(()=>m.selection.size),3);pass('Double-click selects connected component across SVG redraws');
    await page.mouse.move(...await coord(1));assert.equal(await page.evaluate(()=>hover.atomId),await page.evaluate(()=>m.project.atoms[1].id));pass('Hover callback exposes canonical atom ID for element/bond hotkeys');
    await page.evaluate(()=>{options.lasso=true;choose([]);});await page.mouse.move(200,150);await page.mouse.down();for(const q of [[700,150],[700,550],[200,550],[200,150]])await page.mouse.move(...q,{steps:3});await page.mouse.up();assert.equal(await page.evaluate(()=>m.selection.size),3);pass('Lasso polygon selects contained atoms');
    await page.evaluate(()=>{options.anchors=m.project.atoms.map(a=>a.id);options.movingIds=[m.project.atoms[2].id];options.fixedIds=[m.project.atoms[0].id];v.render(m.project,m.selection);});
    assert.equal(await page.locator('.measurement-label').count(),3);assert.equal(await page.locator('.molecule-atom.fixed').count(),1);assert.equal(await page.locator('.molecule-atom.moving').count(),1);pass('A/B/C measurement labels and actual fixed/moving highlights');
    await page.evaluate(async()=>{const data=await(await fetch('/editor/fragments.json')).json();const list=Array.isArray(data)?data:data.fragments;config={template:list[0],portIndex:0,order:1};v.setFragmentPreview(config);});
    const revision=await page.evaluate(()=>m.revision);await page.mouse.move(800,400);assert.equal(await page.evaluate(()=>m.revision),revision);assert(await page.locator('.molecule-atom.ghost').count()>0);
    const count=(await project()).atoms.length;await page.mouse.click(800,400);assert.equal((await project()).atoms.length,count+6);assert.equal(await page.evaluate(()=>m.revision),revision+1);pass('Fragment ghost is transient; click places real graph in one transaction');
    await page.evaluate(async()=>{const data=await(await fetch('/editor/fragments.json')).json(),list=Array.isArray(data)?data:data.fragments;config={template:list[0],portIndex:0,order:1};v.setFragmentPreview(config);});
    const attachmentBefore=await project();await page.mouse.move(...await coord(0));assert.equal(await page.evaluate(()=>v.fragmentVisual.error),null);await page.mouse.click(...await coord(0));
    p=await project();assert.equal(p.atoms.length,attachmentBefore.atoms.length+6);assert.equal(p.bonds.length,attachmentBefore.bonds.length+7);assert.deepEqual(p.atoms.slice(0,attachmentBefore.atoms.length),attachmentBefore.atoms);pass('Fragment ghost aligns to target and creates valid attachment without moving existing atoms');
    await page.evaluate(async()=>{const data=await(await fetch('/editor/fragments.json')).json(),list=Array.isArray(data)?data:data.fragments;const c60=list.find(t=>t.name.startsWith('C60')),base=m.placeFragment(c60,{pos2:[4,4],pos3:[4,4,0],dimension:2});window.c60Target=base.atomIds[0];config={template:list[0],portIndex:0,order:1};v.setFragmentPreview(config);});
    const c60Point=await page.evaluate(()=>v.screen2(m.project.atoms.find(a=>a.id===c60Target).pos2)),invalidBefore=await project();await page.mouse.move(...c60Point);
    assert.match(await page.evaluate(()=>v.fragmentVisual.error),/价态|连接/);assert(await page.locator('.molecule-atom.invalid').count()>0);await page.mouse.click(...c60Point);assert.deepEqual(await project(),invalidBefore);pass('Saturated C60 target shows invalid ghost and rejects placement without edits');
    await page.evaluate(()=>{config=null;v.setFragmentPreview(null);v.cancelGesture();});assert.equal(await page.locator('.molecule-atom.ghost,.molecule-atom.invalid').count(),0);pass('Cancel removes fragment preview without history changes');
    await page.evaluate(()=>{options={numbers:false,anchors:[],movingIds:[],fixedIds:[]};m.commit('Enable3D',p=>p.has3d=true);mode='bond';choose([]);v.setDimension(3);v.reset();});
    await page.evaluate(()=>choose(m.project.fragments.find(f=>f.name.startsWith('C60')).atomIds));assert.equal(await page.locator('.labels-3d span').count(),0);
    await page.evaluate(()=>{options.numbers=true;v.render(m.project,m.selection);});assert((await page.locator('.labels-3d span').count())>=60);
    await page.evaluate(()=>{options.numbers=false;choose([]);});pass('Selecting all 60 C60 atoms shows no element labels by default; numbering explicitly restores labels');
    const cameraBefore=await page.evaluate(()=>v.camera.position.toArray()),moleculeBefore=JSON.stringify(await project());await drag([700,250],[740,300],'right');
    assert.notDeepEqual(await page.evaluate(()=>v.camera.position.toArray()),cameraBefore);assert.equal(JSON.stringify(await project()),moleculeBefore);
    const targetBefore=await page.evaluate(()=>v.controls.target.toArray());await drag([700,250],[740,270],'middle');assert.notDeepEqual(await page.evaluate(()=>v.controls.target.toArray()),targetBefore);await page.mouse.wheel(0,100);
    assert.equal(await page.evaluate(()=>mode),'bond');assert.equal(JSON.stringify(await project()),moleculeBefore);pass('Right orbit + middle pan + wheel work in bond tool without model changes');
    await page.evaluate(()=>{mode='move';choose([m.project.atoms[0].id]);v.updateMode();});assert.equal(await page.evaluate(()=>v.transform.object===v.pivot),true);pass('Official r180 translation gizmo attaches to selected atoms');
    // Pick a visible screen position belonging to one actual gizmo axis.
    const gizmoPoint=async wanted=>page.evaluate(wanted=>{const c=v.screen3(v.pivot.position.toArray());for(let dy=-100;dy<=100;dy+=4)for(let dx=-100;dx<=100;dx+=4){if(Math.hypot(dx,dy)<20)continue;const x=c[0]+dx,y=c[1]+dy;if(x<2||y<2||x>v.width-2||y>v.height-2)continue;v.transform.pointerHover({x:x/v.width*2-1,y:-(y/v.height)*2+1,button:0});if(v.transform.axis===wanted)return [x,y];}throw Error('No visible '+wanted+' gizmo axis');},wanted);
    let point=await gizmoPoint('X'),gizmoBefore=await project(),historyBefore=await page.evaluate(()=>m._undo.length);await drag(point,[point[0]+30,point[1]+15]);
    p=await project();assert.notDeepEqual(p.atoms[0].pos3,gizmoBefore.atoms[0].pos3);assert.deepEqual(p.bonds,gizmoBefore.bonds);assert.deepEqual(p.atoms.map(a=>a.pos2),gizmoBefore.atoms.map(a=>a.pos2));assert.equal(await page.evaluate(()=>m._undo.length),historyBefore+1);
    pass('Real gizmo pointer drag changes 3D coordinates in one transaction, topology/2D preserved');
    await page.evaluate(()=>{mode='rotate';choose(m.project.atoms.slice(0,3).map(a=>a.id));v.updateMode();});assert.equal(await page.evaluate(()=>v.transform.getMode()),'rotate');point=await gizmoPoint('Z');gizmoBefore=await project();await drag(point,[point[0]+23,point[1]+17]);
    p=await project();assert.notDeepEqual(p.atoms.slice(0,3).map(a=>a.pos3),gizmoBefore.atoms.slice(0,3).map(a=>a.pos3));assert.deepEqual(p.bonds,gizmoBefore.bonds);const lastDrag=await page.evaluate(()=>drags.at(-1));assert.equal(lastDrag.axis.length,3);assert.equal(lastDrag.origin.length,3);
    pass('Real rotation ring drag outputs world axis/pivot and edits the selected 3D fragment');
    await setMode('move');point=await gizmoPoint('X');gizmoBefore=await project();historyBefore=await page.evaluate(()=>m._undo.length);
    await page.mouse.move(...point);await page.mouse.down();await page.mouse.move(point[0]+18,point[1]+13,{steps:4});assert.deepEqual(await project(),gizmoBefore);
    await page.evaluate(()=>v.cancelGesture());await page.mouse.up();assert.deepEqual(await project(),gizmoBefore);assert.equal(await page.evaluate(()=>m._undo.length),historyBefore);pass('Cancelling an in-progress real gizmo drag discards preview and does not commit');
    assert.deepEqual(pageErrors,[]);const report={date:new Date().toISOString(),passed:checks.length,pageErrors,checks};const out=path.join(root,'outputs/editor-ux-check');fs.mkdirSync(out,{recursive:true});fs.writeFileSync(path.join(out,'views-report.json'),JSON.stringify(report,null,2));console.log(JSON.stringify({passed:checks.length,pageErrors},null,2));
  }finally{await browser.close();}
})().catch(error=>{console.error(error.stack);process.exitCode=1;});
