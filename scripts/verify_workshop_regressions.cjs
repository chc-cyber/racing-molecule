// Isolated static preview: no race preparation, no existing service modified.
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),{spawn}=require('node:child_process');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root=path.resolve(__dirname,'..'),out=path.join(root,process.env.WEB_TEST_URL?'outputs/workshop-bugfix-public':'outputs/workshop-bugfix');fs.mkdirSync(out,{recursive:true});
const serverCode=`from http.server import ThreadingHTTPServer
from functools import partial
from wasm.tools.serve import StaticHandler, ROOT
class PrefixHandler(StaticHandler):
 def do_GET(self):
  if not self.path.startswith('/racing-molecule/'):
   self.send_error(404);return
  self.path=self.path[len('/racing-molecule'):]
  super().do_GET()
s=ThreadingHTTPServer(('127.0.0.1',0),partial(PrefixHandler,directory=str(ROOT/'wasm/dist')))
print(f'http://127.0.0.1:{s.server_port}/racing-molecule',flush=True)
s.serve_forever()`;
const server=process.env.WEB_TEST_URL?null:spawn(path.join(root,'.venv/Scripts/python.exe'),['-u','-c',serverCode],{cwd:root,windowsHide:true});
const report={checks:[],errors:[],apiCalls:[]};let browser;
async function main(){
 const base=process.env.WEB_TEST_URL?.replace(/\/$/,'')||await new Promise((resolve,reject)=>{server.stdout.on('data',d=>{const m=String(d).match(/http:\/\/127.0.0.1:\d+\/racing-molecule/);if(m)resolve(m[0]);});server.once('error',reject);});
 server?.stderr.on('data',()=>{});
 browser=await chromium.launch({headless:true,executablePath:process.env.CHROME_PATH||'C:/Program Files/Google/Chrome/Application/chrome.exe'});
 report.url=base;
 const page=await browser.newPage({viewport:{width:1440,height:1000}});page.setDefaultTimeout(60000);page.on('pageerror',e=>report.errors.push(e.message));page.on('request',r=>{if(new URL(r.url()).pathname.startsWith('/racing-molecule/api/'))report.apiCalls.push(r.url());});
 const pass=(name,data={})=>{report.checks.push({name,...data});console.log('PASS',name,JSON.stringify(data));};
 const manifest=await page.request.get(base+'/site-manifest.json');assert(manifest.ok());
 const remote=await manifest.json(),local=JSON.parse(fs.readFileSync(path.join(root,'wasm/dist/site-manifest.json')));
 for(const name of ['browser/runtime/molecular-car.wasm','editor/vendor/ketcher-c60-manifest.json','editor/vendor/ketcher-c60/index.html'])assert.equal(remote.files.find(f=>f.path===name)?.sha256,local.files.find(f=>f.path===name)?.sha256,name+' deployed hash');
 report.deployedAssetsVerified=true;
 await page.goto(base+'/editor/index.html');
 const settled=async()=>{await page.waitForFunction(()=>window.workshop?.mature?.ready&&workshop.views?.ready&&!workshop.pendingGeometry&&!workshop.mature.suspended&&!workshop.views.suspended);await page.evaluate(async()=>{await workshop.mature.settled();await workshop.views.settled();});};
 await settled();
 const snapshot=()=>page.evaluate(()=>structuredClone(workshop.model.project));
 const load=async p=>{await page.evaluate(p=>workshop.model.replace(p),p);await settled();};
 const dimensional=async d=>{await page.locator('#view-'+d+'d').click();await settled();assert.equal(await page.evaluate(()=>workshop.dimension),d,await page.locator('#feedback').textContent());};
 await page.evaluate(content=>workshop.importContent(content,'hexahelicene.mol'),fs.readFileSync(path.join(root,'tests/data/workshop/hexahelicene.mol'),'utf8'));await settled();
 const before=await snapshot(),start=performance.now();await dimensional(3);const built=await snapshot();assert.equal(built.atoms.length,42);assert.deepEqual(built.atoms.slice(0,26).map(a=>[a.id,a.pos2]),before.atoms.map(a=>[a.id,a.pos2]));pass('[6] helicene through actual browser Worker and JSmol',{milliseconds:Math.round(performance.now()-start),atoms:built.atoms.length});
 await page.screenshot({path:path.join(out,'hexahelicene.png')});
 await dimensional(2);await dimensional(3);assert.deepEqual(await snapshot(),built);await page.locator('#undo').click();await settled();assert.deepEqual(await snapshot(),before);await page.locator('#redo').click();await settled();assert.deepEqual(await snapshot(),built);pass('helicene switching, undo/redo retain full project and IDs');
 await dimensional(2);await load({...before,atoms:[],bonds:[],fragments:[],metadata:{},has3d:false});
 const frame=await page.locator('#ketcher-frame').elementHandle().then(e=>e.contentFrame());
 await frame.locator('[data-testid="template-8"]').click();const box=await page.locator('#ketcher-frame').boundingBox();await page.mouse.click(box.x+450,box.y+280);await page.mouse.move(box.x+5,box.y+5);await page.keyboard.press('Escape');await settled();
 let cage=await snapshot();assert.equal(cage.atoms.length,60);assert.equal(cage.metadata.ketcherGroups[0].kind,'c60');const t=performance.now();await dimensional(3);let ready=await snapshot();assert.deepEqual(ready.atoms.map(a=>[a.id,a.pos3]),cage.atoms.map(a=>[a.id,a.pos3]));pass('native C60 toolbar placement to spherical JSmol',{milliseconds:Math.round(performance.now()-t)});await page.screenshot({path:path.join(out,'c60.png')});
 await dimensional(2);
 const imine=await page.evaluate(async()=>{const {createProject,validateProject}=await import('./model.js');return validateProject({...createProject('imine'),atoms:[{id:'n',element:'N',pos2:[0,0],pos3:[0,0,0]},{id:'c1',element:'C',pos2:[-1.3,1],pos3:[-1.3,1,0]},{id:'c2',element:'C',pos2:[1.3,1],pos3:[1.3,1,0]}],bonds:[{id:'b1',a:'n',b:'c1',order:1},{id:'b2',a:'n',b:'c2',order:2}]});});await load(imine);
 const point=await frame.evaluate(()=>{const e=ketcher.editor,r=e.render,a=[...e.struct().atoms.values()].find(a=>a.label==='N'),p=a.pp,q=r.page2obj({clientX:0,clientY:0}),x=r.page2obj({clientX:100,clientY:0}),y=r.page2obj({clientX:0,clientY:100}),aa=(x.x-q.x)/100,b=(y.x-q.x)/100,c=(x.y-q.y)/100,d=(y.y-q.y)/100,det=aa*d-b*c;return{x:(d*(p.x-q.x)-b*(p.y-q.y))/det,y:(-c*(p.x-q.x)+aa*(p.y-q.y))/det};});
 await frame.locator('[data-testid="bonds"]').last().click();const fb=await page.locator('#ketcher-frame').boundingBox();await page.mouse.move(fb.x+point.x,fb.y+point.y);await page.mouse.down();await page.mouse.move(fb.x+point.x,fb.y+point.y+90,{steps:12});await page.mouse.up();await page.mouse.move(fb.x+5,fb.y+5);await page.keyboard.press('Escape');await settled();
 const nitrogen=()=>frame.evaluate(()=>{const a=[...ketcher.editor.struct().atoms.values()].find(a=>a.label==='N');return{badConn:a.badConn,implicitH:a.implicitH,charge:a.charge};});
 let n=await nitrogen();assert.equal(n.badConn,true);assert.equal(n.implicitH,0);assert.equal((await snapshot()).bonds.length,3);await page.screenshot({path:path.join(out,'imine-invalid.png')});pass('native bond drag: neutral imine shows valence error without phantom H',n);
 await page.locator('#undo').click();await settled();assert.equal((await snapshot()).bonds.length,2);assert(!(await nitrogen()).badConn);await page.locator('#redo').click();await settled();assert.equal((await nitrogen()).badConn,true);pass('nitrogen error decoration follows shared undo and redo');
 await page.evaluate(()=>workshop.model.commit('charge',p=>p.atoms.find(a=>a.element==='N').charge=1));await settled();n=await nitrogen();assert(!n.badConn);assert.equal(n.implicitH,0);await dimensional(3);pass('positive iminium charge restores valid nitrogen and permits 3D');
 assert.deepEqual(report.errors,[]);assert.deepEqual(report.apiCalls,[]);report.passed=true;
}
main().catch(e=>{report.failure=e.stack;process.exitCode=1;console.error(e);}).finally(async()=>{fs.writeFileSync(path.join(out,'browser-report.json'),JSON.stringify(report,null,2));await browser?.close();server?.kill();});


