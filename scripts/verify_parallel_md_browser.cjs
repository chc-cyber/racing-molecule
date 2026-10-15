/* Isolated static server and real Chromium Worker; never touches user sessions. */
const fs = require('node:fs');
const path = require('node:path');
const {spawn} = require('node:child_process');
const assert = require('node:assert/strict');
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root = path.resolve(__dirname, '..');
const out = path.join(root, 'outputs/parallel-md-validation');
(async () => {
  const server = spawn(path.join(root, '.venv/Scripts/python.exe'), ['wasm/tools/serve.py', '--port', '0'], {cwd:root, windowsHide:true});
  let browser;
  const report = {passed:false, scope:'Chromium module Worker, four full 20 nm races, one exploratory 1 fs case; serial browser Worker timing includes message scheduling', cases:[], errors:[], request_failures:[]};
  let serverLog='';server.stderr.on('data',data=>{serverLog+=String(data);});
  try {
    const base = await new Promise((resolve, reject) => {
      const timeout = setTimeout(() => reject(new Error('Static server startup timed out')), 15000);
      server.once('error', reject);
      server.stdout.on('data', data => {
        const match = String(data).match(/http:\/\/127\.0\.0\.1:\d+\//);
        if(match) {clearTimeout(timeout);resolve(match[0]);}
      });
    });
    browser = await chromium.launch({executablePath:process.env.CHROME_PATH || 'C:/Program Files/Google/Chrome/Application/chrome.exe', headless:true});
    report.browser = browser.version();
    const page = await browser.newPage();const requests=[];page.on('request',r=>requests.push(r.url()));
    // A same-origin inert document avoids starting a second, unrelated MD worker.
    await page.route('**/__optimization_probe', route => route.fulfill({contentType:'text/html', body:'<!doctype html><title>WASM verification</title>'}));
    page.on('pageerror', e => report.errors.push(e.message));
    // Worker termination can emit late ERR_ABORTED events for already consumed
    // streams. Actual load errors reject via worker messages, and ready requires
    // all three decompressed asset hashes to pass inside the production Worker.
    page.on('requestfailed',r=>report.request_failures.push({url:r.url(),error:r.failure()?.errorText}));
    report.base_url=base;
    const response=await page.request.get(base+'browser/assets/model.json');
    assert(response.ok(),`Static model request: ${response.status()}`);
    await page.goto(base+'__optimization_probe');
    for(const {car,scale,dt} of [{car:'fullerene-c60',scale:1.5,dt:.5},{car:'fullerene-c60',scale:2,dt:.5},{car:'fullerene-c60',scale:2,dt:1},{car:'pyridine-11',scale:2,dt:.5}]) {
      const content=fs.readFileSync(path.join(root,'outputs/crystal-contrast-guarded/inputs',car+'.mol'),'utf8');
      const result = await page.evaluate(async ({content,scale,dt,car}) => {
        const model = await (await fetch('./browser/assets/model.json')).json();
        return new Promise((resolve,reject) => {
          const worker = new Worker('./browser/worker.mjs', {type:'module'});
          const timeout=setTimeout(() => {worker.terminate();reject(new Error('Worker timed out'));}, 180000);
          let ready,frames=0,firstTime,last;
          worker.onerror=e=>{clearTimeout(timeout);worker.terminate();reject(new Error(e.message));};
          worker.onmessage=({data:e})=>{
            if(e.kind==='error'){clearTimeout(timeout);worker.terminate();reject(new Error(e.error));}
            if(e.kind==='ready'){ready=e.data;firstTime=performance.now();worker.postMessage({kind:'action',action:'start'});}
            if(e.kind==='frame') {
              last=e.frame;frames++;
              if(!Number.isFinite(last.energy)||!last.positions.every(Number.isFinite)) {
                clearTimeout(timeout);worker.terminate();reject(new Error('Nonfinite Worker frame'));return;
              }
              if(last.step>=800000){worker.postMessage({kind:'action',action:'pause'});}
            }
            if(e.kind==='status'&&e.status==='finished'&&last?.center[0]>=200) {
              clearTimeout(timeout);worker.terminate();
              resolve({car,dt,scale,frames,step:last.step,center:last.center,wall_seconds:(performance.now()-firstTime)/1000,runtime:ready.runtime});
            }
          };
          worker.postMessage({kind:'load',content,filename:car+'.mol',start_pose:{mode:'auto'},
            settings:{...model.settings,surface:'fourier-au',temperature:150,lattice_scale:scale,timestep_fs:dt},drive:{mode:'gradient',strength:.6,growth:1}});
        });
      }, {content,scale,dt,car});
      assert(result.runtime.preparation.accepted);
      assert.equal(result.runtime.surface.optimization_level,3);assert.equal(result.runtime.surface.fourier_enabled,true);
      assert.equal(result.runtime.timestep_fs,dt);assert.equal(result.runtime.surface.proxy_charts,3);
      assert(result.step>4000);assert(result.center[0]>=200);
      report.cases.push(result);
    }
    assert.equal(report.errors.length,0);
    assert(report.request_failures.every(r=>r.error==='net::ERR_ABORTED'));
    assert(requests.some(r=>r.endsWith('fourier-atlas-0.bin.gz')));assert(!requests.some(r=>r.endsWith('fourier-1.bin.gz')));report.passed=true;
  } catch(error) {report.errors.push(error.stack||String(error));process.exitCode=1;}
  finally {
    if(browser)await browser.close();server.kill();
    report.server_log=serverLog;
    fs.writeFileSync(path.join(out,'browser-v2-races.json'),JSON.stringify(report,null,2)+'\n');
    console.log(JSON.stringify({passed:report.passed,cases:report.cases.length,errors:report.errors}));
  }
})();
