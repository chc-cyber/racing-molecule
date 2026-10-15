/* Functional integration only: timings are deliberately not benchmark results. */
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const {spawn}=require('node:child_process');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const root=path.resolve(__dirname,'..'),out=path.join(root,'outputs/wasm-fourier-validation');
(async()=>{
  const server=spawn(path.join(root,'.venv/Scripts/python.exe'),['wasm/tools/serve.py','--port','0'],{cwd:root,windowsHide:true});
  server.stderr.on('data',()=>{});
  let browser;const report={passed:false,errors:[],states:[]};
  try {
    const base=await new Promise((resolve,reject)=>{const t=setTimeout(()=>reject(Error('server timeout')),15000);server.once('error',reject);server.stdout.on('data',d=>{const m=String(d).match(/http:\/\/127\.0\.0\.1:\d+\//);if(m){clearTimeout(t);resolve(m[0]);}});});
    browser=await chromium.launch({executablePath:process.env.CHROME_PATH||'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-unsafe-swiftshader']});
    report.browser=browser.version();const page=await browser.newPage({viewport:{width:1440,height:1000}}),requests=[];
    page.on('pageerror',e=>report.errors.push(e.message));page.on('request',r=>requests.push(r.url()));
    await page.goto(base);
    const ready=async mode=>{await page.waitForFunction(mode=>{const s=window.__molecularCarDebug?.getSnapshot();return s?.connected&&s.status==='paused'&&s.surface===mode&&s.atomCount>0;},mode,{timeout:180000});const s=await page.evaluate(()=>window.__molecularCarDebug.getSnapshot());assert(s.renderer.triangles>0);report.states.push(s);};
    await ready('periodic-au');assert(!requests.some(u=>u.endsWith('fourier-0.bin.gz')));
    await page.locator('.advanced-section').evaluate(e=>e.open=true);
    await page.selectOption('#surface-mode','fourier-au');await ready('fourier-au');
    assert(report.states.at(-1).surfaceMetadata.fourier_enabled);
    await page.screenshot({path:path.join(out,'browser-ui.png'),fullPage:true});
    await page.selectOption('#surface-mode','periodic-au');await ready('periodic-au');
    assert(!report.states.at(-1).surfaceMetadata.fourier_enabled);
    assert(requests.some(u=>u.endsWith('fourier-0.bin.gz')));assert(!requests.some(u=>u.endsWith('fourier-1.bin.gz')));
    assert.equal(report.errors.length,0);report.passed=true;
  }catch(e){report.errors.push(e.stack||String(e));process.exitCode=1;}
  finally{if(browser)await browser.close();server.kill();fs.writeFileSync(path.join(out,'browser-ui.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify({passed:report.passed,errors:report.errors}));}
})();
