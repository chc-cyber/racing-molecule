/* Capture the owned 8770 preview after the finite weak-drive scan.
   This changes only this new service; source defaults and older races stay untouched. */
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const out = path.resolve(__dirname, '../outputs/20nm-game-check');
const apiReport = JSON.parse(fs.readFileSync(path.join(out, 'api-report.json'), 'utf8'));
const base = 'http://127.0.0.1:8770';
assert(apiReport.passed && apiReport.base_url === base && apiReport.server_pid === 20096);
const report = {base_url: base, server_pid: apiReport.server_pid, created_utc: new Date().toISOString(), passed: false,
  protocol_note: 'The scientific parameter scan includes 5 ps zero-drive thermal pre-equilibration. This preview uses the normal game preparation and reset workflow and has no additional 5 ps pre-equilibration.',
  scope: 'Only the owned new 8770 preview is changed; no source default settings or older services are changed.'};
async function request(endpoint, data) {
  const response = await fetch(base + endpoint, data === undefined ? {} : {method: 'POST', headers: {'Content-Type':'application/json'}, body: JSON.stringify(data)});
  const body = await response.json(); assert(response.ok, JSON.stringify(body)); return body;
}
async function wait(predicate) {
  const end = Date.now() + 240000;
  while (Date.now() < end) {
    const s = await request('/api/state'); assert.notEqual(s.status, 'error', s.error);
    if (predicate(s)) return s;
    await new Promise(r => setTimeout(r, 50));
  }
  throw new Error('Owned preview did not become ready');
}
(async () => {
  let browser;
  try {
    const before = await request('/api/state');
    report.before = {generation: before.generation, status: before.status, step: before.frame?.step, settings: before.settings, control: before.control, molecule: before.molecule?.name, atom_count: before.molecule?.atoms.length};
    const unchangedOwnedState = before.status === 'paused' && before.generation === 11 && before.frame?.step === 0 && before.molecule?.atoms.length === 52 && before.molecule?.name === '1,2,4,5-Tetraphenylbenzene test car' && before.settings.surface === 'continuum' && before.settings.temperature === 10 && before.settings.damping_fs === 350 && before.control.drive.mode === 'gradient' && before.control.drive.strength === .75;
    if (!unchangedOwnedState) {report.preserved = true; report.reason = 'The preview differs from the final owned acceptance state; preserved without changes.'; console.log(report.reason); return;}
    await request('/api/control', {temperature: 600, damping_fs: 350, drive: {mode: 'gradient', strength: .25}});
    await wait(s => s.status === 'paused' && s.settings.temperature === 600 && s.control.drive.strength === .25);
    await request('/api/control', {surface: 'periodic-au'});
    await wait(s => s.status === 'paused' && s.settings.surface === 'periodic-au');
    await request('/api/action', {action:'reset'});
    const actual = await wait(s => s.status === 'paused' && s.frame?.step === 0 && s.settings.surface === 'periodic-au');
    assert.equal(actual.molecule.atoms.length, 52); assert.equal(actual.settings.temperature, 600);
    assert.equal(actual.settings.damping_fs, 350); assert.deepEqual(actual.control.drive, {mode:'gradient', strength:.25});
    assert.deepEqual(actual.track, {start_x:0, finish_x:200, half_width:100, surface_z:0});
    assert.equal(actual.runtime.timestep_fs, .5); assert.equal(actual.runtime.surface.mode, 'periodic-au');
    assert.equal(actual.runtime.surface.accelerator, 'numba'); assert.equal(actual.runtime.external_field_kernel.accelerator, 'numba');
    assert.equal(actual.runtime.external_field_kernel.surface_enabled, false); assert.equal(actual.runtime.mobile_surface_atoms, 0);
    assert.equal(actual.runtime.integrated_atoms, 52);
    report.actual = {generation: actual.generation, status: actual.status, step: actual.frame.step, settings:actual.settings, control:actual.control, track:actual.track, atom_count:52, runtime:actual.runtime};
    browser = await chromium.launch({executablePath: process.env.CHROME_PATH || 'C:/Program Files/Google/Chrome/Application/chrome.exe', headless:true, args:['--enable-unsafe-swiftshader']});
    const page = await browser.newPage({viewport:{width:1480,height:1200}});
    await page.route('**/favicon.ico', r => r.fulfill({status:204}));
    report.page_errors=[]; report.console_errors=[];
    page.on('pageerror', e=>report.page_errors.push(e.message));
    page.on('console', m=>{if(m.type()==='error')report.console_errors.push(m.text());});
    await page.goto(base+'/', {waitUntil:'domcontentloaded'});
    await page.waitForFunction(() => {const s=window.__molecularCarDebug?.getSnapshot();return s?.status==='paused' && s.atomCount===52 && s.surface==='periodic-au' && s.drive.mode==='gradient' && s.drive.strength===.25 && s.renderedStep===0;});
    await page.locator('.advanced-section > summary').click();
    assert.equal(await page.locator('#surface-mode').inputValue(),'periodic-au');
    assert.equal(await page.locator('#temperature').inputValue(),'600');
    assert.equal(await page.locator('#strength').inputValue(),'0.25');
    assert.equal(await page.locator('#damping_fs').inputValue(),'350');
    assert.equal(await page.locator('#temperature-value').textContent(),'600 K');
    assert.equal(await page.locator('#strength-value').textContent(),'0.25');
    assert.equal(await page.locator('#damping_fs-value').textContent(),'350 fs');
    assert.equal(await page.locator('#track-size-label').textContent(),'20 nm × 20 nm');
    const label = await page.locator('#surface-mode option:checked').textContent();
    assert(label.includes('冻结 Au') && label.includes('周期势'));
    const stable=()=>page.evaluate(()=>new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r))));
    await page.locator('#home-camera').click(); await stable();
    report.near = await page.evaluate(()=>window.__molecularCarDebug.getSnapshot());
    assert(report.near.goldSurface.matchedPhysicalLattice && report.near.goldSurface.triangles>0);
    await page.screenshot({path:path.join(out,'weak-drive-ready.png'),fullPage:true});
    await page.locator('#overview-camera').click(); await stable();
    report.overview=await page.evaluate(()=>window.__molecularCarDebug.getSnapshot());
    assert.equal(report.overview.cameraMode,'overview'); assert.equal(report.overview.finishX,200);
    assert.equal(report.overview.goldSurface.triangles,0);
    await page.screenshot({path:path.join(out,'weak-drive-overview.png'),fullPage:true});
    await page.locator('#home-camera').click(); await stable();
    const final=await request('/api/state');
    assert.equal(final.generation,actual.generation); assert.equal(final.status,'paused'); assert.equal(final.frame.step,0);
    assert.deepEqual(final.frame.positions,actual.frame.positions);
    assert.deepEqual(report.page_errors,[]); assert.deepEqual(report.console_errors,[]);
    report.ui={surface_label:label,temperature_label:'600 K',strength_label:'0.25',damping_label:'350 fs',track_label:'20 nm × 20 nm'};
    report.passed=true;
    console.log(JSON.stringify({passed:true,base_url:base,server_pid:report.server_pid,generation:final.generation,status:final.status,step:final.frame.step,settings:final.settings,drive:final.control.drive}));
  } catch(error) {report.error=error.stack; throw error;}
  finally {fs.writeFileSync(path.join(out,'weak-drive-preview.json'),JSON.stringify(report,null,2)+'\n'); if(browser)await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
