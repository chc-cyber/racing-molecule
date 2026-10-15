/* Browser acceptance on the isolated service created by verify_fourier_game.py.
   Never addresses the user's 8765, 8766 or 8777 race. */
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');

const out = path.resolve(__dirname, '../outputs/fourier-game-check');
const apiReport = JSON.parse(fs.readFileSync(path.join(out, 'api-report.json'), 'utf8'));
const base = process.argv[2] || apiReport.base_url;
const url = new URL(base);
assert(apiReport.passed && base === apiReport.base_url, 'Use the successfully verified isolated service');
assert(url.hostname === '127.0.0.1' && !['8765', '8766', '8777'].includes(url.port), 'Do not touch an existing user race');

(async () => {
  const browser = await chromium.launch({
    executablePath: process.env.CHROME_PATH || 'C:/Program Files/Google/Chrome/Application/chrome.exe',
    headless: true, args: ['--enable-unsafe-swiftshader'],
  });
  const page = await browser.newPage({viewport: {width: 1480, height: 1080}, acceptDownloads: true});
  // The prototype has no favicon; suppress only that unrelated automatic fetch.
  await page.route('**/favicon.ico', route => route.fulfill({status: 204}));
  const report = {base_url: base, created_utc: new Date().toISOString(), checks: [], pageErrors: [], consoleErrors: [], consoleWarnings: [], externalRequests: []};
  page.on('pageerror', e => report.pageErrors.push(e.message));
  page.on('console', m => {if (m.type() === 'error') report.consoleErrors.push(m.text()); if (m.type() === 'warning') report.consoleWarnings.push(m.text());});
  page.on('request', r => {if (!r.url().startsWith(base) && !r.url().startsWith('blob:') && !r.url().startsWith('data:')) report.externalRequests.push(r.url());});
  const record = (name, evidence = {}) => {report.checks.push({name, passed: true, evidence}); console.log('PASS', name);};
  const state = async () => (await page.request.get(base + '/api/state')).json();
  const snapshot = () => page.evaluate(() => window.__molecularCarDebug.getSnapshot());
  const ready = (surface, atoms) => page.waitForFunction(({surface, atoms}) => {
    const s = window.__molecularCarDebug?.getSnapshot();
    return s?.status === 'paused' && s.surface === surface && s.atomCount === atoms && s.positions?.length === atoms * 3;
  }, {surface, atoms}, {timeout: 240000});
  const stable = async () => {await page.evaluate(() => new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r))));};
  const comparePaused = (a, b) => {
    assert.equal(b.status, 'paused'); assert.equal(b.generation, a.generation);
    assert.equal(b.frame.step, a.frame.step); assert.deepEqual(b.frame.positions, a.frame.positions);
  };
  try {
    // This UI polls real MD continuously, so network-idle is never a load signal.
    await page.goto(base + '/', {waitUntil: 'domcontentloaded'});
    await ready('fourier-au', 52);
    await page.locator('.advanced-section > summary').click();
    await page.locator('#surface-mode').selectOption('continuum');
    await ready('continuum', 52);
    assert.equal(await page.locator('#ball-stick-model').getAttribute('aria-pressed'), 'true');
    const initial = await snapshot();
    assert.equal(initial.goldSurface.spacing, 2.88); assert.equal(initial.goldSurface.matchedPhysicalLattice, false);
    record('Original continuum default and ball-stick visual remain available');

    const oldGeneration = (await state()).generation;
    await page.locator('#surface-mode').selectOption('fourier-au');
    await ready('fourier-au', 52);
    await page.waitForFunction(() => window.__molecularCarDebug.getSnapshot().goldSurface?.matchedPhysicalLattice);
    const periodic = await state(), near = await snapshot(), meta = periodic.runtime.surface;
    assert(periodic.generation > oldGeneration); assert.equal(periodic.frame.step, 0);
    assert.equal(periodic.runtime.lammps_atoms, 52); assert.equal(periodic.runtime.mobile_surface_atoms, 0);
    assert.equal(meta.accelerator, 'numba'); assert.equal(meta.numba_version, '0.61.2');
    assert.equal(near.goldSurface.spacing, meta.surface_period_xy_angstrom[0]);
    assert.equal(near.goldSurface.radius, 2.1); // Au space-filling display radius; centres retain the physical lattice.
    assert.deepEqual(near.goldSurface.phaseXY, [meta.surface_period_xy_angstrom[0], meta.surface_period_xy_angstrom[1] / 3]);
    assert.deepEqual(near.goldSurface.centreHeights, meta.layer_heights_angstrom.slice(-3).reverse());
    assert.equal(near.goldSurface.stackDirection, -1); assert(near.goldSurface.opaque);
    assert(await page.locator('#corrugation').isDisabled());
    assert((await page.locator('#surface-description').textContent()).includes('仅支持 C、H'));
    assert.equal(await page.locator('#status-label').textContent(), '准备就绪');
    report.surface_metadata = meta;
    assert.equal(periodic.runtime.model, 'uff-fourier-au-gradient-soft-tip-v1');
    assert.equal(meta.mode, 'fourier-au');
    record('Default 52-atom car prepares with compiled Fourier Au and physically matched visual lattice', {generation: periodic.generation, accelerator: meta.accelerator, atoms: 52, gold: near.goldSurface});

    // Controls persist across surface replacement. Set the chosen drive
    // explicitly, including when a previously interrupted test is retried.
    await page.locator('#gradient-drive').click();
    await page.waitForFunction(() => window.__molecularCarDebug.getSnapshot().drive.mode === 'gradient');
    await page.waitForFunction(() => window.__molecularCarDebug.getSnapshot().gradientVisible);

    await page.locator('#start-button').click();
    await page.waitForFunction(() => window.__molecularCarDebug.getSnapshot().renderedStep >= 100);
    await page.waitForFunction(() => document.querySelector('#start-label').textContent === '暂停模拟');
    await page.locator('#start-button').click();
    await ready('fourier-au', 52);
    const paused = await state();
    assert(paused.frame.step >= 100); assert(paused.frame.positions.every(Number.isFinite));
    await page.waitForFunction(step => window.__molecularCarDebug.getSnapshot().renderedStep === step, paused.frame.step);
    await stable();
    const rendered = await snapshot();
    assert.deepEqual(rendered.positions, paused.frame.positions);
    assert.equal(await page.locator('#status-label').textContent(), '已暂停');
    record('Default car executes real Fourier-Au MD then pauses with exact displayed coordinates', {step: paused.frame.step, time_fs: paused.frame.time_fs, atomCount: 52});

    await page.locator('#home-camera').click(); await stable();
    await page.screenshot({path: path.join(out, 'fourier-au-ball-stick.png'), fullPage: true});
    const ball = await snapshot(); assert(ball.bondsVisible && ball.representation === 'ball-stick');
    await page.locator('#space-filling-model').click(); await stable();
    const filling = await snapshot(); assert(!filling.bondsVisible && filling.representation === 'space-filling');
    for (const a of filling.atomVisuals) assert(Math.abs(a.radius - (a.element === 'C' ? 1.7 : 1.2)) < 1e-6);
    assert.deepEqual(filling.positions, ball.positions);
    await page.screenshot({path: path.join(out, 'fourier-au-space-filling.png'), fullPage: true});
    await page.locator('#ball-stick-model').click(); await stable();
    const restored = await snapshot(); assert(restored.bondsVisible && restored.representation === 'ball-stick');
    assert.deepEqual(restored.atomVisuals, ball.atomVisuals); comparePaused(paused, await state());
    record('Actual atom matrices use element van der Waals radii; representation toggles preserve paused MD', {radiiSource: filling.radiiSource, elements: [...new Set(filling.atomVisuals.map(a => a.element))]});

    const canvas = await page.locator('#scene canvas').boundingBox();
    await page.mouse.move(canvas.x + canvas.width * .57, canvas.y + canvas.height * .58);
    await page.mouse.down(); await page.mouse.move(canvas.x + canvas.width * .66, canvas.y + canvas.height * .61, {steps: 8}); await page.mouse.up(); await stable();
    const free = await snapshot(); assert.equal(free.cameraMode, 'free'); assert(free.goldSurface.triangles > 0);
    await page.screenshot({path: path.join(out, 'fourier-au-free-view.png'), fullPage: true});
    await page.locator('#overview-camera').click(); await stable();
    const overview = await snapshot(); assert.equal(overview.cameraMode, 'overview');
    assert.equal(overview.goldSurface.triangles, 0); assert(overview.goldSurface.fullFieldRelief);
    assert.equal(overview.track.finish_x, 1000); assert.equal(overview.track.half_width, 250);
    await page.screenshot({path: path.join(out, 'fourier-au-overview.png'), fullPage: true});
    comparePaused(paused, await state());
    record('Free and whole 100 nm track views preserve coordinates and bounded Au LOD', {nearTriangles: ball.goldSurface.triangles, freeTriangles: free.goldSurface.triangles, overviewTriangles: overview.goldSurface.triangles});

    await page.locator('#tip-drive').click();
    await page.waitForFunction(() => window.__molecularCarDebug.getSnapshot().drive.mode === 'tip');
    await page.waitForFunction(() => {
      const s = window.__molecularCarDebug.getSnapshot(); return s.tipVisible && !s.gradientVisible;
    });
    assert(await page.locator('#tip-controls').isVisible());
    const tipMode = await snapshot(); assert(tipMode.tipVisible && !tipMode.gradientVisible);
    await page.locator('#autopilot').uncheck();
    await page.locator('#tip-target-x').fill('2');
    await page.locator('#tip-target-x').dispatchEvent('change');
    await page.waitForFunction(() => Math.abs(window.__molecularCarDebug.getSnapshot().tip.x - 20) < 1e-6);
    await page.locator('#gradient-drive').click();
    await page.waitForFunction(() => window.__molecularCarDebug.getSnapshot().drive.mode === 'gradient');
    const gradientMode = await snapshot(); assert(!gradientMode.tipVisible && gradientMode.gradientVisible);
    comparePaused(paused, await state());
    record('Needle controls and fixed gradient remain selectable without moving paused atoms');

    await page.locator('#surface-mode').selectOption('continuum'); await ready('continuum', 52);
    const continuum = await state(), continuumView = await snapshot();
    assert(continuum.generation > paused.generation && continuum.frame.step === 0);
    assert.equal(continuumView.goldSurface.spacing, 2.88); assert.equal(continuumView.goldSurface.matchedPhysicalLattice, false);
    assert(!(await page.locator('#corrugation').isDisabled()));
    assert.equal((await (await page.request.get(base + '/api/replay')).json()).surface.mode, 'continuum');
    record('Returning to continuum resets safely and restores the original visual and surface control');
    await page.locator('#surface-mode').selectOption('periodic-au'); await ready('periodic-au', 52);
    const existing = await state();
    assert.equal(existing.runtime.surface.mode, 'periodic-au');
    assert.equal(existing.runtime.model, 'uff-periodic-au-gradient-soft-tip-v1');
    assert.equal((await snapshot()).goldSurface.matchedPhysicalLattice, true);
    record('Existing cached periodic surface remains a separate usable mode');
    await page.locator('#surface-mode').selectOption('fourier-au'); await ready('fourier-au', 52);
    await page.locator('#home-camera').click(); await stable();
    await page.screenshot({path: path.join(out, 'fourier-au-ready.png'), fullPage: true});
    assert.deepEqual(report.pageErrors, []); assert.deepEqual(report.consoleErrors, []); assert.deepEqual(report.externalRequests, []);
    report.passed = true;
  } catch (error) {report.passed = false; report.failure = error.stack; throw error;}
  finally {
    fs.writeFileSync(path.join(out, 'browser-report.json'), JSON.stringify(report, null, 2) + '\n');
    await page.screenshot({path: path.join(out, 'browser-last-state.png'), fullPage: true}).catch(() => {});
    await browser.close();
  }
  console.log(JSON.stringify({passed: report.passed, checks: report.checks.length, pageErrors: report.pageErrors, consoleErrors: report.consoleErrors, consoleWarnings: report.consoleWarnings}));
})().catch(error => {console.error(error); process.exitCode = 1;});
