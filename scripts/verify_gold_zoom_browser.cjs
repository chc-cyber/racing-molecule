/* Paused rendering acceptance on the isolated 8790 visual preview.
 * Camera instrumentation is inserted in the browser response only; product
 * source, atom positions and MD controls are never changed by this script. */
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');
const assert = require('node:assert/strict');
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'C:/Users/22845/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const out = path.resolve(__dirname, '../outputs/gold-zoom-check');
const stage = process.argv[2] || 'before';
assert(['before', 'after'].includes(stage));
const base = 'http://127.0.0.1:8790';
const digest = value => crypto.createHash('sha256').update(JSON.stringify(value)).digest('hex');
const sourcePaths = ['frontend/app.js', 'frontend/gold-surface.js'];
const report = {stage, base_url: base, created_utc: new Date().toISOString(), shots: [], errors: [], warnings: [], source_hashes: {}};
for (const f of sourcePaths) report.source_hashes[f] = crypto.createHash('sha256').update(fs.readFileSync(path.resolve(__dirname, '..', f))).digest('hex');
const injection = `
window.__goldZoomQA = {
  setCameraPose(position, target) {
    cameraMode = 'free'; needsCameraFocus = false; orbit.enableDamping = false;
    camera.position.fromArray(position); orbit.target.fromArray(target); orbit.update(); updateSceneScale();
  },
  setCameraHeight(height) {
    cameraMode = 'free'; needsCameraFocus = false; orbit.enableDamping = false;
    const delta = camera.position.clone().sub(orbit.target);
    delta.multiplyScalar((height - orbit.target.y) / delta.y);
    camera.position.copy(orbit.target).add(delta); orbit.update(); updateSceneScale();
  },
  setLowAngle() {
    cameraMode = 'free'; needsCameraFocus = false; orbit.enableDamping = false;
    orbit.target.set(15, 2.2, 0); camera.position.set(25, 15, 120); orbit.update(); updateSceneScale();
  },
  inspect() {
    const casters = [], receivers = [];
    scene.traverse(o => {if (o.castShadow) casters.push(o.name || o.type); if (o.receiveShadow) receivers.push(o.name || o.type);});
    return {shadowMapEnabled: renderer.shadowMap.enabled, casters, receivers,
      atomCastShadow: atomsMesh?.castShadow, bondCastShadow: bondsMesh?.castShadow,
      goldMeshes: goldSurface.meshes.map(m => ({name: m.name, visible: m.visible, count: m.count, receiveShadow: m.receiveShadow})),
      farPlane: {height: goldSurface.farPlane.position.y, visible: goldSurface.farPlane.visible, receiveShadow: goldSurface.farPlane.receiveShadow},
      backingHeight: goldSurface.fallbackHeight ?? goldSurface.backing.position.y};
  },
};
`;
(async () => {
  fs.mkdirSync(out, {recursive: true});
  const browser = await chromium.launch({executablePath: process.env.CHROME_PATH || 'C:/Program Files/Google/Chrome/Application/chrome.exe', headless: true, args: ['--enable-unsafe-swiftshader']});
  const page = await browser.newPage({viewport: {width: 1480, height: 1080}, deviceScaleFactor: 1});
  await page.route('**/favicon.ico', route => route.fulfill({status: 204}));
  await page.route('**/app.js', async route => {
    const original = await route.fetch();
    await route.fulfill({response: original, body: (await original.text()).replace(/\nawait pollState\(\);/, '\n' + injection + '\nawait pollState();')});
  });
  page.on('pageerror', e => report.errors.push(e.message));
  page.on('console', m => {if (m.type() === 'error') report.errors.push(m.text()); if (m.type() === 'warning') report.warnings.push(m.text());});
  const apiState = async () => (await page.request.get(base + '/api/state')).json();
  const snapshot = () => page.evaluate(() => window.__molecularCarDebug.getSnapshot());
  const stable = () => page.evaluate(() => new Promise(r => requestAnimationFrame(() => requestAnimationFrame(() => requestAnimationFrame(r)))));
  const capture = async name => {
    await stable();
    const s = await snapshot(), q = await page.evaluate(() => window.__goldZoomQA.inspect());
    const filename = `${stage}-${name}.png`;
    await page.screenshot({path: path.join(out, filename), fullPage: true});
    report.shots.push({name, filename, camera: s.camera, cameraTarget: s.cameraTarget, cameraMode: s.cameraMode,
      renderedStep: s.renderedStep, position_hash: digest(s.positions), goldSurface: s.goldSurface, clipping: s.clipping, renderer: s.renderer, inspection: q});
    console.log(stage, name, 'height', s.camera[1], 'Au', s.goldSurface.triangles);
  };
  try {
    await page.goto(base + '/', {waitUntil: 'domcontentloaded'});
    await page.waitForFunction(() => {
      const s = window.__molecularCarDebug?.getSnapshot(); return s?.status === 'paused' && s.atomCount === 52 && s.positions?.length === 156;
    }, undefined, {timeout: 240000});
    const initial = await apiState(); assert.equal(initial.frame.step, 0);
    report.initial_state = {status: initial.status, step: initial.frame.step, generation: initial.generation, surface: initial.settings.surface, position_hash: digest(initial.frame.positions)};
    await page.locator('#home-camera').click(); await capture('home');
    // Real mouse-wheel events exercise OrbitControls before the reproducible
    // exact-height views at either side of the former 100 Å detail cutoff.
    const box = await page.locator('#scene canvas').boundingBox();
    await page.mouse.move(box.x + box.width * .6, box.y + box.height * .55);
    report.wheel_sequence = [];
    for (let i = 0; i < 18; i++) {
      await page.mouse.wheel(0, 250); await stable();
      const s = await snapshot(); report.wheel_sequence.push({height: s.camera[1], radius: s.goldSurface.patch?.radius, triangles: s.goldSurface.triangles});
      if (s.camera[1] >= 200) break;
    }
    if (stage === 'after') {
      const previous = JSON.parse(fs.readFileSync(path.join(out, 'before-report.json'), 'utf8')).shots.find(s => s.name === 'wheel-far');
      await page.evaluate(previous => window.__goldZoomQA.setCameraPose(previous.camera, previous.cameraTarget), previous);
    }
    await capture('wheel-far');
    for (const h of [70, 95, 99, 101, 130, 200]) {
      await page.locator('#home-camera').click(); await stable();
      await page.evaluate(h => window.__goldZoomQA.setCameraHeight(h), h);
      await capture('height-' + h);
    }
    await page.evaluate(() => window.__goldZoomQA.setLowAngle()); await capture('low-angle');
    await page.locator('#overview-camera').click(); await capture('overview');
    await page.locator('#home-camera').click(); await page.locator('#space-filling-model').click(); await capture('space-filling');
    await page.locator('#ball-stick-model').click(); await page.locator('#home-camera').click();
    const end = await apiState();
    assert.equal(end.status, 'paused'); assert.equal(end.generation, initial.generation); assert.equal(end.frame.step, 0); assert.deepEqual(end.frame.positions, initial.frame.positions);
    for (const s of report.shots) assert.equal(s.position_hash, report.initial_state.position_hash, `${s.name}: camera changed molecular coordinates`);
    assert.deepEqual(report.errors, []);
    if (stage === 'after') {
      const before = JSON.parse(fs.readFileSync(path.join(out, 'before-report.json'), 'utf8'));
      assert.equal(report.initial_state.position_hash, before.initial_state.position_hash);
      for (const s of report.shots) {
        const previous = before.shots.find(p => p.name === s.name);
        assert(previous); for (let i = 0; i < 3; i++) assert(Math.abs(s.camera[i] - previous.camera[i]) < 1e-8, `${s.name}: mismatched camera`);
        assert.equal(s.inspection.atomCastShadow, false); assert.equal(s.inspection.bondCastShadow, false);
        assert.equal(s.inspection.shadowMapEnabled, false); assert.deepEqual(s.inspection.casters, []); assert.deepEqual(s.inspection.receivers, []);
        assert(s.inspection.farPlane.height >= s.goldSurface.radius, `${s.name}: far carrier sits below the top layer`);
        assert(s.inspection.backingHeight > 0, `${s.name}: fallback exposes a lower Au layer`);
      }
    }
    report.passed = true;
  } catch (error) {report.passed = false; report.failure = error.stack; throw error;}
  finally {
    fs.writeFileSync(path.join(out, `${stage}-report.json`), JSON.stringify(report, null, 2) + '\n');
    await browser.close();
  }
  console.log(JSON.stringify({stage, passed: report.passed, shots: report.shots.length, errors: report.errors, warnings: report.warnings}));
})().catch(e => {console.error(e); process.exitCode = 1;});
