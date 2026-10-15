// Standalone: bundled Node frontend/gold-surface.test.mjs. No server or MD is touched.
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
const threeUrl = new URL('./vendor/three.module.js', import.meta.url).href;
const THREE = await import(threeUrl);
const source = (await readFile(new URL('./gold-surface.js', import.meta.url), 'utf8'))
  .replace("from 'three'", `from '${threeUrl}'`);
const {GoldSurface, AU_LATTICE, latticePosition, sampleGoldCell, latticeForSurface} =
  await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
const track = {start_x: 0, finish_x: 200, half_width: 100, surface_z: 0};
const close = (actual, expected, tolerance = 1e-6) => assert.ok(Math.abs(actual - expected) <= tolerance,
  `${actual} differs from ${expected}`);
const pointDistance = (a, b) => Math.hypot(...a.map((v, i) => v - b[i]));
const top = latticePosition(0, 0), next = latticePosition(1, 0), row = latticePosition(0, 1);
close(pointDistance(top, next), 2.88); close(pointDistance(top, row), 2.88);
const b = latticePosition(0, 0, 1), c = latticePosition(0, 0, 2);
close(b[1], -2.35); close(c[1], -4.7);
close(pointDistance(top, b), 2.88, .002);
close(sampleGoldCell(0, 0).height, 2.1);
const filledHollowHeight = Math.sqrt(2.1 ** 2 - 2.88 ** 2 / 3);
close(sampleGoldCell(b[0], b[2]).height, filledHollowHeight);
close(sampleGoldCell(c[0], c[2]).height, filledHollowHeight);
let hollowSamples = 0;
for (let iz = 0; iz < 120; iz++) for (let ix = 0; ix < 120; ix++) {
  const x = (ix + .137) / 120 * AU_LATTICE.spacing, z = (iz + .317) / 120 * AU_LATTICE.rowHeight * 2;
  const sample = sampleGoldCell(x, z);
  assert.ok(Number.isFinite(sample.height) && sample.layer >= 0, 'Three ABC layers must cover every cell pixel');
  close(Math.hypot(...sample.normal), 1);
  const shifted = sampleGoldCell(x + 2.88, z + AU_LATTICE.rowHeight * 2);
  close(sample.height, shifted.height);
  for (let i = 0; i < 3; i++) close(sample.normal[i], shifted.normal[i]);
  if (sample.layer > 0) hollowSamples++;
}
assert.equal(hollowSamples, 0, 'Overlapping space-filling atoms close the top-layer gaps');
const gold = new GoldSurface(track), camera = new THREE.PerspectiveCamera(39, 750 / 570, .1, 5000);
const matrix = new THREE.Matrix4();
function instanceSites() {
  const sites = new Set();
  gold.meshes.forEach(mesh => {
    const radius = mesh.geometry.parameters.radius;
    close(radius, 2.1); assert.ok(mesh.boundingSphere?.radius > 0 || !mesh.count);
    for (let i = 0; i < mesh.count; i++) {
      mesh.getMatrixAt(i, matrix); close(matrix.getMaxScaleOnAxis(), 1);
      const x = matrix.elements[12], y = matrix.elements[13], z = matrix.elements[14];
      const layer = Math.round(-y / 2.35), row = Math.round(z / AU_LATTICE.rowHeight - layer / 3);
      const column = Math.round(x / 2.88 - row / 2 - layer / 2);
      const expected = latticePosition(column, row, layer);
      [x, y, z].forEach((v, j) => close(v, expected[j], 1e-4));
      sites.add(`${column}/${row}/${layer}`);
    }
  });
  return sites;
}
const cases = [];
for (const [label, x, y, z] of [['start follow', 7, 25, 30], ['middle free', 100, 25, 0],
  ['finish follow', 207, 25, 30], ['side free', 100, 13, 95], ['oblique free', 100, 9, 40],
  ['zoom transition', 100, 88, 0], ['overview', 100, 400, 500]]) {
  camera.position.set(x, y, z); gold.update(camera);
  const snapshot = gold.getSnapshot();
  instanceSites();
  const count = snapshot.meshes.reduce((sum, mesh) => sum + mesh.count, 0);
  assert.ok(count < 6200, `${label}: bounded local instances`);
  assert.ok(snapshot.detailTriangles < 850000 && snapshot.triangles < 1200000, `${label}: bounded interior and edge triangle budget`);
  assert.equal(snapshot.opaque, true); assert.equal(snapshot.fullFieldRelief, true);
  assert.deepEqual(snapshot.bounds, {minX: -12, maxX: 212, minZ: -100, maxZ: 100});
  if (label === 'overview') assert.equal(count, 0);
  if (label !== 'overview') assert.ok(count > 0);
  cases.push({label, instances: count, triangles: snapshot.triangles, patch: snapshot.patch});
}
camera.position.set(100, 25, 0); gold.update(camera);
// Side geometry must survive after interior spheres are removed by far LOD.
camera.position.set(100, 400, 500); gold.update(camera);
assert.equal(gold.meshes.reduce((n, mesh) => n + mesh.count, 0), 0);
const edgeCount = gold.edgeMeshes.reduce((n, mesh) => n + mesh.count, 0);
assert.ok(edgeCount > 0 && edgeCount < 3500, 'bounded permanent rim survives overview');
gold.group.updateMatrixWorld(true);
const ray = new THREE.Raycaster();
let sideRays = 0;
for (const height of [-1, -3]) for (let i = 0; i < 16; i++) {
  const t = (i + .5) / 16;
  const x = gold.bounds.minX + 3 + t * (gold.bounds.maxX - gold.bounds.minX - 6);
  const z = gold.bounds.minZ + 3 + t * (gold.bounds.maxZ - gold.bounds.minZ - 6);
  for (const [origin, direction] of [
    [[x, height, gold.bounds.maxZ + 15], [0, 0, -1]],
    [[x, height, gold.bounds.minZ - 15], [0, 0, 1]],
    [[gold.bounds.maxX + 15, height, z], [-1, 0, 0]],
    [[gold.bounds.minX - 15, height, z], [1, 0, 0]],
  ]) {
    ray.set(new THREE.Vector3(...origin), new THREE.Vector3(...direction));
    const hit = ray.intersectObjects(gold.edgeMeshes, false).find(hit => {
      const p = hit.point;
      return Math.min(p.x - gold.bounds.minX, gold.bounds.maxX - p.x,
        p.z - gold.bounds.minZ, gold.bounds.maxZ - p.z) <= gold.getSnapshot().edgeWidth;
    });
    assert.ok(hit && hit.distance < 22, `exposed side must contain real spheres: ${origin}`);
    sideRays++;
  }
}
camera.position.set(100, 20, 180); gold.update(camera, new THREE.Vector3(100, 2, 50));
assert.ok(Math.abs(gold.patch.z) < 65, 'off-track camera detail follows its viewing ray');
assert.ok(gold.meshes.some(mesh => mesh.count > 0));
close(gold.backing.position.y + gold.backing.geometry.parameters.height / 2, gold.fallbackHeight);
close(gold.backing.position.y - gold.backing.geometry.parameters.height / 2, gold.bottomHeight);
camera.position.set(100, 25, 0); gold.update(camera);
const fixedSites = instanceSites(), rebuilds = gold.rebuildCount;
camera.position.x += .01; gold.update(camera);
assert.equal(gold.rebuildCount, rebuilds, 'sub-cell camera motion does not rebuild or slide atoms');
camera.position.x += 18; gold.update(camera);
const movedSites = instanceSites();
assert.ok([...fixedSites].filter(site => movedSites.has(site)).length > 3500, 'allocation changes preserve shared absolute lattice sites');
const shader = {uniforms: {}, vertexShader: THREE.ShaderLib.standard.vertexShader,
  fragmentShader: THREE.ShaderLib.standard.fragmentShader};
gold.farMaterial.onBeforeCompile(shader);
assert.ok(shader.vertexShader.includes('goldWorldPosition = (modelMatrix * goldVertex).xyz'));
assert.ok(shader.fragmentShader.includes('textureGrad(goldRelief'));
assert.ok(shader.fragmentShader.includes('gl_FragDepth = goldClip.z'));
assert.ok(shader.fragmentShader.includes('normal = normalize(mat3(viewMatrix)'));
assert.ok(shader.fragmentShader.includes('goldDetail * goldDetail'));
assert.equal(gold.farMaterial.roughness, gold.atomMaterial.roughness);
assert.ok(gold.reliefTexture.userData.averageHeight > filledHollowHeight && gold.reliefTexture.userData.averageHeight < 2.1);
assert.equal(shader.uniforms.goldProjection, gold.uniforms.goldProjection);
const shifted = new GoldSurface({...track, surface_z: 8});
camera.position.y = 33; shifted.update(camera);
assert.deepEqual(shifted.getSnapshot().centreHeights, [8, 5.65, 3.3]);
close(shifted.getSnapshot().topHeight, 10.1);
for (const mesh of shifted.getSnapshot().meshes) for (const sample of mesh.instances) {
  const layer = Math.round((8 - sample.centre[2]) / 2.35);
  close(sample.centre[2], 8 - layer * 2.35, 1e-5);
}
shifted.dispose(); gold.dispose();

// Independent check against the physics bottom-to-top six-layer rectangular
// cell: shifts=(k%3)*(d/2,h/6), basis=(0,0)/(d/2,h/2).
const physicalD = 4.08 / Math.sqrt(2), physicalH = Math.sqrt(3) * physicalD;
const physicalHeights = [-11.75, -9.39, -7.04, -4.697, -2.348, 0];
const physicalLattice = latticeForSurface({mode: 'periodic-au',
  surface_period_xy_angstrom: [physicalD, physicalH], layer_heights_angstrom: physicalHeights});
const fourierLattice = latticeForSurface({mode: 'fourier-au',
  surface_period_xy_angstrom: [physicalD, physicalH], layer_heights_angstrom: physicalHeights});
assert.deepEqual(fourierLattice, physicalLattice);
assert.equal(latticeForSurface({mode: 'continuum'}), AU_LATTICE);
assert.equal(physicalLattice.stackDirection, -1);
close(physicalLattice.phaseZ, physicalH / 3);
const physicalGold = new GoldSurface(track, null, physicalLattice);
camera.position.set(8, 25, 30); physicalGold.update(camera);
let physicalSites = 0;
for (const mesh of physicalGold.meshes) {
  close(mesh.geometry.parameters.radius, 2.1);
  for (let i = 0; i < mesh.count; i++) {
    mesh.getMatrixAt(i, matrix);
    const [x, y, z] = [matrix.elements[12], matrix.elements[13], matrix.elements[14]];
    const layer = physicalLattice.layerHeights.findIndex(height => Math.abs(y - height) < 1e-5);
    assert.ok(layer >= 0);
    const k = 5 - layer, phase = k % 3;
    const rowFraction = (z - phase * physicalH / 6) / physicalH;
    const nearestHalfRow = Math.round(rowFraction * 2);
    close(rowFraction * 2, nearestHalfRow, 3e-5);
    const basis = ((nearestHalfRow % 2) + 2) % 2;
    const column = (x - phase * physicalD / 2 - basis * physicalD / 2) / physicalD;
    close(column, Math.round(column), 3e-5);
    physicalSites++;
  }
}
const physicalSnapshot = physicalGold.getSnapshot();
assert.equal(physicalSnapshot.matchedPhysicalLattice, true);
assert.deepEqual(physicalSnapshot.centreHeights, physicalHeights.slice(-3).reverse());
close(sampleGoldCell(physicalD, physicalH / 3, physicalLattice).height, 2.1);
let physicalHollows = 0;
for (let iz = 0; iz < 80; iz++) for (let ix = 0; ix < 80; ix++) {
  const x = (ix + .3) / 80 * physicalD, z = (iz + .1) / 80 * physicalH;
  const first = sampleGoldCell(x, z, physicalLattice), repeated = sampleGoldCell(x + physicalD, z + physicalH, physicalLattice);
  assert.ok(Number.isFinite(first.height));
  close(first.height, repeated.height);
  if (first.layer > 0) physicalHollows++;
}
assert.equal(physicalHollows, 0);
// Far texture samples the same shifted lattice as the real sphere instances.
const texture = physicalGold.reliefTexture.image.data;
for (const [ix, iz] of [[20, 15], [90, 77], [155, 230], [235, 120]]) {
  const sample = sampleGoldCell((ix+.5)/256*physicalD, (iz+.5)/256*physicalH, physicalLattice);
  const decoded = physicalLattice.layerHeights[2] + texture[(iz*256+ix)*4+3] / 255 *
    (physicalLattice.radius - physicalLattice.layerHeights[2]);
  close(decoded, sample.height, .013);
}
physicalGold.dispose();
console.log(JSON.stringify({passed: true, lattice: AU_LATTICE, samples: 14400, hollowSamples, edgeCount, sideRays,
  cases, unchangedSharedSites: [...fixedSites].filter(site => movedSites.has(site)).length,
  physicalLattice: {spacing: physicalD, phaseXY: physicalSnapshot.phaseXY, physicalSites,
    centreHeights: physicalSnapshot.centreHeights, physicalHollows},
  checks: ['real instance matrices and sphere radii', 'ABC centre planes', 'whole-cell normal/height coverage',
    'absolute lattice across patch changes', 'follow/free/overview bounded LOD', 'stable sub-cell camera motion',
    'full-field relief shader wiring', 'nonzero surface_z convention',
    'real geometry on all four sides after far LOD', 'off-track viewing-ray detail', 'closed core top and bottom']}, null, 2));
