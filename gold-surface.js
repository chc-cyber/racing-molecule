import * as THREE from 'three';

// Au(111) visual lattice, in Å. surface_z denotes the TOP-LAYER ATOM CENTRES.
// This display is independent of the continuous MD wall and never edits car coordinates.
// Match the car's space-filling Au radius: RDKit GetRvdw(79) = 2.1 Å.
// Overlap closes the triangular gaps while preserving the real lattice sites.
export const AU_LATTICE = Object.freeze({
  spacing: 2.88, radius: 2.1, rowHeight: Math.sqrt(3) * 1.44,
  layerSpacing: 2.35, layers: 3, patchRadius: 64, maxDetailHeight: 100,
});
const TEXTURE_SIZE = 256;
const CAPACITY = 3200;
const EDGE_WIDTH = 6;

export function latticeForSurface(metadata) {
  if (!['periodic-au', 'fourier-au'].includes(metadata?.mode)) return AU_LATTICE;
  const period = metadata.surface_period_xy_angstrom, heights = metadata.layer_heights_angstrom;
  if (!Array.isArray(period) || period.length !== 2 || !period.every(v => Number.isFinite(v) && v > 0) ||
      !Array.isArray(heights) || heights.length < 3 || !heights.every(Number.isFinite) ||
      Math.abs(period[1] / period[0] - Math.sqrt(3)) > 1e-6) return AU_LATTICE;
  const top = heights.at(-1), layerHeights = heights.slice(-3).reverse().map(v => v - top);
  // Physics enumerates six layers bottom-to-top: top k=5 has shift (d,h/3).
  // Rendering enumerates top-to-bottom, so the ABC shift runs in reverse.
  return Object.freeze({...AU_LATTICE, spacing: period[0],
    rowHeight: period[1] / 2, layerSpacing: -layerHeights[1], layerHeights,
    phaseX: period[0], phaseZ: period[1] / 3, stackDirection: -1,
    matchedPhysicalLattice: true});
}

export function latticePosition(column, row, layer = 0, lattice = AU_LATTICE) {
  const {spacing, rowHeight, layerSpacing} = lattice;
  const shift = (lattice.stackDirection ?? 1) * layer;
  // a1=(d,0), a2=(d/2,sqrt(3)d/2), and ABC shifts=(a1+a2)/3.
  return [spacing * (column + row / 2 + shift / 2) + (lattice.phaseX || 0),
    lattice.layerHeights?.[layer] ?? -layer * layerSpacing,
    rowHeight * (row + shift / 3) + (lattice.phaseZ || 0)];
}

// Upper envelope of three opaque sphere layers, including both triangular hollows.
export function sampleGoldCell(x, z, lattice = AU_LATTICE) {
  const {spacing, rowHeight, radius} = lattice;
  let height = -Infinity, normal = [0, 1, 0], layerFound = -1;
  for (let layer = 0; layer < 3; layer++) {
    const shift = (lattice.stackDirection ?? 1) * layer;
    const rowNear = Math.round((z - (lattice.phaseZ || 0)) / rowHeight - shift / 3);
    for (let row = rowNear - 1; row <= rowNear + 1; row++) {
      const columnNear = Math.round((x - (lattice.phaseX || 0)) / spacing - row / 2 - shift / 2);
      for (let column = columnNear - 1; column <= columnNear + 1; column++) {
        const centre = latticePosition(column, row, layer, lattice);
        const dx = x - centre[0], dz = z - centre[2];
        const squared = radius * radius - dx * dx - dz * dz;
        if (squared < 0) continue;
        const cap = Math.sqrt(squared), candidate = centre[1] + cap;
        if (candidate > height) {
          height = candidate; normal = [dx / radius, cap / radius, dz / radius]; layerFound = layer;
        }
      }
    }
  }
  return {height, normal, layer: layerFound};
}

function makeReliefTexture(renderer, lattice) {
  const cellHeight = lattice.rowHeight * 2;
  const minHeight = lattice.layerHeights?.[2] ?? -lattice.layerSpacing * 2;
  const maxHeight = lattice.radius;
  const pixels = new Uint8Array(TEXTURE_SIZE * TEXTURE_SIZE * 4);
  let sumHeight = 0, minimumHeight = Infinity;
  for (let z = 0; z < TEXTURE_SIZE; z++) for (let x = 0; x < TEXTURE_SIZE; x++) {
    const sample = sampleGoldCell((x + .5) / TEXTURE_SIZE * lattice.spacing,
      (z + .5) / TEXTURE_SIZE * cellHeight, lattice);
    sumHeight += sample.height;
    minimumHeight = Math.min(minimumHeight, sample.height);
    const offset = (z * TEXTURE_SIZE + x) * 4;
    for (let axis = 0; axis < 3; axis++) pixels[offset + axis] = Math.round((sample.normal[axis] * .5 + .5) * 255);
    pixels[offset + 3] = Math.round((sample.height - minHeight) / (maxHeight - minHeight) * 255);
  }
  const texture = new THREE.DataTexture(pixels, TEXTURE_SIZE, TEXTURE_SIZE, THREE.RGBAFormat);
  texture.wrapS = texture.wrapT = THREE.RepeatWrapping;
  texture.generateMipmaps = true;
  texture.minFilter = THREE.LinearMipmapLinearFilter;
  texture.magFilter = THREE.LinearFilter;
  texture.anisotropy = Math.min(8, renderer?.capabilities?.getMaxAnisotropy?.() || 1);
  texture.needsUpdate = true;
  texture.userData.averageHeight = sumHeight / (TEXTURE_SIZE * TEXTURE_SIZE);
  texture.userData.minimumHeight = minimumHeight;
  return texture;
}

const worldVertex = `
  vec4 goldVertex = vec4(transformed, 1.0);
  #ifdef USE_INSTANCING
    goldVertex = instanceMatrix * goldVertex;
  #endif
  goldWorldPosition = (modelMatrix * goldVertex).xyz;
`;

function worldShader(shader) {
  shader.vertexShader = 'varying vec3 goldWorldPosition;\n' + shader.vertexShader;
  shader.vertexShader = shader.vertexShader.replace('#include <project_vertex>', '#include <project_vertex>\n' + worldVertex);
  shader.fragmentShader = 'varying vec3 goldWorldPosition;\n' + shader.fragmentShader;
}

export class GoldSurface {
  constructor(track, renderer = null, lattice = AU_LATTICE) {
    this.lattice = lattice;
    const CELL_HEIGHT = lattice.rowHeight * 2;
    const MIN_HEIGHT = lattice.layerHeights?.[2] ?? -lattice.layerSpacing * 2;
    const MAX_HEIGHT = lattice.radius;
    const CARRIER_HEIGHT = MAX_HEIGHT + .04;
    this.track = {...track};
    this.bounds = {minX: track.start_x - 12, maxX: track.finish_x + 12,
      minZ: -track.half_width, maxZ: track.half_width};
    this.group = new THREE.Group(); this.group.name = 'Au(111) visual crystal';
    this.reliefTexture = makeReliefTexture(renderer, lattice);
    // A real top-level plane supports the far view. It must not fall back to
    // the bottom of the slab when detailed sphere instances are culled.
    // Leave room below the lowest cap for the medium spheres' polygon edges.
    const FALLBACK_HEIGHT = this.reliefTexture.userData.minimumHeight - .3;
    this.uniforms = {
      goldRelief: {value: this.reliefTexture},
      goldPatch: {value: new THREE.Vector3(0, 0, 0)},
      goldSurfaceZ: {value: track.surface_z},
      goldAverageHeight: {value: this.reliefTexture.userData.averageHeight},
      goldProjection: {value: new THREE.Matrix4()},
      goldBounds: {value: new THREE.Vector4(this.bounds.minX, this.bounds.maxX, this.bounds.minZ, this.bounds.maxZ)},
    };
    this.patch = null; this.rebuildCount = 0;
    // A quiet, matte warm gold keeps the molecule in the foreground.
    // Use the same material at every distance, without studio reflection panels.
    const properties = {color: 0xab945d, metalness: .05, roughness: .95};
    const edgeShader = `
      uniform vec4 goldBounds;
      float goldEdgeDistance(vec2 p) {
        return min(min(p.x - goldBounds.x, goldBounds.y - p.x),
                   min(p.y - goldBounds.z, goldBounds.w - p.y));
      }
    `;
    this.atomMaterial = new THREE.MeshStandardMaterial(properties);
    this.atomMaterial.onBeforeCompile = shader => {
      worldShader(shader);
      shader.uniforms.goldPatch = this.uniforms.goldPatch;
      shader.uniforms.goldBounds = this.uniforms.goldBounds;
      shader.fragmentShader = edgeShader + 'uniform vec3 goldPatch;\n' + shader.fragmentShader;
      shader.fragmentShader = shader.fragmentShader.replace('#include <clipping_planes_fragment>', `
        #include <clipping_planes_fragment>
        if (distance(goldWorldPosition.xz, goldPatch.xy) > goldPatch.z) discard;
        if (goldEdgeDistance(goldWorldPosition.xz) <= ${EDGE_WIDTH.toFixed(8)}) discard;
      `);
    };
    this.atomMaterial.customProgramCacheKey = () => 'au111-interior-space-filling-v3';
    this.geometries = [new THREE.SphereGeometry(lattice.radius, 18, 12),
      new THREE.SphereGeometry(lattice.radius, 12, 8), new THREE.SphereGeometry(lattice.radius, 8, 6)];
    this.meshes = this.geometries.map((geometry, index) => {
      const mesh = new THREE.InstancedMesh(geometry, this.atomMaterial, index === 2 ? CAPACITY * 2 : CAPACITY);
      mesh.name = ['Au top detailed', 'Au top medium', 'Au B+C layers'][index];
      mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
      mesh.count = 0; mesh.visible = false; mesh.receiveShadow = false;
      this.group.add(mesh); return mesh;
    });
    // All four borders retain actual ABC spheres at every camera distance.
    // Complementary fragment clipping joins this strip to interior geometry
    // without duplicate surfaces or isolated atoms outside the slab.
    this.edgeMaterial = new THREE.MeshStandardMaterial(properties);
    this.edgeMaterial.onBeforeCompile = shader => {
      worldShader(shader);
      shader.uniforms.goldBounds = this.uniforms.goldBounds;
      shader.fragmentShader = edgeShader + shader.fragmentShader;
      shader.fragmentShader = shader.fragmentShader.replace('#include <clipping_planes_fragment>', `
        #include <clipping_planes_fragment>
        if (goldEdgeDistance(goldWorldPosition.xz) > ${EDGE_WIDTH.toFixed(8)}) discard;
      `);
    };
    this.edgeMaterial.customProgramCacheKey = () => 'au111-permanent-edges-v1';
    const edgeSites = [[], []];
    const edgeDistance = (x, z) => Math.min(x - this.bounds.minX, this.bounds.maxX - x,
      z - this.bounds.minZ, this.bounds.maxZ - z);
    for (let layer = 0; layer < 3; layer++) {
      const shift = (lattice.stackDirection ?? 1) * layer;
      const minRow = Math.ceil((this.bounds.minZ - (lattice.phaseZ || 0)) / lattice.rowHeight - shift / 3);
      const maxRow = Math.floor((this.bounds.maxZ - (lattice.phaseZ || 0)) / lattice.rowHeight - shift / 3);
      for (let row = minRow; row <= maxRow; row++) {
        const minCol = Math.ceil((this.bounds.minX - (lattice.phaseX || 0)) / lattice.spacing - row / 2 - shift / 2);
        const maxCol = Math.floor((this.bounds.maxX - (lattice.phaseX || 0)) / lattice.spacing - row / 2 - shift / 2);
        for (let column = minCol; column <= maxCol; column++) {
          const p = latticePosition(column, row, layer, lattice);
          if (edgeDistance(p[0], p[2]) <= EDGE_WIDTH + lattice.radius + .05) edgeSites[layer ? 1 : 0].push(p);
        }
      }
    }
    this.edgeMeshes = edgeSites.map((sites, index) => {
      const mesh = new THREE.InstancedMesh(this.geometries[index + 1], this.edgeMaterial, sites.length);
      mesh.name = index ? 'Au permanent ABC sides' : 'Au permanent top rim';
      const matrix = new THREE.Matrix4();
      sites.forEach((p, i) => { matrix.makeTranslation(...p); mesh.setMatrixAt(i, matrix); });
      mesh.instanceMatrix.needsUpdate = true;
      mesh.computeBoundingBox(); mesh.computeBoundingSphere(); this.group.add(mesh);
      return mesh;
    });
    this.farMaterial = new THREE.MeshStandardMaterial(properties);
    this.farMaterial.onBeforeCompile = shader => {
      worldShader(shader); Object.assign(shader.uniforms, this.uniforms);
      shader.fragmentShader = `
        uniform sampler2D goldRelief;
        uniform vec3 goldPatch;
        uniform float goldSurfaceZ;
        uniform float goldAverageHeight;
        uniform mat4 goldProjection;
        vec2 goldUv(vec2 p) { return p / vec2(${lattice.spacing.toFixed(8)}, ${CELL_HEIGHT.toFixed(8)}); }
        float goldHeight(vec4 sampleValue) { return ${MIN_HEIGHT.toFixed(8)} + sampleValue.a * ${(MAX_HEIGHT - MIN_HEIGHT).toFixed(8)}; }
      ` + edgeShader + shader.fragmentShader;
      shader.fragmentShader = shader.fragmentShader.replace('#include <map_fragment>', `
        #include <map_fragment>
        // Parallax relief is a sphere-normal/height crystal cell, never a dot pattern.
        // Explicit texture gradients provide stable mip filtering even inside the search.
        vec2 goldDx = dFdx(goldUv(goldWorldPosition.xz));
        vec2 goldDz = dFdy(goldUv(goldWorldPosition.xz));
        // Below a few pixels per atom, suppress nonlinear specular/height aliasing.
        // Mip-filtered normals alone still produce large concentric moiré in overview.
        float goldFootprint = max(length(goldDx * vec2(1.0, 2.0)), length(goldDz * vec2(1.0, 2.0)));
        float goldDetail = 1.0 - smoothstep(.18, .50, goldFootprint);
        vec3 goldView = normalize(cameraPosition - goldWorldPosition);
        // Start at the actual top layer, then march down to its envelope.
        // Both geometry and fragment depth stay near the visible surface.
        float goldTravel = ${(CARRIER_HEIGHT - FALLBACK_HEIGHT).toFixed(8)} / max(.045, goldView.y);
        float goldUpper = 0.0, goldLower = goldTravel;
        vec3 goldPoint = goldWorldPosition;
        bool goldFound = false;
        for (int goldStep = 0; goldStep <= 18; goldStep++) {
          float travel = goldTravel * float(goldStep) / 18.0;
          vec3 p = goldWorldPosition - goldView * travel;
          float h = mix(goldAverageHeight, goldHeight(textureGrad(goldRelief, goldUv(p.xz), goldDx, goldDz)), goldDetail) + goldSurfaceZ;
          if (!goldFound && p.y <= h) { goldLower = travel; goldFound = true; }
          if (!goldFound) goldUpper = travel;
        }
        for (int refine = 0; refine < 5; refine++) {
          float travel = (goldUpper + goldLower) * .5;
          vec3 p = goldWorldPosition - goldView * travel;
          float h = mix(goldAverageHeight, goldHeight(textureGrad(goldRelief, goldUv(p.xz), goldDx, goldDz)), goldDetail) + goldSurfaceZ;
          if (p.y > h) goldUpper = travel; else goldLower = travel;
        }
        goldPoint = goldWorldPosition - goldView * ((goldUpper + goldLower) * .5);
        if (goldPoint.x < goldBounds.x || goldPoint.x > goldBounds.y ||
            goldPoint.z < goldBounds.z || goldPoint.z > goldBounds.w) discard;
        if (goldEdgeDistance(goldPoint.xz) <= ${EDGE_WIDTH.toFixed(8)}) discard;
        if (goldPatch.z > 0.0 && distance(goldPoint.xz, goldPatch.xy) <= goldPatch.z) discard;
        vec4 goldSample = textureGrad(goldRelief, goldUv(goldPoint.xz), goldDx, goldDz);
        float goldAverageAlpha = (goldAverageHeight - ${MIN_HEIGHT.toFixed(8)}) / ${(MAX_HEIGHT - MIN_HEIGHT).toFixed(8)};
        diffuseColor.rgb *= mix(.76, 1.0, mix(goldAverageAlpha, goldSample.a, goldDetail));
        vec4 goldClip = goldProjection * viewMatrix * vec4(goldPoint, 1.0);
        gl_FragDepth = goldClip.z / goldClip.w * .5 + .5;
      `);
      shader.fragmentShader = shader.fragmentShader.replace('#include <normal_fragment_maps>', `
        #include <normal_fragment_maps>
        vec3 goldNormal = normalize(mix(vec3(0.0, 1.0, 0.0), normalize(goldSample.rgb * 2.0 - 1.0), goldDetail * goldDetail));
        normal = normalize(mat3(viewMatrix) * goldNormal);
      `);
    };
    this.farMaterial.customProgramCacheKey = () => `au111-top-relief-v4/${lattice.spacing}/${CELL_HEIGHT}/${MIN_HEIGHT}/${MAX_HEIGHT}/${FALLBACK_HEIGHT}`;
    const width = this.bounds.maxX - this.bounds.minX, depth = this.bounds.maxZ - this.bounds.minZ;
    // Overscan permits relief rays at oblique viewing angles to reach the field edge.
    this.farPlane = new THREE.Mesh(new THREE.PlaneGeometry(width + 180, depth + 180), this.farMaterial);
    this.farPlane.rotation.x = -Math.PI / 2;
    this.farPlane.position.set((this.bounds.minX + this.bounds.maxX) / 2, CARRIER_HEIGHT, 0);
    // Shader relief changes the depth slightly: do not cull it using a flat
    // carrier's bounds. The full field has only two triangles.
    this.farPlane.frustumCulled = false;
    this.farPlane.receiveShadow = false; this.group.add(this.farPlane);
    // A closed interior core supports the top at every LOD. Inset its sides
    // behind the permanent rim, so it cannot cut off the exposed ABC layers.
    const bottom = MIN_HEIGHT - lattice.radius;
    this.backing = new THREE.Mesh(new THREE.BoxGeometry(Math.max(.1, width - EDGE_WIDTH * 2),
      FALLBACK_HEIGHT - bottom, Math.max(.1, depth - EDGE_WIDTH * 2)), new THREE.MeshStandardMaterial(properties));
    this.backing.name = 'Au closed interior core';
    this.backing.position.set(this.farPlane.position.x, (FALLBACK_HEIGHT + bottom) / 2, 0);
    this.fallbackHeight = FALLBACK_HEIGHT;
    this.bottomHeight = bottom;
    this.group.add(this.backing);
  }

  update(camera, target = null) {
    const lattice = this.lattice;
    this.uniforms.goldProjection.value.copy(camera.projectionMatrix);
    const height = Math.max(0, camera.position.y - this.track.surface_z);
    const rawRadius = height >= lattice.maxDetailHeight ? 0 :
      Math.min(lattice.patchRadius, Math.sqrt(lattice.maxDetailHeight ** 2 - height ** 2));
    const radius = rawRadius < 8 ? 0 : Math.min(lattice.patchRadius, Math.ceil(rawRadius / 4) * 4);
    const tileX = lattice.spacing * 6, tileZ = lattice.rowHeight * 6;
    // Snap only the allocation window. Every matrix stays at its absolute lattice site.
    // Cover the visible surface around the viewing ray, rather than a patch
    // beneath an off-track camera. Keep the allocation on the actual slab.
    const focus = camera.position.clone();
    if (target && camera.position.y > target.y + .001) {
      const t = (camera.position.y - this.track.surface_z - this.reliefTexture.userData.averageHeight) /
        (camera.position.y - target.y);
      focus.lerp(target, Math.max(0, t));
    }
    const x = THREE.MathUtils.clamp(Math.round(focus.x / tileX) * tileX, this.bounds.minX, this.bounds.maxX);
    const z = THREE.MathUtils.clamp(Math.round(focus.z / tileZ) * tileZ, this.bounds.minZ, this.bounds.maxZ);
    const detailBand = Math.round(height / 8);
    const key = `${x}/${z}/${radius}/${detailBand}`;
    if (this.patch?.key === key) return;
    this.patch = {x, z, radius, key}; this.rebuildCount++;
    this.uniforms.goldPatch.value.set(x, z, radius);
    const counts = [0, 0, 0], matrix = new THREE.Matrix4();
    const overscan = radius + lattice.radius + .05;
    if (radius > 0) {
      for (let layer = 0; layer < 3; layer++) {
        const shift = (lattice.stackDirection ?? 1) * layer;
        const minRow = Math.floor((z - overscan - (lattice.phaseZ || 0)) / lattice.rowHeight - shift / 3);
        const maxRow = Math.ceil((z + overscan - (lattice.phaseZ || 0)) / lattice.rowHeight - shift / 3);
        for (let row = minRow; row <= maxRow; row++) {
          const minCol = Math.floor((x - overscan - (lattice.phaseX || 0)) / lattice.spacing - row / 2 - shift / 2);
          const maxCol = Math.ceil((x + overscan - (lattice.phaseX || 0)) / lattice.spacing - row / 2 - shift / 2);
          for (let col = minCol; col <= maxCol; col++) {
            const p = latticePosition(col, row, layer, lattice);
            if (Math.hypot(p[0] - x, p[2] - z) > overscan ||
              p[0] < this.bounds.minX || p[0] > this.bounds.maxX ||
              p[2] < this.bounds.minZ || p[2] > this.bounds.maxZ ||
              Math.min(p[0] - this.bounds.minX, this.bounds.maxX - p[0],
                p[2] - this.bounds.minZ, this.bounds.maxZ - p[2]) < EDGE_WIDTH - lattice.radius) continue;
            const distance = camera.position.distanceTo(new THREE.Vector3(p[0], p[1] + this.track.surface_z, p[2]));
            const level = layer > 0 ? 2 : distance < 45 ? 0 : 1;
            if (counts[level] >= this.meshes[level].instanceMatrix.count) throw new Error('Au visual LOD capacity exceeded');
            matrix.makeTranslation(...p); this.meshes[level].setMatrixAt(counts[level]++, matrix);
          }
        }
      }
    }
    this.meshes.forEach((mesh, i) => {
      mesh.count = counts[i]; mesh.visible = counts[i] > 0;
      mesh.instanceMatrix.needsUpdate = true;
      if (mesh.visible) { mesh.computeBoundingBox(); mesh.computeBoundingSphere(); }
    });
  }

  getSnapshot() {
    const lattice = this.lattice;
    const matrix = new THREE.Matrix4();
    const examples = this.meshes.map(mesh => {
      const sample = [];
      for (let i = 0; i < Math.min(5, mesh.count); i++) {
        mesh.getMatrixAt(i, matrix);
        sample.push({centre: [matrix.elements[12], matrix.elements[14], matrix.elements[13] + this.track.surface_z],
          scale: matrix.getMaxScaleOnAxis()});
      }
      return {name: mesh.name, count: mesh.count, trianglesPerAtom: mesh.geometry.index.count / 3, instances: sample};
    });
    return {planeConvention: 'surface_z is top-layer Au atom-centre plane (Å)',
      representation: 'space-filling', radiusSource: 'RDKit GetRvdw(Au)',
      appearance: {finish: 'matte', color: `#${this.atomMaterial.color.getHexString()}`,
        roughness: this.atomMaterial.roughness, metalness: this.atomMaterial.metalness},
      spacing: lattice.spacing, radius: lattice.radius, layerSpacing: lattice.layerSpacing,
      phaseXY: [lattice.phaseX || 0, lattice.phaseZ || 0], stackDirection: lattice.stackDirection ?? 1,
      matchedPhysicalLattice: Boolean(lattice.matchedPhysicalLattice),
      centreHeights: (lattice.layerHeights || [0, -2.35, -4.7]).map(v => v + this.track.surface_z),
      topHeight: this.track.surface_z + lattice.radius, bounds: {...this.bounds},
      farCarrierHeight: this.track.surface_z + this.farPlane.position.y,
      fallbackHeight: this.track.surface_z + this.fallbackHeight,
      bottomHeight: this.track.surface_z + this.bottomHeight,
      edgeWidth: EDGE_WIDTH,
      permanentEdges: this.edgeMeshes.map(mesh => ({name: mesh.name, count: mesh.count,
        triangles: mesh.count * mesh.geometry.index.count / 3})),
      patch: this.patch ? {...this.patch} : null, rebuildCount: this.rebuildCount,
      fullFieldRelief: true, opaque: !this.atomMaterial.transparent && this.atomMaterial.opacity === 1,
      reliefTextureSize: [TEXTURE_SIZE, TEXTURE_SIZE], meshes: examples,
      detailTriangles: examples.reduce((sum, mesh) => sum + mesh.count * mesh.trianglesPerAtom, 0),
      edgeTriangles: this.edgeMeshes.reduce((sum, mesh) => sum + mesh.count * mesh.geometry.index.count / 3, 0),
      triangles: [...this.meshes, ...this.edgeMeshes].reduce((sum, mesh) => sum + mesh.count * mesh.geometry.index.count / 3, 0),
    };
  }

  dispose() {
    this.group.removeFromParent();
    [...this.meshes, ...this.edgeMeshes].forEach(mesh => mesh.dispose());
    [...this.geometries, this.farPlane.geometry, this.backing.geometry].forEach(geometry => geometry.dispose());
    [this.atomMaterial, this.edgeMaterial, this.farMaterial, this.backing.material].forEach(material => material.dispose());
    this.reliefTexture.dispose(); this.group.clear();
  }
}
