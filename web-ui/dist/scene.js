import * as THREE from './assets/three.module.js';

const MODEL_URL = new URL('./assets/molecule.json', import.meta.url);
const FALLBACK_URL = new URL('./assets/molecule-fallback.svg', import.meta.url);

/** A presentation scene built from an actual 52-atom molecular structure.
 * This scene never advances molecular dynamics or changes the molecule.
 */
export function mountMoleculeScene(container, { variant = 'hero', reducedMotion = false } = {}) {
  if (!(container instanceof HTMLElement)) throw new TypeError('A scene container is required.');
  const media = window.matchMedia('(prefers-reduced-motion: reduce)');
  const abort = new AbortController();
  const disposables = new Set();
  let disposed = false;
  let representation = 'ball-stick';
  let motionEnabled = !reducedMotion;
  let inViewport = true;
  let raf = 0;
  let activeSeconds = 0;
  let lastTime = 0;
  let renderer;
  let scene;
  let camera;
  let world;
  let molecule;
  let bondMesh;
  let carbonMesh;
  let hydrogenMesh;
  let atomRecords;
  let environmentTarget;
  let shadowTexture;
  let ready = false;

  const fallback = document.createElement('img');
  fallback.src = FALLBACK_URL.href;
  fallback.alt = '四苯基苯分子赛车，置于金原子密堆表面上方';
  fallback.className = 'molecule-scene-fallback';
  Object.assign(fallback.style, {
    width: '100%', height: '100%', objectFit: 'contain', display: 'block',
    position: 'absolute', inset: '0', pointerEvents: 'none',
  });
  container.append(fallback);

  function own(resource) {
    disposables.add(resource);
    return resource;
  }

  function resize() {
    if (!ready || disposed) return;
    const width = Math.max(1, container.clientWidth);
    const height = Math.max(1, container.clientHeight);
    renderer.setSize(width, height, false);
    camera.aspect = width / height;
    // Keep the car readable in narrow layouts rather than clipping its wheels.
    camera.fov = camera.aspect < 1 ? 36 / Math.max(.62, camera.aspect) : 36;
    camera.updateProjectionMatrix();
    renderFrame();
  }

  function renderFrame() {
    if (!ready || disposed || !renderer) return;
    renderer.render(scene, camera);
  }

  function wantsAnimation() {
    return ready && !disposed && motionEnabled && !reducedMotion && !media.matches && inViewport && !document.hidden;
  }

  function synchronizeAnimation() {
    if (!wantsAnimation()) {
      if (raf) cancelAnimationFrame(raf);
      raf = 0;
      lastTime = 0;
      return;
    }
    if (!raf) raf = requestAnimationFrame(animate);
  }

  function animate(time) {
    raf = 0;
    if (!wantsAnimation()) return;
    if (lastTime) activeSeconds += Math.min(.05, (time - lastTime) / 1000);
    lastTime = time;
    // A slow showroom movement, not a fabricated molecular trajectory.
    world.rotation.y = -.36 + Math.sin(activeSeconds * .13) * (variant === 'garage' ? .18 : .12);
    camera.position.x = .65 + Math.sin(activeSeconds * .09) * .35;
    camera.lookAt(0, 3.25, 0);
    renderFrame();
    raf = requestAnimationFrame(animate);
  }

  function createEnvironment() {
    const room = new THREE.Scene();
    room.background = new THREE.Color(variant === 'hero' ? 0xc8c2b8 : 0x333638);
    const temporary = [];
    const panel = (color, intensity, width, height, position, target) => {
      const geometry = new THREE.PlaneGeometry(width, height);
      const material = new THREE.MeshBasicMaterial({ color, side: THREE.DoubleSide });
      material.color.multiplyScalar(intensity);
      const mesh = new THREE.Mesh(geometry, material);
      mesh.position.set(...position);
      mesh.lookAt(...target);
      room.add(mesh);
      temporary.push(geometry, material);
    };
    panel(0xf3f4ed, 3.3, 18, 13, [0, 14, 2], [0, 0, 0]);
    panel(0xffffff, 2.2, 4, 15, [-11, 6, 5], [0, 3, 0]);
    panel(0xf1ca80, 2.1, 7, 13, [11, 6, -4], [0, 3, 0]);
    panel(0xd6e2ea, 1.7, 10, 7, [0, 7, -12], [0, 3, 0]);
    const generator = new THREE.PMREMGenerator(renderer);
    environmentTarget = generator.fromScene(room, .03);
    generator.dispose();
    temporary.forEach(resource => resource.dispose());
    scene.environment = environmentTarget.texture;
  }

  function createSurface() {
    const radius = 1.405;
    const spacing = 2.885;
    const rowSpacing = spacing * Math.sqrt(3) / 2;
    const positions = [];
    for (let row = -6; row <= 6; row++) {
      for (let column = -7; column <= 7; column++) {
        const x = column * spacing + (Math.abs(row) % 2 ? spacing / 2 : 0);
        const z = row * rowSpacing;
        positions.push([x, -.7, z]);
      }
    }
    const geometry = own(new THREE.SphereGeometry(radius, 22, 14));
    const material = own(new THREE.MeshStandardMaterial({
      color: variant === 'hero' ? 0xb88c47 : 0x927240,
      metalness: .94, roughness: .38, envMapIntensity: variant === 'hero' ? .92 : .66,
    }));
    const surface = new THREE.InstancedMesh(geometry, material, positions.length);
    const transform = new THREE.Object3D();
    const color = new THREE.Color();
    positions.forEach((position, index) => {
      transform.position.set(...position);
      transform.updateMatrix();
      surface.setMatrixAt(index, transform.matrix);
      const distance = Math.hypot(position[0] * .75, position[2]);
      const lightness = Math.max(variant === 'hero' ? .65 : .35, 1 - distance / 42);
      color.setRGB(lightness, lightness, lightness);
      surface.setColorAt(index, color);
    });
    surface.instanceMatrix.needsUpdate = true;
    surface.instanceColor.needsUpdate = true;
    world.add(surface);

    const canvas = document.createElement('canvas');
    canvas.width = canvas.height = 128;
    const ctx = canvas.getContext('2d');
    const gradient = ctx.createRadialGradient(64, 64, 5, 64, 64, 64);
    gradient.addColorStop(0, 'rgba(0,0,0,.56)');
    gradient.addColorStop(.45, 'rgba(0,0,0,.22)');
    gradient.addColorStop(1, 'rgba(0,0,0,0)');
    ctx.fillStyle = gradient;
    ctx.fillRect(0, 0, 128, 128);
    shadowTexture = own(new THREE.CanvasTexture(canvas));
    const shadow = new THREE.Mesh(
      own(new THREE.PlaneGeometry(14, 11)),
      own(new THREE.MeshBasicMaterial({ map: shadowTexture, transparent: true, depthWrite: false, opacity: .66 })),
    );
    shadow.rotation.x = -Math.PI / 2;
    shadow.position.y = .82;
    world.add(shadow);
  }

  function createMolecule(data) {
    const center = data.atoms.reduce((sum, atom) => sum.add(new THREE.Vector3(...atom.position)), new THREE.Vector3()).divideScalar(data.atoms.length);
    atomRecords = data.atoms.map(atom => ({
      element: atom.element,
      position: new THREE.Vector3(atom.position[0] - center.x, atom.position[2] - center.z, atom.position[1] - center.y),
    }));
    molecule = new THREE.Group();
    molecule.position.y = 4.95;
    world.add(molecule);
    const sphere = own(new THREE.SphereGeometry(1, 30, 20));
    const carbonMaterial = own(new THREE.MeshPhysicalMaterial({
      color: 0x333b41, metalness: .62, roughness: .26,
      clearcoat: .5, clearcoatRoughness: .19, envMapIntensity: 1.32,
    }));
    const hydrogenMaterial = own(new THREE.MeshPhysicalMaterial({
      color: 0xe8e2d8, metalness: .18, roughness: .3,
      clearcoat: .35, envMapIntensity: .8,
    }));
    const carbons = atomRecords.filter(atom => atom.element === 'C');
    const hydrogens = atomRecords.filter(atom => atom.element === 'H');
    carbonMesh = new THREE.InstancedMesh(sphere, carbonMaterial, carbons.length);
    hydrogenMesh = new THREE.InstancedMesh(sphere, hydrogenMaterial, hydrogens.length);
    carbonMesh.userData.atoms = carbons;
    hydrogenMesh.userData.atoms = hydrogens;
    molecule.add(carbonMesh, hydrogenMesh);

    const bondGeometry = own(new THREE.CylinderGeometry(1, 1, 1, 14, 1));
    const bondMaterial = own(new THREE.MeshStandardMaterial({
      color: 0x9da09c, metalness: .78, roughness: .27, envMapIntensity: 1.15,
    }));
    bondMesh = new THREE.InstancedMesh(bondGeometry, bondMaterial, data.bonds.length);
    const transform = new THREE.Object3D();
    const up = new THREE.Vector3(0, 1, 0);
    data.bonds.forEach((bond, index) => {
      const a = atomRecords[bond.a].position;
      const b = atomRecords[bond.b].position;
      const delta = new THREE.Vector3().subVectors(b, a);
      transform.position.copy(a).add(b).multiplyScalar(.5);
      transform.quaternion.setFromUnitVectors(up, delta.clone().normalize());
      transform.scale.set(.105, delta.length(), .105);
      transform.updateMatrix();
      bondMesh.setMatrixAt(index, transform.matrix);
    });
    bondMesh.instanceMatrix.needsUpdate = true;
    molecule.add(bondMesh);
    applyRepresentation();
  }

  function applyRepresentation() {
    if (!carbonMesh || !hydrogenMesh) return;
    const transform = new THREE.Object3D();
    [carbonMesh, hydrogenMesh].forEach(mesh => {
      const isCarbon = mesh === carbonMesh;
      const radius = representation === 'space-fill' ? (isCarbon ? 1.7 : 1.2) : (isCarbon ? .49 : .285);
      mesh.userData.atoms.forEach((atom, index) => {
        transform.position.copy(atom.position);
        transform.scale.setScalar(radius);
        transform.updateMatrix();
        mesh.setMatrixAt(index, transform.matrix);
      });
      mesh.instanceMatrix.needsUpdate = true;
      mesh.computeBoundingSphere();
    });
    bondMesh.visible = representation === 'ball-stick';
    renderFrame();
  }

  async function initialize() {
    try {
      const response = await fetch(MODEL_URL, { signal: abort.signal });
      if (!response.ok) throw new Error('Molecule asset could not be loaded.');
      const data = await response.json();
      if (disposed) return;
      renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, powerPreference: 'low-power' });
      renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 1.7));
      renderer.setClearColor(0x000000, 0);
      renderer.outputColorSpace = THREE.SRGBColorSpace;
      renderer.toneMapping = THREE.ACESFilmicToneMapping;
      renderer.toneMappingExposure = variant === 'hero' ? 1.12 : 1.05;
      renderer.domElement.className = 'molecule-scene-canvas';
      renderer.domElement.setAttribute('role', 'img');
      renderer.domElement.setAttribute('aria-label', '四苯基苯分子赛车与 Au(111) 金表面，三维展台预览');
      Object.assign(renderer.domElement.style, { width: '100%', height: '100%', display: 'block', position: 'absolute', inset: '0', pointerEvents: 'none' });
      scene = new THREE.Scene();
      scene.fog = new THREE.FogExp2(variant === 'hero' ? 0xf3f1eb : 0x101214, variant === 'hero' ? .022 : .034);
      camera = new THREE.PerspectiveCamera(36, 1, .1, 150);
      camera.position.set(.65, 16.5, 22.5);
      camera.lookAt(0, 3.25, 0);
      world = new THREE.Group();
      world.rotation.y = -.36;
      scene.add(world);
      createEnvironment();
      scene.add(new THREE.HemisphereLight(0xf0f3ed, 0x2e2115, variant === 'hero' ? 2.2 : 1.55));
      const key = new THREE.DirectionalLight(0xf3f2e9, variant === 'hero' ? 3.3 : 2.7);
      key.position.set(-8, 14, 10);
      scene.add(key);
      const edge = new THREE.DirectionalLight(0xcda665, 2.3);
      edge.position.set(10, 7, -7);
      scene.add(edge);
      const fill = new THREE.DirectionalLight(0xc1ccd6, variant === 'hero' ? 1.2 : .85);
      fill.position.set(-4, 5, -12);
      scene.add(fill);
      createSurface();
      createMolecule(data);
      renderer.domElement.addEventListener('webglcontextlost', onContextLost);
      container.append(renderer.domElement);
      ready = true;
      resize();
      fallback.hidden = true;
      fallback.style.display = 'none';
      container.dataset.sceneState = 'webgl';
      synchronizeAnimation();
    } catch (error) {
      if (disposed || error.name === 'AbortError') return;
      container.dataset.sceneState = 'fallback';
      fallback.hidden = false;
      fallback.style.display = 'block';
      if (renderer) {
        renderer.domElement.remove();
        renderer.dispose();
      }
      ready = false;
    }
  }

  function onContextLost(event) {
    event.preventDefault();
    ready = false;
    synchronizeAnimation();
    renderer.domElement.style.display = 'none';
    fallback.hidden = false;
    fallback.style.display = 'block';
    container.dataset.sceneState = 'fallback';
  }

  const resizeObserver = new ResizeObserver(resize);
  resizeObserver.observe(container);
  const intersectionObserver = new IntersectionObserver(entries => {
    inViewport = entries[0]?.isIntersecting ?? true;
    synchronizeAnimation();
  }, { rootMargin: '30px' });
  intersectionObserver.observe(container);
  document.addEventListener('visibilitychange', synchronizeAnimation);
  media.addEventListener('change', synchronizeAnimation);
  initialize();

  return {
    setRepresentation(value) {
      if (!['ball-stick', 'space-fill'].includes(value)) return;
      representation = value;
      applyRepresentation();
    },
    setMotion(value) {
      motionEnabled = Boolean(value);
      synchronizeAnimation();
    },
    dispose() {
      if (disposed) return;
      disposed = true;
      abort.abort();
      if (raf) cancelAnimationFrame(raf);
      resizeObserver.disconnect();
      intersectionObserver.disconnect();
      document.removeEventListener('visibilitychange', synchronizeAnimation);
      media.removeEventListener('change', synchronizeAnimation);
      if (renderer) {
        renderer.domElement.removeEventListener('webglcontextlost', onContextLost);
        renderer.domElement.remove();
        renderer.dispose();
      }
      environmentTarget?.dispose();
      disposables.forEach(resource => resource.dispose());
      fallback.remove();
      delete container.dataset.sceneState;
    },
  };
}
