import {installRaceBridge} from './browser/embedded.mjs';
import {readPreferences} from './app-state.mjs';
import * as THREE from 'three';
import { OrbitControls } from './vendor/OrbitControls.js';
import { GoldSurface, AU_LATTICE, latticeForSurface } from './gold-surface.js';
import { request, downloadReplay, execution } from './transport.js';
import { raceCourse } from './race-course.mjs';
import {RaceResultDialog} from './race-result.mjs';
import {PlaybackClock} from './playback.mjs';
import { validateStartPose } from './start-pose.js';
import {RaceSetup, describeCompetition, isRaceInProgress} from './competition.mjs';

const $ = (id) => document.getElementById(id);
const raceResults = new RaceResultDialog({again:()=>action('reset'),save:async()=>{
  try{await saveReplay();}catch(error){feedback(error.message,true);}
}});
const playback = new PlaybackClock();
const sceneHost = $('scene');
const statuses = {idle: '等待分子', preparing: '准备分子…', paused: '准备就绪', running: '比赛中', finished: '已完赛', error: '出错'};
const elementStyle = {C: [0x333a32, .52], H: [0xf1f1df, .30], B: [0xb89578, .51], N: [0x6685ae, .49], O: [0xc77357, .49], S: [0xd2b95b, .61], Si: [0x9c9b82, .64], P: [0xb99c63, .58], F: [0x83ab7a, .45], Cl: [0x819e67, .58], Br: [0xa87b53, .63], I: [0x9382ad, .69], Au: [0xb99b45, .64]};
// Element van der Waals radii in Å, queried from RDKit 2025.09.3 GetRvdw.
// Visual appearance only: these are not the MD force field's atom-type parameters.
// API: https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.PeriodicTable.GetRvdw
// RDKit includes 2.0 Å defaults for some heavy elements lacking tabulated data.
const vdwRadii = Object.freeze({
  H: 1.2, He: 1.4, Li: 2.2, Be: 1.9, B: 1.8, C: 1.7, N: 1.6, O: 1.55, F: 1.5, Ne: 1.54,
  Na: 2.4, Mg: 2.2, Al: 2.1, Si: 2.1, P: 1.95, S: 1.8, Cl: 1.8, Ar: 1.88,
  K: 2.8, Ca: 2.4, Sc: 2.3, Ti: 2.15, V: 2.05, Cr: 2.05, Mn: 2.05, Fe: 2.05, Co: 2, Ni: 2, Cu: 2, Zn: 2.1, Ga: 2.1, Ge: 2.1, As: 2.05, Se: 1.9, Br: 1.9, Kr: 2.02,
  Rb: 2.9, Sr: 2.55, Y: 2.4, Zr: 2.3, Nb: 2.15, Mo: 2.1, Tc: 2.05, Ru: 2.05, Rh: 2, Pd: 2.05, Ag: 2.1, Cd: 2.2, In: 2.2, Sn: 2.25, Sb: 2.2, Te: 2.1, I: 2.1, Xe: 2.16,
  Cs: 3, Ba: 2.7, La: 2.5, Ce: 2.48, Pr: 2.47, Nd: 2.45, Pm: 2.43, Sm: 2.42, Eu: 2.4, Gd: 2.38, Tb: 2.37, Dy: 2.35, Ho: 2.33, Er: 2.32, Tm: 2.3, Yb: 2.28, Lu: 2.27,
  Hf: 2.25, Ta: 2.2, W: 2.1, Re: 2.05, Os: 2, Ir: 2, Pt: 2.05, Au: 2.1, Hg: 2.05, Tl: 2.2, Pb: 2.3, Bi: 2.3, Po: 2, At: 2, Rn: 2,
  Fr: 2, Ra: 2, Ac: 2, Th: 2.4, Pa: 2, U: 2.3, Np: 2, Pu: 2, Am: 2, Cm: 2, Bk: 2, Cf: 2, Es: 2, Fm: 2, Md: 2, No: 2, Lr: 2,
  Rf: 2, Db: 2, Sg: 2, Bh: 2, Hs: 2, Mt: 2, Ds: 2, Rg: 2, Cn: 2, Nh: 2, Fl: 2, Mc: 2, Lv: 2, Ts: 2, Og: 2,
});
let savedReplayGeneration = null;
async function saveReplay(){await downloadReplay();savedReplayGeneration=generation;dispatchEvent(new Event('race-state-change'));}
function unsavedReplay(){return Boolean(state?.frame?.step>0&&savedReplayGeneration!==generation);}
addEventListener('beforeunload',event=>{if(parent===window&&(raceLocked()||unsavedReplay())){event.preventDefault();event.returnValue='';}});
let representation = 'ball-stick', moleculeFormula = '';
let state = null, molecule = null, moleculeKey = '', latestSeq = -1, frames = [], generation = null;
let currentFrame = null, latestFrame = null, busy = false, connected = false, epoch = 0, loading = false;
let statePolling = false, framePolling = false, cameraMode = 'follow', needsCameraFocus = true;

let atomsMesh = null, bondsMesh = null, atomCount = 0, displayedPositions = null;
let chartSamples = [], lastChartSeq = -1, feedbackTimer = null;
let poseDraftDirty = false, poseFormGeneration = null, poseAtomSignature = '', posePicking = null, posePointerDown = null;
const poseControlIds = ['start-pose-mode', 'start-pose-heading', 'start-pose-rear', 'start-pose-front',
  'start-pose-yaw', 'start-pose-pitch', 'start-pose-roll', 'apply-start-pose', 'pick-pose-rear', 'pick-pose-front'];
const posePickBanner = document.createElement('div');
posePickBanner.className = 'start-pose-pick-banner'; posePickBanner.hidden = true; sceneHost.append(posePickBanner);
// Initial geometry before the authoritative server track arrives (all values Å).
const track = {start_x: 0, finish_x: 200, half_width: 100, surface_z: 0};
let driveControl = {mode:'gradient',strength:.6,growth:1};
let draftInput = {example:'tetraphenylbenzene.mol'}, draftChanged=false, ruleGeneration=null, preparedGeneration=null, raceStarted=false;
const raceSetup = new RaceSetup($('race-setup'), () => {draftChanged=true;updateRaceDraft();});
function updateRaceDraft(){
  $('race-draft-note').textContent = '';
  $('race-draft-note').hidden = true;
  $('start-pose-panel').hidden = raceSetup.mode !== 'custom';
}
let trackKey = '', runway = null, ground = null, goldSurface = null;

const scene = new THREE.Scene();
scene.background = new THREE.Color(0xe9e9dd);
scene.fog = new THREE.Fog(0xe9e9dd, 65, 1800);
const camera = new THREE.PerspectiveCamera(39, 1, .1, 5000);
camera.up.set(0, 1, 0);
const renderer = new THREE.WebGLRenderer({antialias: true, alpha: false});
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.shadowMap.enabled = false;
renderer.shadowMap.type = THREE.PCFSoftShadowMap;
renderer.outputColorSpace = THREE.SRGBColorSpace;
renderer.toneMapping = THREE.ACESFilmicToneMapping;
renderer.toneMappingExposure = 1.15;
sceneHost.appendChild(renderer.domElement);
const orbit = new OrbitControls(camera, renderer.domElement);
orbit.enableDamping = true;
orbit.dampingFactor = .075;
orbit.maxPolarAngle = Math.PI * .47;
orbit.minDistance = 13;
orbit.maxDistance = 7000;
orbit.enablePan = true;

scene.add(new THREE.HemisphereLight(0xfffff3, 0x718066, 2.5));
const sun = new THREE.DirectionalLight(0xfff7dc, 3.5);
sun.position.set(-15, 45, 24);
sun.castShadow = false;
sun.shadow.mapSize.set(2048, 2048);
sun.shadow.camera.left = -28;
sun.shadow.camera.right = 28;
sun.shadow.camera.top = 28;
sun.shadow.camera.bottom = -28;
sun.shadow.camera.near = 1;
sun.shadow.camera.far = 100;
sun.shadow.normalBias = .07;
scene.add(sun);
scene.add(sun.target);
const fill = new THREE.DirectionalLight(0xeef4e2, 1.5);
fill.position.set(50, 15, -20);
scene.add(fill);

const surface = new THREE.Group();
scene.add(surface);
const finish = new THREE.Group();
const gradientField = new THREE.Group();
scene.add(gradientField);
let gradientMarkers = [], gradientMaterial = null;

function disposeGroup(group) {
  const geometries = new Set(), materials = new Set();
  group.traverse(object => {
    if (object.isInstancedMesh) object.dispose();
    if (object.geometry) geometries.add(object.geometry);
    for (const material of Array.isArray(object.material) ? object.material : object.material ? [object.material] : []) materials.add(material);
  });
  group.clear();
  geometries.forEach(geometry => geometry.dispose());
  materials.forEach(material => material.dispose());
}

function rebuildTrack() {
  const course = raceCourse(track, state?.frame);
  const length = track.finish_x - track.start_x, width = track.half_width * 2;
  const centerX = (track.start_x + track.finish_x) / 2, margin = 12;
  const lattice = latticeForSurface(state?.runtime?.surface);
  const signature = JSON.stringify([track.start_x, track.finish_x, track.half_width, track.surface_z,
    lattice.spacing, lattice.radius, lattice.layerHeights, lattice.phaseX, lattice.phaseZ, lattice.stackDirection,
    course.start_x, course.finish_x]);
  if (signature === trackKey || !(length > 0 && width > 0)) return;
  trackKey = signature;
  goldSurface?.dispose();
  disposeGroup(surface);
  finish.clear();
  disposeGroup(gradientField);
  surface.position.y = track.surface_z;
  ground = new THREE.Mesh(new THREE.PlaneGeometry(length + width, width * 2), new THREE.MeshStandardMaterial({color: 0xdfe1cd, roughness: .95}));
  ground.rotation.x = -Math.PI / 2;
  ground.position.set(centerX, -7.6, 0);
  ground.receiveShadow = false;
  surface.add(ground);
  goldSurface = new GoldSurface(track, renderer, lattice);
  runway = goldSurface.backing;
  surface.add(goldSurface.group);
  // Track lines sit above the visual Au spheres; MD surface_z remains the centre plane.
  const markingHeight = lattice.radius + .025;
  const borderPoints = [];
  for (const side of [-track.half_width, track.half_width]) borderPoints.push(track.start_x - margin, markingHeight, side, track.finish_x + margin, markingHeight, side);
  borderPoints.push(course.start_x, markingHeight, -track.half_width, course.start_x, markingHeight, track.half_width);
  const borderGeometry = new THREE.BufferGeometry();
  borderGeometry.setAttribute('position', new THREE.Float32BufferAttribute(borderPoints, 3));
  surface.add(new THREE.LineSegments(borderGeometry, new THREE.LineBasicMaterial({color: 0xf6ecd0, transparent: true, opacity: .7})));
  const lane = new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(course.start_x, markingHeight, 0), new THREE.Vector3(course.finish_x, markingHeight, 0)]), new THREE.LineDashedMaterial({color: 0xf6e9be, dashSize: 4, gapSize: 8, transparent: true, opacity: .35}));
  lane.computeLineDistances();
  surface.add(lane);
  const stripeWidth = Math.max(1.3, Math.min(8, length * .008));
  const rows = Math.max(2, Math.min(64, Math.ceil(width / stripeWidth)));
  const checker = new THREE.InstancedMesh(new THREE.PlaneGeometry(1, 1), new THREE.MeshBasicMaterial({color: 0xffffff, transparent: true, opacity: .8}), rows * 2);
  const tile = new THREE.Object3D();
  tile.rotation.x = -Math.PI / 2;
  for (let col = 0; col < 2; col++) for (let row = 0; row < rows; row++) {
    tile.position.set((col - .5) * stripeWidth / 2, markingHeight + .01, -track.half_width + (row + .5) * width / rows);
    tile.scale.set(stripeWidth / 2, width / rows, 1);
    tile.updateMatrix();
    checker.setMatrixAt(col * rows + row, tile.matrix);
    checker.setColorAt(col * rows + row, new THREE.Color((col + row) % 2 ? 0x637052 : 0xf2eccf));
  }
  finish.position.set(course.finish_x, 0, 0);
  finish.add(checker);
  surface.add(finish);
  const arrowPoints = [], arrowLength = Math.max(4, Math.min(12, length * .015));
  gradientMarkers = [];
  for (let i = 0; i < 6; i++) for (const side of [-1, 1]) {
    const x = track.start_x + 14 + i * Math.max(0, length - 28) / 5;
    const z = side * Math.min(16, track.half_width * .55), y = track.surface_z + markingHeight + .02;
    gradientMarkers.push([x, z]);
    arrowPoints.push(x, y, z, x + arrowLength, y, z,
      x + arrowLength, y, z, x + arrowLength * .72, y, z - arrowLength * .22,
      x + arrowLength, y, z, x + arrowLength * .72, y, z + arrowLength * .22);
  }
  const arrowGeometry = new THREE.BufferGeometry();
  arrowGeometry.setAttribute('position', new THREE.Float32BufferAttribute(arrowPoints, 3));
  gradientMaterial = new THREE.LineBasicMaterial({color: 0xe7efcf, transparent: true, opacity: .32});
  gradientField.add(new THREE.LineSegments(arrowGeometry, gradientMaterial));
  gradientField.visible = driveControl.mode === 'gradient';
  orbit.maxDistance = Math.max(150, Math.hypot(length, width) * 6);
  $('track-size-label').textContent = `${length / 10} nm × ${width / 10} nm`;
  $('finish-note').querySelector('span').textContent = `${length / 10} nm`;
}
rebuildTrack();

const moleculeGroup = new THREE.Group();
scene.add(moleculeGroup);
const atomGeometry = new THREE.SphereGeometry(1, 22, 16);
const atomMaterial = new THREE.MeshStandardMaterial({color: 0xffffff, roughness: .31, metalness: .16});
const bondGeometry = new THREE.CylinderGeometry(.115, .115, 1, 12);
const bondMaterial = new THREE.MeshStandardMaterial({color: 0x66715a, roughness: .5, metalness: .1});
const dummy = new THREE.Object3D();
const bondAxis = new THREE.Vector3(0, 1, 0);
const start = new THREE.Vector3(), end = new THREE.Vector3(), direction = new THREE.Vector3();

const cameraAnchor = new THREE.Vector3(), cameraOffset = new THREE.Vector3(12, 22, 28);
function focusCenter(frame = currentFrame || latestFrame) {
  const center = frame?.center || [track.start_x, 0, 3.5];
  // Leave the left HUD clear while the molecule stays in the main viewing area.
  const clearance = sceneHost.clientWidth < 500 ? 4.5 : 5;
  return new THREE.Vector3(center[0] - clearance, Math.max(1.5, center[2] - .7), center[1] + clearance * .4);
}
function homeCamera() {
  cameraMode = 'follow';
  cameraAnchor.copy(focusCenter());
  orbit.target.copy(cameraAnchor);
  camera.position.copy(cameraAnchor).add(cameraOffset);
  needsCameraFocus = false;
  orbit.update();
}
function overviewCamera() {
  cameraMode = 'overview';
  const center = new THREE.Vector3((track.start_x + track.finish_x) / 2, track.surface_z + 8, 0);
  const normal = new THREE.Vector3(.12, .74, 1).normalize();
  const right = new THREE.Vector3().crossVectors(camera.up, normal).normalize();
  const up = new THREE.Vector3().crossVectors(normal, right).normalize();
  const tanV = Math.tan(THREE.MathUtils.degToRad(camera.fov / 2)), tanH = tanV * camera.aspect;
  let distance = 1;
  for (const x of [track.start_x - 12, track.finish_x + 12]) for (const z of [-track.half_width - 10, track.half_width + 10]) for (const y of [track.surface_z, track.surface_z + 16]) {
    const relative = new THREE.Vector3(x, y, z).sub(center);
    distance = Math.max(distance, relative.dot(normal) + 1.14 * Math.max(Math.abs(relative.dot(right)) / tanH, Math.abs(relative.dot(up)) / tanV));
  }
  orbit.target.copy(center);
  camera.position.copy(center).addScaledVector(normal, distance);
  orbit.update();
  updateSceneScale();
}
function updateSceneScale() {
  const center = new THREE.Vector3((track.start_x + track.finish_x) / 2, track.surface_z, 0);
  const radius = Math.hypot((track.finish_x - track.start_x) / 2 + 12, track.half_width + 10, 16);
  const distance = camera.position.distanceTo(center), viewDistance = camera.position.distanceTo(orbit.target);
  const near = cameraMode === 'overview' ? Math.max(.1, distance - radius - 100) : .1;
  const far = Math.max(220, distance + radius * 2 + 100);
  if (camera.near !== near || Math.abs(camera.far - far) > 5) {
    camera.near = near;
    camera.far = far;
    camera.updateProjectionMatrix();
  }
  scene.fog.near = Math.max(65, viewDistance - radius * .6);
  scene.fog.far = Math.max(150, viewDistance + radius * 3);
}
orbit.addEventListener('start', () => { cameraMode = 'free'; });
homeCamera();
needsCameraFocus = true;
function resize() {
  const width = sceneHost.clientWidth, height = sceneHost.clientHeight;
  if (!width || !height) return;
  renderer.setSize(width, height, false);
  camera.aspect = width / height;
  camera.updateProjectionMatrix();
  if (cameraMode === 'overview') overviewCamera();
  else updateSceneScale();
}
new ResizeObserver(resize).observe(sceneHost);
resize();

function feedback(message, isError = false) {
  clearTimeout(feedbackTimer);
  $('feedback').hidden = !message;
  $('feedback').textContent = message || '';
  $('feedback').classList.toggle('error', isError);
  if (message && !isError) feedbackTimer = setTimeout(() => { $('feedback').hidden = true; }, 5500);
}

async function api(path, data) {
  return request(path, data);
}

function clearPlayback() {
  frames = []; currentFrame = latestFrame = null; latestSeq = -1; playback.reset();
  chartSamples = []; lastChartSeq = -1; displayedPositions = null;
  if (atomsMesh) atomsMesh.visible = false;
  if (bondsMesh) bondsMesh.visible = false;
  gradientField.visible = false;
  $('distance-chart-line').setAttribute('d', '');
  for (const id of ['distance-value', 'speed-value', 'time-value']) $(id).textContent = '—';
  $('step-value').textContent = '— 步';
  $('progress-fill').style.width = '0%';
  $('progress-value').textContent = '0%';
  $('download-button').disabled = true;
  needsCameraFocus = true;
}

function setControlsDisabled(disabled) {
  for(const id of ['start-button','reset-button'])$(id).disabled=disabled;
  for(const id of poseControlIds)$(id).disabled=disabled;
}

function readStartPose() {
  const atoms = $('start-pose-mode').value === 'auto' && $('start-pose-heading').value === 'atoms' ? [Number($('start-pose-rear').value), Number($('start-pose-front').value)] : [];
  return validateStartPose({mode: $('start-pose-mode').value, yaw_deg: Number($('start-pose-yaw').value),
    pitch_deg: Number($('start-pose-pitch').value), roll_deg: Number($('start-pose-roll').value), axis_atoms: atoms}, molecule?.atoms.length);
}

function setPosePicking(kind) {
  posePicking = kind;
  posePickBanner.hidden = !kind;
  posePickBanner.textContent = kind === 'rear' ? '点击一个原子作为车尾 · Esc 取消' : '点击一个原子作为车头 · Esc 取消';
  for (const id of ['rear', 'front']) $(`pick-pose-${id}`).setAttribute('aria-pressed', String(kind === id));
}

function updateStartPoseUI(data = state) {
  if (!data) return;
  const available = data.competition?.mode === 'custom' && Boolean(data.start_pose) && !busy && !raceLocked(data) && !['idle', 'preparing', 'running'].includes(data.status);
  if (posePicking && !available) setPosePicking(null);
  const signature = (molecule?.atoms || []).map(atom => `${atom.id}:${atom.element}`).join(',');
  if (signature && signature !== poseAtomSignature) {
    poseAtomSignature = signature;
    for (const id of ['rear', 'front']) {
      const select = $(`start-pose-${id}`);select.replaceChildren();
      for (const atom of molecule.atoms) {const option=document.createElement('option');option.value=atom.id;option.textContent=`${atom.element} · ${atom.id}`;select.append(option);}
    }
    $('start-pose-rear').value = molecule.atoms[0].id;
    $('start-pose-front').value = molecule.atoms.at(-1).id;
  }
  if (data.start_pose && (data.generation !== poseFormGeneration || !poseDraftDirty)) {
    const pose = validateStartPose(data.start_pose);poseFormGeneration = data.generation;poseDraftDirty = false;
    $('start-pose-mode').value = pose.mode;
    $('start-pose-heading').value = pose.axis_atoms.length ? 'atoms' : 'geometry';
    if (signature && pose.axis_atoms.length) {
      $('start-pose-rear').value = pose.axis_atoms[0];$('start-pose-front').value = pose.axis_atoms[1];
    }
    for (const [id,key] of [['yaw','yaw_deg'],['pitch','pitch_deg'],['roll','roll_deg']]) $(`start-pose-${id}`).value=pose[key];
  }
  const mode = $('start-pose-mode').value;
  $('start-pose-auto').hidden = mode !== 'auto';$('start-pose-manual').hidden = mode !== 'manual';
  $('start-pose-atoms').hidden = $('start-pose-heading').value !== 'atoms';
  for (const id of poseControlIds) $(id).disabled = !available;
  for (const id of ['pick-pose-rear','pick-pose-front','start-pose-rear','start-pose-front']) $(id).disabled = !available || !atomCount;
  let valid = true;try {readStartPose();}catch {valid=false;}
  $('apply-start-pose').disabled = !available || !valid;
  $('start-pose-description').textContent = !data.start_pose ? '起跑摆放不可用，请重启游戏服务。' :
    mode === 'manual' ? '按文件坐标轴依次调整侧倾、俯仰和平面旋转。应用后会清空轨迹并重新准备赛车。' :
    mode === 'preserve' ? '采用文件朝向。应用后会清空轨迹并重新准备赛车。' :
    '自动识别贴面姿态与车头参考。应用后会清空轨迹并重新准备赛车。';
  const report = data.runtime?.start_pose;
  const strategy = {'automatic-plane-and-heading':'已自动摆放。',
    'automatic-heading-preserve-roll':'已调整车头朝向。',
    'automatic-plane-preserve-yaw':'已调整贴面姿态。',
    'ambiguous-geometry-input-preserved':'未识别到车头，已采用文件朝向。',
    'manual-file-axes-roll-pitch-yaw':'已应用手动摆放。','input-preserved':'已采用文件朝向。'};
  const after = report?.after_relaxation;
  $('start-pose-result').textContent = report ? (strategy[report.strategy] || '起跑姿态已就绪。') +
    (after?.heading_defined ? ` 朝向 ${after.yaw_deg.toFixed(1)}° · 倾斜 ${after.tilt_deg.toFixed(1)}°。` : '') : '';
}

async function applyStartPose() {
  let pose;try {pose=readStartPose();}catch(error){feedback(error.message,true);return;}
  if (busy) return;
  busy = loading = true;epoch++;clearPlayback();setMolecule(null);setControlsDisabled(true);
  setPosePicking(null);showPreparing();
  try {await api('/api/start-pose',pose);feedback('起跑姿态已提交，正在重新准备赛车。');}
  catch(error){feedback(error.message,true);}
  finally {busy=loading=false;await pollState();}
}

function showPreparing() {
  state = {...state, status: 'preparing', frame: null, finish_reason: null};
  $('status-label').textContent = '准备分子…';
  $('live-dot').classList.remove('running', 'error');
  $('start-label').textContent = '准备分子…';
  $('start-icon').textContent = '▶';
  $('start-button').classList.remove('pausing');
  $('progress-label').textContent = '准备分子…';
  $('scene-message').hidden = false;
  $('scene-message').querySelector('.preparing-spinner').hidden = false;
  $('scene-message').querySelector('b').textContent = '准备分子…';
  $('scene-message').querySelector('small').hidden = true;
  $('scene-message').querySelector('small').textContent = '';
}

function updateRepresentationUI() {
  const spaceFilling = representation === 'space-filling';
  $('ball-stick-model').setAttribute('aria-pressed', String(!spaceFilling));
  $('space-filling-model').setAttribute('aria-pressed', String(spaceFilling));
  if (molecule) $('molecule-info').textContent = moleculeFormula;
}

function setRepresentation(next) {
  if (!['ball-stick', 'space-filling'].includes(next)) return;
  representation = next;
  updateRepresentationUI();
  if (displayedPositions) renderAtoms(displayedPositions);
}

function atomDisplayRadius(element) {
  return representation === 'space-filling' ? vdwRadii[element] : (elementStyle[element] || [0xa79887, .55])[1];
}

function setMolecule(data) {
  const signature = data ? JSON.stringify([data.name, data.atoms, data.bonds]) : '';
  if (signature === moleculeKey) return;
  moleculeKey = signature;
  molecule = data;
  atomCount = data?.atoms?.length || 0;
  for (const mesh of [atomsMesh, bondsMesh]) {
    if (mesh) { moleculeGroup.remove(mesh); mesh.dispose(); }
  }
  atomsMesh = bondsMesh = null;
  clearPlayback();
  if (!atomCount) {
    moleculeFormula = '';
    $('molecule-name').textContent = '准备中';
    $('molecule-info').textContent = '';
    $('atom-chip').textContent = '— 原子';
    $('bond-chip').textContent = '— 键';
    return;
  }
  atomsMesh = new THREE.InstancedMesh(atomGeometry, atomMaterial, atomCount);
  atomsMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
  atomsMesh.castShadow = false;
  atomsMesh.receiveShadow = false;
  atomsMesh.frustumCulled = false;
  for (let i = 0; i < atomCount; i++) {
    const [color] = elementStyle[data.atoms[i].element] || [0xa79887, .55];
    atomsMesh.setColorAt(i, new THREE.Color(color));
  }
  atomsMesh.visible = false;
  moleculeGroup.add(atomsMesh);
  if (data.bonds.length) {
    bondsMesh = new THREE.InstancedMesh(bondGeometry, bondMaterial, data.bonds.length);
    bondsMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    bondsMesh.castShadow = false;
    bondsMesh.frustumCulled = false;
    bondsMesh.visible = false;
    moleculeGroup.add(bondsMesh);
  }
  const knownNames = {'1,2,4,5-Tetraphenylbenzene test car': '四苯基苯', 'Benzene test molecule': '苯'};
  $('molecule-name').textContent = knownNames[data.name] || data.name || '已载入分子';
  $('molecule-name').title = knownNames[data.name] || data.name || '';
  const composition = {};
  for (const atom of data.atoms) composition[atom.element] = (composition[atom.element] || 0) + 1;
  moleculeFormula = Object.entries(composition).sort(([a], [b]) => (a === 'C' ? -1 : b === 'C' ? 1 : a === 'H' ? -1 : b === 'H' ? 1 : a.localeCompare(b)))
    .map(([element, count]) => element + (count > 1 ? String(count).replace(/\d/g, x => '₀₁₂₃₄₅₆₇₈₉'[Number(x)]) : '')).join('');
  updateRepresentationUI();
  $('atom-chip').textContent = `${atomCount} 原子`;
  $('bond-chip').textContent = `${data.bonds.length} 键`;
}

function validFrame(frame) {
  return frame && Array.isArray(frame.positions) && frame.positions.length === atomCount * 3 && frame.positions.every(Number.isFinite);
}

function acceptFrames(incoming, dropped = false) {
  const newFrames = incoming.filter(frame => validFrame(frame) && frame.seq > latestSeq).sort((a,b)=>a.seq-b.seq);
  if (!newFrames.length) return;
  // Bootstrap from retained samples when reconnecting; state polls must never
  // advance the stream cursor past coordinates that have not been fetched.
  if (!frames.length) playback.timeFs = newFrames[0].time_fs;
  if (dropped && frames.length && newFrames[0].seq > latestSeq + 1)
    feedback('部分历史轨迹已不可用，正在接续可用坐标。', true);
  frames.push(...newFrames.map(frame=>({frame})));
  latestFrame = newFrames.at(-1); latestSeq = latestFrame.seq;
  $('scene-message').hidden = true;
  $('download-button').disabled = false;
}

function playbackPending(data = state) {
  return data?.status === 'finished' && (!currentFrame || playback.timeFs < data.frame?.time_fs);
}
function displayStatus(data = state) {
  if (playbackPending(data) || data?.status === 'running') return playback.playing ? 'running' : 'paused';
  return data?.status || 'idle';
}
function raceLocked(data = state) {
  return isRaceInProgress(data, raceStarted) || playbackPending(data);
}
function frameStreamPath() {
  return `/api/frames?after=${latestSeq}&played=${playback.timeFs}&generation=${generation}`;
}

function rangePaint(input) {
  const fraction = (Number(input.value) - Number(input.min)) / (Number(input.max) - Number(input.min));
  input.style.setProperty('--range-fill', `${Math.max(0, Math.min(100, fraction * 100))}%`);
}

function updateState(data) {
  if (generation !== null && Number.isFinite(data.generation) && data.generation < generation) return;
  if (data.generation !== generation) {
    generation = data.generation;
    raceStarted = false;
    clearPlayback();
  }
  const previousStatus = state?.status;
  state = data;
  if (data.status === 'running' && previousStatus !== 'running') {playback.playing = true; playback.lastWall = null;}
  if (data.status === 'paused' && previousStatus === 'running') playback.playing = false;
  raceResults.update(data);
  if (data.status === 'running' || Number(data.frame?.step) > 0) raceStarted = true;
  connected = true;
  setMolecule(data.molecule);
  const previousTrack = trackKey;
  if (data.track) Object.assign(track, data.track);
  rebuildTrack();
  if (trackKey !== previousTrack && cameraMode === 'overview') overviewCamera();
  if (latestSeq < 0 && data.frame?.step === 0 && validFrame(data.frame)) acceptFrames([data.frame]);
  if(data.frame && preparedGeneration!==data.generation){preparedGeneration=data.generation;$('race-setup-panel').open=false;}
  driveControl = {...driveControl,...data.control?.drive};
  gradientField.visible=Boolean(latestFrame);
  if(gradientMaterial)gradientMaterial.opacity=.1+driveControl.strength*.25;
  if(data.competition){
    if(ruleGeneration!==data.generation){raceSetup.set(data.competition);ruleGeneration=data.generation;draftChanged=false;if(data.status!=='idle'||data.selected_input)draftInput={reuse:true};updateRaceDraft();if(data.selected_input){$('race-draft-note').textContent=data.selected_input;$('race-draft-note').hidden=false;}}
    $('race-rules').closest('details').hidden = data.competition.mode === 'standard';
    $('race-rules').textContent = data.competition.mode === 'standard' ? '' : describeCompetition(data.competition);
  }
  const status = displayStatus(data);
  const setupLocked = busy || status === 'preparing' || raceLocked(data);
  const previousDownload = $('download-previous');
  if (previousDownload) {
    previousDownload.hidden = execution === 'browser-wasm' || !data.previousReplay;
    previousDownload.title = data.previousReplay ? data.previousReplay.name : '';
  }
  const timeLimit = status === 'finished' && data.finish_reason === 'time_limit';
  const unavailable = busy || !atomCount || !['paused', 'running'].includes(status);
  setControlsDisabled(unavailable);
  $('load-example-button').disabled = $('molecule-file').disabled = $('example-select').disabled = setupLocked;
  raceSetup.disable(setupLocked);
  $('prepare-race').disabled=setupLocked;
  $('status-label').textContent = status === 'paused' ? ((raceStarted || (data.frame?.step || 0) > 0) ? '已暂停' : '准备就绪') : timeLimit ? '本轮结束' : statuses[status] || status;
  $('live-dot').classList.toggle('running', status === 'running');
  $('live-dot').classList.toggle('error', status === 'error');
  const runnable = atomCount > 0 && !['preparing', 'idle', 'error', 'finished'].includes(status);
  $('start-button').disabled = busy || !runnable;
  $('reset-button').disabled = busy || !atomCount || status === 'preparing';
  updateStartPoseUI(data);
  $('start-label').textContent = status === 'preparing' ? '准备分子…' : status === 'idle' ? '等待分子' : status === 'error' ? '无法试跑' : status === 'running' ? '暂停' : status === 'finished' ? (timeLimit ? '本轮结束' : '已完赛') : status === 'paused' && (raceStarted || (data.frame?.step || 0) > 0) ? '继续试跑' : '开始试跑';
  $('start-icon').textContent = status === 'running' ? 'Ⅱ' : '▶';
  $('start-button').classList.toggle('pausing', status === 'running');
  $('progress-label').textContent = status === 'finished' ? (timeLimit ? '达到模拟时长上限' : '已完赛') : status === 'running' ? '比赛中' : status === 'preparing' ? '准备分子…' : `${(track.finish_x - track.start_x) / 10} nm 跑道`;
  $('engine-badge').classList.remove('disconnected');
  $('engine-badge').querySelector('span').textContent = execution === 'browser-wasm' ? '本机计算' : '已连接';
  $('engine-version').textContent = data.runtime?.lammps_version ? `LAMMPS ${data.runtime.lammps_version}` : 'LAMMPS';
  if (data.error) {
    feedback(data.error, true);
    if (!latestFrame) {
      $('scene-message').hidden = false;
      $('scene-message').querySelector('.preparing-spinner').hidden = true;
      $('scene-message').querySelector('b').textContent = '分子准备失败';
      $('scene-message').querySelector('small').hidden = false;
      $('scene-message').querySelector('small').textContent = '查看错误提示，修改分子后重试。';
    }
  }
  else if (status === 'idle') {
    $('scene-message').hidden=false;
    $('scene-message').querySelector('.preparing-spinner').hidden=true;
    $('scene-message').querySelector('b').textContent='选择赛事，准备赛车';
    $('scene-message').querySelector('small').hidden=true;
  }
  else if (status === 'preparing' && !latestFrame) {
    $('scene-message').hidden = false;
    $('scene-message').querySelector('.preparing-spinner').hidden = false;
    $('scene-message').querySelector('b').textContent = '准备分子…';
    $('scene-message').querySelector('small').hidden = true;
    $('scene-message').querySelector('small').textContent = '';
  }
}

async function pollState() {
  if (statePolling || loading) return;
  statePolling = true;
  const requestEpoch = epoch;
  try {
    const data = await api('/api/state');
    if (requestEpoch === epoch && !loading) {updateState(data);dispatchEvent(new Event('race-state-change'));}
  } catch (error) {
    if (requestEpoch !== epoch || loading) return;
    connected = false;
    $('engine-badge').classList.add('disconnected');
    $('engine-badge').querySelector('span').textContent = execution === 'browser-wasm' ? '计算模块出错' : '未连接';
    $('status-label').textContent = '连接已中断';
    setControlsDisabled(true);
    feedback(execution === 'browser-wasm' ? `计算模块不可用：${error.message} 请刷新网页重试。` : `无法连接游戏服务：${error.message} 请确认服务已启动后重试。`, true);
  } finally { statePolling = false; }
}
async function pollFrames() {
  if (framePolling || loading || !molecule || generation === null) return;
  framePolling = true;
  const requestEpoch = epoch, requestGeneration = generation;
  try {
    const data = await api(frameStreamPath());
    if (requestEpoch !== epoch || loading || requestGeneration !== generation || !molecule) return;
    if (data.generation !== generation) { await pollState(); return; }
    acceptFrames(data.frames || [], Boolean(data.dropped));
  } catch { /* Connection feedback is handled by the less frequent state poll. */ }
  finally { framePolling = false; }
}

async function action(kind) {
  if (busy) return;
  busy = true;
  raceSetup.disable(true);
  $('prepare-race').disabled=true;
  $('start-button').disabled = $('reset-button').disabled = true;
  try {
    if (kind === 'reset') {
      epoch++;
      loading = true;
      clearPlayback();

      setControlsDisabled(true);
      showPreparing();
    }
    if (kind === 'start' || kind === 'pause') {
      playback.playing = kind === 'start'; playback.lastWall = null;
      if (kind === 'start') {
        raceStarted = true;
        // Register the compute window before allowing the first MD batch.
        const buffered = await api(frameStreamPath());
        if (buffered.generation === generation) acceptFrames(buffered.frames || [], buffered.dropped);
      }
      const fresh = await api('/api/state');
      if (fresh.status !== 'finished') {
        try {await api('/api/action', {action: kind});}
        catch (error) {if ((await api('/api/state')).status !== 'finished') throw error;}
      }
    } else await api('/api/action', {action: kind});
  } catch (error) { feedback(error.message, true); }
  finally { busy = loading = false; await pollState(); }
}

const raycaster = new THREE.Raycaster(), pointer = new THREE.Vector2();

async function loadStructure(data) {
  if (busy) return;
  if (raceLocked()) {feedback('请先结束或重置当前比赛。', true);return;}
  busy = loading = true;
  epoch++;
  clearPlayback();
  setMolecule(null);
  setControlsDisabled(true);

  showPreparing();
  feedback('正在导入分子…');
  try {
    await api('/api/load', {...data, competition:raceSetup.read()});
    feedback('分子已导入。');
  } catch (error) { feedback(error.message, true); }
  finally {
    busy = loading = false;
    $('load-example-button').disabled = false;
    await pollState();
  }
}

$('start-button').addEventListener('click', () => action(playback.playing ? 'pause' : 'start'));
$('reset-button').addEventListener('click', () => action('reset'));
for (const id of poseControlIds.filter(id=>!id.startsWith('pick-')&&id!=='apply-start-pose')) $(id).addEventListener('change',()=>{
  poseDraftDirty=true;setPosePicking(null);updateStartPoseUI();
});
$('apply-start-pose').addEventListener('click',applyStartPose);
for (const kind of ['rear','front']) $(`pick-pose-${kind}`).addEventListener('click',()=>setPosePicking(posePicking===kind?null:kind));
renderer.domElement.addEventListener('click',event=>{
  if(!posePicking||busy||!atomsMesh)return;
  if(!posePointerDown||Math.hypot(event.clientX-posePointerDown.x,event.clientY-posePointerDown.y)>5)return;
  posePointerDown=null;
  const rect=renderer.domElement.getBoundingClientRect();
  pointer.set((event.clientX-rect.left)/rect.width*2-1,-(event.clientY-rect.top)/rect.height*2+1);
  // Instances move during playback; refresh the picking bound only on demand.
  atomsMesh.computeBoundingSphere();atomsMesh.updateWorldMatrix(true,false);
  raycaster.setFromCamera(pointer,camera);const hit=raycaster.intersectObject(atomsMesh,false)[0];
  if(!hit||!Number.isInteger(hit.instanceId))return;
  $(`start-pose-${posePicking}`).value=molecule.atoms[hit.instanceId].id;
  poseDraftDirty=true;setPosePicking(null);updateStartPoseUI();
});
renderer.domElement.addEventListener('pointerdown',event=>{posePointerDown=posePicking?{x:event.clientX,y:event.clientY}:null;});
renderer.domElement.addEventListener('pointercancel',()=>{posePointerDown=null;});
document.addEventListener('keydown',event=>{if(event.key==='Escape'&&posePicking)setPosePicking(null);});
$('home-camera').addEventListener('click', homeCamera);
$('overview-camera').addEventListener('click', overviewCamera);
$('ball-stick-model').addEventListener('click', () => setRepresentation('ball-stick'));
$('space-filling-model').addEventListener('click', () => setRepresentation('space-filling'));
for (const rate of [1,2]) $(`playback-${rate}x`).addEventListener('click',()=>{
  playback.setRate(rate);
  for (const value of [1,2]) $(`playback-${value}x`).setAttribute('aria-pressed',String(value===rate));
  $('playback-note').textContent = `100 ps / ${50/rate} 秒`;
});
document.addEventListener('visibilitychange',()=>{
  playback.lastWall=null;
  if(document.hidden && playback.playing) action('pause');
});
$('load-example-button').addEventListener('click', () => {
  draftInput={example:$('example-select').value};draftChanged=true;
  feedback('已选择 '+$('example-select').selectedOptions[0].textContent+'，点击准备赛车。');updateRaceDraft();
});
$('prepare-race').addEventListener('click', async()=>{
  try {raceSetup.read();await loadStructure(draftInput);}catch(error){feedback(error.message,true);}
});
$('molecule-file').addEventListener('change', async event => {
  const file = event.target.files?.[0];
  if (!file) return;
  if (file.size > 2 * 1024 * 1024) { feedback('请导入小于 2 MB 的分子结构文件。', true); event.target.value = ''; return; }
  try { draftInput={filename:file.name,content:await file.text()};draftChanged=true;updateRaceDraft();feedback('已选择 '+file.name+'，点击准备赛车。'); }
  catch (error) { feedback(error.message, true); }
  event.target.value = '';
});
$('download-button').addEventListener('click', async () => {
  try {
    await saveReplay();
    const replayMenu = $('replay-menu');
    if (replayMenu) replayMenu.open = false;
    feedback('已发起轨迹下载。');
  } catch (error) { feedback(`轨迹导出失败：${error.message} 请重试。`, true); }
});
$('download-previous')?.addEventListener('click', () => {
  const link = document.createElement('a');
  link.href = '/api/replay/previous';
  link.download = 'molecular-car-previous-replay.json';
  link.click();
  $('replay-menu').open = false;
  feedback('已发起上一场轨迹下载。');
});

function renderAtoms(positions) {
  if (!atomsMesh || positions.length !== atomCount * 3) return;
  displayedPositions = positions;
  for (let i = 0; i < atomCount; i++) {
    const radius = atomDisplayRadius(molecule.atoms[i].element);
    dummy.position.set(positions[i * 3], positions[i * 3 + 2], positions[i * 3 + 1]);
    dummy.quaternion.identity();
    dummy.scale.setScalar(radius);
    dummy.updateMatrix();
    atomsMesh.setMatrixAt(i, dummy.matrix);
  }
  atomsMesh.instanceMatrix.needsUpdate = true;
  atomsMesh.visible = true;
  if (bondsMesh) {
    bondsMesh.visible = representation === 'ball-stick';
    if (!bondsMesh.visible) return;
    const atomIndex = new Map(molecule.atoms.map((atom, index) => [atom.id, index]));
    molecule.bonds.forEach((bond, index) => {
      const a = atomIndex.get(bond.a), b = atomIndex.get(bond.b);
      if (a === undefined || b === undefined) { dummy.scale.setScalar(0); dummy.updateMatrix(); bondsMesh.setMatrixAt(index, dummy.matrix); return; }
      start.set(positions[a * 3], positions[a * 3 + 2], positions[a * 3 + 1]);
      end.set(positions[b * 3], positions[b * 3 + 2], positions[b * 3 + 1]);
      direction.subVectors(end, start);
      dummy.position.copy(start).add(end).multiplyScalar(.5);
      dummy.quaternion.setFromUnitVectors(bondAxis, direction.clone().normalize());
      dummy.scale.set(1, direction.length(), 1);
      dummy.updateMatrix();
      bondsMesh.setMatrixAt(index, dummy.matrix);
    });
    bondsMesh.instanceMatrix.needsUpdate = true;
    bondsMesh.visible = true;
  }
}

function updateMetrics(frame) {
  if (!frame) return;
  const distance = Number(frame.distance_angstrom ?? 0) / 10;
  const speed = Number(frame.speed_angstrom_ps ?? 0) / 10;
  const time = Number(frame.time_fs ?? 0);
  $('distance-value').textContent = Number.isFinite(distance) ? distance.toFixed(2) : '—';
  $('speed-value').textContent = Number.isFinite(speed) ? speed.toFixed(3) : '—';
  $('time-value').textContent = Number.isFinite(time) ? (time < 1000 ? time.toFixed(0) : (time / 1000).toFixed(2)) : '—';
  $('time-unit').textContent = time < 1000 ? ' fs' : ' ps';
  $('step-value').textContent = `${Number(frame.step || 0).toLocaleString('zh-CN')} 步`;
  const progress = Math.max(0, Math.min(1, Number(frame.progress || 0)));
  $('progress-fill').style.width = `${progress * 100}%`;
  $('progress-value').textContent = `${(progress * 100).toFixed(0)}%`;
  if (frame.seq !== lastChartSeq) {
    chartSamples.push(distance);
    if (chartSamples.length > 60) chartSamples.shift();
    lastChartSeq = frame.seq;
    const low = Math.min(0, ...chartSamples), high = Math.max(.2, ...chartSamples);
    const d = chartSamples.map((value, i) => `${i ? 'L' : 'M'}${(i / Math.max(1, chartSamples.length - 1) * 180).toFixed(2)} ${(25 - (value - low) / (high - low) * 23).toFixed(2)}`).join(' ');
    $('distance-chart-line').setAttribute('d', d);
  }
}

let metricsTick = 0;
function animate(now) {
  requestAnimationFrame(animate);
  const visible = !document.hidden && document.documentElement.dataset.appActive !== 'false';
  const targetTime = playback.advance(now, latestFrame?.time_fs ?? playback.timeFs, visible && connected && !loading);
  if (!visible) return;
  if (frames.length) {
    const previousDisplayStatus = displayStatus();
    while (frames.length > 1 && frames[1].frame.time_fs <= targetTime) frames.shift();
    const a = frames[0], b = frames[1] || a;
    const alpha = a === b ? 0 : Math.max(0, Math.min(1, (targetTime - a.frame.time_fs) / (b.frame.time_fs - a.frame.time_fs)));
    const interpolated = a.frame.positions.map((value, index) => value + (b.frame.positions[index] - value) * alpha);
    renderAtoms(interpolated);
    currentFrame = {...a.frame};
    for (const key of ['time_fs','distance_angstrom','speed_angstrom_ps','progress'])
      currentFrame[key] = a.frame[key] + (b.frame[key] - a.frame[key]) * alpha;
    currentFrame.center = a.frame.center.map((value,index)=>value+(b.frame.center[index]-value)*alpha);
    if (needsCameraFocus) homeCamera();
    if (cameraMode === 'follow') {
      const nextAnchor = focusCenter(currentFrame);
      const shift = nextAnchor.sub(orbit.target).multiplyScalar(.1);
      camera.position.add(shift);
      orbit.target.add(shift);
    }
    const center = currentFrame.center;
    if (Array.isArray(center) && center.length === 3 && center.every(Number.isFinite)) {
      sun.position.set(center[0] - 15, center[2] + 45, center[1] + 24);
      sun.target.position.set(center[0], track.surface_z, center[1]);
    }
    gradientField.visible = driveControl.mode === 'gradient';
    if (now - metricsTick > 100) { updateMetrics(currentFrame); metricsTick = now; }
    raceResults.present(a.frame.seq);
    if (state?.status === 'finished' && !playbackPending()) {
      playback.playing = false; updateMetrics(currentFrame);
    }
    if (previousDisplayStatus !== displayStatus() || (state?.status === 'finished' && !playbackPending() && !$('start-button').disabled)) updateState(state);
    if (displayStatus() === 'running' && playback.buffering) $('status-label').textContent='正在缓冲';
    else if (displayStatus() === 'running') $('status-label').textContent='比赛中';
  }
  orbit.update();
  updateSceneScale();
  goldSurface?.update(camera, orbit.target);
  renderer.render(scene, camera);
}
requestAnimationFrame(animate);

window.__molecularCarDebug = {
  getSnapshot: () => ({connected, status: state?.status, atomCount, bondCount: molecule?.bonds.length || 0,
    generation, latestSeq, renderedSeq: currentFrame?.seq, renderedStep: currentFrame?.step, bufferedFrames: frames.length,
    playbackTimeFs: playback.timeFs, playbackRate: playback.rate, playbackPlaying: playback.playing,
    playbackStatus: displayStatus(), buffering: playback.buffering, computedTimeFs: state?.frame?.time_fs,
    positions: displayedPositions ? [...displayedPositions] : null, competition:state?.competition,
    drive: {...driveControl}, gradientVisible: gradientField.visible,
    startPose: state?.start_pose || null, startPoseMetadata: state?.runtime?.start_pose || null, posePicking,
    surface: state?.settings?.surface || 'fourier-au', surfaceMetadata: state?.runtime?.surface || null,
    surfaceEvaluation: currentFrame?.surface_evaluation || null,
    gradientMarkers: gradientMarkers.map(marker => [...marker]),
    track: {...track}, runway: runway ? {position: runway.position.toArray(), size: [runway.geometry.parameters.width, runway.geometry.parameters.height]} : null,
    goldSurface: goldSurface?.getSnapshot() || null,
    finishX: finish.position.x, clipping: {near: camera.near, far: camera.far, maxDistance: orbit.maxDistance, fogNear: scene.fog.near, fogFar: scene.fog.far},
    shadow: {enabled: renderer.shadowMap.enabled, carCastsShadow: Boolean(atomsMesh?.castShadow || bondsMesh?.castShadow),
      light: sun.position.toArray(), target: sun.target.position.toArray(), width: sun.shadow.camera.right - sun.shadow.camera.left},
    representation, atomRadii: molecule?.atoms.map(atom => atomDisplayRadius(atom.element)) || [],
    atomVisuals: displayedPositions && atomsMesh ? molecule.atoms.map((atom, index) => {
      const matrix = new THREE.Matrix4();
      const color = new THREE.Color();
      atomsMesh.getMatrixAt(index, matrix);
      atomsMesh.getColorAt(index, color);
      return {id: atom.id, element: atom.element, radius: matrix.getMaxScaleOnAxis(),
        position: [matrix.elements[12], matrix.elements[14], matrix.elements[13]], color: color.getHex()};
    }) : [],
    radiiSource: 'RDKit 2025.09.3 GetRvdw (Å)', bondsVisible: Boolean(bondsMesh?.visible),
    camera: camera.position.toArray(), cameraTarget: orbit.target.toArray(), cameraMode, sceneSize: [sceneHost.clientWidth, sceneHost.clientHeight],
    renderer: {calls: renderer.info.render.calls, triangles: renderer.info.render.triangles},
  }),
};

await pollState();
setInterval(pollState, 650);
setInterval(pollFrames, 150);


if(execution==='browser-wasm'){
  const applyPreferences=data=>{
    setRepresentation(data.representation);
    renderer.setPixelRatio(Math.min(devicePixelRatio,data.quality==='low'?1:data.quality==='high'?2:1.5));
    orbit.enableDamping=data.motion;
    document.documentElement.dataset.reducedMotion=String(!data.motion);
  };
  applyPreferences(readPreferences());
  installRaceBridge({request,pause:()=>action('pause'),snapshot:()=>state?{...state,status:displayStatus(),started:raceStarted,unsavedReplay:unsavedReplay()}:null,preferences:applyPreferences,load:async data=>{
    if(busy||state?.status==='preparing')throw new Error('赛道正在准备，请稍后再试。');
    if(raceLocked())throw new Error('请先重置本场比赛。');
    const previous=state?.generation;
    if(data.competition)raceSetup.set(data.competition);draftInput={filename:data.filename,content:data.content};draftChanged=true;
    await loadStructure(draftInput);
    const loaded=await request('/api/state');
    if(loaded.generation===previous)throw new Error(loaded.error||'赛车未能载入，请重新准备。');
  }});
}
