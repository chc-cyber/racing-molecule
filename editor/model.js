// Preserve useful errors when the JSmol Java runtime exports its Error class.
const NativeError=globalThis.Clazz?._Error??globalThis.Error;
/**
 * Molecular Car workshop's canonical graph and edit history.
 * No coordinate operation infers connectivity. pos2 and pos3 are independent.
 * Units: angstroms and degrees. The state accepts incomplete chemical designs;
 * chemistryWarnings is advisory, except the explicit fragment connect command.
 */
export const PROJECT_FORMAT = 'molecular-car-workshop-v1';
const ELEMENTS = new Set('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og'.split(' '));
const CHIRAL = new Set(['CHI_UNSPECIFIED', 'CHI_TETRAHEDRAL_CW', 'CHI_TETRAHEDRAL_CCW']);
const STEREO = new Set(['STEREONONE', 'STEREOE', 'STEREOZ', 'STEREOCIS', 'STEREOTRANS', 'STEREOANY']);
const DIRECTION = new Set(['NONE', 'BEGINWEDGE', 'BEGINDASH', 'EITHERDOUBLE', 'UNKNOWN', 'ENDUPRIGHT', 'ENDDOWNRIGHT']);
const ORDERS = new Set([1, 2, 3, 1.5]);
let sequence = 0;

export function clone(value) { return structuredClone(value); }
function id(prefix) {
  return `${prefix}-${globalThis.crypto?.randomUUID?.() ?? `${Date.now().toString(36)}-${(++sequence).toString(36)}`}`;
}
function ensure(condition, message) { if (!condition) throw new NativeError(message); }
function finite(value, name) { ensure(typeof value === 'number' && Number.isFinite(value), `${name} 必须是有限数值。`); return value; }
function point(value, dimension, name) {
  ensure(Array.isArray(value) && value.length === dimension && value.every(v => typeof v === 'number' && Number.isFinite(v)), `${name} 必须包含 ${dimension} 个有限坐标。`);
  return [...value];
}
function identifier(value, name) { ensure(typeof value === 'string' && value.length > 0 && value.length <= 200, `${name} 不能为空或过长。`); }
function enumeration(value, values, name) { ensure(values.has(value), `${name} 不受支持：${value}`); }
const add = (a, b) => a.map((v, i) => v + b[i]);
const sub = (a, b) => a.map((v, i) => v - b[i]);
const scale = (a, s) => a.map(v => v * s);
const dot = (a, b) => a.reduce((sum, v, i) => sum + v * b[i], 0);
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = a => Math.hypot(...a);
function unit(a) { const length = norm(a); ensure(length > 1e-10, '几何方向未定义：所选原子的坐标重合或轴长度为零。'); return scale(a, 1 / length); }
const radians = degrees => degrees * Math.PI / 180;
const degrees = angle => angle * 180 / Math.PI;
function rotateVector(vector, axis, value) {
  const c = Math.cos(value), s = Math.sin(value);
  return add(add(scale(vector, c), scale(cross(axis, vector), s)), scale(axis, dot(axis, vector) * (1 - c)));
}

export function createProject(name = '未命名赛车') {
  return {format: PROJECT_FORMAT, name, atoms: [], bonds: [], fragments: [], has3d: false};
}

/** Structural validation only: do not mistake simulator admission for editing validity. */
export function validateProject(input) {
  ensure(input && typeof input === 'object' && !Array.isArray(input), '项目必须是一个对象。');
  const p = clone(input);
  ensure(p.format === PROJECT_FORMAT, '这不是受支持的分子赛车工坊项目格式。');
  ensure(typeof p.name === 'string', '项目名称必须是文本。');
  ensure(Array.isArray(p.atoms) && Array.isArray(p.bonds) && Array.isArray(p.fragments), '项目必须包含原子、键和片段列表。');
  ensure(typeof p.has3d === 'boolean', '项目必须标记是否已有三维构象。');
  const atomIds = new Set(), bondIds = new Set(), fragmentIds = new Set(), pairs = new Set();
  for (const atom of p.atoms) {
    identifier(atom.id, '原子 ID');
    ensure(!atomIds.has(atom.id), `原子 ID 重复：${atom.id}`); atomIds.add(atom.id);
    ensure(ELEMENTS.has(atom.element), `不支持的元素符号：${atom.element}`);
    atom.charge ??= 0;
    ensure(Number.isSafeInteger(atom.charge), `原子 ${atom.id} 的形式电荷必须是整数。`);
    point(atom.pos2, 2, `${atom.id} 的二维坐标`); point(atom.pos3, 3, `${atom.id} 的三维坐标`);
    if (atom.isotope !== undefined) ensure(Number.isSafeInteger(atom.isotope) && atom.isotope >= 0, '同位素质量数必须是非负整数。');
    if (atom.chiralTag !== undefined) enumeration(atom.chiralTag, CHIRAL, '原子手性标记');
    for (const key of ['explicitHs', 'radicalElectrons']) if (atom[key] !== undefined) ensure(Number.isSafeInteger(atom[key]) && atom[key] >= 0, `${key} 必须是非负整数。`);
    if (atom.noImplicit !== undefined) ensure(typeof atom.noImplicit === 'boolean', 'noImplicit 必须是布尔值。');
    if (atom.generatedHydrogen !== undefined) ensure(typeof atom.generatedHydrogen === 'boolean' && atom.element === 'H', '自动补氢标记只能用于氢原子。');
  }
  for (const bond of p.bonds) {
    identifier(bond.id, '键 ID');
    ensure(!bondIds.has(bond.id), `键 ID 重复：${bond.id}`); bondIds.add(bond.id);
    ensure(atomIds.has(bond.a) && atomIds.has(bond.b), `键 ${bond.id} 引用了不存在的原子。`);
    ensure(bond.a !== bond.b, '一个原子不能与自己成键。');
    const pair = JSON.stringify([bond.a, bond.b].sort());
    ensure(!pairs.has(pair), '同一对原子只能保存一条键；请修改已有键的键级。'); pairs.add(pair);
    enumeration(bond.order, ORDERS, '键级');
    if (bond.stereo !== undefined) enumeration(bond.stereo, STEREO, '键立体标记');
    if (bond.direction !== undefined) enumeration(bond.direction, DIRECTION, '键方向');
    if (bond.stereoAtoms !== undefined) {
      ensure(Array.isArray(bond.stereoAtoms) && [0, 2].includes(bond.stereoAtoms.length) && bond.stereoAtoms.every(i => atomIds.has(i)), '键立体参考原子必须为空列表或两个有效的稳定原子 ID。');
      if (bond.stereoAtoms.length) {
        ensure(bond.stereoAtoms[0] !== bond.stereoAtoms[1], '键立体参考原子不能相同。');
        ensure(validStereoReferences(p, bond), '双键立体参考必须分别与双键两端直接成键；请按 a–b–c–d 顺序选择四个连续原子。');
      }
    }
  }
  for (const fragment of p.fragments) {
    identifier(fragment.id, '片段 ID');
    ensure(!fragmentIds.has(fragment.id), `片段 ID 重复：${fragment.id}`); fragmentIds.add(fragment.id);
    ensure(typeof fragment.name === 'string', '片段必须有名称。');
    ensure(Array.isArray(fragment.atomIds) && fragment.atomIds.length > 0 && fragment.atomIds.every(i => atomIds.has(i)), '片段包含不存在的原子或没有原子。');
    ensure(new Set(fragment.atomIds).size === fragment.atomIds.length, '片段中的原子 ID 不能重复。');
    ensure(typeof fragment.rigid === 'boolean', '片段必须标记是否默认整体移动。');
    ensure(Array.isArray(fragment.ports), '片段必须有连接位点列表。');
    for (const port of fragment.ports) {
      ensure(fragment.atomIds.includes(port.atomId), '连接位点必须属于其片段。');
      if (port.maxOrder !== undefined && port.maxOrder !== 0) enumeration(port.maxOrder, ORDERS, '连接位点允许键级');
    }
  }
  const groups = p.metadata?.ketcherGroups;
  if (groups !== undefined) {
    ensure(Array.isArray(groups), '原生折叠基团必须是列表。');
    const groupIds = new Set();
    for (const group of groups) {
      identifier(group.id, '基团 ID');
      ensure(!groupIds.has(group.id), '基团 ID 重复。'); groupIds.add(group.id);
      ensure(group.type === 'SUP' && typeof group.name === 'string' && group.name.length <= 200, '当前仅支持有名称的原生 Superatom 基团。');
      ensure(typeof group.expanded === 'boolean', '基团必须保存展开状态。');
      ensure(Array.isArray(group.atomIds) && group.atomIds.length > 0 && new Set(group.atomIds).size === group.atomIds.length && group.atomIds.every(i => atomIds.has(i)), '折叠基团引用的原子无效。');
      ensure(Array.isArray(group.attachmentPoints ?? []), '基团连接点必须是列表。');
      for (const port of group.attachmentPoints ?? []) ensure(group.atomIds.includes(port.atomId), '基团连接点必须属于其内部真实原子。');
      if (group.template?.indexToAtomId) ensure(group.template.indexToAtomId.length === group.atomIds.length && group.template.indexToAtomId.every(i => group.atomIds.includes(i)), '模板原子来源映射无效。');
    }
  }
  const pending = p.metadata?.geometryPending;
  if (pending !== undefined) for (const key of ['atomIds', 'changedAtomIds', ...(pending.editedAtomIds ? ['editedAtomIds'] : [])]) ensure(Array.isArray(pending[key]) && new Set(pending[key]).size === pending[key].length && pending[key].every(i => atomIds.has(i)), '待准备三维的原子引用无效。');
  return p;
}

/** A substitution can consume only H generated by this editor. Imported or
 * deliberately drawn hydrogens remain explicit and produce a valence warning.
 * The removal is part of the user's bond transaction and its shared undo. */
export function reconcileGeneratedHydrogens(previous, next) {
  const oldAtoms = new Map(previous.atoms.map(a => [a.id, a]));
  const oldBonds = new Map(previous.bonds.map(b => [b.id, b]));
  const changed = new Set();
  for (const b of next.bonds) {
    const old = oldBonds.get(b.id);
    if (old && old.order >= b.order && ((old.a === b.a && old.b === b.b) || (old.a === b.b && old.b === b.a))) continue;
    const ends = [b.a, b.b].map(i => next.atoms.find(a => a.id === i));
    if (ends.every(a => a && a.element !== 'H')) for (const a of ends) changed.add(a.id);
  }
  const removed = new Set();
  for (const atomId of changed) {
    const atom = next.atoms.find(a => a.id === atomId);
    const maximum = maximumValence(atom);
    if (maximum === null || !oldAtoms.has(atomId)) continue;
    let total = (atom.explicitHs ?? 0) + next.bonds.filter(b => !removed.has(b.a) && !removed.has(b.b) && (b.a === atomId || b.b === atomId)).reduce((sum, b) => sum + b.order, 0);
    for (const b of next.bonds) {
      if (total <= maximum + 1e-8) break;
      const other = b.a === atomId ? b.b : b.b === atomId ? b.a : null;
      const h = other && next.atoms.find(a => a.id === other);
      if (b.order !== 1 || h?.element !== 'H' || !h.generatedHydrogen || !oldAtoms.has(h.id) || next.bonds.filter(x => x.a === h.id || x.b === h.id).length !== 1) continue;
      removed.add(h.id); total--;
    }
  }
  if (removed.size) {
    next.atoms = next.atoms.filter(a => !removed.has(a.id));
    next.bonds = next.bonds.filter(b => !removed.has(b.a) && !removed.has(b.b));
    next.metadata ??= {};
    next.metadata.hydrogenNotice = `新增连接替代了 ${removed.size} 个自动补入的氢。`;
  }
  for (const atom of next.atoms) if (atom.element !== 'H') delete atom.generatedHydrogen;
  return [...removed];
}

function pruneEditorMetadata(project) {
  const available = new Set(project.atoms.map(a => a.id));
  project.fragments = project.fragments.map(f => ({...f, atomIds: f.atomIds.filter(i => available.has(i)), ports: f.ports.filter(p => available.has(p.atomId))})).filter(f => f.atomIds.length);
  if (project.metadata?.ketcherGroups) project.metadata.ketcherGroups = project.metadata.ketcherGroups.map(group => {
    const atomIds = group.atomIds.filter(i => available.has(i));
    const next = {...group, atomIds, attachmentPoints: (group.attachmentPoints ?? []).filter(p => atomIds.includes(p.atomId))};
    if (atomIds.length !== group.atomIds.length) { delete next.template; delete next.kind; }
    return next;
  }).filter(g => g.atomIds.length);
  if (project.metadata?.geometryPending) for (const key of ['atomIds', 'changedAtomIds', ...(project.metadata.geometryPending.editedAtomIds ? ['editedAtomIds'] : [])]) project.metadata.geometryPending[key] = project.metadata.geometryPending[key].filter(i => available.has(i));
}

/** Track chemical edits, independently of layout, selection and 3D movement.
 * Incoming imports / version-checked geometry results already carry state. */
export function reconcileC60Addition(previous,next){
  const oldBonds=new Map(previous.bonds.map(b=>[b.id,b]));let changed=false;
  for(const group of next.metadata?.ketcherGroups??[]){
    if(group.kind!=='c60')continue;const members=new Set(group.atomIds);
    for(const bond of next.bonds){
      if(members.has(bond.a)===members.has(bond.b))continue;
      const old=oldBonds.get(bond.id);if(old&&old.a===bond.a&&old.b===bond.b&&old.order===bond.order)continue;
      const inside=members.has(bond.a)?bond.a:bond.b,atom=atomById(next,inside);
      ensure(atom.element==='C'&&atom.charge===0&&bond.order===1,'C60 加成连接目前支持中性笼碳的外接单键；此连接未应用。');
      if(valence(next,inside)>4){
        ensure(Math.abs(valence(next,inside)-5)<1e-8,'C60 连接位点存在多个过价连接。');
        const double=next.bonds.find(b=>b.order===2&&((b.a===inside&&members.has(b.b))||(b.b===inside&&members.has(b.a))));
        ensure(double,'C60 连接位点已满价，没有可用于加成的笼内双键。');
        Object.assign(double,{order:1,direction:'NONE',stereo:'STEREONONE',stereoAtoms:[]});changed=true;
      }
    }
  }
  return changed;
}
export function validateBondEndpoints(project,a,b){
  for(const id of [a,b]){const atom=atomById(project,id),maximum=maximumValence(atom);ensure(maximum===null||valence(project,id)<=maximum+1e-8,`${atom.element} 的连接位点已满价；请先处理该位点的显式氢或键级。`);}
}
export function validateNativeGroupConnections(previous, next) {
  const oldGroups = new Set((previous.metadata?.ketcherGroups ?? []).map(g => g.id));
  const oldBonds = new Map(previous.bonds.map(b => [b.id, b]));
  const atoms = new Map(next.atoms.map(a => [a.id, a]));
  for (const group of next.metadata?.ketcherGroups ?? []) {
    if (group.kind !== 'c60') continue;
    const members = new Set(group.atomIds);
    for (const bond of next.bonds) {
      if (members.has(bond.a) === members.has(bond.b)) continue;
      const old = oldBonds.get(bond.id);
      if (oldGroups.has(group.id) && old && old.order >= bond.order && ((old.a === bond.a && old.b === bond.b) || (old.a === bond.b && old.b === bond.a))) continue;
      for (const id of [bond.a, bond.b]) {
        const atom = atoms.get(id), maximum = maximumValence(atom);
        ensure(maximum !== null && valence(next, id) <= maximum + 1e-8,
          `C60 连接未应用：${atom.element} 的连接后价态为 ${valence(next, id)}。请先降低接点相邻双键级或移除待替代的显式氢，再连接。`);
      }
    }
  }
}

function trackGeometryPending(previous, next) {
  if (!previous.has3d && !next.has3d) return;
  const oldAtoms = new Map(previous.atoms.map(a => [a.id, a])), atoms = new Map(next.atoms.map(a => [a.id, a]));
  const pending = next.metadata?.geometryPending ?? {atomIds: [], changedAtomIds: []};
  const added = new Set(pending.atomIds), changed = new Set(pending.changedAtomIds), edited = new Set(pending.editedAtomIds ?? []);
  const chemistry = a => JSON.stringify(['element', 'charge', 'isotope', 'explicitHs', 'noImplicit', 'radicalElectrons', 'chiralTag'].map(k => a[k] ?? null));
  const geometryChemistry = a => JSON.stringify(['element', 'charge', 'isotope', 'explicitHs', 'noImplicit', 'radicalElectrons'].map(k => a[k] ?? null));
  for (const atom of next.atoms) { const old = oldAtoms.get(atom.id); if (!old) added.add(atom.id); else if (chemistry(old) !== chemistry(atom)) { changed.add(atom.id); if(geometryChemistry(old)!==geometryChemistry(atom))edited.add(atom.id); } }
  const oldBonds = new Map(previous.bonds.map(b => [b.id, b])), bonds = new Map(next.bonds.map(b => [b.id, b]));
  const chemicalBond = b => JSON.stringify([[b.a, b.b].sort(), b.order, b.stereo ?? 'STEREONONE', b.direction ?? 'NONE', b.stereoAtoms ?? []]);
  const geometryBond = b => JSON.stringify([[b.a, b.b].sort(), b.order]);
  for (const [id, b] of oldBonds) {
    const current = bonds.get(id);
    if (!current || chemicalBond(b) !== chemicalBond(current)) for (const a of [b.a, b.b]) if (atoms.has(a)) {
      changed.add(a);
      if ((!current||geometryBond(b)!==geometryBond(current))&&[b.a, b.b].every(i => !added.has(i) && oldAtoms.get(i)?.element !== 'H')) edited.add(a);
    }
  }
  for (const [id, b] of bonds) {
    const old = oldBonds.get(id);
    if (!old || chemicalBond(old) !== chemicalBond(b)) for (const a of [b.a, b.b]) if (oldAtoms.has(a)) {
      changed.add(a);
      if (old&&geometryBond(old)!==geometryBond(b)&&[b.a, b.b].every(i => !added.has(i) && oldAtoms.has(i) && oldAtoms.get(i)?.element !== 'H')) edited.add(a);
    }
  }
  for (const i of [...added]) if (!atoms.has(i)) added.delete(i);
  for (const i of [...changed]) if (!atoms.has(i) || added.has(i)) changed.delete(i);
  for (const i of [...edited]) if (!atoms.has(i) || added.has(i)) edited.delete(i);
  if (added.size || changed.size) { next.metadata ??= {}; next.metadata.geometryPending = {atomIds: [...added], changedAtomIds: [...changed], ...(edited.size ? {editedAtomIds: [...edited]} : {})}; }
}

export function serializeProject(project) { return JSON.stringify(validateProject(project), null, 2); }
export function parseProject(text) {
  let parsed;
  try { parsed = JSON.parse(text); } catch { throw new NativeError('无法读取项目 JSON，请检查文件是否完整。'); }
  return validateProject(parsed);
}

export function atomById(project, atomId) {
  const atom = project.atoms.find(atom => atom.id === atomId);
  ensure(atom, `找不到原子：${atomId}`); return atom;
}
function bondBetween(project, a, b) { return project.bonds.find(bond => (bond.a === a && bond.b === b) || (bond.a === b && bond.b === a)); }
export function valence(project, atomId) {
  const atom = atomById(project, atomId);
  return (atom.explicitHs ?? 0) + project.bonds.filter(b => b.a === atomId || b.b === atomId).reduce((sum, b) => sum + b.order, 0);
}
function maximumValence(atom) {
  const common = {H: 1, B: atom.charge < 0 ? 4 : 3, C: atom.charge === 0 ? 4 : 3, N: atom.charge > 0 ? 4 : (atom.charge < 0 ? 2 : 3), O: atom.charge > 0 ? 3 : (atom.charge < 0 ? 1 : 2), F: 1, Si: 4, P: 5, S: 6, Cl: 1, Br: 1, I: 1};
  return common[atom.element] ?? null;
}

/** Advisory checks only; this is deliberately not a full chemical perception engine. */
export function chemistryWarnings(project) {
  const warnings = [];
  if (project.metadata?.stereoNotice) warnings.push(project.metadata.stereoNotice);
  for (const atom of project.atoms) {
    const maximum = maximumValence(atom), actual = valence(project, atom.id);
    if (maximum !== null && actual > maximum + 1e-8) warnings.push(`${atom.id} 的键级与氢数之和为 ${actual}，超过常见价态 ${maximum}。请检查化学键和氢原子。`);
    if (atom.radicalElectrons) warnings.push(`${atom.element} (${atom.id}) 含自由基标记；当前赛道不接受自由基。`);
  }
  return warnings;
}

function validStereoReferences(project, bond) {
  const refs = bond.stereoAtoms;
  return refs?.length === 2 && !refs.includes(bond.a) && !refs.includes(bond.b)
    && !!bondBetween(project, bond.a, refs[0]) && !!bondBetween(project, bond.b, refs[1]);
}
function neighborIds(project, atomId) {
  return project.bonds.filter(b => b.a === atomId || b.b === atomId).map(b => b.a === atomId ? b.b : b.a).sort();
}
function stereoNeighborhood(project, centerId) {
  const properties = atom => [atom.element, atom.charge, atom.isotope ?? 0, atom.explicitHs ?? 0, atom.noImplicit ?? false];
  // CW/CCW tags are tied to the order in which a toolkit receives incident
  // bonds. IDs are stable but a delete/re-add changes that neighbor order.
  // This prototype clears the designation rather than guessing new parity.
  return JSON.stringify([properties(atomById(project, centerId)), ...project.bonds.filter(b => b.a === centerId || b.b === centerId).map(b => {
    const neighbor = b.a === centerId ? b.b : b.a;
    return [neighbor, b.order, properties(atomById(project, neighbor))];
  })]);
}
function chiralityVolume(project, centerId, neighbors) {
  if (neighbors.length !== 3 && neighbors.length !== 4) return null;
  const positions = neighbors.map(i => atomById(project, i).pos3);
  const origin = neighbors.length === 4 ? positions[3] : atomById(project, centerId).pos3;
  const vectors = positions.slice(0, 3).map(p => sub(p, origin));
  const denominator = vectors.reduce((value, v) => value * norm(v), 1);
  return denominator < 1e-12 ? 0 : dot(vectors[0], cross(vectors[1], vectors[2])) / denominator;
}
/** Preserve stereo tags; announce explicit topology changes or geometric inversion. */
function reconcileStereo(previous, next) {
  if (!Array.isArray(next.atoms) || !Array.isArray(next.bonds)) return;
  const notices = [];
  for (const bond of next.bonds) {
    const old = previous.bonds.find(b => b.id === bond.id);
    // Repair only formerly valid data made stale by a topology edit. A new,
    // incorrectly specified stereo relationship is rejected by validation.
    if (old && old.stereoAtoms?.length === 2 && JSON.stringify(old.stereoAtoms) === JSON.stringify(bond.stereoAtoms)
      && validStereoReferences(previous, old) && (!validStereoReferences(next, bond) || (old.order === 2 && bond.order !== 2))) {
      bond.stereoAtoms = []; bond.stereo = 'STEREONONE';
      notices.push('拓扑修改使一处双键立体参考失效，已清除该关系；请重新确认。');
    } else if (old && previous.has3d && next.has3d && old.stereo === bond.stereo && bond.stereo && !['STEREONONE', 'STEREOANY'].includes(bond.stereo)
      && validStereoReferences(previous, old) && validStereoReferences(next, bond)) {
      try {
        const before = Math.cos(radians(dihedral(previous, old.stereoAtoms[0], old.a, old.b, old.stereoAtoms[1])));
        const after = Math.cos(radians(dihedral(next, bond.stereoAtoms[0], bond.a, bond.b, bond.stereoAtoms[1])));
        if (before * after < -0.0025) notices.push('三维移动使双键两侧构象翻转；E/Z 标记保持不变，请确认立体关系或重新生成三维。');
      } catch { /* Degenerate geometry remains editable and is checked on export. */ }
    }
  }
  for (const atom of next.atoms) {
    const old = previous.atoms.find(a => a.id === atom.id);
    if (!old || !old.chiralTag || old.chiralTag === 'CHI_UNSPECIFIED' || atom.chiralTag !== old.chiralTag) continue;
    const before = neighborIds(previous, atom.id), after = neighborIds(next, atom.id);
    if (stereoNeighborhood(previous, atom.id) !== stereoNeighborhood(next, atom.id)) {
      atom.chiralTag = 'CHI_UNSPECIFIED';
      // Wedge/dash bonds could otherwise silently reconstruct the old or
      // opposite assignment when exporting the newly ordered graph as 2D.
      for (const bond of next.bonds) if ((bond.a === atom.id || bond.b === atom.id) && ['BEGINWEDGE', 'BEGINDASH'].includes(bond.direction)) bond.direction = 'NONE';
      notices.push(`原子 ${atom.id} 的手性中心邻接关系或邻居顺序已改变，旧手性标记及相关楔形键已清除；请重新指定手性。`);
    } else if (previous.has3d && next.has3d) {
      const v0 = chiralityVolume(previous, atom.id, before), v1 = chiralityVolume(next, atom.id, after);
      if (v0 !== null && Math.abs(v0) > 1e-5 && v1 !== null && (v0 * v1 < 0 || Math.abs(v1) <= 1e-5)) {
        notices.push(`原子 ${atom.id} 的三维手性方向已翻转或变为共面；请重新确认手性或重建三维构象。`);
      }
    }
  }
  if (notices.length) { next.metadata ??= {}; next.metadata.stereoNotice = [...new Set(notices)].join(' '); }
}

/** Expand whole rigid fragments until a fixed point, including overlapping custom fragments. */
export function selectedIds(project, ids, rigid = true) {
  const available = new Set(project.atoms.map(atom => atom.id));
  const selected = new Set([...ids].filter(atomId => available.has(atomId)));
  if (rigid) {
    let changed;
    do {
      changed = false;
      for (const fragment of project.fragments) if (fragment.rigid && fragment.atomIds.some(atomId => selected.has(atomId))) {
        for (const atomId of fragment.atomIds) if (!selected.has(atomId)) { selected.add(atomId); changed = true; }
      }
      for (const group of project.metadata?.ketcherGroups ?? []) if (!group.expanded && group.atomIds.some(atomId => selected.has(atomId))) {
        for (const atomId of group.atomIds) if (!selected.has(atomId)) { selected.add(atomId); changed = true; }
      }
    } while (changed);
  }
  return selected;
}

export function distance(project, a, b) { return norm(sub(atomById(project, a).pos3, atomById(project, b).pos3)); }
export function angle(project, a, b, c) {
  const center = atomById(project, b).pos3;
  return degrees(Math.acos(Math.max(-1, Math.min(1, dot(unit(sub(atomById(project, a).pos3, center)), unit(sub(atomById(project, c).pos3, center)))))));
}
export function dihedral(project, a, b, c, d) {
  const axis = unit(sub(atomById(project, c).pos3, atomById(project, b).pos3));
  const ab = sub(atomById(project, a).pos3, atomById(project, b).pos3), cd = sub(atomById(project, d).pos3, atomById(project, c).pos3);
  const v = unit(sub(ab, scale(axis, dot(ab, axis)))), w = unit(sub(cd, scale(axis, dot(cd, axis))));
  return degrees(Math.atan2(dot(cross(axis, v), w), dot(v, w)));
}
export const measures = {distance, angle, dihedral};

/** Find the moving side after cutting one existing bond; ring cuts are ambiguous and rejected. */
function cutComponent(project, fixed, moving) {
  ensure(bondBetween(project, fixed, moving), '所选原子之间没有化学键。');
  const neighbors = new Map(project.atoms.map(a => [a.id, []]));
  for (const bond of project.bonds) {
    if ((bond.a === fixed && bond.b === moving) || (bond.a === moving && bond.b === fixed)) continue;
    neighbors.get(bond.a).push(bond.b); neighbors.get(bond.b).push(bond.a);
  }
  const result = new Set([moving]), queue = [moving];
  for (let i = 0; i < queue.length; ++i) for (const next of neighbors.get(queue[i])) if (!result.has(next)) { result.add(next); queue.push(next); }
  ensure(!result.has(fixed), '这条键位于环中，无法唯一分离可移动一侧。请直接选择原子移动，或选择连接两个片段的非环键。');
  return result;
}
/** Pure ring preflight/highlight helper. No coordinates, IDs or history change. */
export function movableSide(project, fixed, moving) {
  atomById(project, fixed); atomById(project, moving);
  return cutComponent(project, fixed, moving);
}
function rotateAtoms(project, ids, axis, value, origin) {
  for (const atom of project.atoms) if (ids.has(atom.id)) atom.pos3 = add(origin, rotateVector(sub(atom.pos3, origin), axis, value));
}

/** Preview a constrained geometric edit from the original drag snapshot.
 * reverse moves the opposite connected side. Nothing is committed here;
 * callers can discard this object on cancel or apply it once on release.
 */
export function geometryPreview(project, type, orderedIds, value, {reverse = false} = {}) {
  const arity = {distance: 2, angle: 3, dihedral: 4}[type];
  ensure(arity, '未知的几何编辑类型。');
  ensure(Array.isArray(orderedIds) && orderedIds.length === arity && new Set(orderedIds).size === arity,
    `请按顺序选择 ${arity} 个不同、连续成键的原子。`);
  finite(value, {distance: '键长', angle: '键角', dihedral: '二面角'}[type]);
  const source = validateProject(project), draft = clone(source), ids = reverse ? [...orderedIds].reverse() : [...orderedIds];
  ids.forEach(atomId => atomById(source, atomId));
  let moving;
  if (type === 'distance') {
    ensure(value > 0, '键长必须大于零。');
    const [a, b] = ids;
    moving = movableSide(source, a, b);
    const vector = sub(atomById(source, b).pos3, atomById(source, a).pos3);
    const delta = scale(norm(vector) > 1e-10 ? unit(vector) : [1, 0, 0], value - norm(vector));
    for (const atom of draft.atoms) if (moving.has(atom.id)) atom.pos3 = add(atom.pos3, delta);
  } else if (type === 'angle') {
    ensure(value >= 0 && value <= 180, '键角应位于 0–180°。');
    const [a, b, c] = ids;
    ensure(bondBetween(source, a, b), '请按顺序选择三个不同、连续成键的原子。');
    moving = movableSide(source, b, c);
    ensure(!moving.has(a), '所选键角的固定原子处于可移动一侧。');
    const origin = atomById(source, b).pos3, va = unit(sub(atomById(source, a).pos3, origin)), vc = unit(sub(atomById(source, c).pos3, origin));
    let axis = cross(va, vc);
    if (norm(axis) < 1e-10) axis = cross(va, Math.abs(va[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]);
    rotateAtoms(draft, moving, unit(axis), radians(value - angle(source, a, b, c)), origin);
  } else {
    const [a, b, c, d] = ids;
    ensure(bondBetween(source, a, b) && bondBetween(source, c, d), '请按顺序选择四个不同、连续成键的原子。');
    moving = movableSide(source, b, c);
    ensure(!moving.has(a) && moving.has(d), '二面角两端必须位于中间键的不同侧。');
    const origin = atomById(source, b).pos3, axis = unit(sub(atomById(source, c).pos3, origin));
    rotateAtoms(draft, moving, axis, radians(value - dihedral(source, a, b, c, d)), origin);
  }
  draft.has3d = true;
  reconcileStereo(source, draft);
  return {project: validateProject(draft), movingIds: [...moving], fixedIds: source.atoms.filter(a => !moving.has(a.id)).map(a => a.id)};
}

function assertConnection(project, a, b, order) {
  ensure(a !== b, '请选择两个不同原子的连接位点。'); enumeration(order, ORDERS, '连接键级');
  ensure(!bondBetween(project, a, b), '这两个原子已经成键；请使用修改键级。');
  for (const atomId of [a, b]) {
    const atom = atomById(project, atomId), maximum = maximumValence(atom);
    ensure(maximum !== null, `${atom.element} 的连接价态暂不支持。`);
    ensure(valence(project, atomId) + order <= maximum + 1e-8, `${atom.element} 连接位点没有足够价态（现有 ${valence(project, atomId)}，添加 ${order}，常见上限 ${maximum}）。请先移除显式氢、调整键级或选择其他位点。`);
    for (const fragment of project.fragments.filter(f => f.atomIds.includes(atomId))) {
      const port = fragment.ports.find(p => p.atomId === atomId);
      ensure(port, `所选原子不是“${fragment.name}”声明的连接位点。`);
      ensure(port.maxOrder === undefined || order <= port.maxOrder, `“${fragment.name}”此连接位点最多允许 ${port.maxOrder} 级键。`);
    }
  }
}

const COVALENT_RADII = {H: .31, B: .84, C: .76, N: .71, O: .66, F: .57, Si: 1.11, P: 1.07, S: 1.05, Cl: 1.02, Br: 1.20, I: 1.39};
function estimatedBondLength(a, b, order) {
  return ((COVALENT_RADII[a.element] ?? .77) + (COVALENT_RADII[b.element] ?? .77)) * ({1: 1, 1.5: .93, 2: .90, 3: .82}[order]);
}
function centroid(atoms, key) {
  return scale(atoms.reduce((sum, a) => add(sum, a[key]), atoms[0][key].map(() => 0)), 1 / atoms.length);
}
function outward(project, atomId, key) {
  const atom = atomById(project, atomId), neighbors = neighborIds(project, atomId).map(i => atomById(project, i));
  let vector = neighbors.length ? sub(atom[key], centroid(neighbors, key)) : atom[key].map(() => 0);
  if (norm(vector) < 1e-8) vector = sub(atom[key], centroid(project.atoms, key));
  return norm(vector) < 1e-8 ? atom[key].map((_, i) => i === 0 ? 1 : 0) : unit(vector);
}
function alignedVector(vector, from, to, dimension) {
  from = unit(from); to = unit(to);
  if (dimension === 2) {
    const turn = Math.atan2(to[1], to[0]) - Math.atan2(from[1], from[0]), c = Math.cos(turn), s = Math.sin(turn);
    return [c * vector[0] - s * vector[1], s * vector[0] + c * vector[1]];
  }
  const cosine = Math.max(-1, Math.min(1, dot(from, to)));
  if (cosine > 1 - 1e-12) return [...vector];
  const axis = cosine < -1 + 1e-12 ? unit(cross(from, Math.abs(from[0]) < .9 ? [1, 0, 0] : [0, 1, 0])) : unit(cross(from, to));
  return rotateVector(vector, axis, Math.acos(cosine));
}
function placementNamespace(project, requested) {
  const allIds = [...project.atoms, ...project.bonds, ...project.fragments].map(item => item.id);
  const free = namespace => !allIds.some(existing => existing === namespace || existing.startsWith(`${namespace}:`));
  if (requested !== undefined) {
    identifier(requested, '片段放置 ID'); ensure(requested.length <= 150, '片段放置 ID 过长。');
    ensure(free(requested), '这次片段放置已经应用，或其 ID 已被使用；请开始新的放置操作。');
    return requested;
  }
  let suffix = 1;
  while (!free(`placement-${suffix}`)) ++suffix;
  return `placement-${suffix}`;
}

/** Build an atomic insertion/connection preview without mutating the template.
 * Repeated calls against the same base graph use identical IDs. For a drag,
 * retain placementId (or this returned value) until apply/cancel. pos2/pos3
 * locate the chosen port; portless free fragments use their centroid instead.
 * Only the active layout rotates. The inactive layout is translated only when
 * explicitly supplied, and no existing atom coordinates are ever overwritten.
 */
export function fragmentPlacementPreview(project, template, options = {}) {
  const {targetId, portIndex = 0, pos2, pos3, dimension = 3, order = 1, rigid = template.rigid ?? true} = options;
  ensure(dimension === 2 || dimension === 3, '片段只能放置到二维或三维视图。');
  if (pos2 !== undefined) point(pos2, 2, '片段二维位置');
  if (pos3 !== undefined) point(pos3, 3, '片段三维位置');
  enumeration(order, ORDERS, '连接键级');
  const base = validateProject(project), draft = clone(base);
  const source = validateProject({...createProject(template.name ?? '自定义片段'), ...clone(template), format: PROJECT_FORMAT, fragments: [], has3d: template.has3d ?? true});
  ensure(source.atoms.length > 0, '不能插入空片段。');
  const sourcePorts = template.ports ?? template.fragments?.flatMap(f => f.ports) ?? [];
  ensure(Array.isArray(sourcePorts), '片段连接位点必须是列表。');
  ensure(Number.isSafeInteger(portIndex) && portIndex >= 0 && (sourcePorts.length === 0 ? portIndex === 0 : portIndex < sourcePorts.length), '请选择存在的片段连接位点。');
  ensure(!targetId || sourcePorts.length > 0, '此片段没有声明连接位点；可自由放置，但不能直接连接。');
  const sourcePort = sourcePorts[portIndex], sourceAnchor = sourcePort ? atomById(source, sourcePort.atomId) : null;
  const placementId = placementNamespace(base, options.placementId);
  const idMap = Object.fromEntries(source.atoms.map((a, index) => [a.id, `${placementId}:a:${index}`]));
  const atoms = source.atoms.map(a => ({...clone(a), id: idMap[a.id]}));
  const bonds = source.bonds.map((b, index) => ({...clone(b), id: `${placementId}:b:${index}`, a: idMap[b.a], b: idMap[b.b], ...(b.stereoAtoms ? {stereoAtoms: b.stereoAtoms.map(a => idMap[a])} : {})}));
  const ports = sourcePorts.map(p => ({...clone(p), atomId: idMap[p.atomId]}));
  const fragment = {id: `${placementId}:fragment`, name: source.name, atomIds: atoms.map(a => a.id), ports, rigid};
  draft.atoms.push(...atoms); draft.bonds.push(...bonds); draft.fragments.push(fragment);
  let bondId;
  if (targetId) {
    const newPortId = idMap[sourcePort.atomId];
    assertConnection(draft, targetId, newPortId, order);
    bondId = `${placementId}:attachment`;
    draft.bonds.push({id: bondId, a: targetId, b: newPortId, order});
  }
  for (const currentDimension of [2, 3]) {
    const key = currentDimension === 2 ? 'pos2' : 'pos3', requested = currentDimension === 2 ? pos2 : pos3;
    const anchor = sourceAnchor ? sourceAnchor[key] : centroid(source.atoms, key);
    if (currentDimension === dimension && targetId) {
      const target = atomById(base, targetId);
      const proposed = requested ? sub(requested, target[key]) : outward(base, targetId, key);
      const direction = norm(proposed) > 1e-8 ? unit(proposed) : outward(base, targetId, key);
      const at = add(target[key], scale(direction, estimatedBondLength(target, sourceAnchor, order)));
      const from = outward(source, sourcePort.atomId, key), to = scale(direction, -1);
      for (const atom of atoms) atom[key] = add(at, alignedVector(sub(atom[key], anchor), from, to, currentDimension));
    } else if (requested) {
      for (const atom of atoms) atom[key] = add(atom[key], sub(requested, anchor));
    }
  }
  if (source.has3d) draft.has3d = true;
  reconcileStereo(base, draft);
  trackGeometryPending(base, draft);
  return {project: validateProject(draft), fragmentId: fragment.id, atomIds: [...fragment.atomIds], ports: clone(ports), idMap, bondId, placementId};
}

export class WorkshopModel extends EventTarget {
  constructor(project = createProject(), {historyLimit = 150} = {}) {
    super(); this.project = validateProject(project); this.selection = new Set(); this.revision = 0;
    this.historyLimit = historyLimit; this._undo = []; this._redo = [];
  }
  get canUndo() { return this._undo.length > 0; }
  get canRedo() { return this._redo.length > 0; }
  get undoLabel() { return this._undo.at(-1)?.label ?? ''; }
  get redoLabel() { return this._redo.at(-1)?.label ?? ''; }
  expandedSelection(ids = this.selection) { return selectedIds(this.project, ids, true); }
  _snapshot(label) { return {project: clone(this.project), selection: [...this.selection], label}; }
  _emit(label, kind) {
    const event = new Event('change');
    Object.defineProperty(event, 'detail', {value: {label, revision: this.revision, kind}});
    this.dispatchEvent(event);
  }
  commit(label, mutator, {incomingProject = false, geometryConfirmedIds = []} = {}) {
    const draft = clone(this.project), result = mutator(draft);
    if (!incomingProject) {
      reconcileGeneratedHydrogens(this.project, draft);
      validateNativeGroupConnections(this.project, draft);
      reconcileStereo(this.project, draft);
      trackGeometryPending(this.project, draft);
      // Explicit 3D construction / hand placement confirms those coordinates
      // in this same history transaction. A later view switch must not rebuild
      // them merely because they originally entered the graph as new atoms.
      ensure(Array.isArray(geometryConfirmedIds), '已确认构象的原子 ID 必须为数组。');
      if (geometryConfirmedIds.length && draft.has3d) {
        const available = new Set(draft.atoms.map(atom => atom.id));
        ensure(geometryConfirmedIds.every(id => available.has(id)), '已确认构象引用了不存在的原子。');
        const pending = draft.metadata?.geometryPending;
        if (pending) {
          const confirmed = new Set(geometryConfirmedIds);
          pending.atomIds = pending.atomIds.filter(id => !confirmed.has(id));
          if (!pending.atomIds.length && !pending.changedAtomIds.length && !pending.editedAtomIds?.length) delete draft.metadata.geometryPending;
        }
      }
    }
    pruneEditorMetadata(draft);
    const next = validateProject(draft);
    if (JSON.stringify(next) === JSON.stringify(this.project)) return result;
    this._undo.push(this._snapshot(label));
    if (this._undo.length > this.historyLimit) this._undo.shift();
    this._redo = []; this.project = next;
    const available = new Set(next.atoms.map(atom => atom.id));
    this.selection = new Set([...this.selection].filter(atomId => available.has(atomId)));
    ++this.revision; this._emit(label, 'commit'); return result;
  }
  replace(project, label = '载入项目') {
    const next = validateProject(project);
    return this.commit(label, draft => { for (const key of Object.keys(draft)) delete draft[key]; Object.assign(draft, next); }, {incomingProject: true});
  }
  undo() {
    if (!this.canUndo) return false;
    const previous = this._undo.pop(); this._redo.push(this._snapshot(previous.label));
    this.project = previous.project; this.selection = new Set(previous.selection); ++this.revision;
    this._emit(previous.label, 'undo'); return true;
  }
  redo() {
    if (!this.canRedo) return false;
    const next = this._redo.pop(); this._undo.push(this._snapshot(next.label));
    this.project = next.project; this.selection = new Set(next.selection); ++this.revision;
    this._emit(next.label, 'redo'); return true;
  }
  applyAsync(baseRevision, project, label = '应用几何结果') {
    if (baseRevision !== this.revision) return false;
    this.replace(project, label); return true;
  }
  addAtom(options = {}) {
    const atom = {id: id('a'), element: 'C', charge: 0, pos2: [0, 0], pos3: [0, 0, 0], ...clone(options)};
    this.commit('添加原子', draft => draft.atoms.push(atom)); return atom.id;
  }
  addBond(a, b, order = 1, extra = {}) {
    const bond = {id: id('b'), a, b, order, ...clone(extra)};
    this.commit('添加化学键', draft => draft.bonds.push(bond)); return bond.id;
  }
  /** A canvas drag adds zero, one or two atoms and one bond in one history step.
   * Existing endpoints retain element, stable ID and both layouts. Re-drawing
   * an existing bond changes its order/direction instead of duplicating it.
   */
  drawBondGesture({startId, endId, startPoint, endPoint, dimension = 2, element = 'C', order = 1, direction = 'NONE'} = {}) {
    ensure(dimension === 2 || dimension === 3, '只能在二维或三维视图绘制化学键。');
    ensure(ELEMENTS.has(element), `不支持的元素符号：${element}`);
    enumeration(order, ORDERS, '键级'); enumeration(direction, DIRECTION, '键方向');
    const key = dimension === 2 ? 'pos2' : 'pos3', otherKey = dimension === 2 ? 'pos3' : 'pos2';
    const oldStart = startId ? atomById(this.project, startId) : null, oldEnd = endId ? atomById(this.project, endId) : null;
    const start = oldStart?.[key] ?? point(startPoint, dimension, '起点位置');
    const end = oldEnd?.[key] ?? point(endPoint, dimension, '终点位置');
    ensure(!startId || !endId || startId !== endId, '一个原子不能与自己成键。');
    const created = [];
    const make = (position, opposite) => {
      const atom = {id: id('a'), element, charge: 0, pos2: position.slice(0, 2), pos3: [...position.slice(0, 2), dimension === 3 ? position[2] : 0]};
      if (opposite) {
        const delta = sub(position, opposite[key]);
        atom[otherKey] = add(opposite[otherKey], dimension === 2 ? [...delta, 0] : delta.slice(0, 2));
      }
      created.push(atom); return atom.id;
    };
    startId ??= make(start, oldEnd);
    endId ??= make(end, oldStart);
    const existing = bondBetween(this.project, startId, endId), bondId = existing?.id ?? id('b');
    this.commit('绘制化学键', draft => {
      draft.atoms.push(...created);
      if (existing) {
        const bond = draft.bonds.find(b => b.id === bondId);
        bond.order = order; bond.direction = direction;
        if (['BEGINWEDGE', 'BEGINDASH'].includes(direction)) {
          if (bond.a !== startId && bond.stereoAtoms?.length === 2) bond.stereoAtoms.reverse();
          bond.a = startId; bond.b = endId;
        }
      } else draft.bonds.push({id: bondId, a: startId, b: endId, order, direction});
      if (dimension === 3) draft.has3d = true;
    });
    return {startId, endId, bondId, createdAtomIds: created.map(a => a.id)};
  }
  updateAtom(atomId, patch) {
    ensure(patch.id === undefined || patch.id === atomId, '稳定原子 ID 不能修改。');
    this.commit('修改原子', draft => Object.assign(atomById(draft, atomId), clone(patch)));
  }
  updateBond(bondId, patch) {
    ensure(patch.id === undefined || patch.id === bondId, '稳定键 ID 不能修改。');
    this.commit('修改化学键', draft => {
      const bond = draft.bonds.find(b => b.id === bondId); ensure(bond, '找不到待修改的化学键。'); Object.assign(bond, clone(patch));
    });
  }
  deleteBond(bondId) { this.commit('删除化学键', draft => { draft.bonds = draft.bonds.filter(b => b.id !== bondId); }); }
  deleteSelection({rigid = false, ids = this.selection} = {}) {
    const removing = selectedIds(this.project, ids, rigid);
    this.commit('删除所选原子', draft => {
      draft.atoms = draft.atoms.filter(a => !removing.has(a.id));
      draft.bonds = draft.bonds.filter(b => !removing.has(b.a) && !removing.has(b.b));
      for (const bond of draft.bonds) if (bond.stereoAtoms?.some(a => removing.has(a))) {
        // Deleting an explicit stereo reference makes that relationship undefined.
        // Keep a visible audit note instead of retaining a dangling atom mapping.
        delete bond.stereoAtoms; bond.stereo = 'STEREONONE';
        draft.metadata ??= {}; draft.metadata.stereoNotice = '删除立体参考原子后，相应双键立体关系已清除；请重新确认。';
      }
      draft.fragments = draft.fragments.map(f => ({...f, atomIds: f.atomIds.filter(a => !removing.has(a)), ports: f.ports.filter(p => !removing.has(p.atomId))})).filter(f => f.atomIds.length);
    });
  }
  addFragment(template, {pos2 = [0, 0], pos3 = [0, 0, 0], rigid = template.rigid ?? true} = {}) {
    point(pos2, 2, '片段二维位置'); point(pos3, 3, '片段三维位置');
    const source = validateProject({...createProject(template.name ?? '自定义片段'), ...clone(template), format: PROJECT_FORMAT, fragments: [], has3d: template.has3d ?? true});
    ensure(source.atoms.length > 0, '不能插入空片段。');
    const idMap = Object.fromEntries(source.atoms.map(a => [a.id, id('a')]));
    const atoms = source.atoms.map(a => ({...a, id: idMap[a.id], pos2: add(a.pos2, pos2), pos3: add(a.pos3, pos3)}));
    const bonds = source.bonds.map(b => ({...b, id: id('b'), a: idMap[b.a], b: idMap[b.b], ...(b.stereoAtoms ? {stereoAtoms: b.stereoAtoms.map(a => idMap[a])} : {})}));
    const ports = (template.ports ?? template.fragments?.flatMap(f => f.ports) ?? []).map(p => ({...p, atomId: idMap[p.atomId]}));
    const fragment = {id: id('f'), name: source.name, atomIds: atoms.map(a => a.id), ports, rigid};
    this.commit(`添加${source.name}`, draft => { draft.atoms.push(...atoms); draft.bonds.push(...bonds); draft.fragments.push(fragment); if (source.has3d) draft.has3d = true; });
    return {fragmentId: fragment.id, atomIds: [...fragment.atomIds], ports: clone(ports), idMap};
  }
  placeFragment(template, options = {}) {
    const preview = fragmentPlacementPreview(this.project, template, options);
    this.commit(options.targetId ? `放置并连接${template.name ?? '片段'}` : `放置${template.name ?? '片段'}`, draft => Object.assign(draft, clone(preview.project)));
    return preview;
  }
  connect(a, b, order = 1) {
    assertConnection(this.project, a, b, order);
    return this.addBond(a, b, order);
  }
  moveSelection(delta, {dimension = 3, rigid = true, ids = this.selection} = {}) {
    ensure(dimension === 2 || dimension === 3, '只能移动二维或三维坐标。'); point(delta, dimension, '移动距离');
    const selected = selectedIds(this.project, ids, rigid), key = dimension === 2 ? 'pos2' : 'pos3';
    this.commit(dimension === 2 ? '调整二维排版' : '移动三维原子/片段', draft => {
      for (const atom of draft.atoms) if (selected.has(atom.id)) atom[key] = add(atom[key], delta);
      if (dimension === 3 && selected.size) draft.has3d = true;
    }, {geometryConfirmedIds: dimension === 3 ? [...selected] : []});
  }
  rotateSelection(axis, value, {origin, rigid = true, ids = this.selection} = {}) {
    finite(value, '旋转角度');
    if (typeof axis === 'string') axis = {x: [1, 0, 0], y: [0, 1, 0], z: [0, 0, 1]}[axis];
    const normalized = unit(point(axis, 3, '旋转轴')), selected = selectedIds(this.project, ids, rigid);
    ensure(selected.size > 0, '请先选择要旋转的原子或片段。');
    const center = origin ? point(origin, 3, '旋转中心') : scale([...selected].reduce((sum, a) => add(sum, atomById(this.project, a).pos3), [0, 0, 0]), 1 / selected.size);
    this.commit('旋转三维原子/片段', draft => { rotateAtoms(draft, selected, normalized, radians(value), center); draft.has3d = true; }, {geometryConfirmedIds: [...selected]});
  }
  setDistance(a, b, value, options = {}) {
    const preview = geometryPreview(this.project, 'distance', [a, b], value, options);
    this.commit('调整键长', draft => Object.assign(draft, clone(preview.project)), {geometryConfirmedIds: preview.movingIds});
  }
  setAngle(a, b, c, value, options = {}) {
    const preview = geometryPreview(this.project, 'angle', [a, b, c], value, options);
    this.commit('调整键角', draft => Object.assign(draft, clone(preview.project)), {geometryConfirmedIds: preview.movingIds});
  }
  setDihedral(a, b, c, d, value, options = {}) {
    const preview = geometryPreview(this.project, 'dihedral', [a, b, c, d], value, options);
    this.commit('调整二面角', draft => Object.assign(draft, clone(preview.project)), {geometryConfirmedIds: preview.movingIds});
  }
}
