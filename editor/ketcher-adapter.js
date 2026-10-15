/** Ketcher 3.18 live Struct ↔ canonical workshop graph.
 * Identity comes exclusively from an explicitly bound live Pool session.
 * KET is a transport for deliberate rebuilds, never an identity matching tool.
 * No user atom-map (aam) field is repurposed. Coordinates are in angstroms.
 */
import {clone, validateProject, valence, reconcileGeneratedHydrogens} from './model.js';
import {C60_POS3, C60_BONDS, C60_TEMPLATE_VERSION} from './native-c60-data.js';
import {C60_ADDITION_NOTICE,stripCompletedC60Notice} from './editor-notices.js';
const NativeError=globalThis.Clazz?._Error??globalThis.Error;

export const KETCHER_SIDECAR_FORMAT = 'molecular-car-ketcher-sidecar-v1';
const MARKER = '__molecularCarEditorUUID';
const DIRECTIONS = {0: 'NONE', 1: 'BEGINWEDGE', 4: 'UNKNOWN', 6: 'BEGINDASH', 3: 'EITHERDOUBLE'};
const TO_STEREO = {NONE: 0, BEGINWEDGE: 1, BEGINDASH: 6, UNKNOWN: 4, EITHERDOUBLE: 3, ENDUPRIGHT: 0, ENDDOWNRIGHT: 0};
const ORDERS = {1: 1, 2: 2, 3: 3, 4: 1.5};
let sequence = 0;
const uniqueId = prefix => `${prefix}-${globalThis.crypto?.randomUUID?.() ?? `${Date.now().toString(36)}-${++sequence}`}`;
const assert = (condition, message) => { if (!condition) throw new NativeError(message); };
const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
function finite(value, name) { assert(typeof value === 'number' && Number.isFinite(value), `${name} 不是有限数值。`); return value; }
function entries(pool, name) {
  assert(pool && typeof pool.entries === 'function', `Ketcher ${name} 不提供已验证的 live Pool API。`);
  const result = [...pool.entries()];
  assert(result.every(([key, value]) => Number.isSafeInteger(key) && value && typeof value === 'object'), `Ketcher ${name} 池记录不受支持。`);
  return result;
}
function stamp(object, stableId) {
  // Convenience for diagnostics only. Atom.clone/new Atom intentionally drops
  // this field; the pool-session table, not the property, owns identity.
  try { Object.defineProperty(object, MARKER, {value: stableId, writable: true, configurable: true, enumerable: false}); } catch { /* A sealed toolkit object still has the external table. */ }
}
function atomAttributes(atom) {
  return {
    label: atom.label, charge: atom.charge ?? 0, isotope: atom.isotope ?? 0,
    radical: atom.radical ?? 0, implicitHCount: atom.implicitHCount ?? null,
    aam: atom.aam ?? 0, alias: atom.alias ?? null,
    stereoParity: atom.stereoParity ?? 0, stereoLabel: atom.stereoLabel ?? null,
  };
}
function bondAttributes(bond) { return {begin: bond.begin, end: bond.end, type: bond.type, stereo: bond.stereo ?? 0}; }
function rejectUnsupportedStereoLabel(label) {
  // Pinned Ketcher 3.18 uses 'abs', '&<group number>' and 'or<group number>'.
  // Storing the latter in a UI sidecar cannot preserve their group semantics
  // in the canonical graph or MOL export. Reject before accepting any change.
  assert(label === undefined || label === null || label === '' || label === 'abs',
    '当前共享分子图与 MOL 导出不能无损保存 AND/OR 增强立体组或未知立体标签；此编辑未应用。请使用 ABS（绝对立体）或清除标签。');
}
function rejectUnsupported(struct) {
  for (const name of ['rgroups', 'rxnArrows', 'rxnPluses', 'simpleObjects', 'texts', 'images', 'multitailArrows']) {
    assert(!struct[name]?.size, `当前共享分子图不支持 Ketcher ${name}；此操作尚未应用，请回到普通原子和化学键编辑。`);
  }
  for (const [,group] of entries(struct.sgroups,'sgroups')) {
    assert(group.type==='SUP'&&!group.isMonomer&&!group.data?.class,'当前共享分子图只支持普通 SUP 缩写基团；聚合物、数据组和宏分子未应用。');
    assert(Array.isArray(group.atoms)&&group.atoms.length>0&&new Set(group.atoms).size===group.atoms.length,'SUP 缩写基团的原子引用无效；此编辑未应用。');
  }
  for (const [, atom] of entries(struct.atoms, 'atoms')) {
    rejectUnsupportedStereoLabel(atom.stereoLabel);
    const query = atom.queryProperties && Object.values(atom.queryProperties).some(v => v !== null && v !== undefined);
    assert(!atom.atomList && !query && !atom.hCount && !atom.ringBondCount && !atom.substitutionCount && !atom.unsaturatedAtom
      && !(atom.explicitValence >= 0) && !atom.rglabel && !atom.attachmentPoints,
    '当前共享分子图不支持查询原子、指定查询价态或 R 基团；此编辑未应用。');
    assert(!atom.radical || [2, 3].includes(atom.radical), '当前分子图只支持自由基电子数，不能无损表示 Ketcher 单重态自由基。');
  }
  for (const [, bond] of entries(struct.bonds, 'bonds')) {
    assert(Object.hasOwn(ORDERS, bond.type) && !bond.customQuery && !bond.topology && !bond.reactingCenterStatus,
      '当前共享分子图支持单键、双键、三键和芳香键；查询键、配位键或反应中心未应用。');
    assert(Object.hasOwn(DIRECTIONS, bond.stereo ?? 0), 'Ketcher 键立体编码不受支持；当前设计未被覆盖。');
  }
}
function canonicalComparable(project) {
  const p = clone(project);
  if (p.metadata) { delete p.metadata.ketcherAdapter; if (!Object.keys(p.metadata).length) delete p.metadata; }
  return p;
}
function maxValence(atom) {
  return {H: 1, B: atom.charge < 0 ? 4 : 3, C: atom.charge === 0 ? 4 : 3,
    N: atom.charge > 0 ? 4 : atom.charge < 0 ? 2 : 3, O: atom.charge > 0 ? 3 : atom.charge < 0 ? 1 : 2,
    F: 1, Si: 4, P: 5, S: 6, Cl: 1, Br: 1, I: 1}[atom.element];
}

/** Native Ketcher drawing must not bypass ports on workshop-managed fragments. */
export function validateFragmentConnections(previous, next) {
  const oldBonds = new Map(previous.bonds.map(b => [b.id, b]));
  for (const bond of next.bonds) {
    const old = oldBonds.get(bond.id);
    if (old && old.a === bond.a && old.b === bond.b && old.order === bond.order) continue;
    const restrictions = [];
    for (const fragment of next.fragments) for (const [inside, outside] of [[bond.a, bond.b], [bond.b, bond.a]]) {
      if (!fragment.atomIds.includes(inside) || fragment.atomIds.includes(outside)) continue;
      restrictions.push(fragment);
      const port = fragment.ports.find(p => p.atomId === inside);
      assert(port, `“${fragment.name}”所选原子不是声明的连接位点；请选择允许的位点。`);
      assert(port.maxOrder === undefined || bond.order <= port.maxOrder, `“${fragment.name}”此位点最多允许 ${port.maxOrder} 级连接键。`);
    }
    if (!restrictions.length) continue;
    for (const atomId of [bond.a, bond.b]) {
      const atom = next.atoms.find(a => a.id === atomId), maximum = maxValence(atom);
      assert(maximum !== undefined, `${atom.element} 的片段连接价态规则尚未提供；请先在独立结构中检查。`);
      assert(valence(next, atomId) <= maximum + 1e-8,
        `${atom.element} 连接后显式价态为 ${valence(next, atomId)}，超过常见上限 ${maximum}。C60 等饱和位点须先明确降低相邻键级或移除显式氢。`);
    }
  }
}

export class KetcherAdapter {
  constructor({makeId = uniqueId, toolVersion = '3.18.0'} = {}) {
    this.makeId = makeId; this.toolVersion = toolVersion; this.epoch = 0;
    this.atomPool = null; this.bondPool = null; this.groupPool=null; this.atoms = new Map(); this.bonds = new Map(); this.groups=new Map();
    this._pendingPushes = 0; this._commitDepth = 0; this._pushQueue = Promise.resolve();
    this.lastCommit = {changed: false, requiresPush: false, warnings: []};
  }
  get suppressing() { return this._pendingPushes > 0; }
  get applyingChange() { return this._commitDepth > 0; }
  selectionIds(struct, selection = {}) {
    assert(struct.atoms === this.atomPool && struct.bonds === this.bondPool, '选择来自尚未绑定的 Ketcher 池。');
    const result = new Set();
    for (const poolId of selection.atoms ?? []) { const record = this.atoms.get(poolId); if (record && struct.atoms.has(poolId)) result.add(record.id); }
    for (const poolId of selection.bonds ?? []) {
      const bond = struct.bonds.get(poolId); if (!bond) continue;
      for (const endpoint of [bond.begin, bond.end]) { const record = this.atoms.get(endpoint); if (record) result.add(record.id); }
    }
    for(const poolId of [...selection.sgroups??[],...selection.functionalGroups??[]])for(const atomId of this.groups.get(poolId)?.canonical?.atomIds??[])result.add(atomId);
    return result;
  }
  toolSelection(struct, canonicalIds, {bondId} = {}) {
    assert(struct.atoms === this.atomPool && struct.bonds === this.bondPool, '选择目标是尚未绑定的 Ketcher 池。');
    const selected = new Set(canonicalIds);
    return {atoms: [...struct.atoms.keys()].filter(poolId => selected.has(this.atoms.get(poolId)?.id)),
      bonds: bondId ? [...struct.bonds.keys()].filter(poolId => this.bonds.get(poolId)?.id === bondId) : []};
  }

  /** A single ordered KET molecule plus an out-of-band UUID reference table. */
  toKet(project) {
    const p = validateProject(project);
    // Indigo emits disconnected components consecutively. Declare that order
    // before import, using only the canonical edges and original record order.
    // This is an explicit transport table, not matching a returned/symmetric graph.
    const remaining=new Set(p.atoms.map(a=>a.id)),orderedAtoms=[],orderedBonds=[];
    for(const atom of p.atoms){
      if(!remaining.has(atom.id))continue;const members=new Set([atom.id]),todo=[atom.id];
      for(let i=0;i<todo.length;i++)for(const b of p.bonds){const next=b.a===todo[i]?b.b:b.b===todo[i]?b.a:null;if(next&&!members.has(next)){members.add(next);todo.push(next);}}
      orderedAtoms.push(...p.atoms.filter(a=>members.has(a.id)));orderedBonds.push(...p.bonds.filter(b=>members.has(b.a)));
      for(const id of members)remaining.delete(id);
    }
    const indices=new Map(orderedAtoms.map((a,i)=>[a.id,i]));
    const stored = p.metadata?.ketcherAdapter?.atomData ?? {};
    const atoms = orderedAtoms.map(atom => {
      const remembered = stored[atom.id] ?? {};
      rejectUnsupportedStereoLabel(remembered.stereoLabel);
      const output = {label: atom.element, location: [atom.pos2[0], atom.pos2[1], 0], charge: atom.charge};
      // In Ketcher an explicit isotope:0 is drawn as a real "0C" label;
      // natural abundance must use its absent/null default, not numeric zero.
      if (atom.isotope) output.isotope = atom.isotope;
      const atomMap = atom.aam ?? remembered.aam;
      if (atomMap) output.mapping = atomMap;
      const alias = atom.alias ?? remembered.alias;
      if (alias) output.alias = alias;
      // Computed implicitH is never a declaration. Explicit/noImplicit values
      // are kept in canonical state and only mirrored as a Ketcher constraint.
      if (atom.noImplicit || atom.explicitHs) output.implicitHCount = atom.explicitHs ?? 0;
      else if (remembered.implicitHCount !== undefined && remembered.implicitHCount !== null) output.implicitHCount = remembered.implicitHCount;
      const radical = atom.radicalElectrons ?? 0;
      assert(radical <= 2, '二维编辑器支持最多两个未成对电子。请调整自由基标记。');
      if (radical) output.radical = radical === 1 ? 2 : 3;
      // Ketcher parity is not RDKit's neighbor-order-dependent chiralTag.
      // Reuse only a recorded toolkit declaration; never equate the numbers.
      const sameTag = remembered.canonicalChiralTag === (atom.chiralTag ?? 'CHI_UNSPECIFIED');
      if (sameTag && remembered.stereoParity) output.stereoParity = remembered.stereoParity;
      if (sameTag && remembered.stereoLabel) output.stereoLabel = remembered.stereoLabel;
      return output;
    });
    const bonds = orderedBonds.map(bond => ({type: bond.order === 1.5 ? 4 : bond.order,
      atoms: [indices.get(bond.a), indices.get(bond.b)], stereo: TO_STEREO[bond.direction ?? 'NONE'] ?? 0}));
    const ket = p.atoms.length ? {root: {nodes: [{$ref: 'mol0'}]}, mol0: {type: 'molecule', atoms, bonds}} : {root: {nodes: []}};
    const groups=p.metadata?.ketcherGroups??[];
    if(groups.length)ket.mol0.sgroups=groups.map(group=>{
      assert(group.type==='SUP'&&group.atomIds.every(id=>indices.has(id)), '缩写基团引用失效；请恢复或移除缩写。');
      const result={type:'SUP',name:group.name,atoms:group.atomIds.map(id=>indices.get(id)),expanded:group.expanded??false};
      if(group.attachmentPoints?.length)result.attachmentPoints=group.attachmentPoints.map(point=>{
        const output={attachmentAtom:indices.get(point.atomId),attachmentId:point.attachmentId??'1'};
        if(point.leavingAtomId!==undefined)output.leavingAtom=indices.get(point.leavingAtomId);
        return output;
      });
      return result;
    });
    return {ket, reference: {format: KETCHER_SIDECAR_FORMAT, atomIds: orderedAtoms.map(a => a.id), bondIds: orderedBonds.map(b => b.id),
      bondEndpoints: orderedBonds.map(b => [b.a, b.b]),groupIds:groups.map(g=>g.id),groupAtomIds:groups.map(g=>g.atomIds)}};
  }

  /** Call after a known ordered KET import, never after an arbitrary tool load.
   * Optional pool ID arrays are explicit importer output, not guessed mapping.
   */
  bindRebuild(struct, project, reference) {
    rejectUnsupported(struct);
    const p = validateProject(project), atomEntries = entries(struct.atoms, 'atoms'), bondEntries = entries(struct.bonds, 'bonds');
    assert(reference?.format === KETCHER_SIDECAR_FORMAT, '重建 Ketcher 池必须提供此适配器生成的显式 UUID reference。');
    const atomKeys = reference.atomPoolIds ?? atomEntries.map(([key]) => key), bondKeys = reference.bondPoolIds ?? bondEntries.map(([key]) => key);
    assert(reference.atomIds.length === p.atoms.length && reference.atomIds.length === atomEntries.length && new Set(reference.atomIds).size === reference.atomIds.length && atomKeys.length === atomEntries.length && new Set(atomKeys).size === atomKeys.length,
      'Ketcher 原子数量或导入引用与源项目不一致，不能猜测对称分子的身份。');
    assert(reference.bondIds.length === p.bonds.length && reference.bondIds.length === bondEntries.length && new Set(reference.bondIds).size === reference.bondIds.length && bondKeys.length === bondEntries.length && new Set(bondKeys).size === bondKeys.length,
      'Ketcher 化学键数量或导入引用与源项目不一致。');
    const atomLookup = new Map(p.atoms.map(a => [a.id, a])), bondLookup = new Map(p.bonds.map(b => [b.id, b]));
    const newAtoms = new Map(), newBonds = new Map();
    for (let i = 0; i < atomKeys.length; ++i) {
      const key = atomKeys[i], object = struct.atoms.get(key), canonical = atomLookup.get(reference.atomIds[i]);
      assert(object && canonical, 'Ketcher 导入的原子引用不完整。');
      assert(object.label === canonical.element || (canonical.element === 'H' && ['D', 'T'].includes(object.label)), 'Ketcher 导入顺序与显式原子引用不一致。');
      newAtoms.set(key, {id: canonical.id, canonical: clone(canonical), source: atomAttributes(object), acceptedSource: atomAttributes(object)});
    }
    for (let i = 0; i < bondKeys.length; ++i) {
      const key = bondKeys[i], object = struct.bonds.get(key), canonical = bondLookup.get(reference.bondIds[i]);
      assert(object && canonical, 'Ketcher 导入的键引用不完整。');
      assert(newAtoms.get(object.begin)?.id === canonical.a && newAtoms.get(object.end)?.id === canonical.b,
        'Ketcher 导入后的键端点与显式 UUID 引用不一致；拒绝按图相似性重建身份。');
      newBonds.set(key, {id: canonical.id, canonical: clone(canonical), source: bondAttributes(object)});
    }
    const groupEntries=entries(struct.sgroups,'sgroups'), canonicalGroups=p.metadata?.ketcherGroups??[],groupKeys=reference.groupPoolIds??groupEntries.map(([key])=>key);
    assert(groupEntries.length===canonicalGroups.length&&groupKeys.length===canonicalGroups.length,'Ketcher 缩写基团导入引用与共享项目不一致。');
    const newGroups=new Map();
    for(let i=0;i<groupKeys.length;i++){
      const key=groupKeys[i], object=struct.sgroups.get(key),canonical=canonicalGroups.find(g=>g.id===reference.groupIds[i]);
      assert(canonical&&same(object.atoms.map(pool=>newAtoms.get(pool)?.id),canonical.atomIds),'Ketcher SUP 有序原子引用与显式 UUID 引用不一致。');
      newGroups.set(key,{id:canonical.id,canonical:clone(canonical)});
      if(canonical.kind==='c60'&&canonical.template)object.data.molecularCarTemplate={id:'c60',version:canonical.template.version};
    }
    this.atomPool = struct.atoms; this.bondPool = struct.bonds; this.groupPool=struct.sgroups; this.atoms = newAtoms; this.bonds = newBonds;this.groups=newGroups; this.epoch++;
    for (const [key, object] of atomEntries) stamp(object, this.atoms.get(key).id);
    for (const [key, object] of bondEntries) stamp(object, this.bonds.get(key).id);
    return this.snapshotSidecar(struct);
  }

  snapshotSidecar(struct) {
    const atomIdsByPool = {}, bondIdsByPool = {}, groupIdsByPool={},atomData = {};
    for (const [poolId] of entries(struct.atoms, 'atoms')) {
      const record = this.atoms.get(poolId);
      assert(record, '此 Ketcher 原子尚未绑定 UUID。');
      atomIdsByPool[poolId] = record.id;
      Object.defineProperty(atomData, record.id, {value: {...clone(record.source), canonicalChiralTag: record.canonical?.chiralTag ?? 'CHI_UNSPECIFIED'}, enumerable: true, configurable: true, writable: true});
    }
    for (const [poolId] of entries(struct.bonds, 'bonds')) {
      const record = this.bonds.get(poolId); assert(record, '此 Ketcher 键尚未绑定 UUID。'); bondIdsByPool[poolId] = record.id;
    }
    for(const [poolId] of entries(struct.sgroups,'sgroups')){const record=this.groups.get(poolId);assert(record,'此 Ketcher SUP 尚未绑定 UUID。');groupIdsByPool[poolId]=record.id;}
    return {format: KETCHER_SIDECAR_FORMAT, toolVersion: this.toolVersion, epoch: this.epoch, atomIdsByPool, bondIdsByPool, groupIdsByPool,atomData};
  }

  read(struct, project) {
    assert(struct.atoms === this.atomPool && struct.bonds === this.bondPool,
      'Ketcher 的原子池已重建；必须从共享项目重新推送并显式绑定 UUID，不能按编号或对称图猜测身份。');
    rejectUnsupported(struct);
    this._lastReadToolChanges = false;
    this._lastReadNeedsPush=false;
    const p = clone(project), warnings = [], oldAtoms = new Map(p.atoms.map(a => [a.id, a])), oldBonds = new Map(p.bonds.map(b => [b.id, b]));
    const foundAtoms = new Map(), foundBonds = new Map(), newAtomIds = new Set(), changedStereoCenters = new Set();
    const atomSourceUpdates = new Map(), bondSourceUpdates = new Map();
    for (const [poolId, object] of entries(struct.atoms, 'atoms')) {
      let record = this.atoms.get(poolId);
      if (!record) { record = {id: this.makeId('a'), canonical: null, source: null}; this.atoms.set(poolId, record); }
      // New pool IDs always mean new identity, even if paste copied a custom
      // field, aam, identical coordinates or an entire symmetric C60 graph.
      stamp(object, record.id);
      const known = oldAtoms.get(record.id) ?? record.canonical;
      const source = atomAttributes(object), element = source.label === 'D' || source.label === 'T' ? 'H' : source.label;
      if (record.acceptedSource && ['aam', 'alias', 'stereoParity', 'stereoLabel'].some(key => source[key] !== record.acceptedSource[key])) this._lastReadToolChanges = true;
      const atom = known ? clone(known) : {id: record.id, element, charge: 0, pos2: [0, 0], pos3: [0, 0, 0], chiralTag: 'CHI_UNSPECIFIED'};
      atom.element = element; atom.charge = source.charge;
      if (source.aam || known?.aam !== undefined) atom.aam = source.aam;
      if (source.alias || known?.alias !== undefined) atom.alias = source.alias;
      const isotope = source.label === 'D' ? 2 : source.label === 'T' ? 3 : source.isotope;
      if (known?.isotope !== undefined || isotope) atom.isotope = isotope;
      atom.pos2 = [finite(object.pp?.x, 'Ketcher x'), -finite(object.pp?.y, 'Ketcher y')];
      if (!known) newAtomIds.add(record.id);
      if (!known || source.radical !== record.source?.radical) atom.radicalElectrons = source.radical === 2 ? 1 : source.radical === 3 ? 2 : 0;
      if (!known || source.implicitHCount !== record.source?.implicitHCount) {
        assert(source.implicitHCount === null || (Number.isSafeInteger(source.implicitHCount) && source.implicitHCount >= 0), '指定氢数必须是非负整数。');
        if (source.implicitHCount !== null || known) { atom.explicitHs = source.implicitHCount ?? 0; atom.noImplicit = source.implicitHCount !== null; }
      }
      if (record.source && (source.stereoParity !== record.source.stereoParity || source.stereoLabel !== record.source.stereoLabel)) {
        atom.chiralTag = 'CHI_UNSPECIFIED'; changedStereoCenters.add(atom.id);
        warnings.push('Ketcher 的立体声明已改变；已清除旧手性标记，请通过楔形键重新指定并显式生成三维。');
      }
      foundAtoms.set(atom.id, atom); atomSourceUpdates.set(poolId, source);
    }
    for (const [poolId, object] of entries(struct.bonds, 'bonds')) {
      let record = this.bonds.get(poolId);
      if (!record) { record = {id: this.makeId('b'), canonical: null, source: null}; this.bonds.set(poolId, record); }
      stamp(object, record.id);
      const a = this.atoms.get(object.begin)?.id, b = this.atoms.get(object.end)?.id;
      assert(foundAtoms.has(a) && foundAtoms.has(b), 'Ketcher 键引用了当前池之外的原子。');
      const known = oldBonds.get(record.id) ?? record.canonical, source = bondAttributes(object);
      const sameEndpoints = known && ((known.a === a && known.b === b) || (known.a === b && known.b === a));
      const bond = sameEndpoints ? clone(known) : {id: record.id, a, b, order: ORDERS[source.type], stereo: 'STEREONONE', direction: 'NONE', stereoAtoms: []};
      bond.order = ORDERS[source.type];
      const directionChanged = !record.source || source.stereo !== record.source.stereo || source.begin !== record.source.begin || source.end !== record.source.end;
      if (directionChanged) {
        bond.a = a; bond.b = b; bond.direction = DIRECTIONS[source.stereo];
        if (known && known.a !== a && bond.stereoAtoms?.length === 2) bond.stereoAtoms.reverse();
        if ([1, 6].includes(source.stereo) || [1, 6].includes(record.source?.stereo)) {
          changedStereoCenters.add(a); if (known) changedStereoCenters.add(known.a);
          warnings.push('楔形键已修改；旧手性标记已清除，二维声明会在显式生成三维时解析。');
        }
      }
      foundBonds.set(bond.id, bond); bondSourceUpdates.set(poolId, source);
    }
    // Keep established canonical ordering, including RDKit's neighbor order.
    p.atoms = [...project.atoms.filter(a => foundAtoms.has(a.id)).map(a => foundAtoms.get(a.id)), ...[...foundAtoms.values()].filter(a => !oldAtoms.has(a.id))];
    p.bonds = [...project.bonds.filter(b => foundBonds.has(b.id)).map(b => foundBonds.get(b.id)), ...[...foundBonds.values()].filter(b => !oldBonds.has(b.id))];
    for (const atom of p.atoms) if (changedStereoCenters.has(atom.id)) atom.chiralTag = 'CHI_UNSPECIFIED';
    for (const bond of p.bonds) if (bond.stereoAtoms?.length && !bond.stereoAtoms.every(id => foundAtoms.has(id))) {
      bond.stereoAtoms = []; bond.stereo = 'STEREONONE'; warnings.push('删除立体参考原子后已清除失效 E/Z 指定。');
    }
    p.fragments = p.fragments.map(f => ({...f, atomIds: f.atomIds.filter(id => foundAtoms.has(id)), ports: f.ports.filter(port => foundAtoms.has(port.atomId))})).filter(f => f.atomIds.length);
    this._readGroups(struct,p);
    this._initializeNew3d(p, newAtomIds);
    // Run the shared substitution preflight before testing external C60
    // valence. Otherwise a generated H would cause an early false rejection
    // before WorkshopModel can remove it in this same bond transaction.
    if(reconcileGeneratedHydrogens(project,p).length)this._lastReadNeedsPush=true;
    this._checkC60Connections(project,p,warnings);
    validateFragmentConnections(project, p);
    const persistent=warnings.filter(text=>text!==C60_ADDITION_NOTICE);
    if(p.metadata?.stereoNotice){p.metadata.stereoNotice=stripCompletedC60Notice(p.metadata.stereoNotice);if(!p.metadata.stereoNotice)delete p.metadata.stereoNotice;}
    if (persistent.length) { p.metadata ??= {}; p.metadata.stereoNotice = [...new Set(persistent)].join(' '); }
    this._lastReadWarnings=[...new Set(warnings)];
    // Stamp baselines only after a usable graph has been built. Deleted records
    // remain as tombstones for same-pool AtomAdd undo recreation.
    for (const [poolId, source] of atomSourceUpdates) { const record = this.atoms.get(poolId); record.source = source; record.canonical = clone(foundAtoms.get(record.id)); }
    for (const [poolId, source] of bondSourceUpdates) { const record = this.bonds.get(poolId); record.source = source; record.canonical = clone(foundBonds.get(record.id)); }
    p.metadata ??= {}; p.metadata.ketcherAdapter = this.snapshotSidecar(struct);
    return p;
  }

  _readGroups(struct,project){
    assert(struct.sgroups===this.groupPool,'Ketcher SUP 池已重建，必须显式绑定引用，不能猜测基团身份。');
    const result=[];
    for(const [poolId,object] of entries(struct.sgroups,'sgroups')){
      let record=this.groups.get(poolId);if(!record){record={id:this.makeId('g'),canonical:null};this.groups.set(poolId,record);}
      const group=record.canonical?clone(record.canonical):{id:record.id,type:'SUP'};
      group.name=String(object.data.name??'');group.atomIds=object.atoms.map(id=>this.atoms.get(id)?.id);
      assert(group.atomIds.every(id=>project.atoms.some(a=>a.id===id)),'SUP 缩写基团引用了当前池以外的原子。');
      group.expanded=Boolean(object.data.expanded);
      group.attachmentPoints=(object.getAttachmentPoints?.()??[]).map(point=>{
        const result={atomId:this.atoms.get(point.atomId)?.id,attachmentId:String(point.attachmentId??point.attachmentPointNumber??'1')};
        if(point.leaveAtomId!==undefined)result.leavingAtomId=this.atoms.get(point.leaveAtomId)?.id;
        assert(result.atomId&&group.atomIds.includes(result.atomId),'SUP 连接位点不属于基团。');return result;
      });
      const marker=object.data.molecularCarTemplate;
      if(!record.canonical&&marker?.id==='c60'&&marker.version===C60_TEMPLATE_VERSION){
        assert(group.atomIds.length===C60_POS3.length,'原生 C60 模板缺少明确的60原子引用。');
        // Structural assertions verify this source mapping; they never search
        // the symmetric graph for identity or choose an alternate atom order.
        for(const expected of C60_BONDS)assert(project.bonds.some(b=>(b.a===group.atomIds[expected.a]&&b.b===group.atomIds[expected.b])||(b.a===group.atomIds[expected.b]&&b.b===group.atomIds[expected.a])),'原生 C60 模板的显式索引/键引用不一致。');
        group.kind='c60';group.template={id:'c60',version:C60_TEMPLATE_VERSION,indexToAtomId:[...group.atomIds]};
      }
      if(group.template)assert(same(group.template.indexToAtomId,group.atomIds),'C60 模板成员顺序已改变，不能猜测对称原子身份。');
      record.canonical=clone(group);result.push(group);
    }
    project.metadata??={};if(result.length)project.metadata.ketcherGroups=result;else delete project.metadata.ketcherGroups;
  }

  _checkC60Connections(previous,next,warnings){
    const oldBonds=new Map(previous.bonds.map(b=>[b.id,b]));
    for(const group of next.metadata?.ketcherGroups??[]){
      if(group.kind!=='c60')continue;const members=new Set(group.atomIds);
      for(const bond of next.bonds){
        if(members.has(bond.a)===members.has(bond.b))continue;
        const old=oldBonds.get(bond.id);if(old&&old.a===bond.a&&old.b===bond.b&&old.order===bond.order)continue;
        const inside=members.has(bond.a)?bond.a:bond.b,outside=inside===bond.a?bond.b:bond.a;
        const anchor=next.atoms.find(a=>a.id===inside),other=next.atoms.find(a=>a.id===outside);
        assert(anchor.element==='C'&&anchor.charge===0&&bond.order===1,'C60 加成连接目前支持中性笼碳的外接单键；请展开基团后明确编辑其他加成。');
        if(valence(next,inside)>4){
          assert(Math.abs(valence(next,inside)-5)<1e-8,'C60 所选碳有多个过价连接；此连接未应用，请先展开修正。');
          const double=next.bonds.find(b=>b.order===2&&((b.a===inside&&members.has(b.b))||(b.b===inside&&members.has(b.a))));
          assert(double,'C60 位点已经满价且无可用于加成的内部双键；连接未应用。');
          double.order=1;double.direction='NONE';double.stereo='STEREONONE';double.stereoAtoms=[];this._lastReadNeedsPush=true;
          warnings.push(C60_ADDITION_NOTICE);
        }
        const maximum=maxValence(other);assert(maximum===undefined||valence(next,outside)<=maximum+1e-8,`${other.element} 的连接位点已满价；请移除该位点的显式氢或修改键级后再连接。`);
      }
    }
  }

  _initializeNew3d(project, pending) {
    const atoms = new Map(project.atoms.map(a => [a.id, a]));
    const remaining = new Set(pending);
    for(const group of project.metadata?.ketcherGroups??[]){
      if(group.kind!=='c60'||!group.template||!group.atomIds.every(id=>pending.has(id)))continue;
      const port=atoms.get(group.attachmentPoints?.[0]?.atomId??group.atomIds[0]);
      // Native source marker + ordered Pool references select this seed.
      // 2D layout never becomes a physical fullerene cage or rescales it.
      const center=[...port.pos2,0];
      group.template.indexToAtomId.forEach((id,index)=>{atoms.get(id).pos3=C60_POS3[index].map((v,axis)=>v+center[axis]);remaining.delete(id);});
    }
    while (remaining.size) {
      let progress = false;
      for (const atomId of remaining) {
        const atom = atoms.get(atomId);
        const neighbors = project.bonds.filter(b => b.a === atomId || b.b === atomId).map(b => b.a === atomId ? b.b : b.a);
        const anchor = neighbors.map(id => atoms.get(id)).find(a => !remaining.has(a.id));
        if (!anchor) continue;
        atom.pos3 = [anchor.pos3[0] + atom.pos2[0] - anchor.pos2[0], anchor.pos3[1] + atom.pos2[1] - anchor.pos2[1], anchor.pos3[2]];
        remaining.delete(atomId); progress = true;
      }
      if (!progress) { const atomId = remaining.values().next().value, atom = atoms.get(atomId); atom.pos3 = [...atom.pos2, 0]; remaining.delete(atomId); }
    }
    if (pending.size && project.has3d) { project.metadata ??= {}; project.metadata.geometryNotice = '新增原子的三维位置仅按显式邻接初始化；其余玩家构象保持不变，可显式整理几何。'; }
  }

  /** One accepted tool change → one canonical transaction; no-op signals do not create history. */
  commit(struct, model, label = 'Ketcher 二维编辑') {
    if (this.suppressing) return {changed: false, requiresPush: false, suppressed: true, warnings: []};
    const next = this.read(struct, model.project);
    if (!this._lastReadToolChanges && same(canonicalComparable(next), canonicalComparable(model.project))) return {changed: false, requiresPush: false, warnings: []};
    this._commitDepth++;
    try {
      model.commit(label, draft => { for (const key of Object.keys(draft)) delete draft[key]; Object.assign(draft, next); });
    } finally { this._commitDepth--; }
    const requiresPush = this._lastReadNeedsPush||!same(canonicalComparable(next), canonicalComparable(model.project));
    for (const record of this.atoms.values()) { const atom = model.project.atoms.find(a => a.id === record.id); if (atom) { record.canonical = clone(atom); record.acceptedSource = clone(record.source); } }
    for (const record of this.bonds.values()) { const bond = model.project.bonds.find(b => b.id === record.id); if (bond) record.canonical = clone(bond); }
    this.lastCommit = {changed: true, requiresPush, warnings: this._lastReadWarnings??[]};
    return this.lastCommit;
  }

  /** Serialized pushes prevent older asynchronous imports from rebinding a newer pool. */
  push(ketcher, project) {
    const p = validateProject(project), payload = this.toKet(p); this._pendingPushes++;
    const operation = this._pushQueue.catch(() => {}).then(async () => {
      assert(typeof ketcher.setMolecule === 'function' && typeof ketcher.editor?.struct === 'function', 'Ketcher 实例尚未就绪。');
      await ketcher.setMolecule(JSON.stringify(payload.ket), {preserveCanvasPosition: true, needZoom: false});
      const struct = ketcher.editor.struct(), sidecar = this.bindRebuild(struct, p, payload.reference);
      return {...payload, sidecar};
    }).finally(() => { this._pendingPushes--; });
    this._pushQueue = operation; return operation;
  }
}

export default KetcherAdapter;
