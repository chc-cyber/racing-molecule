/**
 * JSmol 16.3.31 bridge, independently browser-tested in tests/jsmol-live.cjs.
 * The host owns canonical topology/stereo/2D coordinates and undo history.
 * This bridge never derives bonds or identity from coordinates or SMILES.
 * Raw ModelKit topology edits are PROPOSALS, not silent canonical commits.
 * Native JSmol undo/redo/state restoration are unsupported: restore the host
 * snapshot with loadProject and its explicit row map instead.
 */
const clone=value=>structuredClone(value);
const AdapterError=globalThis.Clazz?._Error??globalThis.Error;
const fail=(condition,message)=>{if(!condition)throw new AdapterError(message);};
const pair=(a,b)=>JSON.stringify([a,b].sort());
const newId=prefix=>`${prefix}-${globalThis.crypto.randomUUID()}`;
const MAX_TOKEN=16777215; // exactly representable even by a Float32 property store
const number=value=>{fail(Number.isFinite(value),'Non-finite molecular coordinate.');return String(value);};
/** A view transport, not the project export. Never sanitizes or guesses bonds. */
export function projectToJSmolMol(project){
  const index=new Map(project.atoms.map((a,i)=>[a.id,i+1]));
  const lines=['Molecular Car view','  Workshop          3D','','  0  0  0     0  0            999 V3000','M  V30 BEGIN CTAB',`M  V30 COUNTS ${project.atoms.length} ${project.bonds.length} 0 0 0`,'M  V30 BEGIN ATOM'];
  for(const a of project.atoms){fail(/^[A-Z][a-z]?$/.test(a.element),'Invalid element symbol.');lines.push(`M  V30 ${index.get(a.id)} ${a.element} ${a.pos3.map(number).join(' ')} 0${a.charge?` CHG=${a.charge}`:''}${a.isotope?` MASS=${a.isotope}`:''}`);}
  lines.push('M  V30 END ATOM','M  V30 BEGIN BOND');
  for(const [i,b] of project.bonds.entries()){fail(index.has(b.a)&&index.has(b.b),'Missing bond endpoint.');lines.push(`M  V30 ${i+1} ${b.order===1.5?4:b.order} ${index.get(b.a)} ${index.get(b.b)}`);}
  lines.push('M  V30 END BOND','M  V30 END CTAB','M  END','');
  return{content:lines.join('\n'),indexToAtomId:project.atoms.map(a=>a.id)};
}
export class JSmolAdapter {
  constructor(applet,{Jmol=globalThis.Jmol,allocateId=newId}={}){
    fail(Jmol&&applet,'JSmol must be ready before creating the adapter.');
    this.Jmol=Jmol;this.applet=applet;this.allocateId=allocateId;
    this.tokenToId=new Map();this.idToToken=new Map();this.bondIds=new Map();
    this.nextToken=1;this.baseProject=null;this.baseRevision=null;
  }
  script(command){
    fail(!/(?:^|[;\n])\s*(?:undo|redo|restore\s+(?!orientation\b))/i.test(command),'Use canonical project history and loadProject; native JSmol history does not reliably restore atom identity.');
    return this._execute(command);
  }
  _execute(command){
    const output=this.Jmol.scriptWait(this.applet,command);
    this.lastScript={command,output:String(output)};
    fail(!/script ERROR|scriptError|ScriptException|unrecognized command|invalid argument/i.test(output),String(output));
    return output;
  }
  value(expression){return this.Jmol.evaluateVar(this.applet,expression);}
  configure({modelkit=true}={}){
    this.script(`set autobond false; modelkit addHydrogens false; modelkit autobond false; set modelKitMode ${modelkit?'true':'false'}; set picking dragAtom;`);
    const info=this.value('getProperty("modelkitInfo")');
    fail(info?.addHydrogens===false&&info?.autobond===false,'ModelKit automatic hydrogen/bond settings were not disabled.');
    return info;
  }
  tokenFor(id){
    fail(typeof id==='string'&&id.length>0,'Canonical atom ID must be a nonempty string.');
    if(!this.idToToken.has(id)){
      fail(this.nextToken<=MAX_TOKEN,'JSmol token space exhausted; recreate the adapter.');
      this.idToToken.set(id,this.nextToken);this.tokenToId.set(this.nextToken++,id);
    }
    return this.idToToken.get(id);
  }
  /** indexToAtomId MUST be the explicit MOL record order produced by the host. */
  loadProject(project,molText,{indexToAtomId,revision=null}={}){
    fail(Array.isArray(indexToAtomId)&&indexToAtomId.length===project.atoms.length,'Explicit MOL-index to stable-ID mapping is required.');
    fail(new Set(indexToAtomId).size===indexToAtomId.length,'MOL mapping IDs must be unique.');
    const canonical=new Map(project.atoms.map(a=>[a.id,a]));
    fail(indexToAtomId.every(id=>canonical.has(id)),'MOL mapping references an unknown atom.');
    fail(typeof molText==='string'&&!/^\s*END\s+"mcar_model"\s*$/mi.test(molText),'Invalid inline MOL data.');
    if(project.atoms.length)this._execute('set modelKitMode false; set autobond false; load DATA "mcar_model"\n'+molText+'\nEND "mcar_model";');
    else this.script('set modelKitMode false; set autobond false; zap;');
    // Enabling ModelKit on an empty scene creates a default methane molecule.
    this.configure({modelkit:project.atoms.length>0});
    const atoms=this.Jmol.getPropertyAsArray(this.applet,'atomInfo');
    fail(atoms.length===indexToAtomId.length,'JSmol changed the atom count during import.');
    this.script(atoms.map((a,i)=>`{atomindex=${a.atomIndex}}.property_editorid=${this.tokenFor(indexToAtomId[i])};`).join(''));
    this.baseProject=clone(project);this.baseRevision=revision;
    this.bondIds=new Map(project.bonds.map(b=>[pair(b.a,b.b),b.id]));
    // Upstream V2000 stereo handling can label a 3D dashed single bond as
    // aromaticSingle. The canonical graph is authoritative; restore only known
    // existing edges through explicit row identities, never invent missing edges.
    const loaded=this.read(),expectedBonds=new Map(project.bonds.map(b=>[pair(b.a,b.b),b]));
    fail(loaded.bonds.length===project.bonds.length,'JSmol changed the bond count during import.');
    for(const bond of loaded.bonds){
      const expected=expectedBonds.get(pair(bond.a,bond.b));fail(expected,'JSmol changed bond endpoints during import.');
      if(expected.order!==bond.order){
        const order={1:'single',2:'double',3:'triple',1.5:'aromatic'}[expected.order];fail(order,'Unsupported canonical bond order.');
        const ia=loaded.atoms.find(a=>a.id===bond.a).jmolIndex,ib=loaded.atoms.find(a=>a.id===bond.b).jmolIndex;
        this.script(`connect {atomindex=${ia}} {atomindex=${ib}} ${order} modify;`);
      }
    }
    // Identifiers come from the explicit row map, not from this chemistry check.
    this.assertSameChemistry(this.read(),project);
    const result=this.read();this.nativeBaseline=new Map(result.atoms.map(a=>[a.id,[...a.pos3]]));return result;
  }
  loadCanonical(project,{revision=null,preserveView=true}={}){
    if(preserveView)this.script('save orientation "mcar_camera";');
    const transport=projectToJSmolMol(project);
    const snapshot=this.loadProject(project,transport.content,{indexToAtomId:transport.indexToAtomId,revision});
    // Some tool readers round their text input. Restore original doubles through
    // explicit row identity; no geometry generation, optimization or bond search.
    if(snapshot.atoms.length)this.script(snapshot.atoms.map(a=>`{atomindex=${a.jmolIndex}}.xyz={${project.atoms.find(x=>x.id===a.id).pos3.map(number).join(' ')}};`).join(''));
    if(preserveView)this.script('restore orientation "mcar_camera" 0;');
    const result=this.read();this.nativeBaseline=new Map(result.atoms.map(a=>[a.id,[...a.pos3]]));return result;
  }
  /** Geometry-only acceptance keeps native transforms responsive without reload.
   * Exact project doubles remain authoritative; future gestures use displacement.
   */
  acceptCanonical(project,{revision=null}={}){
    const snapshot=this.read();this.assertSameChemistry(snapshot,project);
    this.baseProject=clone(project);this.baseRevision=revision;
    this.nativeBaseline=new Map(snapshot.atoms.map(a=>[a.id,[...a.pos3]]));
    return snapshot;
  }
  /** New ModelKit atoms have token 0. Allocate a new ID; never reuse deleted IDs. */
  read({assignNew=true}={}){
    const raw=this.Jmol.getPropertyAsArray(this.applet,'atomInfo');
    const rawBonds=this.Jmol.getPropertyAsArray(this.applet,'bondInfo');
    const tokens=raw.length?this.value('{*}.property_editorid.all'):[];
    fail(Array.isArray(tokens)&&tokens.length===raw.length,'JSmol property/atom arrays are inconsistent.');
    const indices=new Map(),seen=new Set();
    const atoms=raw.map((a,i)=>{
      let token=tokens[i];
      if((!token||!Number.isFinite(token))&&assignNew){const id=this.allocateId('a');token=this.tokenFor(id);this.script(`{atomindex=${a.atomIndex}}.property_editorid=${token};`);}
      fail(Number.isInteger(token)&&token>0&&this.tokenToId.has(token),'Missing or unknown JSmol atom token; reload through the explicit canonical mapping.');
      fail(!seen.has(token),'Duplicate JSmol atom token detected.');seen.add(token);
      const id=this.tokenToId.get(token);indices.set(a.atomIndex,id);
      const symbol=a.sym==='D'||a.sym==='T'?'H':a.sym.replace(/^\d+/, '');
      return{id,jmolIndex:a.atomIndex,token,element:symbol,atomicNumber:a.elemno,charge:a.formalCharge,pos3:[a.x,a.y,a.z]};
    });
    const bonds=rawBonds.map(b=>{
      const a=indices.get(b.atom1.atomIndex),c=indices.get(b.atom2.atomIndex);
      fail(a&&c,'JSmol bond references an unknown atom.');
      const key=pair(a,c);if(!this.bondIds.has(key))this.bondIds.set(key,this.allocateId('b'));
      // The MOL reader can display an aromatic bond as aromaticSingle/Double.
      // Those explicit native types retain aromatic semantics, not single/double.
      const order=/^aromatic/i.test(b.type)?1.5:b.order;
      return{id:this.bondIds.get(key),a,b:c,order,jmolIndex:b._bpt,type:b.type};
    });
    return{atoms,bonds,mapping:{jmolIndexToAtomId:Object.fromEntries(indices)},baseRevision:this.baseRevision};
  }
  assertSameChemistry(snapshot,project=this.baseProject){
    fail(project,'Load a canonical project before capturing changes.');
    fail(snapshot.atoms.length===project.atoms.length&&snapshot.bonds.length===project.bonds.length,'This operation changed topology. Review an edit proposal instead of committing it as geometry.');
    const expectedAtoms=new Map(project.atoms.map(a=>[a.id,a]));
    for(const atom of snapshot.atoms){const old=expectedAtoms.get(atom.id);fail(old&&old.element===atom.element&&(old.charge??0)===atom.charge,'This operation changed atom identity/element/charge; it is not a geometry-only edit.');}
    const expectedBonds=new Map(project.bonds.map(b=>[pair(b.a,b.b),b.order]));
    for(const bond of snapshot.bonds)fail(expectedBonds.get(pair(bond.a,bond.b))===bond.order,`This operation changed a bond/order (${expectedBonds.get(pair(bond.a,bond.b))} → ${bond.order}, ${bond.type}); it is not a geometry-only edit.`);
  }
  /** Only pos3 changes. Caller must revision-gate and commit one host history step. */
  captureGeometry(){
    const snapshot=this.read();this.assertSameChemistry(snapshot);
    const next=clone(this.baseProject),coordinates=new Map(snapshot.atoms.map(a=>[a.id,a.pos3]));
    for(const atom of next.atoms)atom.pos3=coordinates.get(atom.id);
    next.has3d=true;
    return{project:next,baseRevision:this.baseRevision,mapping:snapshot.mapping,
      warnings:['Host must validate stereo against the canonical graph and reject stale revisions before applying this coordinate snapshot.']};
  }
  /** Explicit diff only: host decides which intended and incidental edits to accept. */
  captureEditProposal(){
    fail(this.baseProject,'Load a canonical project before capturing changes.');
    const snapshot=this.read(),beforeAtoms=new Map(this.baseProject.atoms.map(a=>[a.id,a]));
    const afterAtoms=new Map(snapshot.atoms.map(a=>[a.id,a]));
    const beforeBonds=new Map(this.baseProject.bonds.map(b=>[pair(b.a,b.b),b]));
    const afterBonds=new Map(snapshot.bonds.map(b=>[pair(b.a,b.b),b]));
    return{baseRevision:this.baseRevision,snapshot,
      addedAtoms:snapshot.atoms.filter(a=>!beforeAtoms.has(a.id)),
      removedAtoms:this.baseProject.atoms.filter(a=>!afterAtoms.has(a.id)),
      changedAtoms:snapshot.atoms.filter(a=>{const b=beforeAtoms.get(a.id);return b&&(a.element!==b.element||a.charge!==(b.charge??0));}),
      movedAtoms:snapshot.atoms.filter(a=>{const b=beforeAtoms.get(a.id);return b?.pos3&&a.pos3.some((v,i)=>Math.abs(v-b.pos3[i])>1e-8);}),
      addedBonds:snapshot.bonds.filter(b=>!beforeBonds.has(pair(b.a,b.b))),
      removedBonds:this.baseProject.bonds.filter(b=>!afterBonds.has(pair(b.a,b.b))),
      changedBonds:snapshot.bonds.filter(b=>{const old=beforeBonds.get(pair(b.a,b.b));return old&&old.order!==b.order;}),
      warnings:['ModelKit can remove attached H atoms and adjust bond lengths even with addHydrogens/autobond disabled. Review this complete diff.',
        'This proposal does not infer or overwrite canonical stereo, isotope, 2D layout, fragment metadata or history. Host integration is required.']};
  }
  /** Accept only the user's declared operation, never incidental ModelKit edits.
   * Host commits this snapshot with WorkshopModel.commit so its stereo audit runs.
   * Each newAtomId and endpoint comes from the explicit native command's result.
   */
  captureIntent(intent){
    fail(this.baseProject,'Load a project before editing.');
    const snapshot=this.read(),next=clone(this.baseProject),oldAtoms=new Map(next.atoms.map(a=>[a.id,a])),native=new Map(snapshot.atoms.map(a=>[a.id,a]));
    const pick=id=>{const a=oldAtoms.get(id);fail(a,'The edit target no longer exists.');return a;};
    const requireNative=id=>{const a=native.get(id);fail(a,'The native edit lost its intended atom.');return a;};
    const newAtom=id=>{fail(!oldAtoms.has(id),'New atom identity already exists.');const a=requireNative(id);const result={id:a.id,element:intent.element??a.element,charge:0,pos2:a.pos3.slice(0,2),pos3:[...a.pos3]};next.atoms.push(result);oldAtoms.set(id,result);return result;};
    let label='三维编辑';
    if(intent.type==='move'||intent.type==='rotate'){
      this.assertSameChemistry(snapshot);
      let moved=false;
      for(const id of intent.atomIds??[]){const old=pick(id),a=requireNative(id),baseline=intent.nativeBefore?.[id]??this.nativeBaseline?.get(id)??old.pos3,delta=a.pos3.map((v,i)=>v-baseline[i]);if(Math.hypot(...delta)>1e-8){old.pos3=old.pos3.map((v,i)=>v+delta[i]);moved=true;}}
      if(moved)next.has3d=true;label=intent.type==='rotate'?'旋转三维原子/片段':'移动三维原子/片段';
    }else if(intent.type==='atom'){
      const old=pick(intent.atomId),actual=requireNative(intent.atomId);
      fail(actual.element===intent.element,'ModelKit did not perform the requested element change.');
      old.element=intent.element;label='修改三维原子元素';
    }else if(intent.type==='charge'){
      const old=pick(intent.atomId),actual=requireNative(intent.atomId);
      fail(actual.charge===intent.charge,'ModelKit did not perform the requested charge change.');
      old.charge=intent.charge;label='修改形式电荷';
    }else if(intent.type==='add'){
      newAtom(intent.newAtomId);next.has3d=true;label='添加三维原子';
    }else if(intent.type==='bond'||intent.type==='draw-bond'){
      for(const id of intent.newAtomIds??[])newAtom(id);
      const a=pick(intent.a),b=pick(intent.b);fail(a.id!==b.id,'A bond requires two distinct atoms.');
      const actual=snapshot.bonds.find(x=>pair(x.a,x.b)===pair(a.id,b.id));
      fail(actual&&actual.order===intent.order,'ModelKit did not perform the requested bond edit.');
      const old=next.bonds.find(x=>pair(x.a,x.b)===pair(a.id,b.id));
      if(old)old.order=intent.order;else next.bonds.push({id:actual.id,a:a.id,b:b.id,order:intent.order});
      if(intent.newAtomIds?.length)next.has3d=true;label=old?'修改三维化学键':'添加三维化学键';
    }else if(intent.type==='delete-bond'){
      const b=next.bonds.find(x=>x.id===intent.bondId);fail(b,'Missing bond to delete.');
      fail(!snapshot.bonds.some(x=>pair(x.a,x.b)===pair(b.a,b.b)),'ModelKit did not delete the requested bond.');
      next.bonds=next.bonds.filter(x=>x.id!==b.id);label='删除三维化学键';
    }else if(intent.type==='delete-atom'){
      pick(intent.atomId);fail(!native.has(intent.atomId),'ModelKit did not delete the requested atom.');
      const removed=intent.atomId;next.atoms=next.atoms.filter(a=>a.id!==removed);next.bonds=next.bonds.filter(b=>b.a!==removed&&b.b!==removed);
      next.fragments=next.fragments.map(f=>({...f,atomIds:f.atomIds.filter(id=>id!==removed),ports:f.ports.filter(p=>p.atomId!==removed)})).filter(f=>f.atomIds.length);
      label='删除三维原子';
    }else fail(false,'Unsupported native edit intent.');
    return{project:next,label,baseRevision:this.baseRevision,mapping:snapshot.mapping};
  }
  /** This export is a diagnostic tool output, not the native project format. */
  exportMol(){return this.value('write("MOL")');}
}
