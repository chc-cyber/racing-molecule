import {JSmolAdapter} from './jsmol-adapter.js';
import {clone,selectedIds,reconcileGeneratedHydrogens,reconcileC60Addition,validateBondEndpoints} from './model.js';
import {C60_ADDITION_NOTICE} from './editor-notices.js';

let runtimePromise,serial=0;
const HostError=globalThis.Clazz?._Error??globalThis.Error;
const ensure=(value,message)=>{if(!value)throw new HostError(message);};
const suppress=e=>{e.preventDefault();e.stopImmediatePropagation();};
function loadRuntime(){
  if(globalThis.Jmol?.getApplet)return Promise.resolve(globalThis.Jmol);
  if(!runtimePromise)runtimePromise=new Promise((resolve,reject)=>{
    const script=document.createElement('script');script.src=new URL('./vendor/jsmol/JSmol.min.js',import.meta.url).href;
    script.onload=()=>resolve(globalThis.Jmol);script.onerror=()=>reject(new HostError('本地 JSmol 未加载，请检查 vendor/jsmol。'));
    document.head.append(script);
  });
  return runtimePromise;
}

/** JSmol owns molecular rendering and live transforms. WorkshopModel owns design.
 * Native chemistry commands execute first; only the declared intent is merged.
 * No renderer event infers connectivity, calls /api/load or runs an optimizer.
 */
export class JSmolHost {
  constructor(container,{model,onSelection,onInteraction,onError,onReady,onBlank,onHover,onAction,onMode,onNotice,unifiedBonds=false,nativeMenus=true}={}){
    ensure(container&&model,'JSmolHost requires a container and WorkshopModel.');
    this.container=container;this.host=container;this.model=model;this.callbacks={onSelection,onInteraction,onError,onReady,onBlank,onHover,onAction,onMode,onNotice};
    this.unifiedBonds=unifiedBonds;this.nativeMenus=nativeMenus;
    this.ready=false;this.isApplying=false;this.suspended=false;this.gesture=false;this.mode='select';
    this.options={element:'C',order:1,rigid:false,native:true};this.selected=new Set();this.previewProject=null;this.visibleProject=null;
    this.queue=Promise.resolve();this.pushSerial=0;this.boundRevision=-1;this.currentSnapshot=null;this.pending=null;this.drag=null;this.bondStart=null;this.numbers=false;
    this.id=`workshopJmol${++serial}`;this.initialized=this.initialize();
  }
  get width(){return this.container.getBoundingClientRect().width;}
  get height(){return this.container.getBoundingClientRect().height;}
  get viewer(){return this.applet?._applet?.viewer;}
  async initialize(){
    this.Jmol=await loadRuntime();this.Jmol.setDocument(0);
    let ready;const loaded=new Promise(resolve=>{ready=resolve;});
    // Bootstrap at a concrete size even when the parent initially shows 2D.
    // Percentage dimensions inside display:none otherwise initialize to zero.
    const info={width:Math.round(this.width)||900,height:Math.round(this.height)||620,color:'white',use:'HTML5',j2sPath:new URL('./vendor/jsmol/j2s',import.meta.url).href,
      serverURL:'',allowJavaScript:false,disableJ2SLoadMonitor:true,disableInitialConsole:true,
      readyFunction:()=>ready(),script:'set antialiasDisplay false; set autobond false; set disablePopupMenu true; background white;'};
    this.applet=this.Jmol.getApplet(this.id,info);globalThis[this.id]=this.applet;
    this.container.innerHTML=this.Jmol.getAppletHtml(this.applet);
    await loaded;this.adapter=new JSmolAdapter(this.applet,{Jmol:this.Jmol});
    this.viewer._mcarModelKitBridge={action:name=>this.nativeAction(name),mode:picking=>this.nativeMode(picking),canUndo:()=>!!this.model._undo.length,canRedo:()=>!!this.model._redo.length};
    this.applet._mcarModelKitBridge=this.viewer._mcarModelKitBridge;
    // The main workshop uses its sidebars. The isolated ModelKit probe can
    // still opt into native menus, routed through the same owner/history.
    const popup=this.viewer.popupMenu,menus=this.nativeMenus;
    this.viewer.popupMenu=function(x,y,type){if(menus)return popup.call(this,x,y,type==='j'?'m':type);};
    this.guardNativeMouseEntry();this.bindEvents();this.ready=true;
    this.loadNow(this.model.project,{fit:true});
    this.callbacks.onReady?.(this);return this;
  }
  guardNativeMouseEntry(){
    // JSmol 16.3.31 sends a synthetic 502 (release) on canvas mouseenter.
    // Its Mouse.processEvent may turn that into a click when the pointer
    // returns to the last pressed coordinates, even without a new press.
    // Keep native release recovery, but only a real press may produce a click.
    // This guard belongs to this applet instance; upstream files stay intact.
    const mouse=this.viewer.mouse,processEvent=mouse.processEvent;
    mouse.processEvent=function(id,x,y,modifiers,time){
      if(id===502&&!this.isMouseDown){
        this.released(time,x,y,globalThis.J.awtjs2d.Mouse.applyLeftMouse(modifiers));
        return true;
      }
      return processEvent.call(this,id,x,y,modifiers,time);
    };
  }
  bindEvents(){
    // Keep keyboard shortcuts in the workshop. The upstream canvas key handler
    // otherwise consumes events before global save/undo/tool shortcuts can run.
    this.Jmol.$(this.canvas()).off('keydown keyup');
    this.container.addEventListener('mousedown',e=>this.pointerDown(e),true);
    this.container.ownerDocument.addEventListener('mouseup',e=>this.pointerUp(e),true);
    this.container.addEventListener('mousemove',e=>{if(!this.ready||this.suspended)return;if(this.drag)this.drag.moved ||= Math.hypot(e.clientX-this.drag.x,e.clientY-this.drag.y)>3;this.callbacks.onHover?.(this.hitTest(e),e);},true);
    this.container.addEventListener('contextmenu',suppress,true);
    this.container.addEventListener('dblclick',e=>this.selectConnected(e),true);
    this.container.addEventListener('keydown',e=>this.keydown(e),true);
    this.container.tabIndex=0;
    this.resizeObserver=new ResizeObserver(()=>{if(this.ready&&this.width&&this.height)this.Jmol.resizeApplet(this.applet,['100%','100%']);});
    this.resizeObserver.observe(this.container);
  }
  error(error){this.callbacks.onError?.(typeof error?.message==='string'?error:new HostError(String(error)));}
  script(command){return this.adapter.script(command);}
  setMode(mode,options={}){
    // A control change owns a transaction boundary even when it only changes
    // the element or bond order. Accept completed work first; an unreleased
    // native gesture is restored from the canonical design before rebinding.
    this.flush();
    const changed=mode!==this.mode||options.nativePicking!==this.options.nativePicking||options.element&&options.element!==this.options.element||options.order!==undefined&&options.order!==this.options.order;
    if(changed&&this.gesture)this.cancelGesture();
    this.mode=mode;
    // Selection and geometry shortcuts use JSmol's own mouse bindings. Exact
    // construction shortcuts retain the explicit owner-controlled fallback.
    const native=options.native??['select','move','molecule','rotateBond','measure'].includes(mode);
    this.options={...this.options,...options,native,nativePicking:options.nativePicking??null};
    if(changed)this.bondStart=null;
    if(this.ready&&!this.suspended&&!this.gesture)this.configureMode();
  }
  configureMode(){
    if(!this.adapter)return;
    this.script(`unbind; modelkit addHydrogens false; modelkit autobond false; modelkit hidden ${this.previewProject?'true':'false'}; set modelKitMode ${this.activeProject().atoms.length&&!this.previewProject?'true':'false'}; set autobond false; set disablePopupMenu ${this.previewProject||!this.nativeMenus?'true':'false'}; ${this.nativeMenus?'':'frank off;'} set antialiasDisplay false; set picking identifyBond; set picking atom; set bondPicking false; set atomPicking true; set dragSelected false;`);
    if(this.options.native&&!this.previewProject){
      const picking=this.options.nativePicking??({move:'dragAtom',molecule:'dragMolecule',rotateBond:'rotateBond',select:'atom',measure:'measure',atom:`assignAtom_${this.options.element}`,bond:`assignBond_${this.options.order}`,delete:'assignAtom_X'}[this.mode]);
      if(picking)this.script(`set picking ${picking};`);
      // Delete exposes both native targets; the owner captures only the atom
      // or bond explicitly clicked, never ModelKit's incidental chemistry.
      if(this.mode==='delete')this.script('set picking assignBond_0;');
      this.nativeGeometryBindings();
      this.viewer.acm.clearMouseInfo();
      this.paintSelection();return;
    }
    // Construction fallbacks also use owner-controlled transactions.
    this.script('unbind "LEFT+drag"; unbind "LEFT+double+click"; unbind "LEFT+double+drag";');
    if(!this.previewProject&&['move','rotate'].includes(this.mode)){const action=this.mode==='rotate'?'_rotateSelected':'_dragSelected';this.script(`set modelKitMode false; set dragSelected true; set allowRotateSelected true; bind "LEFT+drag" "${action}"; bind "LEFT+double+drag" "${action}";`);}
    else this.script('bind "LEFT+drag" "_rotate";');
    this.paintSelection();
  }
  nativeMode(picking){
    if(this.suspended||this.previewProject)return;
    this.flush();if(this.gesture)this.cancelGesture();this.bondStart=null;
    let mode,options={native:true,nativePicking:picking,unified:false};
    if(picking==='dragAtom')mode='move';
    else if(picking==='dragMolecule')mode='molecule';
    else if(picking==='rotateBond')mode='rotateBond';
    else if(picking==='atom')mode='select';
    else if(picking==='measure')mode='measure';
    else if(picking.startsWith('assignAtom_')){
      const value=picking.slice(11);
      if(value==='X')mode='delete';
      else if(/^(?:pl|mi)$/i.test(value)){mode='charge';options.delta=value.toLowerCase()==='pl'?1:-1;}
      else if(value==='Xx'){mode='bond';options.order=1;options.unified=this.unifiedBonds;}
      else{mode='atom';options.element=value;}
    }else if(picking.startsWith('assignBond_')){mode='bond';const value=picking.slice(11);options.order=/^[0123]$/.test(value)?Number(value):value;options.unified=this.unifiedBonds&&/^[123]$/.test(value);}
    if(!mode)return;
    this.mode=mode;this.options={...this.options,...options};
    if(!this.viewer.g.modelKitMode)this.script('set modelKitMode true;');
    this.script('unbind; modelkit addHydrogens false; modelkit autobond false; set dragSelected false; set picking identifyBond; set bondPicking false;');
    this.nativeGeometryBindings();this.viewer.acm.clearMouseInfo();this.paintSelection();
    this.callbacks.onMode?.(mode,{...this.options});
  }
  nativeGeometryBindings(){
    // Upstream interprets rapid consecutive drags as double drags. Keep its
    // geometry action in that case instead of falling through to camera motion.
    if(['move','molecule','rotateBond'].includes(this.mode))this.script('bind "LEFT+double+drag" "_dragAtom"; bind "ALT+LEFT+double+drag" "_rotateSelected";');
  }
  nativeAction(name){
    if(this.suspended||this.previewProject)return;
    this.flush();this.cancelGesture(false);this.callbacks.onInteraction?.();
    try{
      if(this.callbacks.onAction)Promise.resolve(this.callbacks.onAction(name)).catch(error=>this.error(error));
      else if(name==='undo')this.model.undo();
      else if(name==='redo')this.model.redo();
    }catch(error){this.error(error);}
  }
  groupIds(seed){
    const group=this.model.project.metadata?.ketcherGroups?.find(g=>(g.kind==='c60'||g.type==='SUP'&&!g.expanded)&&g.atomIds?.includes(seed));
    if(!group)return null;
    const p=this.model.project,ids=new Set(group.atomIds.filter(id=>p.atoms.some(a=>a.id===id)));
    // Terminal explicit H atoms follow the cage but do not become template
    // members; its indexToAtomId remains the original 60 carbon identities.
    for(const a of p.atoms.filter(a=>a.element==='H')){
      const bonds=p.bonds.filter(b=>b.a===a.id||b.b===a.id);
      if(bonds.length===1&&ids.has(bonds[0].a===a.id?bonds[0].b:bonds[0].a))ids.add(a.id);
    }
    return [...ids];
  }
  connectedIds(seed){
    const ids=new Set([seed]),todo=[seed];
    for(let i=0;i<todo.length;i++)for(const bond of this.model.project.bonds){const next=bond.a===todo[i]?bond.b:bond.b===todo[i]?bond.a:null;if(next&&!ids.has(next)){ids.add(next);todo.push(next);}}
    return [...ids];
  }
  loadNow(project,{fit=false,preview=false}={}){
    const hadAtoms=!!this.visibleProject?.atoms.length;
    this.adapter.loadCanonical(project,{revision:this.model.revision,preserveView:hadAtoms&&!fit});
    this.currentSnapshot=this.adapter.read();this.visibleProject=clone(project);this.boundRevision=this.model.revision;
    if(this.bondStart&&!project.atoms.some(a=>a.id===this.bondStart))this.bondStart=null;
    this.script('spacefill 22%; wireframe 0.12; color atoms cpk; background white; set selectionHalos true; color selectionHalos [xc44b3a];');
    if(fit||!hadAtoms)this.resetView();
    this.configureMode();this.setNumbers(this.numbers);
    if(!preview)this.previewProject=null;
  }
  push(project,{fit=false}={}){
    const snapshot=clone(project),n=++this.pushSerial;
    this.suspended=true;this.cancelGesture(false);this.previewProject=null;
    const task=this.queue.catch(()=>{}).then(async()=>{
      await this.initialized;if(n!==this.pushSerial)return;
      try{this.loadNow(snapshot,{fit});}
      finally{if(n===this.pushSerial){this.suspended=false;this.configureMode();}}
    });this.queue=task;return task;
  }
  async settled(){
    await this.initialized;
    this.flush();let task;do{task=this.queue;await task;}while(task!==this.queue);
  }
  setSelection(ids){this.selected=new Set(ids);if(this.ready&&!this.suspended)this.paintSelection();}
  paintSelection(){
    if(!this.currentSnapshot)return;const indices=this.currentSnapshot.atoms.filter(a=>this.selected.has(a.id)).map(a=>a.jmolIndex);
    this.script('select '+(indices.length?`{${indices.map(i=>`atomindex=${i}`).join(' or ')}}`:'none')+'; set selectionHalos true;');
  }
  emitSelection(ids,bondId=null){this.selected=new Set(ids);this.paintSelection();this.callbacks.onSelection?.([...this.selected],bondId);}
  selectConnected(e){
    if(!this.ready||this.suspended||this.previewProject||!['select','move','rotate','molecule'].includes(this.mode))return;
    suppress(e);
    this.flush();const hit=this.hitTest(e),seed=hit.atomId??hit.a;if(!seed)return;
    const connected=new Set(this.connectedIds(seed));
    this.emitSelection(e.shiftKey?[...new Set([...this.selected,...connected])]:[...connected]);
  }
  setNumbers(value){this.numbers=!!value;if(this.ready)this.script(`select {*}; label ${this.numbers?'"%a"':'off'};`);if(this.ready&&!this.suspended)this.paintSelection();}
  resetView(){if(!this.ready)return;this.script('set zoomLarge false; reset; center {*}; zoom 85;');}
  reset(){this.resetView();}
  async onVisible(){await this.initialized;await this.queue;if(this.width&&this.height)this.Jmol.resizeApplet(this.applet,['100%','100%']);}
  showPreview(project){
    if(!this.ready)return;this.cancelGesture(false);this.previewProject=clone(project);this.loadNow(project,{preview:true});this.previewProject=clone(project);this.configureMode();
  }
  clearPreview(){if(!this.previewProject)return;this.previewProject=null;if(this.ready)this.loadNow(this.model.project);}
  cancelGesture(reload=true){this.drag=null;this.pending=null;this.gesture=false;if(reload&&this.ready&&!this.suspended)this.push(this.model.project).catch(e=>this.error(e));}
  activeProject(){return this.previewProject??this.model.project;}
  screen3(position){
    if(!this.viewer)return null;const p=globalThis.JU.P3.new3(...position),q=globalThis.JU.P3.new3(0,0,0);this.viewer.tm.transformPt3f(p,q);
    const canvas=this.canvas();if(!canvas)return null;const factor=this.viewer.gdata.isAntialiased()?2:1;
    return[q.x/factor,q.y/factor,q.z];
  }
  canvas(){return this.container.querySelector('canvas');}
  point3(event,reference){
    if(!this.viewer)return[0,0,0];const rect=this.canvas()?.getBoundingClientRect()??this.container.getBoundingClientRect();
    const center=reference??this.viewer.tm.fixedRotationCenter;const xyz=Array.isArray(center)?center:[center?.x??0,center?.y??0,center?.z??0];
    const depth=this.screen3(xyz)?.[2]??1000,factor=this.viewer.gdata.isAntialiased()?2:1;
    const p=globalThis.JU.P3.new3((event.clientX-rect.left)*factor,(event.clientY-rect.top)*factor,depth),q=globalThis.JU.P3.new3(0,0,0);this.viewer.tm.unTransformPoint(p,q);return[q.x,q.y,q.z];
  }
  hitTest(event){
    if(!this.viewer||!this.currentSnapshot)return{};const rect=this.canvas()?.getBoundingClientRect();if(!rect)return{};
    const x=event.clientX-rect.left,y=event.clientY-rect.top;
    if(x<0||y<0||x>rect.width||y>rect.height)return{};
    const i=this.viewer.findNearestAtomIndex(x,y),atom=this.currentSnapshot.atoms.find(a=>a.jmolIndex===i);
    if(atom)return{atomId:atom.id,jmolIndex:i};
    // Picking a displayed bond uses its projection, never inter-atom distance or
    // connectivity reconstruction. Bond topology always comes from the project.
    let best=null,distance=11;
    const atoms=new Map((this.visibleProject?.atoms??[]).map(a=>[a.id,a]));
    for(const b of this.visibleProject?.bonds??[]){
      const a=atoms.get(b.a),c=atoms.get(b.b);if(!a||!c)continue;
      const q=this.screen3(a.pos3),r=this.screen3(c.pos3);if(!q||!r)continue;
      const dx=r[0]-q[0],dy=r[1]-q[1],length=dx*dx+dy*dy;if(length<1)continue;
      // Exclude the end caps, where atom picking has priority; allow the whole
      // visible interior instead of requiring an eleven-pixel midpoint target.
      const t=Math.max(.12,Math.min(.88,((x-q[0])*dx+(y-q[1])*dy)/length));
      const d=Math.hypot(x-q[0]-t*dx,y-q[1]-t*dy);if(d<distance){best=b;distance=d;}
    }
    return best?{bondId:best.id,a:best.a,b:best.b}:{};
  }
  pointerDown(e){
    if(e.button!==0||!e.target.closest?.('canvas')||!this.ready)return;
    if(this.suspended||this.previewProject){suppress(e);return;}
    const rect=this.canvas().getBoundingClientRect(),x=e.clientX-rect.left,y=e.clientY-rect.top;
    if(this.viewer.frankClickedModelKit(x,y)||this.viewer.frankClicked(x,y))return;
    const hit=this.hitTest(e);this.container.focus({preventScroll:true});
    this.callbacks.onInteraction?.();
    if(this.mode==='bond'&&this.options.unified){
      // One UI tool delegates by its explicit target. JSmol still draws and
      // edits the live model; no distance-based bond discovery is introduced.
      this.script(`set picking identifyBond; set bondPicking false; ${hit.atomId?'set picking assignAtom_'+this.options.element+';':''}set picking ${hit.bondId?'assignBond_'+this.options.order:hit.atomId?'assignAtom_Xx':'atom'};`);
      this.viewer.acm.clearMouseInfo();
    }
    if(['move','rotate'].includes(this.mode)){
      if(hit.atomId&&!this.selected.has(hit.atomId))this.emitSelection(this.groupIds(hit.atomId)??[hit.atomId]);
      if(hit.bondId&&!this.selected.has(hit.a)&&!this.selected.has(hit.b))this.emitSelection([hit.a,hit.b],hit.bondId);
      if(!hit.atomId&&!hit.bondId){if(!this.options.native)suppress(e);return;}
      const ids=selectedIds(this.model.project,this.selected,!!this.options.rigid);this.selected=ids;this.paintSelection();
    }
    const native=this.options.native;
    if(native&&this.mode==='molecule'&&hit.atomId)this.emitSelection(this.groupIds(hit.atomId)??this.connectedIds(hit.atomId));
    const group=native&&hit.atomId&&['move','molecule'].includes(this.mode)?(this.groupIds(hit.atomId)??(this.mode==='move'&&this.selected.size>1?[...this.selected]:null)):null;
    if(group){
      this.emitSelection(group);
      // The group is explicitly stored by Ketcher, so JSmol's selected-atom
      // transform moves the cage without also moving the attached vehicle.
      this.script('set picking dragSelected; set dragSelected true; set allowRotateSelected true; unbind "LEFT+drag"; unbind "ALT+LEFT+drag"; bind "LEFT+drag" "_dragAtom"; bind "ALT+LEFT+drag" "_rotateSelected";');
    }
    this.gesture=true;this.drag={mode:this.mode,options:{...this.options},hit,x:e.clientX,y:e.clientY,point:this.point3(e),baseRevision:this.model.revision,atomIds:[...this.selected],nativeBefore:Object.fromEntries(this.adapter.read().atoms.map(a=>[a.id,a.pos3])),shift:e.shiftKey,alt:e.altKey,native,group:!!group,moved:false};
    if(!native&&['atom','bond','delete','charge','measure','fragment'].includes(this.mode))suppress(e);
  }
  pointerUp(e){
    if(!this.drag||e.button!==0)return;const d=this.drag;this.drag=null;this.gesture=false;
    d.endHit=this.hitTest(e);d.endPoint=this.point3(e);d.event={shiftKey:e.shiftKey,clientX:e.clientX,clientY:e.clientY};
    d.moved ||= Math.hypot(e.clientX-d.x,e.clientY-d.y)>3;
    this.pending=d;
    if(!d.native&&['atom','bond','delete','charge','measure','fragment'].includes(d.mode))suppress(e);
    // Native mouseup finishes its transform before we read the edited coordinates.
    setTimeout(()=>{try{this.flush();}catch(error){this.error(error);}},0);
  }
  nativeIndex(id){const a=this.adapter.read().atoms.find(a=>a.id===id);ensure(a,'所选原子已变化，请重新选择。');return a.jmolIndex;}
  nativeAdd(point,element){
    ensure(/^[A-Z][a-z]?$/.test(element),'元素符号无效。');const before=new Set(this.adapter.read().atoms.map(a=>a.id));
    this.script(`modelkit add "${element}" {${point.join(' ')}};`);
    const added=this.adapter.read().atoms.filter(a=>!before.has(a.id));ensure(added.length,'JSmol 未能添加原子。');
    // On a completely empty native scene ModelKit may create its default atom
    // before the requested one. Only the last explicit addition is accepted.
    return added.at(-1).id;
  }
  nativeBond(a,b,order){
    ensure([1,2,3,1.5].includes(order),'键级无效。');const ia=this.nativeIndex(a),ib=this.nativeIndex(b);
    const existing=this.adapter.read().bonds.some(x=>x.a===a&&x.b===b||x.a===b&&x.b===a);
    if(order===1.5)this.script(`connect {atomindex=${ia}} {atomindex=${ib}} aromatic ${existing?'modify':'create'};`);
    else if(existing)this.script(`modelkit assign bond {atomindex=${ia} or atomindex=${ib}} "${order}";`);
    else this.script(`modelkit connect {atomindex=${ia}} {atomindex=${ib}} "${order}";`);
  }
  commitIntent(intent){
    const result=this.adapter.captureIntent(intent);ensure(result.baseRevision===this.model.revision,'操作期间设计已变化，已取消旧操作。');
    let addition=false;
    if(intent.validateConnection){
      reconcileGeneratedHydrogens(this.model.project,result.project);
      addition=reconcileC60Addition(this.model.project,result.project);
      validateBondEndpoints(result.project,intent.a,intent.b);
    }
    const oldAtoms=new Map(this.model.project.atoms.map(a=>[a.id,a]));
    const geometryConfirmedIds=result.project.atoms.filter(a=>!oldAtoms.has(a.id)).map(a=>a.id);
    if(['move','rotate'].includes(intent.type))for(const atom of result.project.atoms){
      const previous=oldAtoms.get(atom.id);
      if(intent.atomIds.includes(atom.id)&&previous&&atom.pos3.some((value,i)=>value!==previous.pos3[i]))geometryConfirmedIds.push(atom.id);
    }
    this.isApplying=true;
    try{this.model.commit(result.label,p=>{for(const key of Object.keys(p))delete p[key];Object.assign(p,result.project);},{geometryConfirmedIds});}
    finally{this.isApplying=false;}
    if(['move','rotate'].includes(intent.type)){
      // A geometry gesture cannot change chemistry. Retain the live native
      // model/camera instead of reloading after every drag, and update the exact
      // canonical baseline for the next displacement transaction.
      this.adapter.acceptCanonical(this.model.project,{revision:this.model.revision});
      this.currentSnapshot=this.adapter.read();this.visibleProject=clone(this.model.project);this.boundRevision=this.model.revision;this.paintSelection();
    }else this.push(this.model.project).catch(e=>this.error(e));
    if(addition)this.callbacks.onNotice?.(C60_ADDITION_NOTICE);
    return result;
  }
  flush(){try{return this.flushPending();}catch(error){if(this.ready)this.push(this.model.project).catch(e=>this.error(e));throw error;}}
  flushPending(){
    if(!this.ready||this.suspended||this.gesture||!this.pending)return false;
    const d=this.pending;this.pending=null;
    if(d.baseRevision!==this.model.revision){this.push(this.model.project).catch(e=>this.error(e));this.error(new HostError('操作期间设计已变化，旧三维操作已取消。'));return false;}
    if(d.native)return this.flushNative(d);
    const hit=d.hit,mode=d.mode,element=d.options.element,order=Number(d.options.order);
    if(mode==='move'||mode==='rotate'){
      if(d.moved)this.commitIntent({type:mode,atomIds:d.atomIds,nativeBefore:d.nativeBefore});return d.moved;
    }
    if(mode==='select'||mode==='measure'||mode==='fragment'){
      if(d.moved&&mode==='select')return false;
      if(mode==='fragment'){this.callbacks.onBlank?.(d.endPoint,d.event,hit);return false;}
      if(hit.bondId)this.emitSelection([hit.a,hit.b],hit.bondId);
      else if(hit.atomId){let ids=mode==='measure'?[hit.atomId]:(d.shift?[...this.selected]:[]);if(d.shift&&ids.includes(hit.atomId))ids=ids.filter(x=>x!==hit.atomId);else if(!ids.includes(hit.atomId))ids.push(hit.atomId);this.emitSelection(ids);}
      else if(!d.shift)this.emitSelection([]);
      return false;
    }
    if(mode==='atom'){
      if(hit.atomId){const i=this.nativeIndex(hit.atomId);this.script(`modelkit assign atom {atomindex=${i}} "${element}";`);this.commitIntent({type:'atom',atomId:hit.atomId,element});this.emitSelection([hit.atomId]);}
      else{const id=this.nativeAdd(d.endPoint,element);this.commitIntent({type:'add',newAtomId:id,element});this.emitSelection([id]);}return true;
    }
    if(mode==='charge'&&hit.atomId){const atom=this.model.project.atoms.find(a=>a.id===hit.atomId),charge=d.options.charge??atom.charge+1,delta=charge-atom.charge;ensure(Number.isInteger(charge),'电荷必须为整数。');for(let n=0;n<Math.abs(delta);n++)this.script(`modelkit assign atom {atomindex=${this.nativeIndex(hit.atomId)}} "${delta>0?'Pl':'Mi'}";`);this.commitIntent({type:'charge',atomId:hit.atomId,charge});return true;}
    if(mode==='delete'){
      if(hit.atomId){this.script(`modelkit assign atom {atomindex=${this.nativeIndex(hit.atomId)}} "X";`);this.commitIntent({type:'delete-atom',atomId:hit.atomId});}
      else if(hit.bondId){this.script(`modelkit assign bond {atomindex=${this.nativeIndex(hit.a)} or atomindex=${this.nativeIndex(hit.b)}} "0";`);this.commitIntent({type:'delete-bond',bondId:hit.bondId});}
      return !!(hit.atomId||hit.bondId);
    }
    if(mode==='bond'){
      if(hit.bondId&&!d.moved){this.nativeBond(hit.a,hit.b,order);this.commitIntent({type:'bond',a:hit.a,b:hit.b,order});this.emitSelection([hit.a,hit.b],hit.bondId);this.bondStart=null;return true;}
      if(!d.moved&&hit.atomId&&!this.bondStart){this.bondStart=hit.atomId;this.emitSelection([hit.atomId]);return false;}
      if(!d.moved&&!hit.atomId&&!this.bondStart){const id=this.nativeAdd(d.endPoint,element);this.commitIntent({type:'add',newAtomId:id,element});this.bondStart=id;this.emitSelection([id]);return true;}
      const newAtomIds=[];let a=d.moved?hit.atomId:this.bondStart,b=d.endHit.atomId;
      if(a&&b&&a===b){this.bondStart=null;return false;}
      if(!a){a=this.nativeAdd(d.point,element);newAtomIds.push(a);}
      if(!b){b=this.nativeAdd(d.endPoint,element);newAtomIds.push(b);}
      this.nativeBond(a,b,order);this.commitIntent({type:'draw-bond',a,b,order,element,newAtomIds});this.bondStart=null;this.emitSelection([a,b]);return true;
    }
    return false;
  }
  flushNative(d){
    const {hit,mode}=d,snapshot=this.adapter.read(),old=this.model.project;
    if(mode==='bond'&&d.options.unified)return this.flushUnifiedBond(d,snapshot);
    if(['move','molecule','rotateBond'].includes(mode)){
      this.adapter.assertSameChemistry(snapshot,old);
      if(d.moved){
        const ids=d.group?d.atomIds:old.atoms.map(a=>a.id);
        const changed=snapshot.atoms.some(a=>ids.includes(a.id)&&a.pos3.some((v,i)=>Math.abs(v-d.nativeBefore[a.id][i])>1e-8));
        if(changed)this.commitIntent({type:d.alt||mode==='rotateBond'?'rotate':'move',atomIds:ids,nativeBefore:d.nativeBefore});
        if(d.group)this.configureMode();return changed;
      }
      if(hit.atomId)this.emitSelection(this.groupIds(hit.atomId)??(mode==='molecule'?this.connectedIds(hit.atomId):[hit.atomId]));
      return false;
    }
    if(mode==='select'){
      this.adapter.assertSameChemistry(snapshot,old);this.restoreNativeSideEffects(snapshot,d);
      if(hit.bondId&&!d.moved){this.emitSelection([hit.a,hit.b],hit.bondId);return false;}
      const selected=new Set(this.adapter.value('{selected}.atomindex.all'));
      this.emitSelection(snapshot.atoms.filter(a=>selected.has(a.jmolIndex)).map(a=>a.id));return false;
    }
    if(mode==='measure'){
      this.adapter.assertSameChemistry(snapshot,old);this.restoreNativeSideEffects(snapshot,d);
      if(hit.atomId){const ids=[...this.selected];if(!ids.includes(hit.atomId))ids.push(hit.atomId);this.emitSelection(ids.slice(-4));}
      return false;
    }
    if(mode==='charge'&&hit.atomId){
      const atom=old.atoms.find(a=>a.id===hit.atomId),actual=snapshot.atoms.find(a=>a.id===hit.atomId);
      if(actual&&actual.charge!==atom.charge)this.commitIntent({type:'charge',atomId:atom.id,charge:atom.charge+d.options.delta});
      else this.restoreNativeSideEffects(snapshot,d);
      return true;
    }
    if(mode==='delete'&&hit.atomId){
      if(!snapshot.atoms.some(a=>a.id===hit.atomId))this.commitIntent({type:'delete-atom',atomId:hit.atomId});return true;
    }
    if(mode==='delete'&&hit.bondId){
      if(!snapshot.bonds.some(b=>b.id===hit.bondId))this.commitIntent({type:'delete-bond',bondId:hit.bondId});
      else this.restoreNativeSideEffects(snapshot,d);
      return true;
    }
    if(mode==='bond'&&hit.bondId){
      const bond=old.bonds.find(b=>b.id===hit.bondId),actual=snapshot.bonds.find(b=>b.id===hit.bondId);
      if(!actual)this.commitIntent({type:'delete-bond',bondId:hit.bondId});
      else if(actual.order!==bond.order)this.commitIntent({type:'bond',a:bond.a,b:bond.b,order:actual.order});
      else this.restoreNativeSideEffects(snapshot,d);
      return true;
    }
    if(mode==='atom'&&hit.atomId&&!d.moved){
      const atom=old.atoms.find(a=>a.id===hit.atomId),actual=snapshot.atoms.find(a=>a.id===hit.atomId);
      if(actual&&actual.element!==atom.element)this.commitIntent({type:'atom',atomId:hit.atomId,element:d.options.element});
      else this.restoreNativeSideEffects(snapshot,d);
      return true;
    }
    if(mode==='atom'||mode==='bond'){
      const added=snapshot.atoms.filter(a=>!old.atoms.some(b=>b.id===a.id));
      const newBonds=snapshot.bonds.filter(b=>!old.bonds.some(c=>c.id===b.id));
      const target=newBonds.find(b=>!hit.atomId||b.a===hit.atomId||b.b===hit.atomId);
      if(target){
        const ids=added.filter(a=>target.a===a.id||target.b===a.id).map(a=>a.id);
        this.commitIntent({type:'draw-bond',a:target.a,b:target.b,order:target.order,element:d.options.element,newAtomIds:ids});return true;
      }
      if(added.length){this.commitIntent({type:'add',newAtomId:added.at(-1).id,element:d.options.element});return true;}
      // ModelKit has no atom to drag from on an empty canvas. Explicit blank
      // insertion is the only fallback; it neither creates default methane nor
      // guesses a chemical bond from the cursor distance.
      if(mode==='atom'&&!hit.atomId&&!hit.bondId&&!d.moved){const id=this.nativeAdd(d.endPoint,d.options.element);this.commitIntent({type:'add',newAtomId:id,element:d.options.element});this.emitSelection([id]);return true;}
    }
    this.restoreNativeSideEffects(snapshot,d);return false;
  }
  flushUnifiedBond(d,snapshot){
    const {hit}=d,order=Number(d.options.order);
    try{
      if(hit.bondId&&!d.moved){
        this.nativeBond(hit.a,hit.b,order);
        this.commitIntent({type:'bond',a:hit.a,b:hit.b,order});
        this.emitSelection([hit.a,hit.b],hit.bondId);return true;
      }
      if(hit.atomId&&!d.moved){
        this.restoreNativeSideEffects(snapshot,d);this.emitSelection([hit.atomId]);return false;
      }
      if(!hit.atomId||d.endHit.atomId===hit.atomId){this.restoreNativeSideEffects(snapshot,d);return false;}
      const added=snapshot.atoms.filter(a=>!this.model.project.atoms.some(old=>old.id===a.id));
      const nativeEdge=snapshot.bonds.find(b=>(b.a===hit.atomId&&added.some(a=>a.id===b.b))||(b.b===hit.atomId&&added.some(a=>a.id===b.a)));
      let target=d.endHit.atomId??(nativeEdge?(nativeEdge.a===hit.atomId?nativeEdge.b:nativeEdge.a):null);
      if(!target)target=this.nativeAdd(d.endPoint,d.options.element);
      const newAtomIds=this.model.project.atoms.some(a=>a.id===target)?[]:[target];
      this.nativeBond(hit.atomId,target,order);
      this.commitIntent({type:'draw-bond',a:hit.atomId,b:target,order,element:d.options.element,newAtomIds,validateConnection:true});
      this.emitSelection([hit.atomId,target]);return true;
    }finally{if(!this.suspended)this.configureMode();}
  }
  restoreNativeSideEffects(snapshot,d){
    const old=this.model.project,atoms=new Map(old.atoms.map(a=>[a.id,a])),bonds=new Map(old.bonds.map(b=>[b.id,b]));
    const changed=snapshot.atoms.length!==old.atoms.length||snapshot.bonds.length!==old.bonds.length||snapshot.atoms.some(a=>{
      const original=atoms.get(a.id),coords=d.nativeBefore[a.id];
      return !original||a.element!==original.element||a.charge!==original.charge||!coords||a.pos3.some((v,i)=>Math.abs(v-coords[i])>1e-8);
    })||snapshot.bonds.some(b=>{const original=bonds.get(b.id);return !original||b.order!==original.order;});
    if(changed)this.push(old).catch(e=>this.error(e));
  }
  keydown(e){
    if(!this.ready)return;const key=e.key.toLowerCase(),command=e.ctrlKey||e.metaKey;
    if(this.suspended){suppress(e);return;}
    if(command&&['z','y'].includes(key)){suppress(e);const unfinished=this.gesture;this.flush();this.cancelGesture(false);key==='y'||e.shiftKey?this.model.redo():this.model.undo();if(unfinished)this.push(this.model.project).catch(e=>this.error(e));return;}
    if(e.key==='Escape'){suppress(e);this.cancelGesture();this.bondStart=null;return;}
    // Other keys bubble to the workshop; native canvas key handlers are removed.
  }
}
