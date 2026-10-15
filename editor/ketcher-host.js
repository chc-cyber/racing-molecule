// Preserve useful errors when the JSmol Java runtime exports its Error class.
const NativeError=globalThis.Clazz?._Error??globalThis.Error;
import {KetcherAdapter} from './ketcher-adapter.js';
import {clone} from './model.js';

/** Own the integration boundary; the drawing UI is the unmodified Ketcher build.
 * One pointer gesture becomes one canonical transaction. Programmatic loads are
 * explicitly mapped and cannot feed back into the model. */
export class KetcherHost {
  constructor(frame, {model, onSelection, onInteraction, onError, onReady, onOpen, onSave, onClear, onBusy, onNotice} = {}) {
    this.frame=frame;this.model=model;this.adapter=new KetcherAdapter({toolVersion:'3.18.0'});
    this.callbacks={onSelection,onInteraction,onError,onReady,onOpen,onSave,onClear,onBusy,onNotice};
    this.ready=false;this.isApplying=false;this.suspended=false;this.gesture=false;this.dirty=false;
    this.needsFit=false;this.lastVisibleSize=null;
    this.pushSerial=0;this.queue=Promise.resolve();this.lastSelection='';this.timer=null;
    this.initialized=this.initialize();
  }
  async initialize() {
    const deadline=Date.now()+60000;
    while(Date.now()<deadline) {
      try {const k=this.frame.contentWindow?.ketcher;if(k?.editor?.struct){this.ketcher=k;break;}}catch{}
      await new Promise(resolve=>setTimeout(resolve,75));
    }
    if(!this.ketcher)throw NativeError('Ketcher 未完成初始化，请刷新页面并检查本地发行文件。');
    const doc=this.frame.contentDocument,k=this.ketcher;
    this.injectConfiguration(doc);
    this.changeSubscription=k.editor.subscribe('change',()=>this.changed());
    // Ketcher's native info area is tied to its attachment-point debug display.
    // Surface our template's recovery instruction in the workshop status bar.
    this.messageSubscription=k.editor.subscribe('message',message=>{
      if(typeof message?.info==='string'&&message.info.startsWith('请先将 C₆₀ 放在空白处'))this.callbacks.onNotice?.(message.info);
    });
    // The native Subscription chain may consume the callback argument. Read
    // the current selection instead of relying on its dispatched payload.
    this.selectionSubscription=k.editor.subscribe('selectionChange',()=>this.readSelection());
    // These entry points and iframe shortcuts all use the same project history.
    k.editor.undo=()=>this.undo();k.editor.redo=()=>this.redo();
    doc.addEventListener('pointerdown',e=>{
      if(this.suspended)return;
      if(this.isTextInput(e.target))return;
      this.gesture=true;this.callbacks.onInteraction?.();
    },true);
    doc.addEventListener('pointerup',()=>{this.gesture=false;this.schedule();},true);
    doc.addEventListener('pointercancel',()=>{this.gesture=false;this.schedule();},true);
    doc.addEventListener('keydown',e=>this.keydown(e),true);
    doc.addEventListener('click',e=>{
      if(!e.target.closest?.('[data-testid="clear-canvas"]'))return;
      e.preventDefault();e.stopImmediatePropagation();this.clear();
    },true);
    doc.addEventListener('keyup',()=>this.schedule(),true);
    this.frame.contentWindow.addEventListener('blur',()=>{this.gesture=false;this.schedule();});
    let revision;
    do {
      revision=this.model.revision;
      await this.push(clone(this.model.project),{fit:true});
    } while(revision!==this.model.revision);
    this.ready=true;this.setSelection(this.model.selection);this.readSelection();this.callbacks.onReady?.(k);
    this.resizeObserver=new ResizeObserver(()=>{
      if(!this.ready||this.gesture)return;
      cancelAnimationFrame(this.resizeTick);
      this.resizeTick=requestAnimationFrame(()=>this.onVisible().catch(error=>this.callbacks.onError?.(error)));
    });
    this.resizeObserver.observe(this.frame);
    return k;
  }
  injectConfiguration(doc) {
    const theme=doc.createElement('link');theme.rel='stylesheet';theme.href=new URL('./ketcher-theme.css',import.meta.url).href;doc.head.append(theme);
    const style=doc.createElement('style');style.dataset.workshop='integration';
    // File I/O, global history and geometry preparation belong to the workshop.
    // Indigo operations which rebuild pools require a separate identity bridge.
    const hidden=['open','save','open-file-button','save-file-button','undo','redo','polymer-toggler','3D Viewer button','Layout button','Clean Up button','Aromatize button','Dearomatize button','Add/Remove explicit hydrogens button','Calculate CIP button','Check Structure button','Calculated Values button','miew','layout','clean','arom','dearom','reaction','reaction-arrow','reaction-plus','reaction-map','reaction-arrow-open-angle','reaction-arrow-filled-triangle','reaction-arrow-filled-bow','reaction-arrow-dashed-open-angle','reaction-arrow-failed','reaction-arrow-both-ends-filled-triangle','reaction-arrow-equilibrium-filled-half-bow','reaction-arrow-equilibrium-filled-triangle','reaction-arrow-equilibrium-open-angle','reaction-arrow-unbalanced-equilibrium-filled-half-bow','reaction-arrow-unbalanced-equilibrium-open-half-angle','reaction-arrow-unbalanced-equilibrium-large-filled-half-bow','reaction-arrow-unbalanced-equilibrium-filled-half-triangle','reaction-arrow-elliptical-arc-arrow-filled-bow','reaction-arrow-elliptical-arc-arrow-open-angle','reaction-arrow-retrosynthetic','rgroup','rgroup-label','sgroup','enhanced-stereo','create-monomer','shape-ellipse','shape-rectangle','shape-line','text','images','any-atom','extended-table','template-lib','analyse','recognize'];
    style.textContent=hidden.map(id=>`[data-testid="${id}"],[data-testid="${id}-button"],[data-testid="${id}-tool"],[data-testid="${id}-button-wrapper"]`).join(',')+',[data-testid^="reaction-arrow-"],[data-testid^="reaction-map-"],[data-testid^="shape-"]{display:none!important}';
    // In the pinned 3.18 DOM these four groups contain only excluded tools.
    // Hide the group itself so its background and dropdown-expand SVG do not
    // remain as empty cards; selection/bond/clipboard dropdowns stay visible.
    // Native hiddenControls now removes the unsupported groups, including
    // their dropdown containers. These stable IDs cover old/cached builds.
    style.textContent+='[data-testid="rgroup-drop-down-button"],[data-testid="arrows-drop-down-button"],[data-testid="reaction-mapping-tools-drop-down-button"],[data-testid="shapes-drop-down-button"]{display:none!important}';
    doc.head.append(style);
  }
  keydown(e) {
    // pointer-events cannot block keys in an already focused iframe. A known
    // programmatic pool rebuild owns the frame until it has been bound.
    if(this.suspended){e.preventDefault();e.stopImmediatePropagation();return;}
    if(e.defaultPrevented||e.isComposing)return;
    const editable=this.isTextInput(e.target);
    const command=e.ctrlKey||e.metaKey;
    const key=e.key.toLowerCase();
    if(command&&!e.altKey&&['s','o'].includes(key)) {
      e.preventDefault();e.stopImmediatePropagation();this.gesture=false;this.flush();
      key==='s'?this.callbacks.onSave?.():this.callbacks.onOpen?.();return;
    }
    if(!editable&&command&&!e.altKey&&(key==='z'||key==='y')) {
      e.preventDefault();e.stopImmediatePropagation();this.gesture=false;
      key==='y'||e.shiftKey?this.redo():this.undo();return;
    }
    if(!editable&&command&&!e.altKey&&['delete','backspace'].includes(key)) {
      e.preventDefault();e.stopImmediatePropagation();this.clear();return;
    }
    if(editable)return;
    // The hidden native actions must not remain reachable through shortcuts.
    // Ordinary element, bond, clipboard and basic ring keys stay native.
    const rebuild=command&&key==='l'||e.altKey&&key==='a';
    const unsupported=command&&!e.altKey&&['g','r','m','p'].includes(key)
      ||command&&e.altKey&&['r','d','p'].includes(key)
      ||e.altKey&&['e','t','s','c'].includes(key)
      ||e.shiftKey&&!command&&!e.altKey&&key==='t';
    if(rebuild||unsupported) {
      e.preventDefault();e.stopImmediatePropagation();
      this.callbacks.onError?.(NativeError(rebuild?'请使用顶部的“整理二维排版”。':'此工具暂不支持。请使用原子、化学键或环模板构造分子。'));
    }
  }
  isTextInput(target) {
    // Ketcher focuses its invisible clipboard textarea after canvas clicks.
    // It is an editor hotkey target, not an atom-property text input.
    return !!target.closest?.('input,textarea,select,[contenteditable=true]')&&!target.closest?.('[data-cliparea="true"]');
  }
  changed() {if(this.suspended)return;this.callbacks.onInteraction?.();this.dirty=true;this.schedule();}
  schedule() {if(this.timer!==null)return;this.timer=setTimeout(()=>{this.timer=null;this.flush();},0);}
  flush() {
    if(!this.ready||this.suspended||this.gesture)return false;
    if(!this.dirty){this.readSelection();return false;}
    this.dirty=false;const before=this.model.revision;this.isApplying=true;
    try {const result=this.adapter.commit(this.ketcher.editor.struct(),this.model,'二维绘图');if(result.requiresPush)this.push(this.model.project).catch(e=>this.callbacks.onError?.(e));if(result.warnings?.length)this.callbacks.onNotice?.(result.warnings.join('\n'));this.readSelection();return this.model.revision!==before;}
    catch(error) {this.callbacks.onError?.(error);this.push(this.model.project).catch(e=>this.callbacks.onError?.(e));return false;}
    finally {this.isApplying=false;}
  }
  readSelection() {
    if(!this.ready||this.suspended)return;
    try {
      const tool=this.ketcher.editor.selection()??{},ids=[...this.adapter.selectionIds(this.ketcher.editor.struct(),tool)],key=this.selectionKey(ids);
      if(key!==this.lastSelection){this.lastSelection=key;this.callbacks.onSelection?.(ids);}
    }catch(error){this.callbacks.onError?.(error);}
  }
  setSelection(ids) {
    if(!this.ready||this.suspended)return;
    const wanted=[...ids],key=this.selectionKey(wanted);
    if(key===this.lastSelection)return;
    // selectionChange dispatches synchronously and the host callback may call
    // back into setSelection. Stamp first, then invoke the native setter.
    const previous=this.lastSelection;this.lastSelection=key;
    try {this.ketcher.editor.selection(this.adapter.toolSelection(this.ketcher.editor.struct(),wanted));}
    catch(error){this.lastSelection=previous;this.callbacks.onError?.(error);}
  }
  selectionKey(ids) {return JSON.stringify([...new Set(ids)].sort());}
  push(project, {fit=false} = {}) {
    const serial=++this.pushSerial;this.suspended=true;this.dirty=false;this.gesture=false;
    this.frame.style.pointerEvents='none';this.callbacks.onBusy?.(true);
    const task=this.queue.catch(()=>{}).then(async()=>{
      if(serial!==this.pushSerial)return;
      try {
        await this.adapter.push(this.ketcher,project);
        this.ketcher.editor.clearHistory?.();
        if(fit)this.resetView();
      } finally {
        // A failed import cannot leave the frame permanently locked. An older
        // request must never unlock a newer pending push.
        if(serial===this.pushSerial){this.suspended=false;this.frame.style.pointerEvents='';this.callbacks.onBusy?.(false);this.lastSelection='';this.setSelection(this.model.selection);}
      }
    });
    this.queue=task;return task;
  }
  undo() {this.flush();this.model.undo();}
  redo() {this.flush();this.model.redo();}
  clear() {
    if(!this.ready||this.suspended)return;
    this.gesture=false;this.flush();
    if(this.callbacks.onClear)this.callbacks.onClear();
    else this.callbacks.onError?.(NativeError('请使用文件菜单新建工程。'));
  }
  resetView() {
    const editor=this.ketcher?.editor,struct=editor?.struct();
    this.needsFit=true;
    if(!editor||!struct||!this.visibleFrame())return false;
    // Verified 3.18 APIs: fit the viewport rather than rewriting coordinates.
    // zoomAccordingContent only shrinks an oversized drawing. Start from the
    // standard scale so Reset also restores a previously zoomed-out view.
    editor.zoom(1);
    editor.centerViewportAccordingToStruct(struct);
    editor.event?.zoomChanged?.dispatch();
    this.needsFit=false;
    const bounds=this.frame.getBoundingClientRect();this.lastVisibleSize=[Math.round(bounds.width),Math.round(bounds.height)];
    return true;
  }
  visibleFrame() {
    const frameRect=this.frame.getBoundingClientRect(),area=this.ketcher?.editor?.render?.clientArea?.getBoundingClientRect();
    return this.frame.isConnected&&frameRect.width>0&&frameRect.height>0&&area?.width>0&&area?.height>0;
  }
  async onVisible() {
    await this.initialized;
    let pending;
    do {pending=this.queue;await pending;} while(pending!==this.queue);
    const win=this.frame.ownerDocument.defaultView;
    // The iframe resumes after the outer 2D panel is shown. Let the parent
    // paint and the native ResizeObserver/resize handlers update its viewport.
    await new Promise(resolve=>win.requestAnimationFrame(resolve));
    const bounds=this.frame.getBoundingClientRect();
    if(!this.frame.isConnected||!bounds.width||!bounds.height)return false;
    const size=[Math.round(bounds.width),Math.round(bounds.height)];
    if(this.lastVisibleSize&&size.some((value,i)=>Math.abs(value-this.lastVisibleSize[i])>1))this.needsFit=true;
    this.lastVisibleSize=size;
    const native=this.frame.contentWindow;
    native.dispatchEvent(new native.Event('resize'));
    await new Promise(resolve=>win.requestAnimationFrame(()=>win.requestAnimationFrame(resolve)));
    if(pending!==this.queue){await this.queue;return this.onVisible();}
    // SVG text bounds measured while an iframe is display:none are unusable.
    // Resize alone leaves cached H/charge labels detached from valid bond ends.
    // Rebuild visual objects only after native canvas dimensions are available;
    // this never lays out atoms, rebuilds pools or changes either conformer.
    this.ketcher.editor.render.update(true);
    return this.needsFit?this.resetView():false;
  }
  fitAfterVisible() {return this.onVisible();}
  async settled() {await this.initialized;this.flush();await this.queue;}
}
