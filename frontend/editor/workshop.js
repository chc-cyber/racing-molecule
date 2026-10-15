import {playerError} from './player-error.mjs';
import {browserMode,browserRequest,prepareBrowserRace,cancelBrowserRace} from './browser-service.mjs';
import {validateCompetition} from '../competition.mjs';
// Preserve useful errors when the JSmol Java runtime exports its Error class.
const NativeError=globalThis.Clazz?._Error??globalThis.Error;
import {WorkshopModel,createProject,clone,parseProject,serializeProject,chemistryWarnings,distance,angle,dihedral,geometryPreview} from './model.js';
import {mountPeriodicTable,findElement} from './periodic-table.js';
import {KetcherHost} from './ketcher-host.js';
import {JSmolHost} from './jsmol-host.js';
import {renderKetcherIcons} from './ketcher-icons.js';
import {activeChemistryWarnings,isC60AdditionNotice} from './editor-notices.js';
import {RaceCompletion} from './race-completion.js';

const $=id=>document.getElementById(id),model=new WorkshopModel();
let mode='select',nativeToolPicking=null,dimension=2,selectedBond=null,interactionEpoch=0,requestSerial=0,pendingGeometry=0,loadingRace=false;
let projectNameDraft=null;
let numbers=false,geometryState=null,lastMessage='',mature=null,views=null,switching=false,hoverHit=null;
const ids=()=>[...model.selection],axis=()=>({x:[1,0,0],y:[0,1,0],z:[0,0,1]}[$('rotate-axis').value]);
const geometryKinds={distance:{name:'键长',count:2,unit:'Å',min:.3,max:5,step:.01,measure:distance,method:'setDistance'},angle:{name:'键角',count:3,unit:'°',min:0,max:180,step:.1,measure:angle,method:'setAngle'},dihedral:{name:'二面角',count:4,unit:'°',min:-180,max:180,step:.1,measure:dihedral,method:'setDihedral'}};
const modeInfo={select:['选择','Shift 多选 · 空白拖动转视角'],atom:['原子','周期表选元素 · 点击空白添加 / 点击原子替换'],bond:['化学键','点击键改键型 · 从原子拖出连接或延伸'],move:['移动原子','拖动原子'],rotate:['旋转','选择要旋转的部分 · 属性面板设置轴与角度'],measure:['测量','依次选点 · 面板输入精确几何'],delete:['删除','点击原子 / 键删除 · Delete 删除选择']};
modeInfo.molecule=['移动基团','拖动分子 / 折叠基团平移 · Alt 拖动旋转 · 空白转视角'];
modeInfo.rotateBond=['绕键旋转','先点击键，再拖动一侧原子'];
modeInfo.charge=['形式电荷','选择原子 · 使用电荷增减按钮'];
function human(message){let text=playerError(message);for(const [i,a] of model.project.atoms.entries())text=text.split(a.id).join(`原子 ${i+1} ${a.element}`);for(const [i,b] of model.project.bonds.entries())text=text.split(b.id).join(`键 ${i+1}`);return text;}
function messageDetails(){const warnings=activeChemistryWarnings(model.project).map(human);return lastMessage+(warnings.length?'\n\n当前结构提示：\n'+warnings.join('\n'):'');}
function tell(message,error=false){lastMessage=human(message);$('feedback').textContent=isC60AdditionNotice(message)?'已完成 C₆₀ 加成连接。':lastMessage.split('\n')[0]||'已完成';$('feedback').classList.toggle('error',error);$('detail-message').textContent=messageDetails();}
function completionMessage(message,warnings=[]){return warnings.length?[warnings[0],message,...warnings.slice(1)].filter(Boolean).join('\n'):message;}
function safely(fn){return async(...args)=>{try{return await fn(...args);}catch(e){tell(e.message||String(e),true);}};}
function listen(id,fn){$(id).addEventListener('click',safely(fn));}
function refreshKetcherIcons(){try{renderKetcherIcons();}catch(e){tell(`图标载入失败：${e.message}`,true);}}
function icons(){refreshKetcherIcons();}
function closeMenus(){document.querySelectorAll('.menu[open]').forEach(x=>x.open=false);}
async function api(path,data){if(browserMode)return browserRequest(path,data);const r=await fetch(path,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(data)});let value;try{value=await r.json();}catch{throw NativeError('本地服务没有返回 JSON，请使用支持工坊接口的服务。');}if(!r.ok)throw NativeError(value.error||`请求失败 (${r.status})`);return value;}
function download(name,content,type='application/json'){const a=document.createElement('a'),url=URL.createObjectURL(new Blob([content],{type}));a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function flushDrawing(){if(dimension===2)mature?.flush();else views?.flush();}
function commitProjectName(){if(projectNameDraft===null)return;const name=projectNameDraft.slice(0,120)||'未命名赛车';projectNameDraft=null;if(name!==model.project.name)model.commit('修改名称',p=>{p.name=name;});}
async function settleEditors(){commitProjectName();flushDrawing();if(dimension===2)await mature?.settled();else await views?.settled();}
function resetView(){if(dimension===2)mature?.resetView();else views?.resetView();}
function syncDrawing(fit=false){if(mature?.ready)return mature.push(clone(model.project),{fit}).catch(e=>tell(e.message,true));}
function sync3D(fit=false){if(views?.ready)return views.push(clone(model.project),{fit}).catch(e=>tell(e.message,true));}
function currentBond(){const ordered=ids();if(ordered.length===4)return model.project.bonds.find(b=>[b.a,b.b].includes(ordered[1])&&[b.a,b.b].includes(ordered[2]));return model.project.bonds.find(b=>b.id===selectedBond)||model.project.bonds.find(b=>ordered.length===2&&ordered.includes(b.a)&&ordered.includes(b.b));}
function configureMode(){
  const order=Number($('bond-order').value),native=mode!=='rotate';
  let picking=nativeToolPicking;
  if(!['atom','bond','charge'].includes(mode))picking=null;
  if(mode==='atom')picking=`assignAtom_${$('element').value}`;
  if(mode==='bond'&&picking!=='assignAtom_Xx'&&!/^assignBond_[pm0]$/.test(picking||'')){
    // The native ModelKit bond picker supports 1/2/3. Aromatic and stereo
    // fields remain editable in the explicit properties panel and Ketcher.
    if(![1,2,3].includes(order)){$('bond-order').value='1';$('bond-symbol').textContent='—';}
    picking=`assignBond_${$('bond-order').value}`;
  }
  nativeToolPicking=picking;
  views?.setMode(mode,{element:$('element').value,order:Number($('bond-order').value),rigid:false,native,nativePicking:picking,unified:mode==='bond'&&!/^assignBond_[pm0]$/.test(picking||'')});
}
function nativeMode(next,options={}){
  if(!modeInfo[next])return;
  cancelGeometry(true);mode=next;nativeToolPicking=options.nativePicking??null;
  if(options.element&&findElement(options.element)){
    $('element').value=options.element;$('element').textContent=options.element;periodic.setCurrent(options.element);
    $('element').setAttribute('aria-label',`当前元素：${findElement(options.element).name}`);
    document.querySelectorAll('[data-element]').forEach(button=>{const active=button.dataset.element===options.element;button.classList.toggle('active',active);button.setAttribute('aria-pressed',String(active));});
  }
  if([1,2,3,1.5].includes(Number(options.order))){$('bond-order').value=String(options.order);$('bond-symbol').textContent={'1':'—','2':'=','3':'≡','1.5':'⌁'}[options.order];}
  if(next==='measure'){geometryState={type:'distance',reverse:false};configureGeometry();}
  // The native menu has already chosen its picking/drag mode. Rendering only
  // the surrounding labels avoids resetting that mode inside its callback.
  render();
}
async function nativeAction(name){
  await settleEditors();
  if(name==='periodic-table'){cancelGeometry();periodic.open();return;}
  if(name==='generate'||name==='rebuild')return geometry('rebuild');
  if(['hydrogenate','relax','prepare'].includes(name))return geometry(name,{show:name!=='hydrogenate'});
  if(name==='undo'||name==='redo'){cancelGeometry(true);return model[name]();}
  if(name==='save'){requireAppliedGeometry('保存');download((model.project.name||'my-car')+'.mcar',serializeProject(model.project));tell('工程已保存。');return;}
  if(name==='export'){$('export-type').value='3:mol';return exportFile();}
  throw NativeError('此原生操作尚未接入共同工程。');
}
const periodic=mountPeriodicTable({dialog:$('periodic-dialog'),getSelection:()=>({count:ids().length,element:ids().length===1?model.project.atoms.find(a=>a.id===ids()[0])?.element:null}),onChoose:(symbol,{replaceSelection})=>{flushDrawing();cancelGeometry();if(replaceSelection){if(ids().length!==1)throw NativeError('请只选一个原子。');model.updateAtom(ids()[0],{element:symbol});setElement(symbol,false);tell(`已替换为 ${symbol}。`);}else{setElement(symbol,true);tell(`${findElement(symbol).name} ${symbol} · 点击空白添加，点击原子替换。`);}}});
function setElement(symbol,activate=true){if(!findElement(symbol))throw NativeError(`未知元素 ${symbol}`);$('element').value=symbol;$('element').textContent=symbol;$('element').setAttribute('aria-label',`当前元素：${findElement(symbol).name}`);periodic.setCurrent(symbol);document.querySelectorAll('[data-element]').forEach(b=>{const active=b.dataset.element===symbol;b.classList.toggle('active',active);b.setAttribute('aria-pressed',String(active));});if(activate)setMode('atom');else configureMode();render();}
function select(atomIds,bond=null,{from3D=false}={}){interactionEpoch++;let next=[...atomIds];if(from3D&&mode==='measure'&&next.length===1){const n=geometryKinds[geometryState?.type||'distance'].count,old=ids();next=old.includes(next[0])?old.filter(id=>id!==next[0]):[...(old.length>=n?[]:old),next[0]];}model.selection=new Set(next);selectedBond=bond;if(mode==='measure')configureGeometry();if(!mature?.isApplying)mature?.setSelection(model.selection);views?.setSelection(model.selection);render();}
function cancelGeometry(close=false){if(geometryState?.preview)views?.clearPreview();if(geometryState)geometryState.preview=null;if(close)geometryState=null;}
async function setMode(next,{nativePicking=null}={}){if(!modeInfo[next])return;flushDrawing();cancelGeometry(true);if(['rotate','measure','molecule','rotateBond'].includes(next)&&dimension===2){if(!await setDimension(3))return;}mode=next;nativeToolPicking=nativePicking;$('inspector').hidden=true;if(next==='measure'){geometryState={type:ids().length===4?'dihedral':ids().length===3?'angle':'distance',reverse:false};configureGeometry();}if(next==='rotate')$('inspector').hidden=false;configureMode();render();}
function showDimension(d){dimension=d;document.body.dataset.dimension=d;$('canvas-3d').hidden=d!==3;$('jsmol-loading').hidden=d!==3||!!views?.ready;$('view-2d').setAttribute('aria-pressed',String(d===2));$('view-3d').setAttribute('aria-pressed',String(d===3));if(d===2&&['rotate','measure'].includes(mode)){mode='select';cancelGeometry(true);}render();requestAnimationFrame(()=>{const tool=d===2?mature:views;tool?.onVisible?.().catch(e=>tell(e.message,true));});}
async function setDimension(d){if(d===dimension)return true;if(switching)return false;flushDrawing();cancelGeometry(true);if(d===3&&model.project.atoms.length){switching=true;$('view-3d').disabled=true;const had3D=model.project.has3d;try{if(!await geometry('prepare',{show:false}))return false;}catch(e){if(!had3D)throw e;showDimension(3);configureMode();tell(`三维准备未完成：${e.message}\n请调整新增结构后重试。`,true);return true;}finally{switching=false;$('view-3d').disabled=false;}}showDimension(d);configureMode();return true;}
function require3D(){if(!model.project.has3d)throw NativeError('请先切换三维以生成构象并补氢，或导入三维文件。');}
function beginGeometry(type,ordered=null){flushDrawing();require3D();cancelGeometry(true);showDimension(3);mode='measure';geometryState={type,reverse:false};if(ordered)model.selection=new Set(ordered);$('inspector').hidden=true;configureGeometry();configureMode();render();}
function configureGeometry(){
  if(mode!=='measure')return;const type=geometryState?.type||'distance',reverse=geometryState?.reverse||false,kind=geometryKinds[type];if(geometryState?.preview)views?.clearPreview();geometryState={type,reverse,ids:ids(),baseRevision:model.revision,baseProject:clone(model.project),preview:null,valid:false};
  const n=ids().length;$('geometry-title').textContent=kind.name;$('geometry-unit').textContent=kind.unit;const input=$('geometry-value'),slider=$('geometry-slider');input.min=slider.min=kind.min;input.max=slider.max=kind.max;input.step=slider.step=kind.step;document.querySelectorAll('[data-geometry]').forEach(b=>b.classList.toggle('active',b.dataset.geometry===type));$('fixed-side').textContent=reverse?'末端侧固定':'A 侧固定';
  let message=`依次选择 ${kind.count} 个原子 · 已选 ${n}`;
  if(model.project.has3d&&n===kind.count){try{const value=kind.measure(model.project,...ids());geometryState.value=value;input.value=value.toFixed(type==='distance'?3:2);slider.min=Math.min(kind.min,value);slider.max=Math.max(kind.max,value);slider.value=value;geometryPreview(model.project,type,ids(),value,{reverse});geometryState.valid=true;message=`${ids().map((_,i)=>String.fromCharCode(65+i)).join('—')} · ${value.toFixed(type==='distance'?3:2)} ${kind.unit}`;}catch(e){message=human(e.message);}}
  else if(!model.project.has3d)message='先生成三维，再调整几何';
  $('measurement').textContent=message;for(const el of [input,slider,$('geometry-swap'),$('geometry-apply')])el.disabled=!geometryState.valid;
}
function previewGeometry(value){if(!geometryState?.valid)return;const g=geometryState;if(g.baseRevision!==model.revision){configureGeometry();throw NativeError('设计已变化，旧几何预览已取消。');}if(!Number.isFinite(value))return;interactionEpoch++;g.preview=geometryPreview(g.baseProject,g.type,g.ids,value,{reverse:g.reverse});g.value=value;$('geometry-value').value=value;$('geometry-slider').value=value;views?.showPreview(g.preview.project);$('measurement').textContent=`预览 ${value.toFixed(g.type==='distance'?3:2)} ${geometryKinds[g.type].unit}`;render();}
function applyGeometry(){const g=geometryState;if(!g?.valid)throw NativeError('请选择完整、可分离的几何对象。');if(g.baseRevision!==model.revision)throw NativeError('设计已变化，请重新选择几何对象。');const value=Number($('geometry-value').value);views?.clearPreview();model[geometryKinds[g.type].method](...g.ids,value,{reverse:g.reverse});geometryState={type:g.type,reverse:g.reverse};configureGeometry();render();tell(`${geometryKinds[g.type].name}已应用。`);}
function torsionForBond(){require3D();const b=currentBond();if(!b)throw NativeError('请先选择一条键。');const adjacent=id=>model.project.bonds.filter(x=>x.a===id||x.b===id).map(x=>x.a===id?x.b:x.a).filter(x=>x!==b.a&&x!==b.b).sort((a,c)=>(model.project.atoms.find(x=>x.id===a).element==='H')-(model.project.atoms.find(x=>x.id===c).element==='H'));for(const a of adjacent(b.a))for(const d of adjacent(b.b)){try{const ordered=[a,b.a,b.b,d],value=dihedral(model.project,...ordered);geometryPreview(model.project,'dihedral',ordered,value);beginGeometry('dihedral',ordered);return;}catch{}}throw NativeError('此键没有可分离且方向明确的两侧。环内键可用直接原子编辑。');}
function render(){
  const p=model.project;$('counts').textContent=`${p.atoms.length} 原子 · ${p.bonds.length} 键`;$('empty-hint').hidden=!!p.atoms.length;$('undo').disabled=!model.canUndo;$('redo').disabled=!model.canRedo;$('revision-label').textContent=`修订 ${model.revision}`;
  $('send-track').disabled=loadingRace||!!pendingGeometry||switching||!!unappliedGeometry();
  if(projectNameDraft===null&&document.activeElement!==$('project-name'))$('project-name').value=p.name;
  const warnings=activeChemistryWarnings(p);$('warning-count').hidden=!warnings.length;$('warning-count').textContent=warnings.length?`⚠ ${warnings.length}`:'';$('feedback-details').setAttribute('aria-label',warnings.length?`查看提示详情，${warnings.length} 处结构提示`:'查看提示详情');
  $('relax').hidden=dimension!==3;$('layout-2d').hidden=dimension!==2;
  const picked=ids().map(id=>p.atoms.find(a=>a.id===id)).filter(Boolean),a=picked.length===1?picked[0]:null,b=currentBond();$('selection-info').textContent=picked.length?`已选 ${picked.length}：`+picked.slice(0,8).map(a=>`${p.atoms.indexOf(a)+1} ${a.element}`).join(' → ')+(picked.length>8?'…':''):'未选择';$('apply-atom').disabled=!a;$('atom-properties').hidden=!a;$('bond-properties').hidden=!(b||picked.length===2||picked.length===4);$('transform-properties').hidden=dimension!==3||!picked.length;
  if(a){$('atom-element').value=a.element;$('atom-element-picker').textContent=`${a.element} · 周期表`;$('atom-charge').value=a.charge;$('atom-isotope').value=a.isotope||0;$('atom-chiral').value=a.chiralTag||'CHI_UNSPECIFIED';$('atom-periodic').textContent=a.element;$('charge-value').textContent=a.charge>0?`+${a.charge}`:a.charge;}
  $('bond-info').textContent=b?`${p.atoms.findIndex(a=>a.id===b.a)+1} — ${p.atoms.findIndex(a=>a.id===b.b)+1}`:'选两端点可添加键';if(b){$('edit-bond-order').value=b.order;$('bond-direction').value=b.direction||'NONE';$('bond-stereo').value=b.stereo||'STEREONONE';}
  $('context-atom').hidden=!a;$('context-bond').hidden=!b;$('context-toolbar').hidden=dimension!==3||!picked.length||!['select','move','rotate'].includes(mode);document.querySelectorAll('[data-context-order]').forEach(x=>x.classList.toggle('active',b?.order===Number(x.dataset.contextOrder)));
  const connecting=mode==='bond'&&nativeToolPicking==='assignAtom_Xx',deletingBond=mode==='bond'&&nativeToolPicking==='assignBond_0';
  $('bond-tool-symbol').textContent=$('bond-symbol').textContent;
  $('geometry-panel').hidden=dimension!==3||mode!=='measure';$('view-badge').textContent='JSmol · 3D · Å';$('view-2d').classList.toggle('active',dimension===2);$('view-3d').classList.toggle('active',dimension===3);$('mode-name').textContent=dimension===2?'二维绘图':deletingBond?'删除键':mode==='atom'?`原子 · ${$('element').value}`:modeInfo[mode][0];$('mode-help').textContent=picked.length?`已选 ${picked.length} 个原子`:'';document.querySelectorAll('[data-mode]').forEach(x=>{const active=x.dataset.mode==='delete'?mode==='delete'||deletingBond:x.dataset.mode==='bond'?mode==='bond'&&!deletingBond:x.dataset.mode===mode;x.classList.toggle('active',active);x.setAttribute('aria-pressed',String(active));});
  refreshKetcherIcons();
}
model.addEventListener('change',()=>{interactionEpoch++;selectedBond=null;cancelGeometry();if(mode==='measure')configureGeometry();render();$('validation-state').textContent='';try{localStorage.setItem('molecular-car-workshop-draft',serializeProject(model.project));$('save-state').title='已暂存到此浏览器';$('save-state').setAttribute('aria-label','已暂存');}catch{$('save-state').title='暂存失败，请保存工程';}const warnings=activeChemistryWarnings(model.project);tell(warnings.length?`${warnings.length} 处结构提示\n${warnings.join('\n')}`:'已更新。');if(!mature?.isApplying)syncDrawing();if(!views?.isApplying)sync3D();});
async function geometry(operation,{show=true}={}){
  await settleEditors();if(operation==='relax'&&dimension!==3)throw NativeError('UFF 优化在三维视图执行，请先进入三维。');cancelGeometry(true);const baseRevision=model.revision,had3D=model.project.has3d,epoch=interactionEpoch,serial=++requestSerial,fixedIds=operation==='relax'&&$('fixed-selection').checked?ids():[];pendingGeometry++;render();$('busy-state').textContent=operation==='prepare'?'构建与补氢中':'计算中…';tell(operation==='prepare'?'正在准备三维结构…':operation==='rebuild'?'正在重建三维构象…':operation==='hydrogenate'?'正在补氢…':operation==='relax'?'正在进行 UFF 优化…':'正在整理二维排版…');
    try{const result=await api('/api/editor/geometry',{project:clone(model.project),operation,fixedIds,revision:baseRevision});flushDrawing();if(result.revision!==baseRevision||epoch!==interactionEpoch||serial!==requestSerial||!model.applyAsync(baseRevision,result.project,operation==='prepare'?'补氢并准备新增三维结构':operation==='rebuild'?'重新生成三维坐标':operation==='hydrogenate'?'补齐氢':operation==='relax'?'UFF 优化':'整理二维')){tell('设计已更新，旧计算结果已丢弃。');return false;}if(show&&operation!=='layout2d'){showDimension(3);configureMode();}if(operation==='rebuild'||operation==='prepare'&&!had3D)await sync3D(true);tell(completionMessage(result.message||'操作已完成。',result.warnings));return true;}finally{pendingGeometry--;$('busy-state').textContent=pendingGeometry?'计算中…':'';render();}
}
async function importContent(content,filename){await settleEditors();cancelGeometry(true);const revision=model.revision,epoch=interactionEpoch;if(/\.(mcar|json)$/i.test(filename)){model.replace(parseProject(content),'打开工程');select([]);mode='select';showDimension(model.project.has3d?3:2);configureMode();syncDrawing(true);sync3D(true);tell('已打开工程。');return;}const result=await api('/api/editor/import',{content,filename});flushDrawing();if(revision!==model.revision||epoch!==interactionEpoch){tell('读取期间设计已变化，请再次导入。');return;}model.replace(result.project,'导入分子');select([]);mode='select';showDimension(result.project.has3d?3:2);configureMode();syncDrawing(true);sync3D(true);tell(completionMessage(`已导入 ${filename}。`,result.warnings));}
function unappliedGeometry(){return geometryState?.preview&&geometryState.value!==geometryKinds[geometryState.type].measure(geometryState.baseProject,...geometryState.ids);}
function requireAppliedGeometry(action){flushDrawing();if(unappliedGeometry())throw NativeError(`请先应用或取消几何预览，再${action}。`);}
async function exportFile(){await settleEditors();requireAppliedGeometry('导出');const [d,format]=$('export-type').value.split(':');const result=await api('/api/editor/export',{project:clone(model.project),dimension:Number(d),format});download(result.filename,result.content,'chemical/x-mdl-molfile');closeMenus();tell(completionMessage(`已导出 ${result.filename}。`,result.warnings));}
async function validateRace(){await settleEditors();const revision=model.revision,result=await api('/api/editor/validate',{project:clone(model.project),revision});if(revision!==model.revision){tell('设计已变化，旧检查已丢弃。');return false;}$('validation-state').textContent=result.eligible?`试跑检查通过 · 补氢后 ${result.preparedAtomCount} 原子`:`暂不能试跑：${result.errors.join('；')}。`;tell($('validation-state').textContent,!result.eligible);if(!result.eligible)$('message-dialog').showModal();return result.eligible;}
const handoffKey='molecular-car-workshop-handoff';
let raceRevision=null;
async function raceRequest(path,data){
  const controller=new AbortController(),timer=setTimeout(()=>controller.abort(),30000);
  try{
    const response=await fetch(path,{cache:'no-store',signal:controller.signal,...(data===undefined?{}:{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(data)})});
    const result=await response.json();
    if(!response.ok){const error=NativeError(result.error||'准备请求失败。');error.missing=path.includes('/status?')&&response.status===400;error.recoverable=response.status>=500;throw error;}
    return result;
  }catch(error){if(error instanceof TypeError||error.name==='AbortError'||error instanceof SyntaxError){const disconnected=NativeError('连接已中断，请重试。');disconnected.recoverable=true;throw disconnected;}throw error;}
  finally{clearTimeout(timer);}
}
function saveSnapshot(project){
  const content=serializeProject(project);
  localStorage.setItem('molecular-car-workshop-draft',content);
  download(`${project.name||'molecule'}.mcar`,content);
}
const raceCompletion=new RaceCompletion({
  uuid:()=>crypto.randomUUID(),
  competition:()=>trackCompetition,
  settle:async()=>{await settleEditors();requireAppliedGeometry('进入赛道');if(pendingGeometry||switching)throw NativeError('请等待当前计算完成。');raceRevision=model.revision;return clone(model.project);},
  save:project=>{try{saveSnapshot(project);}catch{throw NativeError('浏览器暂存失败，请保存工程后重试。');}},
  export:project=>raceRequest('/api/editor/export',{project,dimension:project.has3d?3:2,format:'mol'}),
  prepare:job=>raceRequest('/api/editor/race/prepare',job),
  status:id=>raceRequest('/api/editor/race/status?id='+encodeURIComponent(id)),
  commit:id=>{if((raceRevision!==null&&raceRevision!==model.revision)||raceCompletion.job?.design!==JSON.stringify(model.project))throw NativeError('设计已更新，请返回编辑后重新准备。');return raceRequest('/api/editor/race/commit',{id});},
  cancel:id=>raceRequest('/api/editor/race/cancel',{id}),
  active:()=>raceRequest('/api/state'),
  delay:()=>new Promise(resolve=>setTimeout(resolve,450)),
  remember:job=>{try{if(job)sessionStorage.setItem(handoffKey,JSON.stringify(job));else sessionStorage.removeItem(handoffKey);}catch{/* The current in-memory job remains usable if session storage is unavailable. */}},
  navigate:generation=>location.assign('/?workshop='+generation),
  notify:({state,message,retry,locked})=>{
    loadingRace=locked;
    for(const element of document.querySelectorAll('.app-header,.topbar,main,footer'))element.inert=locked;
    $('tooltip').hidden=true;
    $('race-title').textContent=state==='failed'?'准备未完成':state==='committing'||state==='success'?'进入赛道':'准备赛车';
    $('race-message').textContent=message;
    $('race-retry').hidden=!retry;
    $('race-cancel').textContent=state==='failed'?'返回编辑':'取消';
    $('race-cancel').disabled=['committing','success','cancelling'].includes(state);
    $('race-message').classList.toggle('error',state==='failed');
    if(locked&&!$('race-dialog').open)$('race-dialog').showModal();
    if(!locked&&$('race-dialog').open)$('race-dialog').close();
    tell(message,state==='failed');render();
  }
});
let trackCompetition=validateCompetition();
let browserPreparing=false,browserSendSerial=0;
async function sendToTrack(){
  if(raceCompletion.locked||browserPreparing)return;
  if(browserMode)return browserSendToTrack();
  trackCompetition=validateCompetition((await raceRequest('/api/state')).competition);
  return raceCompletion.start();
}
async function browserSendToTrack(){
  if(browserPreparing)return;
  const serial=++browserSendSerial;browserPreparing=true;loadingRace=true;
  for(const el of document.querySelectorAll('.app-header,.topbar,main,footer'))el.inert=true;
  $('race-title').textContent='准备赛车';$('race-message').textContent='保存并检查结构…';$('race-retry').hidden=true;$('race-cancel').textContent='返回编辑';$('race-cancel').disabled=false;$('race-dialog').showModal();
  try{
    await settleEditors();if(serial!==browserSendSerial)return;requireAppliedGeometry('进入赛道');
    if(pendingGeometry||switching)throw NativeError('请等待当前计算完成。');
    const project=clone(model.project);
    let savedDraft=true;
    try{localStorage.setItem('molecular-car-workshop-draft',serializeProject(project));}
    catch{savedDraft=false;$('save-state').title='暂存失败，请使用保存工程下载备份';$('save-state').setAttribute('aria-label','暂存失败');}
    const exported=await browserRequest('/api/editor/export',{project,dimension:project.has3d?3:2,format:'mol'});
    if(serial!==browserSendSerial)return;
    await prepareBrowserRace(exported,message=>$('race-message').textContent=message);
    const notice=savedDraft?'赛车已就绪。':'赛车已就绪，但浏览器暂存失败。请返回工坊，点击“保存工程”下载备份。';
    tell(notice);
    if(!savedDraft&&parent!==window)parent.postMessage({channel:'racing-molecule',type:'workshop-notice',message:notice},location.origin);
  }catch(error){if(serial===browserSendSerial){tell(error.message,true);$('race-message').textContent=error.message;}}
  finally{if(serial===browserSendSerial){browserPreparing=false;loadingRace=false;for(const el of document.querySelectorAll('.app-header,.topbar,main,footer'))el.inert=false;if($('race-dialog').open)$('race-dialog').close();}}
}
function returnFromBrowserPreparation(){++browserSendSerial;cancelBrowserRace();browserPreparing=false;loadingRace=false;for(const el of document.querySelectorAll('.app-header,.topbar,main,footer'))el.inert=false;if($('race-dialog').open)$('race-dialog').close();tell('已返回编辑。');}
if(browserMode){
  if(parent!==window)parent.postMessage({channel:'racing-molecule',type:'workshop-ready'},location.origin);
  document.addEventListener('click',event=>{const link=event.target.closest('a[data-app-route]');if(link&&parent!==window){event.preventDefault();parent.postMessage({channel:'racing-molecule',type:'navigate',route:link.dataset.appRoute},location.origin);}});
}

document.querySelectorAll('[data-mode]').forEach(b=>b.addEventListener('click',safely(()=>setMode(b.dataset.mode,{nativePicking:b.dataset.nativePicking??null}))));document.querySelectorAll('[data-element]').forEach(b=>b.addEventListener('click',safely(()=>setElement(b.dataset.element))));
for(const id of ['periodic-table-trigger','element'])listen(id,()=>{cancelGeometry();periodic.open();});for(const id of ['atom-periodic','atom-element-picker'])listen(id,()=>periodic.open({replaceSelection:true}));
document.querySelectorAll('[data-bond-style]').forEach(b=>b.addEventListener('click',safely(()=>{const [order]=b.dataset.bondStyle.split(':');$('bond-order').value=order;$('bond-symbol').textContent={'1':'—','2':'=','3':'≡','1.5':'⌁'}[order];closeMenus();return setMode('bond');})));
document.querySelectorAll('[data-context-order]').forEach(b=>b.addEventListener('click',safely(()=>{const bond=currentBond();if(bond)model.updateBond(bond.id,{order:Number(b.dataset.contextOrder)});})));
for(const [id,direction] of [['context-wedge','BEGINWEDGE'],['context-dash','BEGINDASH']])listen(id,()=>{const b=currentBond();if(b)model.updateBond(b.id,{order:1,direction});});listen('flip-bond',()=>{const b=currentBond();if(b)model.updateBond(b.id,{a:b.b,b:b.a,...(b.stereoAtoms?.length?{stereoAtoms:[...b.stereoAtoms].reverse()}:{})});});
for(const [id,delta] of [['charge-down',-1],['charge-up',1]])listen(id,()=>{if(ids().length===1){const a=model.project.atoms.find(a=>a.id===ids()[0]);model.updateAtom(a.id,{charge:a.charge+delta});}});
listen('view-2d',()=>setDimension(2));listen('view-3d',()=>setDimension(3));for(const action of ['undo','redo'])listen(action,()=>{flushDrawing();cancelGeometry();model[action]();});listen('reset-view',resetView);listen('fit-view',resetView);listen('select-all',()=>select(model.project.atoms.map(a=>a.id)));listen('clear-selection',()=>{cancelGeometry();select([]);});listen('delete-selection',()=>{flushDrawing();model.deleteSelection();});
listen('open-properties',()=>{flushDrawing();$('inspector').hidden=!$('inspector').hidden;});listen('retry-editor',()=>location.reload());listen('show-numbers',()=>{numbers=!numbers;$('show-numbers').setAttribute('aria-pressed',String(numbers));views?.setNumbers(numbers);});listen('toggle-inspector',()=>{$('inspector').hidden=!$('inspector').hidden;});listen('close-inspector',()=>{$('inspector').hidden=true;});
listen('new-project',()=>{flushDrawing();cancelGeometry(true);model.replace(createProject(),'新建工程');select([]);mode='select';showDimension(2);configureMode();syncDrawing(true);sync3D(true);closeMenus();tell('已新建工程。');});$('project-name').addEventListener('input',event=>{projectNameDraft=event.target.value;});$('project-name').addEventListener('change',safely(commitProjectName));
$('open-file').addEventListener('change',safely(async e=>{const file=e.target.files[0];if(!file)return;if(file.size>3000000)throw NativeError('文件超过 3 MB。');await importContent(await file.text(),file.name);e.target.value='';closeMenus();}));listen('save-project',async()=>{await settleEditors();requireAppliedGeometry('保存');download(`${model.project.name||'molecule'}.mcar`,serializeProject(model.project));tell('工程已保存。');});listen('export-file',exportFile);
for(const [id,name] of [['example-tpb','tetraphenylbenzene.mol'],['example-chiral','chiral.mol']])listen(id,async()=>{const r=await fetch('./examples/'+name);if(!r.ok)throw NativeError('示例文件未找到。');await importContent(await r.text(),name);closeMenus();});listen('generate-3d',()=>geometry('rebuild'));listen('add-hydrogens',()=>geometry('hydrogenate',{show:false}));listen('relax',()=>geometry('relax'));listen('layout-2d',()=>geometry('layout2d'));
listen('apply-atom',()=>{if(ids().length!==1)throw NativeError('请只选一个原子。');model.updateAtom(ids()[0],{element:$('atom-element').value,charge:Number($('atom-charge').value),isotope:Number($('atom-isotope').value),chiralTag:$('atom-chiral').value});tell('原子属性已应用。');});
listen('apply-bond',()=>{const ordered=ids();if(![2,4].includes(ordered.length))throw NativeError('选 2 个端点；E/Z 选 4 个连续原子。');const b=currentBond(),stereo=$('bond-stereo').value,patch={order:Number($('edit-bond-order').value),direction:$('bond-direction').value,stereo};if(!['STEREONONE','STEREOANY'].includes(stereo)){if(ordered.length!==4||patch.order!==2)throw NativeError('E/Z 需选 A—B—C—D，中间为双键。');patch.stereoAtoms=b&&b.a===ordered[2]?[ordered[3],ordered[0]]:[ordered[0],ordered[3]];}else patch.stereoAtoms=[];if(b)model.updateBond(b.id,patch);else model.addBond(ordered.length===4?ordered[1]:ordered[0],ordered.length===4?ordered[2]:ordered[1],patch.order,patch);tell('键属性已应用。');});
listen('translate',()=>{flushDrawing();require3D();model.moveSelection(['move-x','move-y','move-z'].map(x=>Number($(x).value)),{dimension:3,rigid:false});});listen('rotate',()=>{flushDrawing();require3D();model.rotateSelection(axis(),Number($('rotate-degrees').value),{rigid:false});});
for(const type of Object.keys(geometryKinds))listen(`set-${type}`,()=>beginGeometry(type));listen('bond-torsion',torsionForBond);$('geometry-value').addEventListener('input',safely(()=>previewGeometry(Number($('geometry-value').value))));$('geometry-slider').addEventListener('input',safely(()=>previewGeometry(Number($('geometry-slider').value))));listen('geometry-apply',applyGeometry);listen('geometry-cancel',()=>{cancelGeometry();configureGeometry();render();});listen('geometry-swap',()=>{const value=geometryState.value;geometryState.reverse=!geometryState.reverse;configureGeometry();previewGeometry(value);});listen('close-geometry',()=>{cancelGeometry(true);return setMode('select');});
listen('send-track',sendToTrack);listen('race-retry',()=>raceCompletion.retry());listen('race-cancel',()=>{if(browserMode){returnFromBrowserPreparation();}else raceCompletion.cancel();});$('race-dialog').addEventListener('cancel',e=>{e.preventDefault();if(browserMode){returnFromBrowserPreparation();}else raceCompletion.cancel();});listen('feedback-details',()=>{$('detail-message').textContent=messageDetails();$('message-dialog').showModal();});listen('warning-count',()=>{$('detail-message').textContent=messageDetails();$('message-dialog').showModal();});listen('help',()=>$('help-dialog').showModal());document.querySelectorAll('[data-close-dialog]').forEach(b=>b.addEventListener('click',()=>b.closest('dialog').close()));document.addEventListener('pointerdown',e=>{if(!e.target.closest('.menu'))closeMenus();});
document.addEventListener('keydown',e=>{if(e.defaultPrevented||e.isComposing||document.querySelector('dialog[open]'))return;if((e.ctrlKey||e.metaKey)&&e.key.toLowerCase()==='s'){e.preventDefault();$('save-project').click();return;}if(['INPUT','TEXTAREA','SELECT'].includes(e.target.tagName))return;if((e.ctrlKey||e.metaKey)&&['z','y'].includes(e.key.toLowerCase())){e.preventDefault();flushDrawing();cancelGeometry();e.key.toLowerCase()==='y'||e.shiftKey?model.redo():model.undo();return;}if((e.ctrlKey||e.metaKey)&&e.key.toLowerCase()==='a'){e.preventDefault();select(model.project.atoms.map(a=>a.id));return;}if(e.ctrlKey||e.metaKey||e.altKey)return;if(e.key==='Delete'){e.preventDefault();flushDrawing();safely(()=>model.deleteSelection())();return;}if(e.key==='Escape'){views?.cancelGesture();cancelGeometry(true);select([]);setMode('select');return;}if(dimension!==3)return;const key=e.key.toLowerCase(),symbol={c:'C',n:'N',o:'O',h:'H',f:'F',p:'P',s:'S',i:'I'}[key];if(symbol&&hoverHit?.atomId){e.preventDefault();safely(()=>{model.updateAtom(hoverHit.atomId,{element:symbol});setElement(symbol,false);})();return;}if(['1','2','3'].includes(key)){e.preventDefault();const b=currentBond();if(b)safely(()=>model.updateBond(b.id,{order:Number(key)}))();else{$('bond-order').value=key;$('bond-symbol').textContent={'1':'—','2':'=','3':'≡'}[key];setMode('bond');}return;}const next={v:'select',a:'atom',b:'bond',m:'move',r:'rotate',g:'measure',d:'delete'}[key];if(next){e.preventDefault();setMode(next);}else if(symbol){e.preventDefault();setElement(symbol);}});
let tooltipTarget=null;function showTooltip(target){const text=target?.getAttribute('title');if(!text)return;tooltipTarget=target;const tip=$('tooltip');tip.textContent=text;tip.hidden=false;const r=target.getBoundingClientRect();tip.style.left=`${Math.max(4,Math.min(innerWidth-tip.offsetWidth-5,r.left))}px`;tip.style.top=`${r.bottom+6+tip.offsetHeight>innerHeight?r.top-tip.offsetHeight-6:r.bottom+6}px`;}
document.addEventListener('pointerover',e=>{const target=e.target.closest('[title]');if(target&&!target.closest('dialog'))showTooltip(target);});document.addEventListener('pointerout',e=>{if(tooltipTarget&&!tooltipTarget.contains(e.relatedTarget)){tooltipTarget=null;$('tooltip').hidden=true;}});document.addEventListener('focusin',e=>showTooltip(e.target.closest('[title]')));document.addEventListener('focusout',()=>{$('tooltip').hidden=true;});
window.workshop={model,select,setDimension,setMode,setElement,periodic,geometry,beginGeometry,previewGeometry,applyGeometry,importContent,exportFile,sendToTrack,validateRace,raceCompletion,get mature(){return mature;},get views(){return views;},get jsmol(){return views;},get dimension(){return dimension;},get pendingGeometry(){return pendingGeometry;},get geometryState(){return geometryState;},get mode(){return mode;}};
try{const draft=localStorage.getItem('molecular-car-workshop-draft');if(draft){model.replace(parseProject(draft),'恢复暂存');tell('已恢复浏览器暂存。');}else tell('编辑器已就绪。');}catch(e){tell(`暂存读取：${e.message}`,true);}render();icons();
mature=new KetcherHost($('ketcher-frame'),{model,onSelection:atomIds=>select(atomIds),onInteraction:()=>{interactionEpoch++;},onError:e=>tell(e.message,true),onNotice:message=>tell(message),onOpen:()=>$('open-file').click(),onSave:()=>$('save-project').click(),onClear:()=>{cancelGeometry(true);model.commit('清空分子',p=>{p.atoms=[];p.bonds=[];p.fragments=[];p.has3d=false;});select([]);},onReady:()=>{$('editor-loading').hidden=true;tell('二维编辑器已就绪。');}});
mature.initialized.catch(e=>{$('editor-loading').classList.add('error');$('editor-loading').querySelector('span').textContent=e.message;$('retry-editor').hidden=false;tell(e.message,true);});
views=new JSmolHost($('canvas-3d'),{model,unifiedBonds:true,nativeMenus:false,onSelection:(atomIds,bond)=>select(atomIds,bond,{from3D:true}),onInteraction:()=>{interactionEpoch++;},onHover:hit=>{hoverHit=hit;},onMode:nativeMode,onNotice:message=>tell(message),onAction:safely(nativeAction),onError:e=>tell(e.message,true),onReady:()=>{$('jsmol-loading').hidden=true;configureMode();}});
views.initialized.catch(e=>{$('jsmol-loading').classList.add('error');$('jsmol-loading').textContent=e.message;tell(e.message,true);});
try{const pending=!browserMode&&sessionStorage.getItem(handoffKey);if(pending){const job=JSON.parse(pending);if(typeof job.id==='string'&&typeof job.content==='string'&&typeof job.filename==='string')raceCompletion.restore(job);else sessionStorage.removeItem(handoffKey);}}catch(e){tell(`准备恢复失败：${e.message}`,true);}
