import {isRaceInProgress} from '../competition.mjs';
import {readPreferences,writePreferences,defaults} from '../app-state.mjs';
const main=document.querySelector('#main'),nav=document.querySelector('#main-nav');
const origin=location.origin,appRoot=new URL('../',import.meta.url);
const frames=new Map(),waiters=new Map();
let preferences=readPreferences(),route='',returnRoute='home',toastTimer,preparing=null,raceState=null,raceStarted=false,trackReadyPromise;
const channel='racing-molecule';
let preloadWorker,preloadTimer,needHetero=false,downloadHide;
const downloadStatus=document.createElement('aside');downloadStatus.id='download-status';downloadStatus.hidden=true;downloadStatus.setAttribute('role','status');
downloadStatus.innerHTML='<span></span><progress max="100" value="0" aria-label="赛道资源下载进度"></progress>';document.body.append(downloadStatus);
function downloadProgress(percent,label='下载赛道资源'){
 clearTimeout(downloadHide);downloadStatus.hidden=percent==null;if(percent==null)return;
 downloadStatus.querySelector('span').textContent=`${label} ${percent}%`;downloadStatus.querySelector('progress').value=percent;
 if(percent===100)downloadHide=setTimeout(()=>downloadStatus.hidden=true,1500);
}
function schedulePreload(){
  clearTimeout(preloadTimer);
  preloadTimer=setTimeout(()=>{
    if(preloadWorker||!navigator.locks||!globalThis.caches)return;
    if(route==='garage'){
      const editor=frames.get('garage')?.contentWindow?.workshop;
      if(!editor?.mature?.ready){schedulePreload();return;}
    }
    const start=()=>{
      if(preloadWorker)return;
      preloadWorker=new Worker(new URL('browser/preload.mjs',appRoot),{type:'module'});
      preloadWorker.onmessage=({data})=>{document.documentElement.dataset.racePreload=data.status;if(!preparing)downloadProgress(data.percent);};
      preloadWorker.postMessage({hetero:needHetero});
      preloadWorker.onerror=()=>{document.documentElement.dataset.racePreload='unavailable';};
    };
    if(globalThis.requestIdleCallback)requestIdleCallback(start,{timeout:5000});else start();
  },2000);
}
const post=(frame,data)=>frame?.contentWindow?.postMessage({channel,...data},origin);
function notify(text){const el=document.querySelector('#toast');el.textContent=text;el.classList.add('visible');clearTimeout(toastTimer);toastTimer=setTimeout(()=>el.classList.remove('visible'),4000);}
function locked(){return isRaceInProgress(raceState,raceStarted);}
function ensureFrame(name){
  if(frames.has(name))return frames.get(name);
  const panel=document.createElement('section');panel.className='app-screen';panel.dataset.screen=name;panel.inert=true;panel.setAttribute('aria-hidden','true');
  const iframe=document.createElement('iframe');iframe.title=name==='garage'?'分子赛车工坊':'分子赛道';iframe.id=name+'-frame';iframe.allow='clipboard-read; clipboard-write';
  iframe.src=new URL(name==='garage'?'editor/index.html':'race.html',appRoot).href;
  panel.append(iframe);document.body.append(panel);frames.set(name,iframe);return iframe;
}
function broadcastPreferences(){for(const frame of frames.values())post(frame,{type:'preferences',data:preferences});}
function showFrame(name){const frame=ensureFrame(name);post(frame,{type:'screen-visibility',active:true});frame.parentElement.classList.add('active');frame.parentElement.inert=false;frame.parentElement.setAttribute('aria-hidden','false');}
function settingsPage(){return `<section class="settings-page"><div class="settings-heading"><h1>设置</h1><a class="button button-neutral" href="#${returnRoute}">返回</a></div><form id="settings-form"><section class="settings-section"><h2>显示</h2><div class="preference-row"><label for="representation-preference">分子模型</label><select id="representation-preference"><option value="ball-stick">球棍模型</option><option value="space-filling">空间填充</option></select></div><div class="preference-row"><label for="quality-preference">画面质量</label><select id="quality-preference"><option value="low">低</option><option value="balanced">均衡</option><option value="high">高</option></select></div><label class="motion-setting"><input id="motion-preference" type="checkbox">动画过渡</label></section><p id="settings-error" class="form-errors" role="status" hidden></p><div class="settings-actions"><button type="button" class="button button-neutral" id="reset-preferences">恢复默认</button><button type="submit" class="button button-primary">保存设置</button></div></form></section>`;}
function bindSettings(){
  document.querySelector('#representation-preference').value=preferences.representation;document.querySelector('#quality-preference').value=preferences.quality;document.querySelector('#motion-preference').checked=preferences.motion;
  document.querySelector('#reset-preferences').addEventListener('click',()=>{document.querySelector('#representation-preference').value=defaults.representation;document.querySelector('#quality-preference').value=defaults.quality;document.querySelector('#motion-preference').checked=defaults.motion;});
  document.querySelector('#settings-form').addEventListener('submit',event=>{event.preventDefault();try{
    preferences=writePreferences({representation:document.querySelector('#representation-preference').value,quality:document.querySelector('#quality-preference').value,motion:document.querySelector('#motion-preference').checked});
    document.querySelector('#settings-error').hidden=true;broadcastPreferences();notify('设置已保存。');
  }catch(error){const el=document.querySelector('#settings-error');el.hidden=false;el.textContent=error.message;}});
}
function render(){
  let next=location.hash.slice(1)||'home';if(next==='main')return;
  if(['standard','custom'].includes(next))next='track';
  if(!['home','garage','track','settings'].includes(next))next='home';
  if(next==='settings'&&route&&route!=='settings')returnRoute=route;
  if(route==='track'&&next!=='track')post(frames.get('track'),{type:'pause-race'});
  route=next;document.body.dataset.page=next;
  for(const frame of frames.values()){post(frame,{type:'screen-visibility',active:false});frame.parentElement.classList.remove('active');frame.parentElement.inert=true;frame.parentElement.setAttribute('aria-hidden','true');}
  main.hidden=['garage','track'].includes(next);document.querySelector('.site-header').hidden=['garage','track'].includes(next);
  main.innerHTML=next==='home'?homePage():next==='settings'?settingsPage():'';
  if(next==='settings')bindSettings();else if(next==='garage'||next==='track')showFrame(next);
  document.title=`Racing Molecule — ${{home:'分子赛车',garage:'分子工坊',track:'分子赛道',settings:'设置'}[next]}`;
  main.focus({preventScroll:true});
}
async function trackReady(){
  const frame=ensureFrame('track');if(frame.dataset.ready==='true')return frame;
  trackReadyPromise??=new Promise((resolve,reject)=>{const timer=setTimeout(()=>{waiters.delete('ready');trackReadyPromise=null;reject(new Error('赛道载入失败，请刷新页面。'));},45000);waiters.set('ready',{resolve:()=>{clearTimeout(timer);resolve();}});});await trackReadyPromise;return frame;
}
function finishPreparation(error){
  if(!preparing)return;downloadProgress(null);const task=preparing;preparing=null;clearTimeout(task.timer);post(frames.get('garage'),{type:'prepare-result',id:task.id,...(error?{error}:{})});
  if(!error){location.hash='track';}
}
addEventListener('message',async event=>{
  if(event.origin!==origin||event.data?.channel!==channel)return;
  const workshop=frames.get('garage'),track=frames.get('track');
  if(event.source!==workshop?.contentWindow&&event.source!==track?.contentWindow)return;
  const {type,id,data}=event.data;
  if(type==='navigate'){const target=event.data.route;if(['home','garage','track','settings'].includes(target))location.hash=target;return;}
  if(type==='workshop-notice'&&event.source===workshop?.contentWindow){notify(String(event.data.message||''));return;}
  if(type==='tool-progress'&&event.source===workshop?.contentWindow){downloadProgress(event.data.percent,'下载计算工具');return;}
  if(type==='resource-needs'&&event.source===workshop?.contentWindow){needHetero=!!event.data.hetero;preloadWorker?.postMessage({hetero:needHetero});return;}
  if(type==='workshop-ready'){workshop.dataset.ready='true';schedulePreload();post(workshop,{type:'preferences',data:preferences});return;}
  if(type==='race-ready'){track.dataset.ready='true';waiters.get('ready')?.resolve();waiters.delete('ready');post(track,{type:'preferences',data:preferences});post(track,{type:'screen-visibility',active:route==='track'});return;}
  if(type==='race-status'){
    const changed=raceState?.generation!==event.data.state.generation;
    if(changed)raceStarted=false;
    raceState=event.data.state;if(raceState.started||raceState.status==='running'||Number(raceState.frame?.step)>0)raceStarted=true;
    return;
  }
  if(type==='prepare'&&event.source===workshop?.contentWindow){
    if(locked()){post(workshop,{type:'prepare-result',id,error:'请先在赛道重置本场比赛。'});return;}
    if(preparing){post(workshop,{type:'prepare-result',id,error:'已有赛车正在准备。'});return;}
    preparing={id,timer:setTimeout(()=>{post(frames.get('track'),{type:'cancel-load'});finishPreparation('准备超时，请返回编辑后重试。');},600000)};
    try{const frame=await trackReady();if(preparing?.id!==id)return;post(workshop,{type:'prepare-progress',id,message:'准备赛车…'});post(frame,{type:'load-race',id,data});}
    catch(error){if(preparing?.id===id)finishPreparation(error.message);}return;
  }
  if(type==='cancel-prepare'&&event.source===workshop?.contentWindow){post(frames.get('track'),{type:'cancel-load'});finishPreparation('已返回编辑。');return;}
  if(type==='race-progress'&&event.source===track?.contentWindow&&preparing?.id===id){downloadProgress(event.data.percent,event.data.message);post(workshop,{type:'prepare-progress',id,message:event.data.message,percent:event.data.percent});}
  if(type==='race-loaded'&&event.source===track?.contentWindow&&preparing?.id===id)finishPreparation(event.data.error);
});
// A full navigation destroys the in-memory worker and trajectory. Internal
// screen changes preserve them and must remain free of leave-page prompts.
addEventListener('beforeunload',event=>{
  if(!preparing&&!raceState?.unsavedReplay&&!locked())return;
  event.preventDefault();event.returnValue='';
});
addEventListener('hashchange',render);
nav.innerHTML=nav.querySelector('.settings-link').outerHTML;
render();
schedulePreload();

function homePage() {
  return `<section class="welcome" aria-label="分子赛车首页"><div class="welcome-lockup"><svg class="welcome-mark" aria-hidden="true" color="#191b1a" xmlns="http://www.w3.org/2000/svg" viewBox="-21 -6 130 84" fill="none" stroke="#191b1a" stroke-width="2.6" stroke-linecap="round" stroke-linejoin="round">
  
  <path fill="#c9c9c6" fill-rule="evenodd" d="M12 37.32H32L42 54.64 32 71.96H12L2 54.64ZM11.3 54.64a10.7 10.7 0 1 0 21.4 0a10.7 10.7 0 1 0-21.4 0Z"/>
  <path fill="#c9c9c6" fill-rule="evenodd" d="M72 37.32H92L102 54.64 92 71.96H72L62 54.64ZM71.3 54.64a10.7 10.7 0 1 0 21.4 0a10.7 10.7 0 1 0-21.4 0Z"/>
  <path d="M32 37.32 42 20H62L72 37.32M42 54.64H62"/>
  <path d="M12 37.32 2 20M2 20-15.32 10 2 0Z"/>
  <circle cx="22" cy="54.64" r="10.7"/>
  <circle cx="82" cy="54.64" r="10.7"/>
</svg><h1 class="welcome-name"><span>RACING</span><span>MOLECULE</span><small><i>/</i>分子赛车</small></h1><a href="#garage" class="button button-primary welcome-join">加入赛事</a></div></section>`;
}

