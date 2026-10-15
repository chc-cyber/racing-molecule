import {validateProject} from './model.js';
const NativeError=globalThis.Clazz?._Error??globalThis.Error;
export const browserMode=document.documentElement.dataset.engine==='wasm';
let worker,serial=0;
const pending=new Map();
export async function browserRequest(path,data) {
  if(data?.project)validateProject(data.project);
  const operation=path.split('/').pop();
  if(operation==='import'&&!/\.(mol|sdf)$/i.test(data.filename||''))throw new NativeError('请导入 MOL 或单分子 SDF 文件。');
  if(!['import','export','geometry'].includes(operation))throw new NativeError('未知工坊操作。');
  if(!worker){
    worker=new Worker(new URL('../browser/editor-worker.mjs',import.meta.url),{type:'module',name:'Racing Molecule workshop'});
    worker.onmessage=({data})=>{const task=pending.get(data.id);if(!task)return;clearTimeout(task.timer);pending.delete(data.id);data.error?task.reject(new NativeError(data.error)):task.resolve(data.result);};
    worker.onerror=()=>resetWorker('构象计算中断，请重试。');
  }
  return new Promise((resolve,reject)=>{
    const id=++serial,timer=setTimeout(()=>resetWorker('构象计算超时，请检查结构后重试。'),90000);
    pending.set(id,{resolve,reject,timer});worker.postMessage({id,operation,payload:data});
  });
}
function resetWorker(message){worker?.terminate();worker=null;for(const task of pending.values()){clearTimeout(task.timer);task.reject(new NativeError(message));}pending.clear();}
addEventListener('pagehide',()=>resetWorker('页面已关闭。'));

// The parent owns the race iframe and worker, so crossing screens never loses
// a prepared molecule or silently substitutes a built-in example.
let loadSerial=0;
let cancelPending;
export function prepareBrowserRace(data,onProgress) {
  if(parent===window)throw new NativeError('请从首页进入工坊。');
  const id=crypto.randomUUID(),sequence=++loadSerial;
  return new Promise((resolve,reject)=>{
    cancelPending=()=>{cleanup();reject(new NativeError('已返回编辑。'));};
    const handler=event=>{
      if(event.source!==parent||event.origin!==location.origin||event.data?.channel!=='racing-molecule'||event.data.id!==id)return;
      if(event.data.type==='prepare-progress'){onProgress(event.data.message);return;}
      if(event.data.type!=='prepare-result')return;
      cleanup();if(sequence!==loadSerial)return;
      event.data.error?reject(new NativeError(event.data.error)):resolve(event.data);
    };
    const timer=setTimeout(()=>{cleanup();reject(new NativeError('准备超时，请返回编辑后重试。'));},610000);
    function cleanup(){clearTimeout(timer);removeEventListener('message',handler);cancelPending=null;}
    addEventListener('message',handler);parent.postMessage({channel:'racing-molecule',type:'prepare',id,data},location.origin);
  });
}
export function cancelBrowserRace(){++loadSerial;cancelPending?.();parent.postMessage({channel:'racing-molecule',type:'cancel-prepare'},location.origin);}
