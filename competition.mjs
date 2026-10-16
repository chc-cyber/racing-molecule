import {validateStartPose} from './start-pose.js';
export const DEFAULT_COMPETITION = Object.freeze({mode:'standard', version:'race-v1', length_nm:20,
  temperature:150, strength:.6, surface:'fourier-au', lattice_scale:1.5, corrugation:.08,
  growth:1, timestep_fs:1, damping_fs:350, random_seed:98245165, start_pose:validateStartPose()});
export function validateCompetition(data = {}) {
  if (!data || typeof data !== 'object' || Array.isArray(data) || Object.keys(data).some(k=>!Object.hasOwn(DEFAULT_COMPETITION,k))) throw new Error('赛事配置包含未知参数。');
  const result = {...structuredClone(DEFAULT_COMPETITION), ...structuredClone(data)};
  if (!['standard','custom'].includes(result.mode)) throw new Error('请选择标准赛事或自定义赛事。');
  for (const key of ['version','length_nm','timestep_fs','damping_fs','random_seed']) if(result[key] !== DEFAULT_COMPETITION[key]) throw new Error('当前赛事版本固定为 20 nm、1 fs 步长与 350 fs 阻尼。');
  if (!['fourier-au','periodic-au','continuum'].includes(result.surface)) throw new Error('未知表面模型。');
  for (const [key,low,high] of [['temperature',0,1000],['strength',0,1],['lattice_scale',0,3],['corrugation',0,.6],['growth',.1,2]])
    if(typeof result[key]!=='number'||!Number.isFinite(result[key])||result[key]<low||result[key]>high) throw new Error(`${key} 必须在 ${low}–${high} 范围内。`);
  result.start_pose = validateStartPose(result.start_pose);
  if(result.mode==='standard' && Object.keys(DEFAULT_COMPETITION).filter(k=>k!=='start_pose').some(k=>JSON.stringify(result[k])!==JSON.stringify(DEFAULT_COMPETITION[k]))) throw new Error('标准赛事参数固定；请切换到自定义赛事。');
  return result;
}
export function engineParameters(c) {
  return {settings:Object.fromEntries(['surface','temperature','damping_fs','corrugation','lattice_scale','timestep_fs'].map(k=>[k,c[k]])), drive:{mode:'gradient',strength:c.strength,growth:c.growth}};
}
export function describeCompetition(c) {
  if (c.mode === 'standard') return `标准赛事 · ${c.length_nm} nm`;
  const surface = {'fourier-au':'代理势 v2','periodic-au':'周期势',continuum:'连续表面'}[c.surface];
  return `${c.mode==='standard'?'标准赛事':'自定义赛事'} · ${c.length_nm} nm · ${c.temperature} K · 驱动 ${c.strength} · ${surface} · ${c.surface==='continuum'?'起伏 '+c.corrugation:'晶格起伏 '+c.lattice_scale+' 倍'} · 沿程增长 ${c.growth} · ${c.timestep_fs} fs`;
}

// A paused race remains active until it ends or is reset to a new generation.
export function isRaceInProgress(state, started = false) {
  return state?.status === 'running' || state?.status === 'paused' && (started || Number(state.frame?.step) > 0);
}

// Same preparation form in the workshop, native race page and static WASM site.
export class RaceSetup {
  constructor(host, onChange=()=>{}) {
    this.host=host;this.onChange=onChange;
    host.innerHTML=`<div class="race-mode-options" role="group" aria-label="赛事类型"><button type="button" data-mode="standard" aria-pressed="true">标准赛事</button><button type="button" data-mode="custom" aria-pressed="false">自定义赛事</button></div>
      <div class="race-custom-fields" hidden>
        <label>环境温度 / K<input data-key="temperature" type="number" min="0" max="1000" step="1" value="150"></label>
        <label>驱动强度<input data-key="strength" type="number" min="0" max="1" step="0.005" value="0.6"></label>
        <label>表面模型<select data-key="surface"><option value="fourier-au">金表面 · 代理势 v2</option><option value="periodic-au">金表面 · 周期势</option><option value="continuum">连续表面</option></select></label>
        <label data-for="lattice_scale">晶格起伏倍率<input data-key="lattice_scale" type="number" min="0" max="3" step="0.25" value="1.5"></label>
        <label data-for="corrugation">连续表面起伏<input data-key="corrugation" type="number" min="0" max="0.6" step="0.01" value="0.08"></label>
        <label>电场沿程增长<input data-key="growth" type="number" min="0.1" max="2" step="0.05" value="1"></label>
        <label>起跑摆放<select data-key="pose_mode"><option value="auto">自动摆放</option><option value="preserve">采用文件朝向</option><option value="manual">手动旋转</option></select></label>
        <div class="race-pose-angles" hidden>${['yaw','pitch','roll'].map((key,i)=>`<label>${['平面旋转','俯仰','侧倾'][i]} / °<input data-key="${key}_deg" type="number" min="-180" max="180" value="0"></label>`).join('')}</div>
        <p>20 nm · 1 fs · 350 fs 阻尼</p>
      </div>`;
    this.mode='standard';this.pose=validateStartPose();this.disabled=false;
    host.addEventListener('click', e=>{const button=e.target.closest('[data-mode]');if(button && !this.disabled && !button.disabled){this.mode=button.dataset.mode;this.refresh();onChange();}});
    host.addEventListener('input',()=>{if(!this.disabled){this.refresh();onChange();}});
    host.addEventListener('change',()=>{if(!this.disabled){this.refresh();onChange();}});this.refresh();
  }
  refresh() {
    const custom=this.mode==='custom';
    this.host.querySelectorAll('[data-mode]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.mode===this.mode)));
    this.host.querySelector('.race-custom-fields').hidden=!custom;
    const continuum=this.host.querySelector('[data-key="surface"]').value==='continuum';
    this.host.querySelector('[data-for="lattice_scale"]').hidden=continuum;
    this.host.querySelector('[data-for="corrugation"]').hidden=!continuum;
    this.host.querySelector('.race-pose-angles').hidden=this.host.querySelector('[data-key="pose_mode"]').value!=='manual';
  }
  read() {
    if(this.mode==='standard')return validateCompetition({start_pose:this.pose});
    const result={mode:'custom'};
    for(const key of ['temperature','strength','surface','lattice_scale','corrugation','growth']){const v=this.host.querySelector(`[data-key="${key}"]`).value;result[key]=key==='surface'?v:(v.trim()?Number(v):NaN);}
    result.start_pose={...this.pose, mode:this.host.querySelector('[data-key="pose_mode"]').value};
    for(const key of ['yaw_deg','pitch_deg','roll_deg'])result.start_pose[key]=Number(this.host.querySelector(`[data-key="${key}"]`).value);
    return validateCompetition(result);
  }
  set(value) {
    const c=validateCompetition(value);this.mode=c.mode;this.pose=c.start_pose;
    for(const key of ['temperature','strength','surface','lattice_scale','corrugation','growth'])this.host.querySelector(`[data-key="${key}"]`).value=c[key];
    this.host.querySelector('[data-key="pose_mode"]').value=c.start_pose.mode;
    for(const key of ['yaw_deg','pitch_deg','roll_deg'])this.host.querySelector(`[data-key="${key}"]`).value=c.start_pose[key];
    this.refresh();
  }
  disable(value){this.disabled=Boolean(value);this.host.querySelectorAll('input,select,button').forEach(e=>e.disabled=this.disabled);}
}
