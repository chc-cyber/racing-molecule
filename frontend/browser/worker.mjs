import createMolecularCar from './runtime/molecular-car.mjs';
import {loadTable,loadJSON} from './assets.mjs';
import {raceFinished} from '../race-course.mjs';

let module, engine, metadata, settings, frameSteps, maxTimeFs, running = false, timer;
let playbackLimitFs = Infinity, computedTimeFs = 0;
const send = event => postMessage(event);

function message(error) {
  if (module?.getExceptionMessage) {
    try { const result = module.getExceptionMessage(error); module.decrementExceptionRefcount(error); return Array.isArray(result) ? result.join(': ') : result.message || JSON.stringify(result); } catch {}
  }
  return error?.message || String(error);
}

async function loadAssets(settings) {
  send({kind:'progress',message:'加载赛道资源…'});
  metadata=await loadJSON(new URL('./assets/model.json',import.meta.url));
  const info=settings.surface==='fourier-au'?await loadJSON(new URL('./assets/fourier-atlas.json',import.meta.url)):null;
  if(info&&info.parameter_sha256!==metadata.surface.parameter_sha256)throw new Error('赛道参数版本不一致，请刷新页面。');
  const arrays=metadata.surface.arrays,descriptors=[arrays.tile,arrays['saip-grid'],arrays['hetero-grid'],...(info?[info.arrays[0]]:[])];
  const counts=descriptors.map(()=>0),total=descriptors.reduce((n,d)=>n+d.download_bytes,0);let last=0;
  const runtime=createMolecularCar({print:()=>{},printErr:value=>console.warn(value)});
  const tables=descriptors.map((descriptor,index)=>loadTable(descriptor,{base:import.meta.url,onProgress:(bytes,cached)=>{
    counts[index]=bytes;const now=performance.now();if(now-last<200&&!cached&&bytes!==0)return;last=now;
    const loaded=counts.reduce((a,b)=>a+b,0);
    send({kind:'progress',message:`加载赛道资源… ${(loaded/1e6).toFixed(1)} / ${(total/1e6).toFixed(1)} MB`});
  }}));
  const [loadedModule,values]=await Promise.all([runtime,Promise.all(tables)]);module=loadedModule;
  send({kind:'progress',message:'初始化模拟…'});
  engine=new module.BrowserEngine(metadata,values[0],values[1],values[2]);
  if(info){engine.installFourier(info);engine.loadFourierSpecies(0,values[3]);}
}

function tick() {
  timer = null;
  if (!running || computedTimeFs >= playbackLimitFs) return;
  const started = performance.now();
  try {
    const frame = engine.advance(frameSteps);
    computedTimeFs = frame.time_fs;
    const duration = (performance.now() - started) / 1000;
    send({kind: 'frame', frame, stats: {steps_per_second: frameSteps / Math.max(duration, 1e-9)}, control: engine.controls()});
    if (raceFinished(frame) || frame.time_fs >= maxTimeFs) {
      running = false;
      send({kind: 'status', status: 'finished', reason: raceFinished(frame) ? 'finish' : 'time_limit'});
    }
    // Yield to controls after each 50 fs sample batch, without the former
    // 40 ms playback cap. Browser timer scheduling still adds wall-clock cost.
    if (running) timer = setTimeout(tick, 0);
  } catch (error) {
    running = false; send({kind: 'error', error: message(error)});
  }
}

onmessage = async ({data: event}) => {
  try {
    if (event.kind === 'load') {
      settings=event.settings;await loadAssets(settings);
      send({kind:'progress',message:'整理起跑构象…'});
      const info = engine.prepareWithPose(event.content, event.filename,
        {...settings, drive_strength:event.drive.strength, field_growth:event.drive.growth}, event.start_pose);
      // Keep sampling, controls and race limits in physical time across the
      // timesteps. Default 1 fs produces 50 steps per frame; preparation stays at 0.5 fs.
      frameSteps = metadata.frame_steps * metadata.timestep_fs / info.runtime.timestep_fs;
      if (!Number.isInteger(frameSteps) || frameSteps < 1 || frameSteps > 1000) throw new Error('时间步与轨迹采样间隔不兼容。');
      maxTimeFs = metadata.max_steps * metadata.timestep_fs;
      send({kind: 'ready', data: {molecule: info.molecule, runtime: info.runtime,
        frame: engine.frame(), control: engine.controls(), settings, track: metadata.track, start_pose: info.runtime.start_pose.requested}});
    } else if (event.kind === 'playback') {
      playbackLimitFs = event.until_fs;
      if (running && timer == null && computedTimeFs < playbackLimitFs) timer = setTimeout(tick, 0);
    } else if (event.kind === 'action') {
      if (event.action === 'pause') {running = false; clearTimeout(timer);send({kind: 'status', status: 'paused'});}
      else if (event.action === 'start') {if (!running) {running = true;send({kind: 'status', status: 'running'});timer = setTimeout(tick, 0);}}
    }
  } catch (error) {running = false;send({kind: 'error', error: message(error)});}
};
