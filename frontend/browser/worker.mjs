import createMolecularCar from './runtime/molecular-car.mjs';
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

async function loadAssets() {
  const response = await fetch(new URL('./assets/model.json', import.meta.url));
  if (!response.ok) throw new Error('模型参数下载失败。');
  metadata = await response.json();
  const load = async descriptor => {
    const response = await fetch(new URL(`./${descriptor.path}`, import.meta.url));
    if (!response.ok) throw new Error('金表面势表下载失败。');
    const bytes = await new Response(response.body.pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
    if (bytes.byteLength !== descriptor.bytes) throw new Error('势表大小不正确。');
    const hash = await crypto.subtle.digest('SHA-256', bytes);
    const hex = Array.from(new Uint8Array(hash), n => n.toString(16).padStart(2, '0')).join('');
    if (hex !== descriptor.sha256) throw new Error('势表校验失败，请重新加载网页。');
    return new Float64Array(bytes);
  };
  const arrays = metadata.surface.arrays;
  const [tile, saip, hetero] = await Promise.all([load(arrays.tile), load(arrays['saip-grid']), load(arrays['hetero-grid'])]);
  module = await createMolecularCar({print: () => {}, printErr: value => console.warn(value)});
  engine = new module.BrowserEngine(metadata, tile, saip, hetero);
}

async function loadFourier() {
  const response=await fetch(new URL('./assets/fourier-atlas.json',import.meta.url));
  if(!response.ok)throw new Error('傅里叶参数下载失败。');
  const info=await response.json();
  if(info.parameter_sha256!==metadata.surface.parameter_sha256)throw new Error('傅里叶与周期势参数版本不一致。');
  engine.installFourier(info);
  // Ordinary bonded H has no direction term; uncommon directional H falls back
  // to the complete lattice sum. Download only the carbon coefficients.
  const descriptor=info.arrays[0];
  const data=await fetch(new URL('./'+descriptor.path,import.meta.url));
  if(!data.ok)throw new Error('傅里叶势表下载失败。');
  const bytes=await new Response(data.body.pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
  const hash=await crypto.subtle.digest('SHA-256',bytes);
  const hex=Array.from(new Uint8Array(hash),n=>n.toString(16).padStart(2,'0')).join('');
  if(bytes.byteLength!==descriptor.bytes||hex!==descriptor.sha256)throw new Error('傅里叶势表校验失败。');
  engine.loadFourierSpecies(0,new Float64Array(bytes));
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
      await loadAssets();settings = event.settings;
      if(settings.surface==='fourier-au')await loadFourier();
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
