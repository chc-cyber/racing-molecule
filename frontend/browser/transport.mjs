// Per-tab state; /api-shaped messages never leave the browser.
import {validateStartPose} from '../start-pose.js';
import {validateCompetition,engineParameters} from '../competition.mjs';
export class BrowserTransport {
  constructor() {
    this.competition=validateCompetition();
    this.generation = 0;
    this.frames = [];
    this.settings = {surface: 'fourier-au', temperature: 150, damping_fs: 350, corrugation: .08, lattice_scale: 1.5, timestep_fs: 1};
    this.drive = {mode: 'gradient', strength: .6, growth: 1};
    this.startPose = validateStartPose();
    this.state = {generation: 0, status: 'idle', frame: null, molecule: null,
      settings: {...this.settings}, start_pose: structuredClone(this.startPose), competition:structuredClone(this.competition), control: {drive: {...this.drive}},
      track: {start_x: 0, finish_x: 200, half_width: 100, surface_z: 0}};
    this.worker = null;
    this.droppedFrames = 0;
    this.initial = Promise.resolve();
    document.addEventListener('visibilitychange', () => {
      if (document.hidden && this.state.status === 'running') this.worker?.postMessage({kind: 'action', action: 'pause'});
    });
    addEventListener('pagehide', () => this.worker?.terminate());
  }

  fail(error) {this.worker?.terminate();this.state = {...this.state, status: 'error', error: error.message || String(error)};}

  appendFrame(frame) {
    this.frames.push(frame);
    // Bounded replay and polling cost for weak-drive, long-running races.
    const limit = Math.max(500, Math.floor(24_000_000 / Math.max(1, frame.positions.length * 8 + 512)));
    if (this.frames.length > limit) this.droppedFrames += this.frames.splice(0, this.frames.length - limit).length;
  }

  async load(data, keepPose = false) {
    const competition=validateCompetition(keepPose?this.competition:data.competition);
    if(data.start_pose)competition.start_pose=validateStartPose(data.start_pose);
    validateCompetition(competition);
    const pose=competition.start_pose;
    if(data.reuse){if(!this.input)throw new Error('请先选择分子。');data={...this.input,competition};}
    this.competition=competition;
    const parameters=engineParameters(competition);this.settings=parameters.settings;this.drive=parameters.drive;
    const generation = ++this.generation;
    try {
    this.worker?.terminate(); this.worker = null;
    this.frames = []; this.droppedFrames = 0;
    this.startPose = structuredClone(pose);
    this.state = {...this.state, generation, status: 'preparing', error: null, preparation:'加载赛道资源…', frame: null, molecule: null,
      runtime: null, stats: null, finish_reason: null, competition:structuredClone(this.competition), settings: {...this.settings}, start_pose: structuredClone(this.startPose)};
    let content = data.content, filename = data.filename;
    if (data.example) {
      if (!['tetraphenylbenzene.mol', 'benzene.mol'].includes(data.example)) throw new Error('未知的示例分子。');
      const response = await fetch(new URL(`./assets/${data.example}`, import.meta.url));
      if (!response.ok) throw new Error('示例结构下载失败。');
      content = await response.text(); filename = data.example;
    }
    if (generation !== this.generation) return this.state;
    if (typeof content !== 'string' || !content.trim()) throw new Error('分子文件不能为空。');
    if (typeof filename !== 'string' || !/\.(mol|sdf)$/i.test(filename)) throw new Error('请导入 MOL 或单记录 SDF 文件。');
    if (new TextEncoder().encode(content).length > 1_000_000) throw new Error('文件超过 1 MB。');
    this.input = {content, filename};
    const worker = new Worker(new URL('./worker.mjs', import.meta.url), {type: 'module', name: 'Racing Molecule MD'});
    this.worker = worker;
    worker.onmessage = ({data: event}) => {
      if (generation !== this.generation) return;
      if(event.kind==='progress'){this.state={...this.state,preparation:event.message,preparationPercent:event.percent};
      } else if (event.kind === 'ready') {
        this.state = {...this.state, ...event.data, status: 'paused'};
        this.appendFrame(event.data.frame);
      } else if (event.kind === 'frame') {
        this.state = {...this.state, frame: event.frame, stats: event.stats, control: event.control};
        this.appendFrame(event.frame);
      } else if (event.kind === 'control') {
        this.state = {...this.state, control: event.control, settings: event.settings};
      } else if (event.kind === 'status') {
        this.state = {...this.state, status: event.status, finish_reason: event.reason || null};
      } else if (event.kind === 'error') this.fail(new Error(event.error));
    };
    worker.onerror = event => {if (generation === this.generation) this.fail(new Error(event.message || '浏览器计算模块启动失败。'));};
    worker.postMessage({kind: 'load', content, filename, settings: {...this.settings}, drive: {...this.drive}, start_pose: structuredClone(this.startPose)});
    return this.state;
    } catch (error) {
      if (generation === this.generation) this.fail(error);
      throw error;
    }
  }

  async request(path, data) {
    if(path==='/api/cancel-load'){
      ++this.generation;this.worker?.terminate();this.worker=null;this.frames=[];
      this.state={...this.state,generation:this.generation,status:'idle',frame:null,molecule:null,error:null};
      return this.state;
    }
    if (path.startsWith('/api/state')) return this.state;
    if (path.startsWith('/api/frames')) {
      const query = new URL(path, 'https://browser.invalid').searchParams;
      const after = Number(query.get('after') ?? -1);
      const played = Number(query.get('played'));
      if (query.has('played') && Number.isFinite(played) && played >= 0 && Number(query.get('generation')) === this.generation)
        this.worker?.postMessage({kind:'playback', until_fs:played + 10000});
      const frames = this.frames.filter(f => f.seq > after).slice(0,250);
      return {generation: this.generation, frames, dropped: frames.length > 0 && frames[0].seq > after + 1};
    }
    if (path === '/api/load') return this.load(data);
    if (path === '/api/start-pose') {
      const pose = validateStartPose(data, this.state.molecule?.atoms.length);
      if (!this.input || ['preparing', 'running', 'idle'].includes(this.state.status))
        throw new Error('请等待准备完成或暂停比赛，再修改起跑姿态。');
      this.competition=validateCompetition({...this.competition,start_pose:pose});
      this.startPose = pose;
      return this.load(this.input, true);
    }
    if (path === '/api/action') {
      if (data.action === 'reset') return this.load(this.input, true);
      if (!['start', 'pause'].includes(data.action)) throw new Error('未知的比赛操作。');
      if (!['paused', 'running'].includes(this.state.status)) throw new Error('分子尚未准备就绪。');
      this.worker.postMessage({kind: 'action', action: data.action});
      this.state.status = data.action === 'start' ? 'running' : 'paused';
      return this.state;
    }
    if(path==='/api/control')throw new Error('赛事参数在准备时固定，请选择赛事配置后重新准备。');
    throw new Error('浏览器计算接口不存在。');
  }

  async replay() {
    return {schema: 'molecular-car-browser-replay-v1', execution: 'browser-wasm',
      competition:structuredClone(this.competition), molecule: this.state.molecule, runtime: this.state.runtime, settings: this.settings,
      track: this.state.track, input: this.input, frames: this.frames, generation: this.generation,
      start_pose: structuredClone(this.startPose),
      finish_reason: this.state.finish_reason,
      replay_complete: this.droppedFrames === 0, dropped_frames: this.droppedFrames,
      first_retained_step: this.frames[0]?.step ?? null};
  }
}
