// One saved snapshot, one candidate, one promotion. Retry recovers the same job.
const NativeError=globalThis.Clazz?._Error??globalThis.Error;
export class RaceCompletion {
  constructor(deps) {
    this.deps = deps;
    this.state = 'idle';
    this.job = null;
    this.cleanupJob = null;
    this.cleanupTask = null;
    this.serial = 0;
    this.inFlight = false;
  }
  get locked() { return this.state !== 'idle'; }
  update(state, message, {retry = false} = {}) {
    this.state = state;
    this.deps.notify({state, message, retry, locked: this.locked});
  }
  remember() { this.deps.remember(this.cleanupJob || this.job); }
  async cleanup() {
    if (!this.cleanupJob) return;
    if (this.cleanupTask) return this.cleanupTask;
    const pending = this.cleanupJob;
    this.cleanupTask = (async()=>{
      await this.deps.cancel(pending.id);
      if (this.cleanupJob === pending) { this.cleanupJob = null; this.remember(); }
    })();
    try { await this.cleanupTask; }
    finally { this.cleanupTask = null; }
  }
  async start() {
    if (this.locked || this.inFlight) return;
    const competition = structuredClone(this.deps.competition?.());
    const serial = ++this.serial;
    this.inFlight = true;
    this.update('saving', '保存工程…');
    try {
      const project = await this.deps.settle();
      if (serial !== this.serial) return;
      await this.deps.save(project);
      if (serial !== this.serial) return;
      this.update('checking', '检查结构…');
      const exported = await this.deps.export(project);
      if (serial !== this.serial) return;
      await this.cleanup();
      if (serial !== this.serial) return;
      this.job = {id: this.deps.uuid(), filename: exported.filename, content: exported.content,
                  design: JSON.stringify(project), ...(competition ? {competition} : {})};
      this.remember();
      await this.run(serial, false);
    } catch (error) {
      if (serial === this.serial) this.failure(error);
    } finally { if (serial === this.serial) this.inFlight = false; }
  }
  failure(error) {
    this.update('failed', error.message || String(error), {retry: Boolean(error.recoverable && this.job)});
  }
  async run(serial, recover) {
    let result;
    if (recover) {
      try { result = await this.deps.status(this.job.id); }
      catch (error) {
        if (!error.missing) throw error;
        result = await this.deps.prepare(this.job);
      }
    } else result = await this.deps.prepare(this.job);
    while (serial === this.serial) {
      if (result.status === 'error') throw new NativeError(result.error || '准备失败，请调整结构后重试。');
      if (result.status === 'cancelled') throw new NativeError('这次准备已取消。');
      if (result.status === 'ready') {
        this.update('committing', '进入赛道…');
        result = await this.deps.commit(this.job.id);
        if (serial !== this.serial) return;
      }
      if (result.status === 'committed') {
        const active = await this.deps.active();
        if (serial !== this.serial) return;
        if (active.generation !== result.generation || !active.molecule || !active.frame ||
            !['paused', 'running', 'finished'].includes(active.status)) {
          throw new NativeError(active.error || '赛道已切换，请返回编辑。');
        }
        this.job = null;
        this.remember();
        this.update('success', '赛车已就绪。');
        this.deps.navigate(result.generation);
        return;
      }
      if (!['checking', 'preparing'].includes(result.status)) throw new NativeError('准备状态无效，请返回编辑。');
      this.update(result.status, result.status === 'checking' ? '检查结构…' : '准备赛车…');
      await this.deps.delay();
      if (serial !== this.serial) return;
      result = await this.deps.status(this.job.id);
    }
  }
  async retry() {
    if (!this.job || this.inFlight || this.state !== 'failed') return;
    const serial = ++this.serial;
    this.inFlight = true;
    this.update('checking', '恢复准备…');
    try { await this.run(serial, true); }
    catch (error) { if (serial === this.serial) this.failure(error); }
    finally { if (serial === this.serial) this.inFlight = false; }
  }
  async restore(job) {
    if (this.locked || this.inFlight) return;
    if (job.cancelRequested) {
      this.cleanupJob = job;
      try { await this.cleanup(); } catch { /* Retry cleanup before the next preparation. */ }
      return;
    }
    this.job = job;
    this.update('failed', '恢复准备…', {retry: true});
    return this.retry();
  }
  async cancel() {
    if (this.state === 'committing' || this.state === 'success' || this.state === 'cancelling' || !this.locked) return;
    ++this.serial;
    if (this.job) this.cleanupJob = {...this.job, cancelRequested: true};
    this.job = null;
    this.remember();
    this.inFlight = false;
    this.update('idle', '已返回编辑。');
    // Local cancellation fences all late responses immediately. An unreachable
    // server must not trap the user in the dialog; cleanup retries on reload or
    // before the next preparation, and the server also expires orphaned jobs.
    try { await this.cleanup(); } catch { /* Keep the cancellation record for retry. */ }
  }
}
