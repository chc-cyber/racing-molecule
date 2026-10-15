// Presentation time only. MD steps, samples and scoring remain untouched.
export class PlaybackClock {
  constructor() { this.rate = 1; this.reset(); }
  reset(timeFs = 0) { this.timeFs = timeFs; this.playing = false; this.lastWall = null; this.buffering = false; }
  setRate(rate) {
    if (rate !== 1 && rate !== 2) throw new RangeError('Playback rate must be 1 or 2');
    this.rate = rate;
    this.lastWall = null;
  }
  advance(now, availableFs, visible = true) {
    const elapsed = this.lastWall === null ? 0 : Math.max(0, now - this.lastWall);
    this.lastWall = now;
    this.buffering = false;
    if (!this.playing || !visible) return this.timeFs;
    // A suspended tab must not jump forward on its first visible frame.
    const desired = this.timeFs + Math.min(elapsed, 250) * 2 * this.rate;
    this.timeFs = Math.min(desired, Math.max(this.timeFs, availableFs));
    this.buffering = desired > availableFs;
    return this.timeFs;
  }
}
