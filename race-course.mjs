// New frames carry their post-relaxation origin, even in truncated replays.
// Older replays retain their original fixed-line rule; never infer a new origin
// from the first retained frame, which may be well into the race.
export function raceCourse(track, frame) {
  if (Number.isFinite(frame?.race_start_x) && Number.isFinite(frame?.race_length_angstrom)
      && frame.race_length_angstrom > 0) {
    return {start_x: frame.race_start_x, finish_x: frame.race_start_x + frame.race_length_angstrom,
      length_angstrom: frame.race_length_angstrom, rule: 'longitudinal-displacement-v1'};
  }
  return {start_x: track.start_x, finish_x: track.finish_x,
    length_angstrom: track.finish_x - track.start_x, rule: 'fixed-finish-line-v1'};
}

export function raceFinished(frame) {
  return frame.progress >= 1;
}
