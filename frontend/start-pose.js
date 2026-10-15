// Start placement is a preparation request, never a live force/control.
export const DEFAULT_START_POSE = Object.freeze({mode: 'auto', yaw_deg: 0, pitch_deg: 0, roll_deg: 0, axis_atoms: []});
export function validateStartPose(value = DEFAULT_START_POSE, atomCount) {
  if (!value || typeof value !== 'object' || Array.isArray(value) ||
      Object.keys(value).some(key => !Object.hasOwn(DEFAULT_START_POSE, key))) throw new Error('起跑姿态参数格式不正确。');
  const result = {...DEFAULT_START_POSE, ...value};
  if (!['auto', 'manual', 'preserve'].includes(result.mode)) throw new Error('未知起跑摆放方式。');
  for (const key of ['yaw_deg', 'pitch_deg', 'roll_deg'])
    if (typeof result[key] !== 'number' || !Number.isFinite(result[key]) || result[key] < -180 || result[key] > 180)
      throw new Error('起跑旋转角度必须在 −180° 到 180° 之间。');
  const ids = result.axis_atoms;
  if (!Array.isArray(ids) || ![0, 2].includes(ids.length) || ids.some(id => !Number.isInteger(id) || id < 1 || id > 256) ||
      ids.length === 2 && ids[0] === ids[1]) throw new Error('请选择两个不同的车尾、车头参考原子。');
  if (atomCount !== undefined && ids.some(id => id > atomCount)) throw new Error('车头参考原子不属于当前分子。');
  if (result.mode !== 'manual') for (const key of ['yaw_deg', 'pitch_deg', 'roll_deg']) result[key] = 0;
  result.axis_atoms = result.mode === 'auto' ? [...ids] : [];
  return result;
}
