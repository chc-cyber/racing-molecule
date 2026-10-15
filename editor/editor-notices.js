import {chemistryWarnings} from './model.js';

export const LEGACY_C60_ADDITION_NOTICE='C60 单键连接按加成处理：连接位点相邻的笼内双键已改为单键，另一笼碳将在生成三维时补氢。此操作可整体撤销。';
export const C60_ADDITION_NOTICE='已完成 C₆₀ 加成连接。\n连接点相邻的笼内双键已改为单键。';
const completedNotes=[C60_ADDITION_NOTICE,LEGACY_C60_ADDITION_NOTICE];
export const isC60AdditionNotice=text=>completedNotes.some(note=>String(text).includes(note));
export const stripCompletedC60Notice=text=>completedNotes.reduce((value,note)=>value.replaceAll(note,''),String(text)).trim();
// Older saved projects mixed this completed operation with unresolved stereo.
// Filter only that exact legacy note; do not hide real valence/stereo warnings.
export function activeChemistryWarnings(project){
  return [...new Set(chemistryWarnings(project).map(stripCompletedC60Notice).filter(Boolean))];
}
