import {raceResult} from './race-result.mjs';
import {validateCompetition} from './competition.mjs';
export const SCORE_SCHEMA='racing-molecule-score-v1';
export const SCORE_DB='racing-molecule-scores';
const hash=async value=>Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',new TextEncoder().encode(JSON.stringify(value)))),b=>b.toString(16).padStart(2,'0')).join('');
export function formula(atoms){const counts={};for(const a of atoms)counts[a.element]=(counts[a.element]||0)+1;return Object.keys(counts).sort((a,b)=>(a==='C'?-2:a==='H'?-1:0)-(b==='C'?-2:b==='H'?-1:0)||a.localeCompare(b)).map(e=>e+(counts[e]>1?counts[e]:'')).join('');}
export function structuralInput(content){const lines=content.replace(/\r/g,'').split('\n');const n=Number(lines[3]?.slice(0,3));if(lines[3]?.includes('V2000')&&n>0)return [...lines.slice(4,4+n).map(l=>l.slice(30)),...lines.slice(4+n)].join('\n');return lines.slice(3).join('\n');}
export async function createScore(state,replay,engine){
 const result=raceResult(state);if(!result?.finished||state.competition?.mode!=='standard')return null;
 const competition=validateCompetition(state.competition);
 if(!replay?.input?.content||!replay.initial_frame||!engine)throw new Error('本轮成绩缺少起始记录。');
 const {start_pose,...rules}=competition;
 const structureHash=await hash(structuralInput(replay.input.content));
 const ruleset=await hash({rules,engine});
 return {schema:SCORE_SCHEMA,id:crypto.randomUUID(),createdAt:new Date().toISOString(),verification:'local',ruleset,engine,structureHash,
  name:result.name,formula:formula(state.molecule.atoms),timeFs:state.frame.time_fs,distanceAngstrom:state.frame.distance_angstrom,
  competition,input:structuredClone(replay.input),molecule:structuredClone(state.molecule),initialFrame:structuredClone(replay.initial_frame),
  finalFrame:structuredClone(state.frame),runtime:structuredClone(state.runtime),startPose:structuredClone(state.start_pose)};
}
export function bestScores(records,score){const key=r=>r.ruleset+'|'+r.structureHash,old=records.find(r=>key(r)===key(score));if(old&&old.timeFs<=score.timeFs)return records;return [...records.filter(r=>key(r)!==key(score)),score];}
function database(){return new Promise((resolve,reject)=>{const r=indexedDB.open(SCORE_DB,1);r.onupgradeneeded=()=>r.result.createObjectStore('scores',{keyPath:'id'});r.onsuccess=()=>resolve(r.result);r.onerror=()=>reject(new Error('无法保存本地成绩，请检查浏览器存储空间。'));});}
export class LocalScores {
 async list(){const db=await database();return new Promise((resolve,reject)=>{const tx=db.transaction('scores'),r=tx.objectStore('scores').getAll();r.onsuccess=()=>resolve(r.result.sort((a,b)=>a.timeFs-b.timeFs));r.onerror=()=>reject(r.error);tx.oncomplete=()=>db.close();});}
 async save(score){const db=await database();return new Promise((resolve,reject)=>{const tx=db.transaction('scores','readwrite'),store=tx.objectStore('scores');store.getAll().onsuccess=event=>{const all=event.target.result,next=bestScores(all,score);if(next===all)return;for(const old of all)if(old.ruleset===score.ruleset&&old.structureHash===score.structureHash)store.delete(old.id);store.put(score);};tx.oncomplete=()=>{db.close();resolve(score);};tx.onerror=()=>{db.close();reject(new Error('本地成绩保存失败，请检查浏览器存储空间。'));};});}
}
// Future upload adapters consume this versioned payload; no network submission here.
export function scorePayload(record){if(record.schema!==SCORE_SCHEMA)throw new Error('成绩格式不支持。');return structuredClone(record);}
