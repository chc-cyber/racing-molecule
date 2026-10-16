import {LocalScores,createScore,scorePayload} from './local-scores.mjs';
export class LocalLeaderboard {
 constructor({replay,rerun,engine,onError}){
  this.store=new LocalScores();this.replay=replay;this.rerun=rerun;this.engine=engine;this.onError=onError;this.recorded=new Set();
  this.dialog=document.createElement('dialog');this.dialog.id='leaderboard-dialog';this.dialog.className='leaderboard-dialog';
  this.dialog.innerHTML='<header><h2>本地排行榜</h2><button type="button" aria-label="关闭排行榜">×</button></header><label>赛事版本<select id="score-ruleset"></select></label><p id="score-empty">完成标准赛事后，成绩会自动记录在这里。</p><ol id="score-list"></ol><section id="score-detail" hidden></section>';
  document.body.append(this.dialog);this.dialog.querySelector('header button').onclick=()=>this.dialog.close();this.dialog.querySelector('select').onchange=()=>this.render();
  const button=document.createElement('button');button.id='show-leaderboard';button.className='secondary-button';button.textContent='排行榜';button.onclick=()=>this.open();document.querySelector('.run-section').append(button);
  this.note=document.createElement('p');this.note.id='local-rank';document.querySelector('#race-result-note').after(this.note);
  const resultButton=document.createElement('button');resultButton.type='button';resultButton.className='secondary-button';resultButton.textContent='排行榜';resultButton.onclick=()=>this.open();document.querySelector('.result-actions').prepend(resultButton);
 }
 async update(state){
  if(state.status!=='finished'||this.recorded.has(state.generation))return;
  this.recorded.add(state.generation);const generation=state.generation;this.currentGeneration=generation;
  if(state.competition?.mode!=='standard'||state.finish_reason==='time_limit'||state.frame?.progress<1){this.note.textContent='';return;}
  try{const snapshot=structuredClone(state),replay=await this.replay();const record=await createScore(snapshot,replay,await this.engine());if(!record)return;await this.store.save(record);const rows=(await this.store.list()).filter(r=>r.ruleset===record.ruleset);if(this.currentGeneration===generation)this.note.textContent=`本地最佳排名 #${rows.findIndex(r=>r.structureHash===record.structureHash)+1}`;}
  catch(e){this.onError(e.message);}
 }
 reset(generation){if(this.currentGeneration!==generation){this.currentGeneration=generation;this.note.textContent='';}}
 async open(){try{this.rows=await this.store.list();const select=this.dialog.querySelector('select');const old=select.value;select.replaceChildren();for(const key of [...new Set(this.rows.map(r=>r.ruleset))]){const r=this.rows.find(r=>r.ruleset===key),o=document.createElement('option');o.value=key;o.textContent=`${r.competition.version} · ${r.engine.slice(0,12)}`;select.append(o);}if([...select.options].some(o=>o.value===old))select.value=old;this.render();if(!this.dialog.open)this.dialog.showModal();}catch(e){this.onError(e.message);}}
 render(){const list=this.dialog.querySelector('ol');list.replaceChildren();this.dialog.querySelector('#score-detail').hidden=true;const rows=this.rows.filter(r=>r.ruleset===this.dialog.querySelector('select').value);this.dialog.querySelector('#score-empty').hidden=!!rows.length;this.dialog.querySelector('label').hidden=!rows.length;rows.forEach((r,i)=>{const li=document.createElement('li'),button=document.createElement('button');button.className='score-row';const title=document.createElement('strong'),time=document.createElement('b'),sub=document.createElement('small');title.textContent=`${i+1}. ${r.name}`;time.textContent=`${(r.timeFs/1000).toFixed(2)} ps`;sub.textContent=`${r.formula} · ${new Date(r.createdAt).toLocaleDateString()}`;button.append(title,time,sub);button.onclick=()=>this.detail(r);li.append(button);list.append(li);});}
 detail(record){const panel=this.dialog.querySelector('#score-detail');panel.hidden=false;panel.replaceChildren();const title=document.createElement('h3');title.textContent=record.name;const p=document.createElement('p'),pose=record.startPose;p.textContent=`${record.formula} · ${record.molecule.atoms.length} 原子 · ${record.molecule.bonds.length} 键\n起跑朝向：${{auto:'自动摆放',preserve:'文件朝向',manual:'手动旋转'}[pose.mode]}${pose.mode==='manual'?`（${pose.yaw_deg}° / ${pose.pitch_deg}° / ${pose.roll_deg}°）`:''}`;panel.append(title,p);
  const download=(label,filename,content,type)=>{const b=document.createElement('button');b.textContent=label;b.className='secondary-button';b.onclick=()=>{const url=URL.createObjectURL(new Blob([content],{type})),a=document.createElement('a');a.href=url;a.download=filename;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};panel.append(b);};
  download('下载结构',record.input.filename,record.input.content,'chemical/x-mdl-molfile');download('导出成绩',record.name+'.score.json',JSON.stringify(scorePayload(record),null,2),'application/json');
  const again=document.createElement('button');again.textContent='再次试跑';again.className='primary-button';again.onclick=async()=>{try{await this.rerun(record);this.dialog.close();document.querySelector('#race-result-dialog').close();}catch(e){this.onError(e.message);}};panel.append(again);
 }
}
