// Scores use MD time and forward displacement, never animation wall time.
export function raceResult(data) {
  const frame=data?.frame;
  if(data?.status!=='finished'||!frame)return null;
  const distance=Number(frame.distance_angstrom), time=Number(frame.time_fs);
  if(!Number.isFinite(distance)||!Number.isFinite(time)||time<=0)return null;
  return {
    finished:data.finish_reason!=='time_limit' && frame.progress>=1,
    distanceNm:distance/10, timePs:time/1000,
    // 1 Å / fs = 100,000 m/s. Signed net forward displacement.
    averageSpeed:distance/time*100000,
    name:data.molecule?.name||'分子赛车',
    mode:data.competition?.mode==='custom'?'自定义赛事':'标准赛事',
  };
}

export class RaceResultDialog {
  constructor({again,save}) {
    this.generation=null;this.shown=false;this.result=null;this.finalSeq=null;
    this.dialog=document.createElement('dialog');
    this.dialog.id='race-result-dialog';this.dialog.className='race-result-dialog';
    this.dialog.setAttribute('aria-labelledby','race-result-title');
    this.dialog.innerHTML=`<div class="result-stripe" aria-hidden="true"></div>
      <button type="button" class="result-close" aria-label="关闭成绩">×</button>
      <div class="result-eyebrow"><svg viewBox="0 0 40 40" aria-hidden="true"><path d="M9 34V6m0 1c8-6 15 6 24 0v16c-9 6-16-6-24 0"/></svg><span id="race-result-mode"></span></div>
      <h2 id="race-result-title">冲线成功！</h2><p id="race-result-name"></p>
      <div class="result-hero"><span>平均速度</span><div><strong id="race-result-speed">—</strong><b>m/s</b></div></div>
      <div class="result-stats"><div><span>仿真用时</span><strong><span id="race-result-time">—</span><small>ps</small></strong></div><div><span>净前进距离</span><strong><span id="race-result-distance">—</span><small>nm</small></strong></div></div>
      <p id="race-result-note">沿赛道方向净位移 ÷ 仿真用时</p>
      <div class="result-actions"><button type="button" class="secondary-button" id="result-save">保存轨迹</button><button type="button" class="primary-button" id="result-again">再跑一次 <span aria-hidden="true">↗</span></button></div>`;
    document.body.append(this.dialog);
    this.reopen=document.createElement('button');this.reopen.id='show-race-result';
    this.reopen.type='button';this.reopen.className='secondary-button result-reopen';
    this.reopen.textContent='查看成绩';this.reopen.hidden=true;
    document.querySelector('.run-section').append(this.reopen);
    this.reopen.addEventListener('click',()=>this.open());
    this.dialog.querySelector('.result-close').addEventListener('click',()=>this.dialog.close());
    this.dialog.querySelector('#result-again').addEventListener('click',()=>{this.dialog.close();again();});
    this.dialog.querySelector('#result-save').addEventListener('click',()=>save());
  }
  reset(generation) {
    this.generation=generation;this.shown=false;this.result=null;this.finalSeq=null;
    this.reopen.hidden=true;if(this.dialog.open)this.dialog.close();
  }
  update(data) {
    if(this.generation!==data.generation)this.reset(data.generation);
    this.result=raceResult(data);this.finalSeq=data.frame?.seq;
    if(!this.result){this.reopen.hidden=true;return;}
    const r=this.result, put=(id,value)=>this.dialog.querySelector('#'+id).textContent=value;
    put('race-result-title',r.finished?'冲线成功！':'本轮结束');
    put('race-result-mode',r.mode);put('race-result-name',r.name);
    put('race-result-speed',r.averageSpeed.toFixed(1));
    put('race-result-time',r.timePs.toFixed(2));put('race-result-distance',r.distanceNm.toFixed(2));
    put('race-result-note',r.finished?'沿赛道方向净位移 ÷ 仿真用时':'达到仿真时长上限，尚未完赛。速度按本轮净前进距离计算。');
  }
  // Wait until the actual final coordinates have been displayed, not just received.
  present(renderedSeq) {
    if(!this.result||renderedSeq<this.finalSeq)return;
    this.reopen.hidden=false;
    if(!this.shown)this.open();
  }
  open() {
    if(!this.result||this.dialog.open)return;
    this.shown=true;this.dialog.showModal();
    this.dialog.querySelector('.result-close').focus({preventScroll:true});
  }
}
