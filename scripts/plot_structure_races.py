"""Summarize measured races and render scientific figures; no inferred finishes."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
for variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[variable]='1'
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/structure-race-scan'
POLAR=ROOT/'outputs/structure-race-polar'
sys.path.insert(0,str(ROOT/'.local/analysis-libs'))
os.environ['MPLCONFIGDIR']=str(OUT/'matplotlib-cache')
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt, font_manager
import numpy as np

NAMES=('benzene-12','coronene-36','ethynyl-rigid','folded-eicosane','test-car-52','fullerene-c60','pyridine-11','benzonitrile-13','acetone-10')
LABELS=('苯\nC6H6','蔻\nC24H12','炔桥\nC22H14','柔性链\nC20H42','四臂车\nC30H22','球笼\nC60','吡啶\nC5H5N','苯甲腈\nC7H5N','丙酮\nC3H6O')
COLORS=('#0072B2','#E69F00','#009E73','#CC79A7','#D55E00','#6A51A3','#2C7FB8','#7F8C3E','#A6611A')
font=Path('C:/Windows/Fonts/msyh.ttc')
if font.exists():font_manager.fontManager.addfont(str(font));plt.rcParams['font.family']=['Microsoft YaHei','DejaVu Sans']
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':10,'svg.fonttype':'none','axes.unicode_minus':False})


def save(fig,name):
    for extension in ('png','svg'):
        fig.savefig(OUT/'figures'/(name+'.'+extension),dpi=180,facecolor='white')
    plt.close(fig)


def groups(rows):
    result=[]
    for temp,strength in sorted({(r['config']['temperature_K'],r['config']['strength']) for r in rows}):
        structure=[]
        all_valid=True
        for name in NAMES:
            trials=sorted([r for r in rows if r['config']['molecule']==name and r['config']['temperature_K']==temp and r['config']['strength']==strength],key=lambda r:r['config']['seed'])
            valid=[r for r in trials if r['status']=='passed' and r['competitive_valid']]
            finished=[r['finish_time_ps'] for r in valid if r['outcome']=='finish']
            complete=len(finished)==len(trials) and len(trials)>=3
            all_valid &= complete
            structure.append({'molecule':name,'n':len(trials),'finished':len(finished),'complete':complete,
                'mean_ps':float(np.mean(finished)) if complete else None,
                'sd_ps':float(np.std(finished,ddof=1)) if complete else None,
                'observed_finish_range_ps':[min(finished),max(finished)] if finished else None,
                'censored':sum(r.get('right_censored',False) for r in valid),
                'failures':len(trials)-len(valid),'seeds':[r['config']['seed'] for r in trials]})
        row={'temperature_K':temp,'strength':strength,'structures':structure,'all_finished_valid':all_valid}
        if all_valid:
            means=np.array([s['mean_ps'] for s in structure]);counts=np.array([s['n'] for s in structure])
            grand=float(np.average(means,weights=counts));ssb=float(np.sum(counts*(means-grand)**2))
            ssw=sum((s['n']-1)*s['sd_ps']**2 for s in structure)
            noise=float(np.sqrt(ssw/(counts.sum()-len(counts))))
            fast=int(np.argmin(means));slow=int(np.argmax(means))
            row.update(slow_fast_mean_ratio=float(means[slow]/means[fast]),mean_gap_ps=float(np.ptp(means)),
                fastest=NAMES[fast],slowest=NAMES[slow],pooled_within_structure_sd_ps=noise,
                gap_to_pooled_noise_ratio=float(np.ptp(means)/noise) if noise else None,
                anova_eta_squared=float(ssb/(ssb+ssw)),mean_within_structure_cv=float(np.mean([s['sd_ps']/s['mean_ps'] for s in structure])))
        result.append(row)
    return result


def condition_rows(rows,c):
    return [r for r in rows if r['config']['temperature_K']==c['temperature_K'] and r['config']['strength']==c['strength']]


def main():
    parser=argparse.ArgumentParser(__doc__);parser.add_argument('--screen-only',action='store_true');args=parser.parse_args()
    (OUT/'figures').mkdir(exist_ok=True)
    screen=json.loads((OUT/'screen-report.json').read_text(encoding='utf8'))
    if screen['status']!='passed' or not screen['source_unchanged']:raise ValueError('Screen is not complete/accepted')
    polar=json.loads((POLAR/'report.json').read_text(encoding='utf8'))
    if polar['status']!='passed' or not polar['source_unchanged']:raise ValueError('Polar study is not complete/accepted')
    screen['trials'] += [r for r in polar['trials'] if r['config']['phase']=='screen']
    screen['fixtures'] += polar['fixtures']
    # This is a one-seed fixed-horizon selection diagnostic, not time variance.
    temps=(150.,300.,600.,900.);strengths=(.45,.60,.75)
    fig,axes=plt.subplots(3,3,figsize=(13.5,11.3),layout='constrained')
    fig.suptitle('初筛：同样 40 ps 后，不同分子前进了多远？',fontsize=16)
    for ax,name,label in zip(axes.flat,NAMES,LABELS):
        values=np.full((len(temps),len(strengths)),np.nan)
        for i,temp in enumerate(temps):
            for j,strength in enumerate(strengths):
                trial=next(r for r in screen['trials'] if r['config']['molecule']==name and r['config']['temperature_K']==temp and r['config']['strength']==strength)
                if trial['status']=='passed' and trial['competitive_valid']:values[i,j]=trial['final_forward_nm']
        image=ax.imshow(values,origin='lower',aspect='auto',vmin=0,vmax=16,cmap='viridis')
        ax.set_xticks(range(3),[f'{s:.2f}' for s in strengths]);ax.set_yticks(range(4),[str(int(t)) for t in temps])
        ax.set(xlabel='驱动强度',ylabel='温度 (K)',title=label.replace('\n',' · '))
        for (i,j),value in np.ndenumerate(values):
            ax.text(j,i,f'{value:.1f}' if np.isfinite(value) else '无效',ha='center',va='center',color='white' if value<10 else 'black',fontsize=10)
    fig.colorbar(image,ax=axes.ravel().tolist(),label='40 ps 前向质心位移 (nm)',shrink=.85)
    fig.supxlabel('每格仅一个随机种子：用于选条件，不用于证明完赛时间或结构排名。周期 Au；0.5 fs；350 fs；无额外热整备。',fontsize=10)
    save(fig,'01-screen-progress')
    if args.screen_only:return
    rows=[]
    for directory in (OUT,POLAR):
        for path in sorted((directory/'trials').glob('confirm-*.json')):
            row=json.loads(path.read_text(encoding='utf8'));row['report_file']=str(path.relative_to(ROOT))
            row['coordinates_path']=str((directory/'trials'/row['coordinates_file']).relative_to(ROOT));rows.append(row)
    if not rows:raise ValueError('No confirmation trials')
    unique={(r['config']['molecule'],r['config']['temperature_K'],r['config']['strength'],r['config']['seed']) for r in rows}
    if len(unique)!=len(rows):raise ValueError('Duplicate seed counted as replication')
    hashes={**screen['source_sha256'],**polar['source_sha256']}
    if any(hashes[k]!=v for r in rows for k,v in r['source_sha256'].items()):raise ValueError('Different frozen source versions')
    for name,expected in hashes.items():
        path=Path(rows[0]['runtime']['surface']['parameter_file']) if name=='CHAu.ILP' else ROOT/name
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise ValueError('Changed numerical source: '+name)
    summaries=groups(rows)
    baseline=next(s for s in summaries if s['temperature_K']==600 and s['strength']==.75)
    accepted=[s for s in summaries if s['all_finished_valid']]
    if not accepted:raise ValueError('No fully finishing repeated condition')
    best=max(accepted,key=lambda s:s['slow_fast_mean_ratio'])
    ordered=[baseline]+[s for s in summaries if s is not baseline]
    fig,axes=plt.subplots(2,2,figsize=(13,9.2),layout='constrained')
    fig.suptitle('20 nm 完赛用时：分子结构的差距与随机波动',fontsize=16)
    for ax,c in zip(axes.flat,ordered):
        trials=condition_rows(rows,c)
        for j,(name,color) in enumerate(zip(NAMES,COLORS)):
            subset=sorted((r for r in trials if r['config']['molecule']==name),key=lambda r:r['config']['seed'])
            for k,r in enumerate(subset):
                if r['status']!='passed':continue
                if not r['competitive_valid']:
                    ax.scatter(j+(k-(len(subset)-1)/2)*.065,r['observed_time_ps'],s=65,color='#D62728',marker='x',linewidths=2,zorder=5)
                    continue
                ax.scatter(j+(k-(len(subset)-1)/2)*.065,r['observed_time_ps'],s=24,color=color,
                    marker='o' if r['outcome']=='finish' else '^',facecolors=color if r['outcome']=='finish' else 'none',alpha=.65)
            s=c['structures'][j]
            if s['complete']:
                ax.errorbar(j,s['mean_ps'],yerr=s['sd_ps'],fmt='_',color='black',ms=17,lw=1.7,capsize=5,zorder=4)
        info=(f"慢 / 快均值 = {c['slow_fast_mean_ratio']:.2f} 倍；同结构 CV ≈ {100*c['mean_within_structure_cv']:.1f}%"
              if c['all_finished_valid'] else '含持续离面或未完赛样例，不给出完整均值排名')
        ax.set_title(f"{c['temperature_K']:g} K · 强度 {c['strength']:.2f}\n{info}",fontsize=11)
        ax.set_xticks(range(len(NAMES)),LABELS,fontsize=8);ax.set_ylabel('仿真完赛时间 (ps)');ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.2)
    for ax in axes.flat[len(ordered):]:ax.set_visible(False)
    fig.supxlabel('圆点 = 有效完赛；黑横线及误差棒 = 均值 ± 样本标准差；红叉 = 持续离面，空心三角 = 未完赛。各面板纵轴量程不同。',fontsize=9.5)
    save(fig,'02-finish-times')
    trials=condition_rows(rows,best)
    fig,axes=plt.subplots(1,2,figsize=(12.8,5.8),layout='constrained')
    fig.suptitle(f"结构差异最大且全部完赛的复测条件：{best['temperature_K']:g} K · 强度 {best['strength']:.2f}",fontsize=14)
    for name,label,color in zip(NAMES,LABELS,COLORS):
        for k,r in enumerate(sorted((r for r in trials if r['config']['molecule']==name),key=lambda r:r['config']['seed'])):
            sample=r['samples'][::10]
            if sample[-1]['step']!=r['samples'][-1]['step']:sample.append(r['samples'][-1])
            t=np.array([s['time_ps'] for s in sample]);p=np.array([s['center_nm'] for s in sample]);rel=p-p[0]
            axes[0].plot(t,p[:,0],color=color,alpha=.65,lw=1,label=label.split('\n')[0] if k==0 else None)
            axes[1].plot(rel[:,0],rel[:,1],color=color,alpha=.65,lw=1)
    axes[0].axhline(20,color='gray',ls=':',lw=.8);axes[0].set(xlabel='仿真时间 (ps)',ylabel='质心 x (nm)',ylim=(-.5,20.5));axes[0].legend(fontsize=9)
    axes[1].set(xlabel='前向质心位移 (nm)',ylabel='横向质心位移 (nm)');axes[1].axhline(0,color='gray',ls=':',lw=.8)
    for ax in axes:ax.grid(alpha=.2)
    fig.supxlabel('每条曲线来自一场真实比赛，终点停止；右图横纵坐标量程不同。温度、驱动以外的模型和时间步均未改变。',fontsize=9.5)
    save(fig,'03-selected-trajectories')
    fig=plt.figure(figsize=(13,10));fig.suptitle('九种代表分子的实际松弛后起跑结构',fontsize=15)
    manifest={f['id']:f for f in screen['fixtures']}
    for j,name in enumerate(NAMES):
        r=next(r for r in trials if r['config']['molecule']==name and r['status']=='passed')
        bundle=np.load(ROOT/r['coordinates_path']);p=bundle['positions_A'][0];p=p-np.average(p,axis=0,weights=bundle['masses_u']);bonds=bundle['bonds']
        heavy=bundle['masses_u']>2;ms=bundle['masses_u'][heavy]
        colors=['#2675bf' if np.isclose(m,14.007) else '#dd4e4e' if np.isclose(m,15.999) else '#37474f' for m in ms]
        ax=fig.add_subplot(3,3,j+1,projection='3d')
        for a,b in bonds:ax.plot(*p[[a,b]].T,color='#9ca3af',lw=.7)
        ax.scatter(*p[heavy].T,color=colors,s=28,depthshade=True)
        ax.scatter(*p[~heavy].T,color='#aaaaaa',s=8,depthshade=True)
        if 'initial_polar_response' in r:
            mu=np.array(r['initial_polar_response']['model_dipole_vector_e_A']);mu=mu/np.linalg.norm(mu)*4
            ax.quiver(0,0,0,*mu,color='black',arrow_length_ratio=.25)
        ax.set(xlim=(-10.5,10.5),ylim=(-10.5,10.5),zlim=(-10.5,10.5),title=f"{LABELS[j].split(chr(10))[0]} · {manifest[name]['formula']} · {manifest[name]['atoms']} 原子")
        ax.set_proj_type('ortho');ax.set_box_aspect((1,1,1),zoom=1.5);ax.view_init(elev=20,azim=-65);ax.set_axis_off()
    fig.subplots_adjust(left=0,right=1,top=.93,bottom=.07,wspace=0,hspace=.03)
    fig.text(.5,.025,'相同空间比例；C 深灰、N 蓝、O 红；极性组黑箭头为模型偶极方向（统一长度）。固定构象，其他构象可能改变成绩。',ha='center',fontsize=9.5)
    save(fig,'00-molecular-lineup')
    invalid=[r for r in rows if r['status']=='passed' and not r['competitive_valid']]
    if invalid:
        bad=invalid[0];c=bad['config']
        same=[r for r in rows if r['config']['molecule']==c['molecule'] and r['config']['temperature_K']==c['temperature_K'] and r['config']['strength']==c['strength']]
        fig,axes=plt.subplots(1,2,figsize=(12.6,5.2),layout='constrained')
        fig.suptitle(f"完整赛程检查：{LABELS[NAMES.index(c['molecule'])].split(chr(10))[0]} · {c['temperature_K']:g} K · 强度 {c['strength']:.2f}",fontsize=14)
        for k,r in enumerate(sorted(same,key=lambda r:r['config']['seed'])):
            t=np.array([s['time_ps'] for s in r['samples']]);z=np.array([s['min_z_A']/10 for s in r['samples']]);x=np.array([s['center_nm'][0] for s in r['samples']])
            color='#D62728' if not r['competitive_valid'] else COLORS[k]
            axes[0].plot(t,z,color=color,lw=1,label=f"种子 {r['config']['seed']}"+(' · 持续离面' if not r['competitive_valid'] else ''))
            axes[1].plot(t,x,color=color,lw=1)
        axes[0].axhline(.8,color='gray',ls='--',lw=1);axes[0].set(xlabel='仿真时间 (ps)',ylabel='最低原子高度 (nm)');axes[0].legend(fontsize=8)
        axes[1].set(xlabel='仿真时间 (ps)',ylabel='质心 x (nm)');axes[1].axhline(20,color='gray',ls=':',lw=.8)
        for ax in axes:ax.grid(alpha=.2)
        fig.supxlabel(f"红线：整车最低原子高于 0.8 nm，持续约 {bad['maximum_lift_duration_ps']:.1f} ps。坐标、能量和键长检查仍通过；到达终点不等于有效的贴面比赛。",fontsize=9.5)
        save(fig,'04-surface-contact-check')
    compact=[{k:v for k,v in r.items() if k not in ('samples','runtime','source_sha256','traceback')} for r in rows]
    report={'status':'passed' if all(r['status']=='passed' for r in rows) else 'completed-with-errors',
            'screen_trials':len(screen['trials']),'confirmation_trials':len(rows),'source_sha256':hashes,
            'protocol':screen['protocol'],'fixtures':screen['fixtures'],'conditions':summaries,
            'best_fully_finishing_condition':{k:v for k,v in best.items() if k!='structures'},'trials':compact,
            'qualification':'Screen used one separate seed. Confirmation uses independent seeds and fixed imported conformers; limited tested conditions, not a global optimum or an experimental nanocar ranking.',
            'valid_surface_confirmation_trials':sum(r['status']=='passed' and r['competitive_valid'] for r in rows),
            'invalid_surface_confirmation_trials':[r['key'] for r in rows if not r.get('competitive_valid',False)],
            'model_dipole_reference':'polarity-diagnostics.json; fixed Gasteiger-charge estimates, not measured dipoles',
            'plot_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT/'final-report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')
    with (OUT/'finish-summary.csv').open('w',newline='',encoding='utf-8-sig') as stream:
        fields=('temperature_K','strength','molecule','n','finished','mean_ps','sd_ps','censored','failures')
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        for c in summaries:
            for s in c['structures']:writer.writerow({k:c.get(k,s.get(k,'')) for k in fields})
    print(json.dumps({'report':str(OUT/'final-report.json'),'confirmation_trials':len(rows),'best':report['best_fully_finishing_condition']},ensure_ascii=False))


if __name__=='__main__':main()
