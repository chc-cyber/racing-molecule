"""Honest finite-time race summaries and scientific figures for seven cars."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from scan_crystal_contrast import CARS
from scan_structure_races import write_json
LABELS=('苯','蔻','炔连接结构','四苯基苯','C60','吡啶','苯甲腈')
CONDITIONS=((1.,1.),(1.,.25),(1.5,1.),(1.5,.25),(2.,1.),(2.,.25),(3.,1.))


def panel_shape(count, columns=3):
    """Allocate an axis for every condition, including future added scans."""
    if count < 1 or columns < 1:
        raise ValueError('Panel count and columns must be positive')
    columns=min(columns,count)
    return ((count+columns-1)//columns,columns)


def trial_identity(config):
    return tuple(config[k] for k in ('phase','molecule','temperature_K','strength','lattice_scale','growth','seed','limit_ps'))


def condition_identity(config):
    return (config['temperature_K'],config['strength'],config.get('lattice_scale',1.),config.get('growth',1.))


def summarize(rows, planned=None):
    planned=[r['config'] for r in rows] if planned is None else planned
    configs=sorted({condition_identity(c) for c in planned},
        key=lambda c:(c[:2],CONDITIONS.index(c[2:]) if c[2:] in CONDITIONS else len(CONDITIONS),c[2:]))
    summaries=[]
    for temp,strength,scale,growth in configs:
        chosen=[r for r in rows if condition_identity(r['config'])==(temp,strength,scale,growth)]
        cars=[]
        for car in CARS:
            trials=[r for r in chosen if r['config']['molecule']==car]
            expected=sum(condition_identity(c)==(temp,strength,scale,growth) and c['molecule']==car for c in planned)
            valid=[r for r in trials if r['status']=='passed' and r.get('competitive_valid',False)]
            finished=[r['finish_time_ps'] for r in valid if r['outcome']=='finish']
            complete=len(finished)==expected and expected>0
            cars.append({'molecule':car,'planned':expected,'observed':len(trials),'missing':expected-len(trials),'valid':len(valid),'finished':len(finished),
                         'failed':sum(r['status']!='passed' for r in trials),
                         'lifted':sum(r.get('sustained_lift',False) for r in trials),
                         'censored':sum(r['status']=='passed' and r['outcome']!='finish' for r in trials),
                         'finish_mean_ps':float(np.mean(finished)) if complete else None,
                         'finish_sd_ps':float(np.std(finished,ddof=1)) if complete and len(finished)>1 else None,
                         'mean_forward_nm':float(np.mean([r['final_forward_nm'] for r in valid])) if valid else None})
        complete=all(c['finished']==c['planned'] and c['planned']>1 for c in cars)
        group={'temperature_K':temp,'strength':strength,'lattice_scale':scale,'growth':growth,'cars':cars,'all_finished_valid':complete}
        if complete:
            means=np.array([c['finish_mean_ps'] for c in cars]);sds=np.array([c['finish_sd_ps'] for c in cars]);n=np.array([c['planned'] for c in cars])
            noise=float(np.sqrt(np.sum((n-1)*sds*sds)/np.sum(n-1)));sorted_mean=np.sort(means)
            effects=[abs(means[i]-means[j])/np.sqrt((sds[i]**2+sds[j]**2)/2) for i in range(len(cars)) for j in range(i+1,len(cars))]
            group.update(slow_fast_ratio=float(means.max()/means.min()),mean_gap_ps=float(np.ptp(means)),
                         pooled_within_sd_ps=noise,gap_to_noise=float(np.ptp(means)/noise),
                         trimmed_slow_fast_ratio=float(sorted_mean[-2]/sorted_mean[1]),
                         mean_IQR_to_noise=float((np.percentile(means,75)-np.percentile(means,25))/noise),
                         pairs_with_effect_size_over_2=int(np.sum(np.array(effects)>2)),total_pairs=len(effects),
                         mean_within_cv_percent=float(np.mean(sds/means)*100.))
        summaries.append(group)
    return summaries


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'outputs/crystal-contrast-guarded');args=p.parse_args();out=args.output
    browser=[json.loads(f.read_text(encoding='utf8')) for f in sorted((out/'browser-trials').glob('*.json'))]
    jobs=json.loads((out/'browser-confirm-jobs.json').read_text(encoding='utf8'))
    expected={trial_identity(c) for c in jobs};observed=[trial_identity(r['config']) for r in browser]
    if len(expected)!=len(jobs) or len(set(observed))!=len(observed) or set(observed)-expected:
        raise ValueError('Duplicate or unexpected race records')
    screen=[]
    for row in browser:
        if row['status']!='passed' or row['config']['lattice_scale']!=1 and row['runtime']['surface'].get('corrugation_control_version')!='six-phase-atomic-crystal-contrast-v3':continue
        sample=next((s for s in row['samples'] if s['step']==80000),None)
        if sample:
            screen.append(dict(row,outcome='time-limit',finish_time_ps=None,
                final_forward_nm=sample['center_nm'][0]-row['initial_center_nm'][0],observed_time_ps=40.))
    progress=json.loads((out/'browser-confirm-report.json').read_text(encoding='utf8')) if (out/'browser-confirm-report.json').exists() else {}
    report={'status':'complete' if progress.get('status')=='complete' and set(observed)==expected else 'partial',
            'planned_trials':len(jobs),'observed_trials':len(browser),'missing_trials':len(expected-set(observed)),
            'screen_summary':summarize(screen,jobs),'race_summary':summarize(browser,jobs),
            'valid_surface_trials':sum(r['status']=='passed' and r.get('competitive_valid',False) for r in browser),
            'valid_finishes':sum(r['status']=='passed' and r.get('competitive_valid',False) and r.get('outcome')=='finish' for r in browser),
            'right_censored_trials':sum(r.get('right_censored',False) for r in browser),
            'baseline_source':'42 archived same-backend WASM scale=1 trials; explicit old/current compatibility gate',
            'metrics':'descriptive finite-seed effects, not seven distinct statistically proven ranks',
            'source_sha256':{str(f.relative_to(out)):hashlib.sha256(f.read_bytes()).hexdigest()
                for f in [out/'browser-confirm-jobs.json',out/'browser-confirm-report.json',*sorted((out/'browser-trials').glob('*.json'))]},
            'trial_files':[f'browser-trials/{r["key"]}.json' for r in browser],
            'protocol':progress.get('protocol'),
            'development_note':'earlier native/global-core-guard exploration is excluded from final comparisons'}
    write_json(out/'analysis.json',report)
    with (out/'finish-summary.csv').open('w',newline='',encoding='utf-8-sig') as f:
        names=['temperature_K','strength','lattice_scale','growth','molecule','planned','observed','missing','valid','finished','failed','lifted','censored','finish_mean_ps','finish_sd_ps']
        writer=csv.DictWriter(f,fieldnames=names);writer.writeheader()
        for g in report['race_summary']:
            for car in g['cars']:writer.writerow({k:(g[k] if k in g else car[k]) for k in names})
    sys.path.insert(0,str(ROOT/'.local/analysis-libs'))
    os.environ['MPLCONFIGDIR']=str(ROOT/'.local/matplotlib')
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    font=FontProperties(fname='C:/Windows/Fonts/msyh.ttc');plt.rcParams.update({'font.family':font.get_name(),'axes.unicode_minus':False,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    figures=out/'figures';figures.mkdir(exist_ok=True)
    conditions=report['screen_summary'];matrix=np.array([[next(c['mean_forward_nm'] for c in g['cars'] if c['molecule']==car) for g in conditions] for car in CARS],dtype=float)
    cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#dddddd')
    fig,ax=plt.subplots(figsize=(11,4.7));im=ax.imshow(matrix,aspect='auto',cmap=cmap,vmin=0)
    ax.set_xticks(range(len(conditions)),[f"起伏 {g['lattice_scale']:g}×\n增长 {g['growth']:g}×" for g in conditions]);ax.set_yticks(range(7),LABELS)
    for i,j in np.ndindex(matrix.shape):
        missing='未测' if conditions[j]['cars'][i]['planned']==0 else '拒绝'
        ax.text(j,i,f'{matrix[i,j]:.2f}' if np.isfinite(matrix[i,j]) else missing,ha='center',va='center',color='white' if np.isfinite(matrix[i,j]) and matrix[i,j]<np.nanmax(matrix)/2 else 'black')
    ax.set_title('40 ps 前进距离（nm）：同一模型、多种子平均\n150 K / 驱动 0.60；此图不能替代完整比赛成绩')
    fig.colorbar(im,ax=ax,label='前进距离 / nm');fig.tight_layout()
    for ext in ('png','svg'):fig.savefig(figures/f'01-screen.{ext}',dpi=200)
    plt.close(fig)
    groups=report['race_summary'];nrows,ncols=panel_shape(len(groups))
    fig,axes=plt.subplots(nrows,ncols,figsize=(4.7*ncols,4.15*nrows),squeeze=False,sharey=True)
    colors=plt.get_cmap('tab10')(np.arange(7))
    for ax,g in zip(axes.ravel(),groups):
        for i,car in enumerate(g['cars']):
            trials=[r for r in browser if condition_identity(r['config'])==condition_identity(g) and r['config']['molecule']==car['molecule']]
            for j,row in enumerate(trials):
                offset=(j-(len(trials)-1)/2)*.11
                if row['status']=='passed' and row.get('competitive_valid'):
                    if row['outcome']=='finish':ax.scatter(i+offset,row['finish_time_ps'],s=16,color=colors[i],alpha=.45)
                    else:ax.scatter(i+offset,row['observed_time_ps'],marker='^',s=45,facecolors='none',edgecolors=colors[i])
            if car['finish_mean_ps'] is not None:
                ax.errorbar(i,car['finish_mean_ps'],yerr=car['finish_sd_ps'] or 0,fmt='o',color=colors[i],capsize=3)
            else:ax.text(i,.92,'拒绝' if car['failed'] else f"{car['censored']}/{car['planned']} 超时" if car['censored'] else '样本未齐',rotation=90,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=8)
        ax.set_xticks(range(7),LABELS,rotation=45,ha='right');ax.set_title(f"起伏 {g['lattice_scale']:g}× / 增长 {g['growth']:g}×\n"+(f"慢/快 {g['slow_fast_ratio']:.2f}×" if g['all_finished_valid'] else '含未通过或未完赛样本'))
        ax.grid(axis='y',alpha=.2)
    for ax in axes.ravel()[len(groups):]:ax.set_visible(False)
    for ax in axes[:,0]:ax.set_ylabel('20 nm 完赛时间 / ps（均值 ± 样本标准差）')
    fig.suptitle(f'七种赛车：{len(groups)} 组条件，每分子每条件计划 3 个种子，已记录 {len(browser)}/{len(jobs)} 场\n实点为完赛；误差棒为样本标准差；空心三角为 400 ps 超时（不计作完赛）',y=1.02);fig.tight_layout()
    for ext in ('png','svg'):fig.savefig(figures/f'02-finish-times.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)
    # Every seed is shown; trajectories stop at their observed finish/timeout.
    selected=((1.,1.),(2.,1.),(2.,.25),(3.,1.))
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharey=True)
    for ax,(scale,growth) in zip(axes.ravel(),selected):
        for i,car in enumerate(CARS):
            trials=[r for r in browser if r['config']['molecule']==car and (r['config']['lattice_scale'],r['config']['growth'])==(scale,growth)]
            for j,row in enumerate(trials):
                if row['status']!='passed':continue
                ax.plot([s['time_ps'] for s in row['samples']],[s['center_nm'][0] for s in row['samples']],
                    color=colors[i],alpha=.65,lw=1.,label=LABELS[i] if j==0 else None)
        ax.axhline(20,color='#444444',ls='--',lw=.8);ax.set_ylim(-.4,21)
        ax.set_title(f'起伏 {scale:g}× / 沿程增长 {growth:g}×');ax.set_xlabel('模拟时间 / ps');ax.set_ylabel('质心 x / nm');ax.grid(alpha=.18)
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=7,bbox_to_anchor=(.5,-.025))
    fig.suptitle('实际比赛轨迹：每种分子展示全部 3 个种子\n20 nm 为终点；各面板时间轴随观测时长变化',y=1.02);fig.tight_layout()
    for ext in ('png','svg'):fig.savefig(figures/f'03-trajectories.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5.5))
    for condition,offset,color,label in (((1.,1.),-.14,'#4477aa','原配置'),((2.,1.),.14,'#cc6677','2 倍起伏，原沿程增长')):
        g=next(g for g in groups if (g['lattice_scale'],g['growth'])==condition)
        for i,car in enumerate(g['cars']):
            if car['finish_mean_ps'] is not None:
                ax.errorbar(car['finish_mean_ps'],i+offset,xerr=car['finish_sd_ps'],fmt='o',color=color,capsize=3,label=label if i==0 else None)
            for r in browser:
                if r['config']['molecule']==car['molecule'] and (r['config']['lattice_scale'],r['config']['growth'])==condition and r.get('competitive_valid') and r.get('finish_time_ps') is not None:
                    ax.scatter(r['finish_time_ps'],i+offset,s=14,color=color,alpha=.4)
    ax.set_yticks(range(7),LABELS);ax.invert_yaxis();ax.set_xlabel('20 nm 完赛时间 / ps（均值 ± 样本标准差）')
    ax.set_title('原先候选对照：1 倍与 2 倍晶格起伏\n150 K / 驱动 0.60 / 沿程增长 1；每种分子每条件 3 个种子')
    ax.legend(loc='lower right');ax.grid(axis='x',alpha=.2);fig.tight_layout()
    for ext in ('png','svg'):fig.savefig(figures/f'04-recommended-comparison.{ext}',dpi=200)
    plt.close(fig)
    # Directly compare the new midpoint with both bracket conditions.  Show
    # every seed alongside sample SD rather than suggesting inferred ranks.
    if any(g['lattice_scale']==1.5 for g in groups):
        fig,axes=plt.subplots(1,2,figsize=(15,6.5),sharex=True,sharey=True)
        for ax,growth in zip(axes,(1.,.25)):
            for scale,offset,color in ((1.,-.23,'#4477aa'),(1.5,0.,'#228833'),(2.,.23,'#cc6677')):
                group=next((g for g in groups if (g['lattice_scale'],g['growth'])==(scale,growth)),None)
                if group is None:continue
                ax.plot([],[],'o',color=color,label=f'起伏 {scale:g}×')
                for i,car in enumerate(group['cars']):
                    if car['finish_mean_ps'] is not None:
                        ax.errorbar(car['finish_mean_ps'],i+offset,xerr=car['finish_sd_ps'],fmt='o',color=color,capsize=3)
                    trials=[r for r in browser if condition_identity(r['config'])==condition_identity(group) and r['config']['molecule']==car['molecule']]
                    for j,row in enumerate(trials):
                        position=i+offset+(j-(len(trials)-1)/2)*.025
                        if row['status']=='passed' and row.get('competitive_valid'):
                            if row['outcome']=='finish':ax.scatter(row['finish_time_ps'],position,s=15,color=color,alpha=.4)
                            else:ax.scatter(row['observed_time_ps'],position,marker='>',s=45,facecolors='none',edgecolors=color)
                    if car['finish_mean_ps'] is None:
                        note='拒绝' if car['failed'] else f"{car['censored']}/{car['planned']} 超时" if car['censored'] else '样本未齐'
                        ax.text(.98,i+offset,note,ha='right',va='center',color=color,fontsize=8,transform=ax.get_yaxis_transform())
            ax.set_yticks(range(len(CARS)),LABELS);ax.set_title(f'沿程增长 {growth:g}×')
            ax.set_xlabel('20 nm 完赛时间 / ps');ax.grid(axis='x',alpha=.2);ax.legend(loc='lower right')
        axes[0].invert_yaxis()
        fig.suptitle('1 / 1.5 / 2 倍晶格起伏直接比较\n150 K / 驱动 0.60；每分子每条件 3 个种子；实点为均值 ± 样本标准差，浅点为各次完赛',y=1.01)
        fig.tight_layout()
        for ext in ('png','svg'):fig.savefig(figures/f'05-midpoint-comparison.{ext}',dpi=200,bbox_inches='tight')
        plt.close(fig)

        fig,axes=plt.subplots(1,2,figsize=(14,5.5),sharex=True,sharey=True)
        for ax,growth in zip(axes,(1.,.25)):
            for i,car in enumerate(CARS):
                trials=[r for r in browser if r['config']['molecule']==car and (r['config']['lattice_scale'],r['config']['growth'])==(1.5,growth)]
                for j,row in enumerate(trials):
                    if row['status']!='passed':continue
                    ax.plot([s['time_ps'] for s in row['samples']],[s['center_nm'][0] for s in row['samples']],
                            color=colors[i],alpha=.7,lw=1.,label=LABELS[i] if j==0 else None)
            ax.axhline(20,color='#444444',ls='--',lw=.8);ax.set_ylim(-.4,21)
            ax.set_title(f'起伏 1.5× / 沿程增长 {growth:g}×');ax.set_xlabel('模拟时间 / ps');ax.set_ylabel('质心 x / nm');ax.grid(alpha=.18)
        handles,labels=axes[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='lower center',ncol=len(CARS),bbox_to_anchor=(.5,0.))
        fig.suptitle('新增 1.5 倍起伏：实际比赛轨迹\n150 K / 驱动 0.60；每种分子展示全部 3 个种子；两面板使用相同坐标范围',y=1.02)
        fig.tight_layout(rect=(0.,.10,1.,1.))
        for ext in ('png','svg'):fig.savefig(figures/f'06-midpoint-trajectories.{ext}',dpi=200,bbox_inches='tight')
        plt.close(fig)
    print(json.dumps([{k:v for k,v in g.items() if k!='cars'} for g in groups],ensure_ascii=False,indent=2))


if __name__=='__main__':main()
