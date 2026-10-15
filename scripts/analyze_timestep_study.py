"""Summarize a timestep study without interpreting finite trajectories as equivalence."""
from pathlib import Path
import json
import statistics as stats
import random
from collections import defaultdict

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/timestep-study'
def load(name):
    return json.loads((OUT/name).read_text(encoding='utf-8'))
def summarize(values):
    values=[x for x in values if x is not None]
    return {'n':len(values),'mean':stats.mean(values) if values else None,'sd':stats.stdev(values) if len(values)>1 else None,'min':min(values) if values else None,'max':max(values) if values else None}

nve=load('nve.json')
groups=defaultdict(list)
for r in nve['rows']:groups[r['dt_fs']].append(r)
summary={'nve_status':nve['status'],'nve_by_dt':{},'race_groups':[],'benchmark_groups':[]}
if (OUT/'sampling.json').exists():summary['peak_sampling_check']=load('sampling.json')
for dt,rows in groups.items():
    summary['nve_by_dt'][dt]={'trials':len(rows),'stability_pass':sum(r.get('stability_guard_pass',False) for r in rows),'precision_pass':sum(r.get('numerical_screen_pass',False) for r in rows),'relative_error':summarize([r.get('relative_to_initial_kinetic') for r in rows]),'failures':[{'car':r['car'],'surface':r['surface'],'temperature_K':r['temperature_K'],'scale':r['lattice_scale'],'error':r.get('error')} for r in rows if r['status']!='complete']}
summary['nve_paired_error_ratios']=[]
for row in nve['rows']:
    base=next((r for r in nve['rows'] if all(r[k]==row[k] for k in ['car','surface','temperature_K','lattice_scale']) and r['dt_fs']==.5),None)
    if base and base.get('max_energy_error') and row.get('max_energy_error') is not None:
        summary['nve_paired_error_ratios'].append({**{k:row[k] for k in ['car','surface','temperature_K','lattice_scale','dt_fs']},'error_vs_half_fs':row['max_energy_error']/base['max_energy_error'],'energy_slope_kcal_mol_ps':row['energy_linear_slope_kcal_mol_ps'],'first_half_error':row['first_half_max_error'],'second_half_error':row['second_half_max_error']})

race_groups=defaultdict(list);all_trials={}
for file in OUT.glob('races*.json'):
    report=json.loads(file.read_text(encoding='utf-8'))
    if not isinstance(report,dict) or 'rows' not in report:continue
    for r in report['rows']:all_trials[r['key']]=r
# Include atomically published helper caches before the owner reaches that job.
for file in OUT.glob('races*/*.json'):
    r=json.loads(file.read_text(encoding='utf-8'))
    if 'config' in r and 'key' in r:all_trials[r['key']]=r
for r in all_trials.values():
    c=r['config'];key=(c['car'],c['surface'],c['temperature_K'],c['lattice_scale'],c['dt_fs']);race_groups[key].append(r)
if (OUT/'race-jobs-all.json').exists():
    planned=load('race-jobs-all.json');completed={json.dumps(r['config'],sort_keys=True) for r in all_trials.values()}
    summary['race_progress']={'planned_base_jobs':len(planned),'recorded_base_jobs':sum(json.dumps(c,sort_keys=True) in completed for c in planned),'unique_results_including_supplements':len(all_trials)}
def aggregate_races(race_groups):
    result=[]
    for key,rows in sorted(race_groups.items()):
        row=dict(zip(['car','surface','temperature_K','lattice_scale','dt_fs'],key))
        row.update(trials=len(rows),finished=sum(r.get('outcome')=='finish' for r in rows),valid_finished=sum(r.get('outcome')=='finish' and r.get('competitive_valid',False) for r in rows),right_censored=sum(r.get('right_censored',False) for r in rows),failed=sum(r['status']!='complete' for r in rows),competitive_valid=sum(r.get('competitive_valid',False) for r in rows),invalid_finishes=sum(r.get('outcome')=='finish' and not r.get('competitive_valid',False) for r in rows))
        row['finish_time_ps']=summarize([r.get('finish_time_ps') for r in rows if r.get('competitive_valid',False)])
        row['observed_final_forward_nm']=summarize([r['samples'][-1]['center_nm'][0]-r['initial_center_nm'][0] for r in rows if r.get('samples')])
        row['bond_ratio_min']=min((r['bond_ratio_range'][0] for r in rows if 'bond_ratio_range' in r),default=None)
        row['bond_ratio_max']=max((r['bond_ratio_range'][1] for r in rows if 'bond_ratio_range' in r),default=None)
        for metric in ['lateral_rms_nm','low_progress_fraction','maximum_lateral_excursion_nm','maximum_lift_duration_ps']:row[metric]=summarize([r.get(metric) for r in rows])
        fixed_metrics=defaultdict(list)
        for trial in rows:
            if trial.get('observed_time_ps',0)<50:continue
            samples=[s for s in trial['samples'] if s['time_ps']<=50+1e-9]
            if len(samples)!=1001:continue
            y0=trial['initial_center_nm'][1]
            fixed_metrics['first_50ps_lateral_rms_nm'].append((sum((s['center_nm'][1]-y0)**2 for s in samples)/len(samples))**.5)
            fixed_metrics['first_50ps_low_progress_fraction'].append(sum(samples[i]['center_nm'][0]-samples[i-20]['center_nm'][0]<.02 for i in range(20,len(samples)))/(len(samples)-20))
            fixed_metrics['first_50ps_forward_nm'].append(samples[-1]['center_nm'][0]-samples[0]['center_nm'][0])
        for metric in ['first_50ps_lateral_rms_nm','first_50ps_low_progress_fraction','first_50ps_forward_nm']:row[metric]=summarize(fixed_metrics[metric])
        row['sensitivity_only']=any(r['config'].get('sensitivity_only',False) for r in rows)
        result.append(row)
    for row in result:
        base=next((r for r in result if all(r[k]==row[k] for k in ['car','surface','temperature_K','lattice_scale']) and r['dt_fs']==.5),None)
        if base and base['finish_time_ps']['mean'] and row['finish_time_ps']['mean']:
            row['mean_finish_change_pct']=100*(row['finish_time_ps']['mean']/base['finish_time_ps']['mean']-1)
    return result
summary['race_groups']=aggregate_races(race_groups)
summary['race_outcomes']=[{'key':r['key'],'status':r['status'],'outcome':r.get('outcome'),'valid_finish':r.get('competitive_valid',False) and r.get('outcome')=='finish','right_censored':r.get('right_censored'),'finish_time_ps':r.get('finish_time_ps'),'observed_time_ps':r.get('observed_time_ps'),'final_center_nm':r.get('samples',[{}])[-1].get('center_nm'),'bond_ratio_range':r.get('bond_ratio_range')} for r in all_trials.values()]
base_keys={json.dumps(c,sort_keys=True) for c in load('race-jobs-all.json')}
base_groups=defaultdict(list)
for r in all_trials.values():
    if json.dumps(r['config'],sort_keys=True) in base_keys:
        c=r['config'];base_groups[(c['car'],c['surface'],c['temperature_K'],c['lattice_scale'],c['dt_fs'])].append(r)
summary['base_race_groups']=aggregate_races(base_groups)
summary['race_progress']['recorded_supplement_jobs']=len(all_trials)-summary['race_progress']['recorded_base_jobs']
summary['race_progress']['planned_supplement_jobs']=28
# Exploratory paired bootstrap, fixed analysis RNG; no new MD and no equivalence claim.
summary['paired_bootstrap']={'analysis_seed':20261014,'replicates':10000,'qualification':'Exploratory percentile interval for mean(1fs - 0.5fs), paired by thermostat seed; not an equivalence test','rows':[]}
rng=random.Random(20261014)
def metric_value(r,metric):
    if metric=='finish_time_ps':return r['finish_time_ps']
    samples=[x for x in r['samples'] if x['time_ps']<=50+1e-9]
    if len(samples)!=1001:return None
    if metric=='first_50ps_lateral_rms_nm':return (sum((x['center_nm'][1]-r['initial_center_nm'][1])**2 for x in samples)/len(samples))**.5
    return sum(samples[i]['center_nm'][0]-samples[i-20]['center_nm'][0]<.02 for i in range(20,len(samples)))/(len(samples)-20)
def percentile(v,q):
    x=(len(v)-1)*q;i=int(x);return v[i]+(v[min(i+1,len(v)-1)]-v[i])*(x-i)
for car in ['benzene-12','pyridine-11']:
    matched={dt:{r['config']['seed']:r for r in all_trials.values() if r['config']['car']==car and r['config']['surface']=='fourier-au' and r['config']['temperature_K']==150 and r['config']['lattice_scale']==2 and r['config']['dt_fs']==dt and r.get('competitive_valid',False) and r.get('outcome')=='finish'} for dt in [.5,1]}
    seeds=sorted(matched[.5].keys()&matched[1].keys())
    if not seeds:continue
    for metric in ['finish_time_ps','first_50ps_lateral_rms_nm','first_50ps_low_progress_fraction']:
        a=[metric_value(matched[.5][seed],metric) for seed in seeds];b=[metric_value(matched[1][seed],metric) for seed in seeds]
        if any(x is None for x in a+b):continue
        diff=[y-x for x,y in zip(a,b)]
        boot=sorted(stats.mean(rng.choices(diff,k=len(diff))) for _ in range(10000))
        summary['paired_bootstrap']['rows'].append({'car':car,'metric':metric,'n':len(seeds),'seeds':seeds,'half_fs':summarize(a),'one_fs':summarize(b),'paired_difference':summarize(diff),'mean_difference_percentile_95':[percentile(boot,.025),percentile(boot,.975)]})
if (OUT/'benchmark.json').exists():
    bg=defaultdict(list)
    for row in load('benchmark.json')['rows']:bg[(row['car'],row['surface'],row['lattice_scale'],row['dt_fs'])].append(row['seconds_per_100ps'])
    for key,v in sorted(bg.items()):summary['benchmark_groups'].append({**dict(zip(['car','surface','lattice_scale','dt_fs'],key)),**summarize(v)})
joint_file=ROOT/'outputs/parallel-md-validation/benchmark.json'
if joint_file.exists():
    joint=json.loads(joint_file.read_text(encoding='utf-8'))
    if joint.get('status')=='complete':
        summary['joint_v2_benchmark']={'source':str(joint_file.relative_to(ROOT)),'qualification':'New v2 atlas; independent serial 2 ps timing extrapolated to 100 ps, not the v1 backend used for timestep behavior races','rows':[]}
        for car,scale in sorted({(r['car'],r['scale']) for r in joint['md'] if r['mode']=='v2'}):
            values={dt:stats.mean(r['seconds_per_100ps'] for r in joint['md'] if r['mode']=='v2' and r['car']==car and r['scale']==scale and r['dt_fs']==dt) for dt in [.5,1]}
            summary['joint_v2_benchmark']['rows'].append({'car':car,'scale':scale,'half_fs_seconds_per_100ps':values[.5],'one_fs_seconds_per_100ps':values[1],'speedup':values[.5]/values[1]})
(OUT/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
lines=['# 时间步扫描验证','',f"NVE 扫描状态：{nve['status']}。初始相空间、准备、势模型和物理观察时长匹配；默认步长仍为 0.5 fs。",'',
'短程稳定性检查与精度目标分开：稳定性只表示本次有限时段无非有限数值及键长越界；精度目标为最大能量偏差低于初始动能的 1%，它不是理论稳定极限或游戏可玩性界限；完整比赛用于另外判断实际玩法敏感性。','',
'| 步长 fs | 完成数 | 稳定性 guard | 10 fs 抽样的 1% 目标 | 最大抽样相对能量偏差 |','|---:|---:|---:|---:|---:|']
for dt,row in summary['nve_by_dt'].items():
    error=row['relative_error']['max'];lines.append(f"| {dt} | {row['trials']} | {row['stability_pass']} | {row['precision_pass']} | {100*error:.3f}% |" if error is not None else f"| {dt} | {row['trials']} | {row['stability_pass']} | {row['precision_pass']} | 无有效数据 |")
ratios={dt:stats.median(r['error_vs_half_fs'] for r in summary['nve_paired_error_ratios'] if r['dt_fs']==dt) for dt in [.25,1,2]}
lines+=['',f"配对最大能量误差相对0.5fs的中位数：0.25fs为{ratios[.25]:.3f}倍，1fs为{ratios[1]:.3f}倍，2fs为{ratios[2]:.3f}倍。缩小步长可系统降低误差；1fs约四倍误差，与二阶积分的预期趋势相符。各条件线性能量斜率及前/后半段误差保存在 summary.json；有限短程有界不等于任意长时间稳定。"]
if 'peak_sampling_check' in summary:
    lines+=['','**逐积分步复核修正：** 另对 C60 150 K / 2 倍起伏及吡啶 600 K / 2 倍起伏做 500 fs 每步观测。C60 的 1 fs 最大误差为初始动能的 **1.01918%**，超过预设预算；同条轨迹每 10 fs 抽样会低估峰值约 16.83%。因此表中 C60 的抽样通过不构成精度验收。0.5 fs 的对应逐步峰值为 0.25056%；吡啶 1 fs 为 2.35297%。这属于有界积分误差的精度问题，未观察到数值爆炸。详见 `sampling.json`。','',
    '![时间步能量误差](nve-timestep.png)']
lines+=['','完整比赛统计：基础40场与额外28场轻分子补样分开保存。下表合并苯/吡啶旧傅里叶 T150K 的10个种子，基础筛选的3种子统计另见下表及 summary.json 的 base_race_groups。相同种子用于控制抽样条件，不要求轨迹重合。均值和标准差为描述统计；少量种子不能证明统计等效。成绩均值只含无持续离面的有效完赛；未完赛及无效完赛不外推成绩，详见 summary.json。','',
'| 分子 | 模型 | T/K | 起伏 | dt/fs | 有效完赛/总数 | 删失/异常 | 有效平均成绩 ps | 相对 0.5 fs |','|---|---|---:|---:|---:|---:|---:|---:|---:|']
for row in summary['race_groups']:
    m=row['finish_time_ps']['mean'];d=row.get('mean_finish_change_pct');prefix=f"| {row['car']} | {row['surface']} | {row['temperature_K']} | {row['lattice_scale']} | {row['dt_fs']} | {row['valid_finished']}/{row['trials']} | {row['right_censored']}/{row['failed']} |"
    lines.append(prefix+(f" {m:.3f} | {d:+.2f}% |" if m is not None and d is not None else (f" {m:.3f} | — |" if m is not None else " 无有效完赛 | — |")))
lines+=['','![比赛成绩与删失](race-timestep.png)','','基础40场中苯/吡啶旧傅里叶的三种子筛选（保留扩样前判断）：','','| 分子 | dt/fs | 完赛/总数 | 平均成绩 ps | 相对0.5fs |','|---|---:|---:|---:|---:|']
for r in summary['base_race_groups']:
    if r['car'] not in ['benzene-12','pyridine-11'] or r['surface']!='fourier-au' or r['temperature_K']!=150:continue
    m=r['finish_time_ps']['mean'];d=r.get('mean_finish_change_pct')
    lines.append(f"| {r['car']} | {r['dt_fs']} | {r['valid_finished']}/{r['trials']} | {m:.3f} | {d:+.2f}% |" if m is not None and d is not None else f"| {r['car']} | {r['dt_fs']} | {r['valid_finished']}/{r['trials']} | — | — |")
lines+=['','柔性链逐种子结果（含未完赛，避免成功者选择偏差）：','','| T/K | dt/fs | 种子 | 结果 | 观察至 ps | 最终 x/nm |','|---:|---:|---:|---|---:|---:|']
for r in sorted(all_trials.values(),key=lambda r:r['key']):
    c=r['config']
    if c['car']!='folded-eicosane' or r['status']!='complete':continue
    outcome='有效完赛' if r['outcome']=='finish' and r['competitive_valid'] else ('无效完赛' if r['outcome']=='finish' else '400 ps右删失')
    lines.append(f"| {c['temperature_K']} | {c['dt_fs']} | {c['seed']} | {outcome} | {r['observed_time_ps']:.2f} | {r['samples'][-1]['center_nm'][0]:.3f} |")
if summary['race_groups']:
    lines+=['','运动行为：横向 RMS 以各自初始质心 y 为零点；低前进比例指任意 1 ps 窗口前进不足 0.02 nm（含后退）的占比。完整赛程指标会受到不同完赛时长影响，因此额外报告前 50 ps 的等物理观察窗口。只有覆盖完整 50 ps 的轨迹进入该统计，未覆盖者不外推。','',
    '| 分子/模型/T | dt fs | 赛程横向RMS nm | 赛程低前进 % | 前50ps横向RMS nm | 前50ps低前进 % | 最大连续离面 ps |','|---|---:|---:|---:|---:|---:|---:|']
    for r in summary['race_groups']:
        fields=[('lateral_rms_nm',1),('low_progress_fraction',100),('first_50ps_lateral_rms_nm',1),('first_50ps_low_progress_fraction',100)]
        values=['—' if r[k]['mean'] is None else f"{r[k]['mean']*factor:.4f}" for k,factor in fields]
        lift=r['maximum_lift_duration_ps']['max']
        lines.append(f"| {r['car']} / {r['surface']} / {r['temperature_K']}K | {r['dt_fs']} | "+' | '.join(values)+f" | {lift if lift is not None else '—'} |")
lines+=['','标记 sensitivity_only 的 1 fs 条件未通过预设 1% 精度目标，比赛结果仅回答玩法敏感性，不能用于自动提高默认步长。','',
'赛车每 50 fs 采样，坐标每 500 fs 保存；阻尼保持 350 fs，脉冲保持 400 fs，比赛上限按物理时间设置。时长差异最多含 0.05 ps 终点检测量化。','',
'数据：`nve.json`、`summary.json`、`races*.json` 及各子目录的完整轨迹。源文件和隔离 WASM 版本散列随结果保存。']
lines+=['','![相同物理窗口下的运动行为](behavior-timestep.png)','','轻分子按种子配对的不确定性参考：固定分析随机种子 20261014、10000 次有放回重采样，以下区间描述平均差（1 fs − 0.5 fs）。这是探索性百分位区间，不是等效检验或多重检验校正结果。','','| 分子 | 指标 | 配对n | 0.5 fs均值±SD | 1 fs均值±SD | 平均差95%探索区间 |','|---|---|---:|---:|---:|---:|']
for r in summary['paired_bootstrap']['rows']:
    a,b,ci=r['half_fs'],r['one_fs'],r['mean_difference_percentile_95']
    if a['sd'] is None or b['sd'] is None:continue
    lines.append(f"| {r['car']} | {r['metric']} | {r['n']} | {a['mean']:.4f} ± {a['sd']:.4f} | {b['mean']:.4f} ± {b['sd']:.4f} | [{ci[0]:+.4f}, {ci[1]:+.4f}] |")
if 'joint_v2_benchmark' in summary:
    lines+=['','时间步的成本收益来自另一路联合串行基准：**新 v2 atlas** 的 0.5 / 1 fs 配对测速。下表是 2 ps 段外推的 100 ps 计算耗时，不含下载、准备和渲染。上文比赛行为使用隔离旧 v1 势，两个比较的模型版本不同，不把其耗时混作同一次实验。','',
    '| 分子 | 起伏 | 0.5 fs 秒/100ps | 1 fs 秒/100ps | 加速 |','|---|---:|---:|---:|---:|']
    for r in summary['joint_v2_benchmark']['rows']:lines.append(f"| {r['car']} | {r['scale']} | {r['half_fs_seconds_per_100ps']:.2f} | {r['one_fs_seconds_per_100ps']:.2f} | {r['speedup']:.2f}x |")
    lines+=['','基准原始数据：[联合 benchmark.json](../parallel-md-validation/benchmark.json)。1 fs 耗时仅供候选评估，不代表该步长已获准作为默认。']
(OUT/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps(summary['nve_by_dt'],indent=2))
