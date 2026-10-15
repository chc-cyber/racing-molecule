"""Combine independent timestep and atlas studies; never infer equivalence from n=3."""
from pathlib import Path
import json, statistics, hashlib, sys, os, gzip
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'.local/analysis-libs'))
os.environ['MPLCONFIGDIR']=str(ROOT/'.local/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'outputs/parallel-md-validation'
LABELS={'benzene-12':'Benzene','pyridine-11':'Pyridine','test-car-52':'Tetraphenylbenzene','fullerene-c60':'C60'}
CHINESE={'benzene-12':'苯','pyridine-11':'吡啶','test-car-52':'四苯基苯（52原子）','fullerene-c60':'C60','folded-eicosane':'折叠正二十烷'}
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def check(s):
    for f,h in s.items():
        p=Path(f);p=p if p.is_absolute() else ROOT/p
        assert hashlib.sha256(p.read_bytes()).hexdigest()==h,f
def main():
    bench=read(OUT/'benchmark.json');assert bench['status']=='complete';check(bench['sources'])
    assert len(bench['md'])==96 and len(bench['nve'])==24
    gate=read(ROOT/'outputs/c60-proxy-study/wasm-validation.json');assert gate['status']=='passed';check(gate['sources'])
    races=read(ROOT/'outputs/c60-proxy-study/races/browser-confirm-report.json')
    assert races['status']=='complete' and len(races['trials'])==24
    for f,h in races['source_sha256'].items():
        p=ROOT/'outputs/c60-proxy-study/races'/f if f.startswith('inputs/') else ROOT/f
        assert hashlib.sha256(p.read_bytes()).hexdigest()==h,f
    performance=[];groups=[]
    for car in LABELS:
        for scale in [1.5,2]:
            row=dict(car=car,scale=scale)
            for mode,dt,label in [('periodic',.5,'periodic'),('v1',.5,'v1'),('v2',.5,'v2'),('v2',1,'v2_1fs')]:
                values=[r['seconds_per_100ps'] for r in bench['md'] if r['car']==car and r['scale']==scale and r['mode']==mode and r['dt_fs']==dt]
                assert len(values)==3;row[label]=statistics.median(values)
            row['v2_speedup_over_periodic']=row['periodic']/row['v2'];row['v2_speedup_over_v1']=row['v1']/row['v2']
            performance.append(row)
            rs=[r for r in races['trials'] if r['config']['molecule']==car and r['config']['lattice_scale']==scale]
            assert len(rs)==3
            g=dict(car=car,scale=scale)
            for mode,directory in [('periodic','outputs/crystal-contrast-guarded/browser-trials'),('v1','outputs/wasm-fourier-validation/races/browser-trials'),('v2','outputs/c60-proxy-study/races/browser-trials')]:
                rows=[read(ROOT/directory/(r['key']+'.json')) for r in rs]
                assert all(r['status']=='passed' and r['outcome']=='finish' and r['competitive_valid'] for r in rows)
                values=[r['finish_time_ps'] for r in rows]
                g[mode+'_mean_ps']=statistics.mean(values);g[mode+'_sd_ps']=statistics.stdev(values)
            g['delta_periodic_percent']=100*(g['v2_mean_ps']/g['periodic_mean_ps']-1)
            g['delta_v1_percent']=100*(g['v2_mean_ps']/g['v1_mean_ps']-1)
            counts=np.zeros(3)
            for r in rs:
                assert r['runtime']['surface']['proxy_charts']==3 and r['runtime']['timestep_fs']==.5
                assert r['source_sha256']==races['source_sha256']
                assert r['samples'][-2]['center_nm'][0]<20<=r['samples'][-1]['center_nm'][0]
                counts+=np.sum([s['proxy_counts'] for s in r['samples']],axis=0)
                data=read(ROOT/'outputs/c60-proxy-study/races/browser-trials'/(r['key']+'.json'))
                coords=json.loads(gzip.decompress((ROOT/'outputs/c60-proxy-study/races/browser-trials'/data['coordinates_file']).read_bytes()))
                x=np.asarray(coords['positions']).reshape(-1,r['atoms'],3);m=np.asarray(coords['masses'])
                com=(x*m[None,:,None]).sum(axis=1)/m.sum()/10
                recorded=np.asarray([r['samples'][step//100]['center_nm'] for step in coords['steps']])
                assert np.max(np.abs(com-recorded))<1e-10
            g['coverage']=(counts/counts.sum()).tolist();groups.append(g)
    browser=read(OUT/'browser-v2-races.json');assert browser['passed'] and len(browser['cases'])==4
    assert read(OUT/'browser-ui.json')['passed']
    ts=read(ROOT/'outputs/timestep-study/summary.json')
    assert sum(g['trials'] for g in ts['race_groups'])==68, 'Incomplete timestep race study'
    assert not any(g['failed'] for g in ts['race_groups'])
    step_rows=[]
    for label in ['races-1','races-2','races-supplement']:
        saved=read(ROOT/f'outputs/timestep-study/{label}.json')
        assert saved['status']=='complete';check(saved['sources'])
        for r in saved['rows']:
            dt=r['config']['dt_fs'];assert r['status']=='complete'
            samples=r['samples'];assert all(abs(s['step']*dt/1000-s['time_ps'])<1e-9 for s in samples)
            xyzfile=ROOT/f"outputs/timestep-study/{label}/{r['key']}.json.coords.gz"
            coords=json.loads(gzip.decompress(xyzfile.read_bytes()))
            m=np.asarray(coords['masses']);x=np.asarray(coords['positions']).reshape(-1,len(m),3)
            com=(x*m[None,:,None]).sum(axis=1)/m.sum()/10
            recorded=np.asarray([samples[round(step*dt/50)]['center_nm'] for step in coords['steps']])
            assert np.max(np.abs(com-recorded))<1e-10,r['key']
            if r['outcome']=='finish':assert samples[-2]['center_nm'][0]<20<=samples[-1]['center_nm'][0]
            step_rows.append(r['key'])
    assert len(step_rows)==len(set(step_rows))==68
    hot=read(ROOT/'outputs/c60-proxy-study/temperature-600-validation.json')
    assert hot['status']=='complete' and len(hot['rows'])==18
    assert all(r['stability_guard_pass'] and r['numerical_screen_pass'] for r in hot['rows'])
    check(hot['sources'])
    review=read(OUT/'race-review.json');assert review['status']=='passed'
    result=dict(status='complete-experimental',performance=performance,race_groups=groups,timestep=ts,browser=browser['cases'],high_temperature=hot['protocol'])
    (OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
    for ax,scale in zip(axes,[1.5,2]):
        rows=[r for r in performance if r['scale']==scale];y=np.arange(4)
        for offset,key,color in [(-.25,'periodic','#9099a1'),(0,'v1','#4486ae'),(.25,'v2','#00978c')]:
            ax.barh(y+offset,[r[key] for r in rows],.23,color=color,label=key)
        ax.set_yticks(y,[LABELS[r['car']] for r in rows]);ax.invert_yaxis();ax.set_xlabel('Seconds per 100 ps (2 ps extrapolation)')
        ax.set_title(f'{scale:g}x corrugation | 0.5 fs');ax.axvline(30,color='#b16939',ls='--');ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
    axes[0].legend();fig.suptitle('Serial WASM MD | 150 K | drive 0.6 | median of 3')
    for ext in ['png','svg']:fig.savefig(OUT/f'performance.{ext}',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
    for ax,scale in zip(axes,[1.5,2]):
        selected=[g for g in groups if g['scale']==scale];x=np.arange(4)
        for offset,mode,color in [(-.24,'periodic','#9099a1'),(0,'v1','#4486ae'),(.24,'v2','#00978c')]:
            ax.errorbar(x+offset,[g[mode+'_mean_ps'] for g in selected],yerr=[g[mode+'_sd_ps'] for g in selected],fmt='o',capsize=4,color=color,label=mode)
        ax.set_xticks(x,[LABELS[g['car']] for g in selected],rotation=18)
        ax.set_ylabel('Finish time (ps), mean ± SD');ax.set_title(f'{scale:g}x corrugation | n = 3 per model')
        ax.grid(axis='y',alpha=.2)
    axes[0].legend();fig.suptitle('Full 20 nm races | 150 K | 0.5 fs | descriptive statistics')
    for ext in ['png','svg']:fig.savefig(OUT/f'race-comparison.{ext}',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,7.5),constrained_layout=True)
    for column,car in enumerate(['fullerene-c60','benzene-12']):
        trial=next(r for r in races['trials'] if r['config']['molecule']==car and r['config']['lattice_scale']==2 and r['config']['seed']==2718281)
        for mode,directory,color in [('periodic','outputs/crystal-contrast-guarded/browser-trials','#9099a1'),('v1','outputs/wasm-fourier-validation/races/browser-trials','#4486ae'),('v2','outputs/c60-proxy-study/races/browser-trials','#00978c')]:
            samples=read(ROOT/directory/(trial['key']+'.json'))['samples']
            t=np.array([s['time_ps'] for s in samples]);xyz=np.array([s['center_nm'] for s in samples])
            axes[0,column].plot(t,xyz[:,0],label=mode,color=color,lw=1.5)
            axes[1,column].plot(xyz[:,0],xyz[:,1]-xyz[0,1],label=mode,color=color,lw=1.3)
        axes[0,column].set_title(LABELS[car]);axes[0,column].set_xlabel('Time (ps)');axes[0,column].set_ylabel('Track position x (nm)')
        axes[1,column].set_xlabel('Track position x (nm)');axes[1,column].set_ylabel('Lateral displacement (nm)')
        for ax in axes[:,column]:ax.grid(alpha=.2)
    axes[0,0].legend();fig.suptitle('Example trajectories | 2x corrugation | seed 2718281\nSame seed does not imply identical stochastic trajectories')
    for ext in ['png','svg']:fig.savefig(OUT/f'trajectories.{ext}',dpi=180)
    plt.close(fig)
    lines=['# 时间步与C60傅里叶代理并行优化','',
      '游戏默认仍为完整周期势、0.5 fs。浏览器可选傅里叶已升级为三方向图代理；1 fs只作为本轮诊断与玩法敏感性探索，尚未提升为默认。','',
      '## C60覆盖与实现','',
      '新增三张方向图及其归一化平滑权重，覆盖球侧法向；高度支持扩为2.4–11.2 Å。保留标量势、全部位置/法向/混合导数和完整周期势回退。C表下载30.58 MB、系数36.99 MB；默认周期不下载。v1资产保留用于复现。','',
      '84组原生/WASM与18组Cartesian有限差分通过。126个既有C60快照全部使用代理，最大力差0.03747 kcal/mol/Å；8×8横向势垒变化误差0.000559 kcal/mol。低高度/近退化压力构型误差更大，详见[实现与数值报告](../c60-proxy-study/IMPLEMENTATION.md)，这些不是所有构型的误差界。','',
      '## 串行计算成本','',
      f"测试设备：{bench['hardware']['cpu']}；Node {bench['hardware']['node']}；真实浏览器 Chrome {browser['browser']}。这些数值不是其他电脑或移动端的等待时间保证。",'',
      '本机同初始相空间，2 ps × 3轮交替计时；下表单位秒/100 ps，为短段外推，排除准备、下载、绘制。1 fs列只是性能潜力，不能视作精度通过。','',
      '| 分子 | 起伏 | 周期0.5fs | v1 0.5fs | v2 0.5fs | v2 1fs探索 | v2相对周期加速 |','|---|---:|---:|---:|---:|---:|---:|']
    for r in performance:lines.append(f"| {CHINESE[r['car']]} | {r['scale']:g} | {r['periodic']:.2f} | {r['v1']:.2f} | {r['v2']:.2f} | {r['v2_1fs']:.2f} | {r['v2_speedup_over_periodic']:.2f}× |")
    lines+=['','![性能](performance.png)','','真实Chrome Worker完整20 nm比赛，含消息与定时器调度，不含下载、准备或Three.js绘制：','']
    for r in browser['cases']:lines.append(f"- {CHINESE[r['car']]}，起伏{r['scale']:g}，{r['dt']:g} fs：模拟{r['step']*r['dt']/1000:.2f} ps，实际{r['wall_seconds']:.2f}秒。"+(' **仅探索**' if r['dt']==1 else ''))
    lines+=['','C60在保持0.5 fs时约需一分钟，尚未达到半分钟目标；1 fs探索运行约35秒，但数值误差和玩法敏感性需要单独评估。以上仍排除了首次下载与准备，不能当作从打开网站到比赛结束的总等待时间。']
    lines+=['','## 新代理完整比赛','',
      '24场新代理比赛均有效完赛；检查松弛、结构、持续离面、终点穿越及保存坐标重算COM。对照为同参数历史周期/v1比赛；每组3种子、不是统计等效证明。','',
      '| 分子 | 起伏 | 周期均值±SD / ps | v1均值±SD / ps | v2均值±SD / ps | v2相对周期 | v2完全/混合/回退 |','|---|---:|---:|---:|---:|---:|---|']
    for g in groups:lines.append(f"| {CHINESE[g['car']]} | {g['scale']:g} | {g['periodic_mean_ps']:.2f}±{g['periodic_sd_ps']:.2f} | {g['v1_mean_ps']:.2f}±{g['v1_sd_ps']:.2f} | {g['v2_mean_ps']:.2f}±{g['v2_sd_ps']:.2f} | {g['delta_periodic_percent']:+.1f}% | "+'/'.join(f'{v:.3%}' for v in g['coverage'])+' |')
    lines+=['','![完整比赛均值与离散程度](race-comparison.png)','',
      '平均成绩对周期势变化为−2.86%至+3.04%，C60约0至+0.033%。但苯2倍起伏的SD由1.83增至6.92 ps，52原子车由3.91增至10.85 ps；1.5倍起伏时52原子车与吡啶的均值名次交换。n=3无法区分真实方差变化与抽样波动，因此不宣称成绩分布、轨迹或排名保持不变。详细横向游走和停滞比较见[RACE-REVIEW.md](RACE-REVIEW.md)。',
      '新C60比赛的完整代理查询约99.980%、混合约0.020%、完整回退0；静态旧快照的100%覆盖不代表任何时刻都无需混合。','',
      '![位移时间与横向轨迹示例](trajectories.png)','',
      '另补600 K、0.5 fs的C60/折叠二十烷/吡啶 × 1.5/2倍起伏：6组配对2 ps NVE（共12条）与6条5 ps受驱动短程检查。全部结构、有限性、离面与暂定能量预算通过；这是高温短程检查，不能替代600 K完整比赛统计。数据见[temperature-600-validation.json](../c60-proxy-study/temperature-600-validation.json)。']
    lines+=['','## 时间步与组合验证','',
      '128条2 ps NVE按相同坐标速度、物理时长比较0.25/0.5/1/2 fs。初筛全部通过有限性/键长检查；0.5 fs全部达到暂定1%初始动能误差预算。1 fs误差普遍增大，逐步采样又发现C60峰值约1.019%，因此不能依稀疏采样认定1 fs通过。1%是本项目保守数值预算，不是普适物理边界。',
      '按成对条件统计，1 fs能量误差约为0.5 fs的4.17倍（中位数），2 fs约19.59倍。本轮未观察到数值爆炸；这反映数值精度损失，不能仅凭结构稳定接受更大步长。',
      '1 fs继续做带热浴的完整比赛敏感性对照，但不以热浴掩盖NVE误差；具体成绩、游走、停滞及删失见[时间步报告](../timestep-study/REPORT.md)。2 fs未进入完整比赛测试。','',
      '改变步长也改变了Langevin随机数采样与物理时间的对应关系，同种子不代表同一条噪声路径。因此单场差异不能全部归因于积分误差；NVE用来隔离数值误差，多种子比赛用于观察实际玩法敏感性。','',
      '![步长误差](../timestep-study/nve-timestep.png)','',
      '最终主构建另做新代理0.25/0.5/1 fs共24段NVE以及96段受驱动计时；组合误差保留在benchmark.json。400 fs脉冲、350 fs阻尼和50 fs采样保持物理时间不变。诊断赛限400 ps，浏览器原默认上限5 ns，未混用。准备阶段仍按0.5 fs松弛。','',
      '## 发布与剩余限制','',
      '默认周期势及0.5 fs未改。新增Python代理是开发参考/生成工具，未替换本地Python游戏求力路径；性能收益属于浏览器WASM。完整包保留两代系数，因此ZIP大小不等于玩家本次按需下载量。未公开部署。移动端、更多浏览器、真实网络首载和更广分子统计仍需验收。','',
      '独立审查见CODE-REVIEW.md、INDEPENDENT-REVIEW.md和RACE-REVIEW.md；worker-scheduling.json为调度mock，browser-ui.json为真实页面切换。不同表面模式与步长的成绩应分开记录。完整试玩包：[molecular-car-web-preview.zip](../wasm-preparation/molecular-car-web-preview.zip)。','']
    step_table=['| 分子 | 模型 | T/K | 有效完赛0.5→1 fs | 均值成绩0.5→1 fs / ps | 成绩变化 | 前50ps横向RMS均值0.5→1 fs / nm |',
      '|---|---|---:|---:|---:|---:|---:|']
    for g in ts['race_groups']:
        if g['dt_fs']!=1:continue
        base=next(r for r in ts['race_groups'] if r['dt_fs']==.5 and all(r[k]==g[k] for k in ['car','surface','temperature_K','lattice_scale']))
        fmt=lambda x:'—' if x is None else f'{x:.3f}'
        change=g.get('mean_finish_change_pct');change='—' if change is None else f'{change:+.2f}%'
        name=CHINESE.get(g['car'],g['car'])
        mode='v1' if g['surface']=='fourier-au' else 'periodic'
        step_table.append(f"| {name} | {mode} | {g['temperature_K']} | {base['finish_time_ps']['n']}/{base['trials']} → {g['finish_time_ps']['n']}/{g['trials']} | {fmt(base['finish_time_ps']['mean'])} → {fmt(g['finish_time_ps']['mean'])} | {change} | {fmt(base['first_50ps_lateral_rms_nm']['mean'])} → {fmt(g['first_50ps_lateral_rms_nm']['mean'])} |")
    step_table+=['','共40场初筛与28场针对性补样，所有条件为2倍起伏。表中150 K的v1苯/吡啶每步长各10种子，其余150 K条件各3种子，600 K各1种子；初筛与合并统计在时间步报告中分别保留。未完赛或持续离面的场次不计入成绩均值，不能把剩余均值当作全体成绩。表中v1是冻结的旧代理，用于单独研究时间步；新版v2的组合短程验证和浏览器1 fs探索另列。','',
      '10种子下，苯/吡啶平均成绩分别变化−2.07%/−4.46%，前50ps横向RMS均值则增加约47%/34%。苯横向RMS为0.387±0.154→0.569±0.238 nm，吡啶0.418±0.174→0.562±0.242 nm。探索性配对bootstrap区间与低前进指标均保留在时间步报告；这些结果不能等同于普遍的扩散系数变化，也不是统计等效验收。','',
      '150 K折叠长链在0.5 fs下三场均在400 ps右删失；1 fs下一场287.85 ps完赛、两场同样删失。287.85 ps仅是唯一成功场次的成绩，不能作为该条件总体平均，不能据此宣称长链稳定提速。','',
      '![固定观察窗口的运动行为](../timestep-study/behavior-timestep.png)','']
    index=lines.index('![步长误差](../timestep-study/nve-timestep.png)')
    lines[index:index]=step_table
    (OUT/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(dict(status=result['status'],performance=performance,groups=groups),indent=2))
if __name__=='__main__':main()
