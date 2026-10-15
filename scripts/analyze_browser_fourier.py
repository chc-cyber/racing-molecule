"""Summarize bounded Fourier accuracy, serial costs and full-race behavior."""
from pathlib import Path
import os,sys,json,gzip,hashlib,statistics
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'.local/analysis-libs'));os.environ['MPLCONFIGDIR']=str(ROOT/'.local/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'outputs/wasm-fourier-validation'
LABELS={'benzene-12':'Benzene','pyridine-11':'Pyridine','test-car-52':'Tetraphenylbenzene','fullerene-c60':'C60'}
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def check(s):
    for f,h in s.items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h,f
def motion(r):
    c=np.array([s['center_nm'] for s in r['samples']]);dy=c[:,1]-c[0,1]
    return dict(lateral_rms_nm=float(np.sqrt(np.mean(dy*dy))),max_lateral_nm=float(np.max(np.abs(dy))),
        low_progress_1ps_fraction=float(np.mean(c[20:,0]-c[:-20,0]<.02)))
def main():
    num=read(OUT/'numerics.json');bench=read(OUT/'benchmark.json');browser=read(OUT/'browser-validation.json')
    for r in [num,bench]:assert r['status']=='passed';check(r['sources'])
    assert browser['passed']
    assert read(OUT/'browser-ui.json')['passed']
    reference=read(OUT/'periodic-regression.json')
    assert reference['status']=='passed' and reference['cases']==108
    check(reference['sources'])
    assert len(num['static'])==108 and len(bench['md'])==48 and len(bench['nve'])==16
    performance=[]
    for car in LABELS:
        for scale in [1.5,2]:
            a=[x for x in bench['md'] if x['car']==car and x['scale']==scale]
            times=[statistics.median(x['seconds_per_100ps'] for x in a if x['proxy']==p) for p in [False,True]]
            counts=np.sum([x['counts'] for x in a if x['proxy']],axis=0)
            performance.append(dict(car=car,scale=scale,periodic_s=times[0],fourier_s=times[1],speedup=times[0]/times[1],coverage=(counts/counts.sum()).tolist()))
            drifts=[next(x['max_total_energy_drift'] for x in bench['nve'] if x['car']==car and x['scale']==scale and x['proxy']==p) for p in [False,True]]
            assert drifts[1]<=drifts[0]*1.1+1e-5,(car,scale,drifts)
    races=read(OUT/'races/browser-confirm-report.json')
    assert races['status']=='complete' and races['planned']==len(races['trials'])==24
    for f,h in races['source_sha256'].items():
        source=OUT/'races'/f if f.startswith('inputs/') else ROOT/f
        assert hashlib.sha256(source.read_bytes()).hexdigest()==h,f
    expected={(car,scale,seed) for car in LABELS for scale in [1.5,2] for seed in [2718281,1234567,31415927]}
    assert {(r['config']['molecule'],r['config']['lattice_scale'],r['config']['seed']) for r in races['trials']}==expected
    comparisons=[]
    for r in races['trials']:
        assert r['status']=='passed' and r['outcome']=='finish' and r['competitive_valid'],r.get('error',r['key'])
        assert r['source_sha256']==races['source_sha256']
        assert r['runtime']['surface']['fourier_enabled']
        assert r['runtime']['preparation']['global_force_norm_kcal_mol_A']<=1e-4
        assert r['samples'][-2]['center_nm'][0]<20<=r['samples'][-1]['center_nm'][0]
        old=read(ROOT/'outputs/crystal-contrast-guarded/browser-trials'/(r['key']+'.json'))
        assert r['config']==old['config']
        coords=json.loads(gzip.decompress((OUT/'races/browser-trials'/r['coordinates_file']).read_bytes()))
        p=np.asarray(coords['positions']).reshape(-1,r['atoms'],3);m=np.asarray(coords['masses'])
        com=np.sum(p*m[None,:,None],axis=1)/sum(m)/10
        comerror=float(np.max(np.abs(com-np.asarray([r['samples'][k//100]['center_nm'] for k in coords['steps']]))))
        assert comerror<1e-10
        counts=np.sum([s['proxy_counts'] for s in r['samples']],axis=0)
        comparisons.append(dict(car=r['config']['molecule'],scale=r['config']['lattice_scale'],seed=r['config']['seed'],
            periodic_finish_ps=old['finish_time_ps'],fourier_finish_ps=r['finish_time_ps'],
            periodic_motion=motion(old),fourier_motion=motion(r),coverage=(counts/counts.sum()).tolist(),
            com_error_nm=comerror,bond_ratio_range=r['bond_ratio_range'],height_range_A=r['atom_height_range_A']))
    groups=[]
    for car in LABELS:
        for scale in [1.5,2]:
            rows=[r for r in comparisons if r['car']==car and r['scale']==scale]
            g=dict(car=car,scale=scale)
            for prefix in ['periodic','fourier']:
                values=[r[prefix+'_finish_ps'] for r in rows]
                g[prefix+'_mean_ps']=statistics.mean(values);g[prefix+'_sd_ps']=statistics.stdev(values)
                g[prefix+'_mean_lateral_rms_nm']=statistics.mean(r[prefix+'_motion']['lateral_rms_nm'] for r in rows)
                g[prefix+'_mean_low_progress_fraction']=statistics.mean(r[prefix+'_motion']['low_progress_1ps_fraction'] for r in rows)
            g['finish_mean_change_percent']=100*(g['fourier_mean_ps']/g['periodic_mean_ps']-1)
            groups.append(g)
    summary=dict(status='passed',performance=performance,groups=groups,races=comparisons,browser=[{k:r[k] for k in ['scale','step','wall_seconds']} for r in browser['cases']],
        max_port_force_error=max(r['port_force_error'] for r in num['static']),
        max_proxy_force_error=max(r['approximation_force_error'] for r in num['static']),
        max_proxy_energy_error=max(r['approximation_energy_error'] for r in num['static']),
        nve=bench['nve'],limits=['Full race comparison uses archived periodic model three-seed runs, not identical compiler arithmetic.',
          'Three seeds give descriptive checks, not a population equivalence proof.',
          'Browser timing is Worker-only in an inert page; excludes download, preparation and Three.js rendering.',
          'Proxy support remains bounded; C60 still largely falls back.'])
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),constrained_layout=True)
    for ax,scale in zip(axes,[1.5,2]):
        rows=[r for r in performance if r['scale']==scale];y=np.arange(4)
        ax.barh(y-.18,[r['periodic_s'] for r in rows],.34,color='#84919c',label='Periodic reference')
        ax.barh(y+.18,[r['fourier_s'] for r in rows],.34,color='#008d8f',label='Fourier proxy')
        ax.set_yticks(y,[LABELS[r['car']] for r in rows]);ax.invert_yaxis();ax.axvline(30,color='#bd6242',ls='--')
        ax.set_title(f'Corrugation {scale:g}x');ax.set_xlabel('Seconds per 100 ps (2 ps extrapolation)');ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
    axes[0].legend();fig.suptitle('Serial WASM MD | 150 K | drive 0.6 | dt 0.5 fs | median of 3')
    for fmt in ['png','svg']:fig.savefig(OUT/f'performance.{fmt}',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for ax,scale in zip(axes,[1.5,2]):
        for i,car in enumerate(LABELS):
            g=next(g for g in groups if g['car']==car and g['scale']==scale)
            for offset,prefix,color in [(-.13,'periodic','#84919c'),(.13,'fourier','#008d8f')]:
                ax.errorbar(i+offset,g[prefix+'_mean_ps'],yerr=g[prefix+'_sd_ps'],fmt='o',color=color,capsize=4,label=prefix if i==0 else None)
                rows=[r for r in comparisons if r['car']==car and r['scale']==scale]
                ax.scatter([i+offset]*3,[r[prefix+'_finish_ps'] for r in rows],s=15,color=color,alpha=.5)
        ax.set_xticks(range(4),list(LABELS.values()),rotation=15);ax.set_ylabel('Finish time / ps');ax.set_title(f'{scale:g}x | mean ± sample SD; n=3');ax.grid(axis='y',alpha=.2)
    axes[0].legend();fig.suptitle('20 nm races: periodic reference vs bounded Fourier approximation')
    for fmt in ['png','svg']:fig.savefig(OUT/f'race-scores.{fmt}',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    colors=['#007c91','#b47700','#8f4ea0','#32834d']
    for col,scale in enumerate([1.5,2]):
        for car,color in zip(LABELS,colors):
            for trial in [r for r in races['trials'] if r['config']['molecule']==car and r['config']['lattice_scale']==scale]:
                old=read(ROOT/'outputs/crystal-contrast-guarded/browser-trials'/(trial['key']+'.json'))
                for data,ls in [(old,'--'),(trial,'-')]:
                    t=[s['time_ps'] for s in data['samples']]
                    c=np.asarray([s['center_nm'] for s in data['samples']])
                    axes[0,col].plot(t,c[:,0],color=color,ls=ls,lw=.8,alpha=.65)
                    axes[1,col].plot(t,c[:,1]-c[0,1],color=color,ls=ls,lw=.8,alpha=.65)
            axes[0,col].plot([],[],color=color,label=LABELS[car])
        axes[0,col].set_title(f'Corrugation {scale:g}x');axes[0,col].set_ylabel('COM x / nm')
        axes[1,col].set_ylabel('Lateral displacement / nm')
        for ax in axes[:,col]:ax.set_xlabel('Time / ps');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Three seeds per car | dashed: periodic reference; solid: Fourier')
    for fmt in ['png','svg']:fig.savefig(OUT/f'trajectories.{fmt}',dpi=180)
    plt.close(fig)
    lines=['# 浏览器傅里叶模式实施与验收','',
      '状态：实验性可选模式已接入；周期势继续默认。0.5 fs时间步、UFF、Langevin、驱动及比赛判定未改变。','',
      '## 实现','',
      '- 以当前浏览器参数重新生成系数；沿用12个复谐波和高度/方向Hermite插值。全部系数与计算保持Float64。',
      '- 六相位在系数层合并：m为偶数且n为3的倍数的谐波保持，其余乘局部起伏系数，并保留高度导数。代理部分只插值一次。',
      '- 高度2.8–8 Å、方向坐标|nx/nz|和|ny/nz|≤2；边缘平滑混合，完整保留位置、法向和混合权重导数。域外及近退化法向回退完整16 Å周期势。',
      '- 浏览器只按需下载碳方向表，额外32.55 MB压缩下载、37.59 MB系数内存；罕见氢方向项使用完整回退。默认周期模式不下载这张表。开发/完整ZIP仍携带可选氢表。',
      '- 移除40 ms/100步的主动等待。每100步仍输出轨迹，并向事件循环让出控制；浏览器计时器调度仍有成本。','',
      '## 性能','',
      '以下同一WASM模块串行交替比较周期/傅里叶，每条件3次2 ps，使用相同初始相空间；100 ps成本为外推，排除准备、下载和显示。','',
      '| 分子 | 起伏 | 周期 / s | 傅里叶 / s | 加速 | 完全代理 / 混合 / 完全回退 |','|---|---:|---:|---:|---:|---|']
    for r in performance:lines.append(f"| {LABELS[r['car']]} | {r['scale']:g} | {r['periodic_s']:.2f} | {r['fourier_s']:.2f} | {r['speedup']:.2f}× | "+' / '.join(f'{v:.1%}' for v in r['coverage'])+' |')
    lines+=['','覆盖率来自每0.05 ps采样的方向相互作用中心计数，不是CPU时间占比。','', '![计算成本](performance.png)','',
      '真实Chrome Worker在无Three.js绘制的测试页运行吡啶20 nm完整比赛；以下实际耗时包含Worker消息和计时器调度，不包含下载与准备：','']
    for r in browser['cases']:lines.append(f"- {r['scale']:g}倍：{r['step']/2000:.2f} ps，实际{r['wall_seconds']:.2f}秒。")
    lines+=['','## 数值与动力学','',
      f"108组与原生六次查询实现对照，最大力差{summary['max_port_force_error']:.3e} kcal/mol/Å；12组全笛卡尔有限差分通过。",
      f"相对完整周期势，采样构型最大能量差{summary['max_proxy_energy_error']:.4g} kcal/mol、力分量差{summary['max_proxy_force_error']:.4g} kcal/mol/Å。这些是近似误差，不是数值移植误差，也不是所有构型的误差上界。",'',
      '四种分子2倍起伏的6×6平移扫描最大势能变化误差不超过0.00162 kcal/mol。16段2 ps NVE的代理能量漂移与各自周期参考接近，均通过参考漂移×1.1+1e-5 kcal/mol的回归门槛；这不表示有限步长能量误差为零。','',
      '24场完整比赛均检查准备收敛、键长、持续离面、终点穿越及保存坐标重算COM。对照为历史同条件周期势三种子比赛；其编译器算术与当前优化路径不保证逐位一致。','',
      '| 分子 | 起伏 | 周期均值±SD / ps | 傅里叶均值±SD / ps | 均值变化 |','|---|---:|---:|---:|---:|']
    for g in groups:lines.append(f"| {LABELS[g['car']]} | {g['scale']:g} | {g['periodic_mean_ps']:.2f}±{g['periodic_sd_ps']:.2f} | {g['fourier_mean_ps']:.2f}±{g['fourier_sd_ps']:.2f} | {g['finish_mean_change_percent']:+.1f}% |")
    lines+=['','各组均值变化约−1.9%到+4.6%。2倍起伏下，苯、吡啶、52原子车的样本标准差均高于历史周期对照；1.5倍下吡啶与52原子车均值顺序交换。当前不能把这些变化全部归因于近似误差，也不能宣称排行榜保持不变，故不提升为默认或混用正式成绩。','',
      '![三种子成绩](race-scores.png)','','![位移时间与横向游走](trajectories.png)','',
      '横向RMS、最大横向游走、1 ps低前进比例见summary.json。低前进定义为任意1 ps窗口前进不足0.02 nm，包含后退，不等同于严格的吸附陷阱判定。每组仅三种子，不能宣称已证明统计等效或任意排行榜名次不变。','',
      '## 当前结论与限制','',
      '小分子已有明显收益；较大平面结构仍可能超过半分钟；C60因多数方向项回退完整计算，仍需数分钟。下一步如果要继续加速C60，应针对法向与高度覆盖研究更合适的代理表示，而不是调大时间步。',
      f"关闭代理的108组固定构型与本轮修改前WASM对照通过，最大力差{reference['max_force_error']:.3e} kcal/mol/Å。真实Three.js页面已验证周期→傅里叶→周期切换、准备及按需加载，没有页面错误。",'',
      '首载额外势表仍较大，移动端、多浏览器、真实网络和完整Three.js页面的性能需要继续验收。本轮没有公开部署。保留周期势默认，傅里叶作为玩家可选的实验模式。','',
      '核心证据：numerics.json、benchmark.json、browser-validation.json、races/、summary.json。更新包：outputs/wasm-preparation/molecular-car-web-preview.zip。','']
    (OUT/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(dict(status='passed',performance=performance,groups=groups),indent=2))
if __name__=='__main__':main()
