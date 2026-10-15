"""Audit the candidate-specific trajectories and render a descriptive acceptance report.

No statistical-equivalence claim: three thermostat seeds share one prepared geometry
and velocity seed. Censored races remain in completion and endpoint summaries.
"""
from pathlib import Path
import gzip
import hashlib
import html
import json
import math
import statistics as st
import sys
import os

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/default-candidate-acceptance'
read = lambda p: json.loads(p.read_text(encoding='utf-8'))
names = {'pyridine-11': 'Pyridine', 'test-car-52': '52-atom car',
         'fullerene-c60': 'C60', 'folded-eicosane': 'Flexible C20H42'}


def describe(values):
    return {'mean': st.mean(values), 'sd': st.stdev(values) if len(values) > 1 else 0,
            'values': values}


def main():
    races, nve = read(OUT / 'races.json'), read(OUT / 'nve.json')
    assert races['status'] == nve['status'] == 'complete'
    assert len(races['rows']) == len(races['jobs']) == 24
    assert len(nve['rows']) == 8
    sources = {**races['sources'], **nve['sources']}
    for file, expected in sources.items():
        assert hashlib.sha256(Path(file).read_bytes()).hexdigest() == expected, file
    audit = {'status': 'running', 'sources_verified': len(sources), 'trials': [],
             'groups': [], 'nve': [], 'paired_initial_state_max_error': 0,
             'coordinate_frames': 0, 'max_com_error_nm': 0}
    keys = set()
    for r in races['rows']:
        key, c, samples = r['key'], r['config'], r['samples']
        assert key not in keys
        keys.add(key)
        assert r['status'] == 'complete', (key, r.get('error'))
        assert c['temperature_K'] == 150 and c['strength'] == .6 and c['lattice_scale'] == 1.5
        assert r['runtime']['preparation']['accepted']
        assert r['runtime']['surface']['proxy_charts'] == 3
        assert r['runtime']['timestep_fs'] == c['dt_fs']
        lift = peak = 0
        for i, s in enumerate(samples):
            assert s['step'] == i * 50 / c['dt_fs']
            assert abs(s['time_ps'] - i * .05) < 1e-9
            assert all(math.isfinite(v) for v in [*s['center_nm'], s['energy_kcal_mol'], s['internal_temperature_K']])
            if i:
                lift = lift + 1 if s['min_z_A'] > 8 else 0
                peak = max(peak, lift)
        assert abs(peak * .05 - r['maximum_lift_duration_ps']) < 1e-9
        assert r['competitive_valid'] == (peak < 100)
        if r['outcome'] == 'finish':
            assert samples[-2]['center_nm'][0] < 20 <= samples[-1]['center_nm'][0]
            assert r['finish_time_ps'] == samples[-1]['time_ps']
            assert not r['right_censored']
        else:
            assert r['right_censored'] and r['finish_time_ps'] is None
            assert samples[-1]['time_ps'] == 400 and samples[-1]['center_nm'][0] < 20
        coords = json.loads(gzip.decompress((OUT / 'races' / (key + '.json.coords.gz')).read_bytes()))
        masses = coords['masses']
        observed = {s['step']: s for s in samples}
        for step, p in zip(coords['steps'], coords['positions']):
            assert len(p) == 3 * len(masses) and step in observed
            com = [sum(m * p[3*i+d] for i, m in enumerate(masses))/sum(masses)/10 for d in range(3)]
            error = max(abs(a-b) for a, b in zip(com, observed[step]['center_nm']))
            audit['max_com_error_nm'] = max(audit['max_com_error_nm'], error)
            z = p[2::3]
            assert abs(min(z) - observed[step]['min_z_A']) < 1e-9
            assert abs(max(z) - observed[step]['max_z_A']) < 1e-9
            audit['coordinate_frames'] += 1
        assert audit['max_com_error_nm'] < 1e-9
        s50 = [s for s in samples if s['time_ps'] <= 50]
        assert s50[-1]['time_ps'] == 50
        y0 = samples[0]['center_nm'][1]
        rms = math.sqrt(st.mean((s['center_nm'][1]-y0)**2 for s in s50))
        slow = st.mean(s50[i]['center_nm'][0]-s50[i-20]['center_nm'][0] < .02 for i in range(20, len(s50)))
        audit['trials'].append({'key': key, 'car': c['car'], 'dt_fs': c['dt_fs'], 'seed': c['seed'],
            'finish_ps': r['finish_time_ps'], 'right_censored': r['right_censored'],
            'observed_ps': r['observed_time_ps'], 'endpoint_x_nm': samples[-1]['center_nm'][0],
            'rms_first50_nm': rms, 'slow_fraction_first50': slow,
            'temperature_after5_K': st.mean(s['internal_temperature_K'] for s in samples if s['time_ps'] >= 5),
            'maximum_lateral_excursion_nm': r['maximum_lateral_excursion_nm'],
            'wall_s': r['propagation_seconds'], 'seconds_per_100ps': r['propagation_seconds']*100/r['observed_time_ps'],
            'bond_range': r['bond_ratio_range'], 'max_lift_ps': r['maximum_lift_duration_ps'],
            'preparation_s': r['preparation_seconds'], 'competitive_valid': r['competitive_valid']})
    for car in names:
        for seed in [2718281, 1234567, 31415927]:
            pair = [r for r in races['rows'] if r['config']['car'] == car and r['config']['seed'] == seed]
            assert len(pair) == 2
            for field in ['positions', 'velocities']:
                err = max(abs(a-b) for a, b in zip(pair[0]['initial_phase_space'][field], pair[1]['initial_phase_space'][field]))
                audit['paired_initial_state_max_error'] = max(err, audit['paired_initial_state_max_error'])
                assert err < 1e-12
        for dt in [.5, 1]:
            rows = [r for r in audit['trials'] if r['car'] == car and r['dt_fs'] == dt]
            finish = [r['finish_ps'] for r in rows if not r['right_censored']]
            g = {'car': car, 'dt_fs': dt, 'n': len(rows), 'finish_count': len(finish),
                 'finish_ps': describe(finish) if len(finish) == len(rows) else None}
            for metric in ['endpoint_x_nm', 'rms_first50_nm', 'slow_fraction_first50', 'temperature_after5_K', 'maximum_lateral_excursion_nm', 'wall_s', 'seconds_per_100ps']:
                g[metric] = describe([r[metric] for r in rows])
            audit['groups'].append(g)
    for row in nve['rows']:
        assert row['status'] == 'complete' and row['stability_guard_pass']
        samples, dt = row['samples'], row['dt_fs']
        assert len(samples) == 2000/dt
        assert all(abs(s['time_ps']-(i+1)*dt/1000) < 1e-10 for i,s in enumerate(samples))
        peak = max(abs(s['energy_delta']) for s in samples)
        assert abs(peak-row['max_energy_error']) < 1e-10
        audit['nve'].append({k: row[k] for k in ['car','dt_fs','relative_to_initial_kinetic','bond_ratio_range',
            'energy_linear_slope_kcal_mol_ps','first_half_max_error','second_half_max_error']})
    audit['status'] = 'passed'
    audit['race_stability_pass'] = all(t['competitive_valid'] and t['bond_range'][0] >= .65 and t['bond_range'][1] <= 1.45 for t in audit['trials'])
    audit['finished_races'] = sum(not t['right_censored'] for t in audit['trials'])
    audit['precision_screen_1percent'] = {str(dt): sum(r['relative_to_initial_kinetic'] < .01 for r in audit['nve'] if r['dt_fs'] == dt) for dt in [.5, 1]}
    audit['comparisons'] = []
    for car in names:
        a, b = [next(g for g in audit['groups'] if g['car'] == car and g['dt_fs'] == dt) for dt in [.5, 1]]
        audit['comparisons'].append({'car': car,
            'finish_mean_delta_pct': 100*(b['finish_ps']['mean']/a['finish_ps']['mean']-1) if a['finish_ps'] and b['finish_ps'] else None,
            'runtime_per100ps_speedup': a['seconds_per_100ps']['mean']/b['seconds_per_100ps']['mean'],
            'rms_first50_delta_pct': 100*(b['rms_first50_nm']['mean']/a['rms_first50_nm']['mean']-1)})
    audit['scope'] = 'Data integrity and sampled stability audit; NOT statistical equivalence or universal molecule approval.'
    (OUT/'audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
    render(audit, races, nve)
    print(json.dumps({k:v for k,v in audit.items() if k not in ['trials','nve']}, ensure_ascii=False))


def render(audit, races, nve):
    sys.path.insert(0, str(ROOT / '.local/analysis-libs'))
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.local/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    colors = {.5:'#4575b4', 1:'#d95f02'}
    fig, axs = plt.subplots(2, 2, figsize=(11,7), constrained_layout=True)
    for ax, car in zip(axs.flat, names):
        for row in races['rows']:
            if row['config']['car'] != car: continue
            dt = row['config']['dt_fs']; s = row['samples']
            ax.plot([x['time_ps'] for x in s], [x['center_nm'][0] for x in s], color=colors[dt], alpha=.8,
                    label=f'{dt:g} fs' if row['config']['seed']==2718281 else None)
        ax.axhline(20,color='gray',ls=':',lw=1)
        ax.set(title=names[car], xlabel='Time (ps)', ylabel='COM x (nm)')
        ax.legend()
    fig.suptitle('v2 / 150 K / drive 0.6 / corrugation 1.5 / 3 thermostat seeds')
    fig.savefig(OUT/'trajectories.png',dpi=160); fig.savefig(OUT/'trajectories.svg'); plt.close(fig)
    fig, axs = plt.subplots(1, 3, figsize=(13,4), constrained_layout=True)
    for ax, metric, title in zip(axs, ['rms_first50_nm','seconds_per_100ps',None],
                                ['Lateral RMS, first 50 ps (nm)','Serial runtime (s / 100 ps)','NVE peak |E-E0| / initial KE (%)']):
        for i,car in enumerate(names):
            for dt,offset in [(.5,-.18),(1,.18)]:
                if metric:
                    g=next(g for g in audit['groups'] if g['car']==car and g['dt_fs']==dt)[metric]
                    val,err=g['mean'],g['sd']
                else:
                    val=100*next(g for g in audit['nve'] if g['car']==car and g['dt_fs']==dt)['relative_to_initial_kinetic']; err=0
                error_options = {'yerr': err, 'capsize': 3} if metric else {}
                ax.bar(i+offset,val,.34,color=colors[dt],label=f'{dt:g} fs' if i==0 else None,**error_options)
        ax.set_title(title);ax.set_xticks(range(4),names.values(),rotation=30,ha='right');ax.legend()
    fig.savefig(OUT/'comparison.png',dpi=160);fig.savefig(OUT/'comparison.svg');plt.close(fig)
    table=[]
    for g in audit['groups']:
        f=g['finish_ps']; ftext=f"{f['mean']:.2f} ± {f['sd']:.2f}" if f else '含 400 ps 截尾，不计算完赛均值'
        table.append(f"<tr><td>{names[g['car']]}</td><td>{g['dt_fs']:g}</td><td>{g['finish_count']}/3</td><td>{ftext}</td><td>{g['rms_first50_nm']['mean']:.3f}</td><td>{g['seconds_per_100ps']['mean']:.2f}</td><td>{g['wall_s']['mean']:.2f}</td></tr>")
    body='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>默认候选参数验收</title>
<style>body{font:16px/1.65 system-ui;max-width:1120px;margin:40px auto;padding:0 25px;color:#203040}table{border-collapse:collapse;width:100%}td,th{padding:9px;border-bottom:1px solid #ccc;text-align:left}img{max-width:100%}h1,h2{color:#174669}</style>
<h1>默认候选参数验收</h1><p>150 K · 驱动 0.6 · 代理势 v2 · 晶格起伏 1.5 · 比较 0.5 / 1 fs。赛道 20 × 20 nm，阻尼 350 fs。</p>
<p>4 种赛车 × 3 个热浴随机种子 × 2 种步长，共 24 场；最长 400 ps。全部采用 0.5 fs 周期势松弛，同组初始位置和速度逐项相同。比赛每 50 fs 采样，每 500 fs 保存坐标。另有每种分子两种步长的 2 ps NVE，每个积分步检查能量。未更改游戏默认配置。</p>
<p>这是有限样本的稳定性与玩法验收，不是统计等效证明。种子改变热浴随机力，初始速度种子固定；不同步长下相同种子不代表相同的连续时间随机力路径。超时按右截尾保留。</p>
<p>沿用诊断界限：采样键长相对初始值应在 0.65–1.45 内；全部原子高于表面 8 Å 持续 5 ps 视为异常离面。NVE 的 1% 初始动能误差是保守精度筛选值，不是崩溃边界。稳定性通过不代表精度筛选通过，也不代表任意外部分子都可靠。</p>
<table><tr><th>赛车</th><th>步长 fs</th><th>完赛</th><th>完赛 ps（均值 ± SD）</th><th>前 50 ps 横向 RMS nm</th><th>串行秒 / 100 ps</th><th>每场计算秒（均值）</th></tr>'''+''.join(table)+'''</table>
<p>耗时来自完整 Node/WASM 比赛，按实际模拟时长归一化；不含下载、准备或绘制。柱图误差线为 3 种子的标准差；NVE 图为单个初始状态的逐步峰值，无误差线。浏览器 Worker 实测另见 browser-v2-races.json。</p>
<p>编译运行时和 v2 势表已核对与此前受力、有限差分验证的版本完全一致，本次不重复扫描未变更的求力核。工作站的后台负载未受控；串行计时不代表其他设备的性能保证。</p>
<h2>位移与数值诊断</h2><img src="trajectories.png"><img src="comparison.png">
<p>原始结果：<a href="races.json">完整比赛</a> · <a href="nve.json">逐步能量</a> · <a href="audit.json">独立坐标重算与完整性检查</a> · <a href="browser-v2-races.json">真实浏览器</a></p></html>'''
    (OUT/'report.html').write_text(body,encoding='utf-8')


if __name__ == '__main__':
    main()
