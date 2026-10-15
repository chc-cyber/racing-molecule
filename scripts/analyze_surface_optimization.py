"""Recheck saved optimization evidence and summarize the serial benchmark."""
from pathlib import Path
import gzip
import hashlib
import json
import os
import statistics
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '.local/analysis-libs'))
os.environ['MPLCONFIGDIR'] = str(ROOT / '.local/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = ROOT / 'outputs/surface-kernel-optimization'
OLD = ROOT / 'outputs/crystal-contrast-guarded'
NAMES = {'benzene-12': 'Benzene', 'coronene-36': 'Coronene',
         'ethynyl-rigid': 'Ethynyl-linked', 'test-car-52': 'Tetraphenylbenzene',
         'fullerene-c60': 'C60', 'pyridine-11': 'Pyridine', 'benzonitrile-13': 'Benzonitrile'}


def read(p):
    return json.loads(p.read_text(encoding='utf-8'))


def check_sources(mapping):
    for f, h in mapping.items():
        assert hashlib.sha256((ROOT / f).read_bytes()).hexdigest() == h, f


def main():
    numerical, benchmark = [read(OUT / f) for f in ['numerical-validation.json', 'benchmark.json']]
    for report, key in [(numerical, 'sources'), (benchmark, 'source_sha256')]:
        assert report['status'] == 'passed' and report['source_integrity'] == 'passed'
        check_sources(report[key])
    element = read(OUT / 'element-validation.json')
    assert element['status'] == 'passed' and len(element['cases']) == 280
    check_sources(element['sources'])
    assert read(OUT / 'browser-validation.json')['passed']
    assert read(ROOT / 'outputs/wasm-preparation/package-validation.json')['passed']
    assert len(benchmark['md']) == 84 and len(benchmark['static_edges']) == 35
    performance = []
    for car in NAMES:
        for scale in [1.5, 2]:
            rows = [r for r in benchmark['md'] if r['car'] == car + '.mol' and r['scale'] == scale]
            times = [statistics.median(r['seconds_per_100ps'] for r in rows if r['level'] == level) for level in [0, 3]]
            performance.append(dict(car=car, scale=scale, old_s_per_100ps=times[0], new_s_per_100ps=times[1],
                                    speedup=times[0]/times[1], reduction_percent=100*(1-times[1]/times[0])))
    races = read(OUT / 'races/browser-confirm-report.json')
    assert races['status'] == 'complete' and races['planned'] == len(races['trials']) == 10
    for f, h in races['source_sha256'].items():
        source = OUT / 'races' / f if f.startswith('inputs/') else ROOT / f
        assert hashlib.sha256(source.read_bytes()).hexdigest() == h, f
    expected = {(c, s) for c in ['benzene-12', 'pyridine-11', 'test-car-52', 'fullerene-c60'] for s in [1.5, 2]}
    assert {(r['config']['molecule'], r['config']['lattice_scale']) for r in races['trials']} == expected
    expected_jobs = {(c, s, 2718281) for c, s in expected} | {('benzene-12', 2, seed) for seed in [1234567, 31415927]}
    assert {(r['config']['molecule'], r['config']['lattice_scale'], r['config']['seed']) for r in races['trials']} == expected_jobs
    comparison = []
    for r in races['trials']:
        assert r['status'] == 'passed' and r['competitive_valid'] and r['outcome'] == 'finish'
        old = read(OLD / 'browser-trials' / (r['key'] + '.json'))
        assert r['config'] == old['config']
        assert r['source_sha256'] == races['source_sha256']
        assert r['samples'][-1]['center_nm'][0] >= 20
        assert r['samples'][-2]['center_nm'][0] < 20
        assert r['finish_time_ps'] == r['samples'][-1]['step']*.0005
        assert all(s['step'] == i*100 for i, s in enumerate(r['samples']))
        assert r['runtime']['preparation']['accepted']
        assert r['runtime']['preparation']['global_force_norm_kcal_mol_A'] <= 1e-4
        samples = np.array([s['center_nm'] for s in r['samples']])
        before = np.array([s['center_nm'] for s in old['samples']])
        length = min(len(samples), len(before))
        trajectory = json.loads(gzip.decompress((OUT / 'races/browser-trials' / r['coordinates_file']).read_bytes()))
        p = np.asarray(trajectory['positions']).reshape(-1, r['atoms'], 3)
        m = np.asarray(trajectory['masses'])
        com = np.sum(p*m[None, :, None], axis=1)/sum(m)/10
        saved = np.asarray([r['samples'][step//100]['center_nm'] for step in trajectory['steps']])
        com_error = float(np.max(np.abs(com-saved)))
        assert com_error < 1e-10
        previous_seeds = [read(p)['finish_time_ps'] for p in (OLD / 'browser-trials').glob(
            f"browser-confirm-{r['config']['molecule']}-T150-S0.6-L{r['config']['lattice_scale']:g}-G1-seed*-400ps.json")]
        assert len(previous_seeds) == 3 and all(v is not None for v in previous_seeds)
        comparison.append(dict(car=r['config']['molecule'], scale=r['config']['lattice_scale'], seed=r['config']['seed'],
            old_finish_ps=old['finish_time_ps'], new_finish_ps=r['finish_time_ps'],
            finish_delta_ps=r['finish_time_ps']-old['finish_time_ps'],
            max_center_difference_nm=float(np.max(np.abs(samples[:length]-before[:length]))),
            center_difference_2ps_nm=float(np.max(np.abs(samples[:41]-before[:41]))),
            old_lateral_nm=old['maximum_lateral_excursion_nm'], new_lateral_nm=r['maximum_lateral_excursion_nm'],
            old_three_seed_finish_range_ps=[min(previous_seeds), max(previous_seeds)],
            com_reconstruction_error_nm=com_error, bond_ratio_range=r['bond_ratio_range'],
            sustained_lift=r['sustained_lift']))
    sensitive = [r for r in comparison if r['car'] == 'benzene-12' and r['scale'] == 2]
    sensitivity = {phase: {'mean_ps': statistics.mean(r[phase+'_finish_ps'] for r in sensitive),
                           'sd_ps': statistics.stdev(r[phase+'_finish_ps'] for r in sensitive)} for phase in ['old', 'new']}
    result = dict(status='passed', hardware=benchmark['hardware'], performance=performance, races=comparison,
        benzene_2x_three_seed=sensitivity,
        static_cases=len(numerical['static']),
        element_static_cases=len(element['cases']), elements=element['elements'],
        max_energy_error_kcal_mol=max(r['energy_error'] for r in numerical['static']),
        max_force_error_kcal_mol_A=max(r['force_error'] for r in numerical['static']),
        finite_difference_cases=len(numerical['finite_difference']), paired_nve=numerical['paired_nve'],
        limits=['100 ps costs extrapolate serial 2 ps intervals; not complete race timings.',
                'Node executes the real WASM; excludes browser download, preparation, rendering and playback pacing.',
                'Eight base races use one Langevin seed; benzene 2x adds two seeds; not a population-level equivalence proof.',
                'Concurrent full-race wall times are not single-player performance benchmarks.',
                'The browser Worker retains a 40 ms target per 100 steps: nominal 80 s pacing per 100 ps.'])
    (OUT / 'summary.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, scale in zip(axes, [1.5, 2]):
        rows = [r for r in performance if r['scale'] == scale]
        y = np.arange(len(rows))
        ax.barh(y-.18, [r['old_s_per_100ps'] for r in rows], .34, label='Original algorithm', color='#7b8794')
        ax.barh(y+.18, [r['new_s_per_100ps'] for r in rows], .34, label='Optimized', color='#007e87')
        ax.set_yticks(y, [NAMES[r['car']] for r in rows]); ax.invert_yaxis()
        ax.axvline(30, color='#bc593d', linestyle='--', linewidth=1, label='30 s target')
        ax.set_title(f'Corrugation {scale:g}x'); ax.set_xlabel('Seconds per 100 ps (2 ps timing extrapolation)')
        ax.grid(axis='x', alpha=.2); ax.set_axisbelow(True)
    axes[0].legend(loc='lower right', fontsize=8)
    fig.suptitle('Serial WASM MD | 150 K | drive 0.6 | dt 0.5 fs | median of 3')
    for suffix in ['png', 'svg']: fig.savefig(OUT / ('performance.'+suffix), dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 4, figsize=(13, 6), constrained_layout=True, sharey=True)
    for column, car in enumerate(['benzene-12', 'pyridine-11', 'test-car-52', 'fullerene-c60']):
        for row, scale in enumerate([1.5, 2]):
            trial = next(r for r in races['trials'] if r['config']['molecule'] == car and r['config']['lattice_scale'] == scale and r['config']['seed'] == 2718281)
            old = read(OLD / 'browser-trials' / (trial['key']+'.json'))
            ax = axes[row, column]
            for data, label, style, color in [(old, 'Original', '-', '#7b8794'), (trial, 'Optimized', '--', '#007e87')]:
                ax.plot([s['time_ps'] for s in data['samples']], [s['center_nm'][0] for s in data['samples']],
                        style, color=color, label=label, linewidth=1.3)
            ax.set_title(f'{NAMES[car]} | {scale:g}x');ax.set_xlabel('Time / ps');ax.grid(alpha=.2)
            if column == 0: ax.set_ylabel('COM x / nm')
    axes[0, 0].legend(fontsize=8)
    fig.suptitle('Complete 20 nm races | same seed 2718281 | no change to timestep or dynamics')
    for suffix in ['png', 'svg']: fig.savefig(OUT / ('race-trajectories.'+suffix), dpi=180)
    nve_position = max(r['position_error_A'] for r in numerical['paired_nve'])
    nve_energy = max(r['energy_error'] for r in numerical['paired_nve'])
    lines = [
        '# 浏览器周期 Au 求力优化验收', '',
        '日期：2026-10-14。三项优化已实现并打包；本机总MD耗时降低约13%–19%，尚未达到所有赛车半分钟内完成计算的目标。', '',
        '## 实现范围', '',
        '1. 六相位合并：整车法向只计算一次，复用工作数组，合并法向能量梯度后只做一次链式求导。碳原子常见接触路径独立编译，减少分支与临时量。',
        '2. 对齐网格预平均：32×60网格上的半周期/三分之一周期移位是整数格点，先平均势能表再求同一插值函数及其导数。近表面保护门函数仍在查询后施加，保留其z导数。方向相关SAIP项仍完整计算六相位。',
        '3. 保守Au候选表：使用4×6横向分箱与1 Å高度分箱；闭合分箱到位点的最小距离可能小于16 Å才保留，求力时仍检查原精确截断。每个原子分别回卷查询，原分子坐标用于法向，长分子跨胞不会被当成单个最小镜像分子。超出候选表高度范围回退原遍历。', '',
        '未改变UFF、原子质量、势能参数、截断、步长0.5 fs、Langevin或赛车驱动。1倍起伏仍走原路径。此次改动只作用于浏览器C++/WASM；Python/Numba参考实现保留。', '',
        '## 串行性能', '',
        f"设备：{benchmark['hardware']['cpu']}；{benchmark['hardware']['node']}。每条件交替运行旧算法(level 0)和新版(level 3)，同初始相空间，3次2 ps样本取中位数，共84次计时。", '',
        '**下表100 ps秒数由2 ps计时线性折算，不是实际跑满100 ps，也不包括准备、下载、渲染或Worker播放节流。** 旧算法保留在新二进制中；数值另与归档的旧WASM独立对照。完整比赛复核使用并行任务，其耗时不用于本表。', '',
        '| 分子 | 1.5倍：旧→新 / s | 耗时下降 | 2倍：旧→新 / s | 耗时下降 |',
        '|---|---:|---:|---:|---:|']
    for car, name in NAMES.items():
        a, b = [next(r for r in performance if r['car'] == car and r['scale'] == s) for s in [1.5, 2]]
        lines.append(f"| {name} | {a['old_s_per_100ps']:.2f} → {a['new_s_per_100ps']:.2f} | {a['reduction_percent']:.1f}% | {b['old_s_per_100ps']:.2f} → {b['new_s_per_100ps']:.2f} | {b['reduction_percent']:.1f}% |")
    lines += ['', '![串行计算成本](performance.png)', '',
        '固定构型分项计时中，新版表面核仍占UFF+表面+驱动三个求力内核时间的约95%–97%（苯/C60）；这不是整个程序的CPU采样占比。六相位的方向势指数计算仍是主要成本，继续单独优化UFF不足以达到大分子半分钟目标。', '',
        '## 数值与比赛复核', '',
        f"- 498组原/新静态检查：表面能最大差 {result['max_energy_error_kcal_mol']:.3e} kcal/mol，表面力最大分量差 {result['max_force_error_kcal_mol_A']:.3e} kcal/mol/Å；UFF能量和力逐值一致。覆盖倾斜结构、保护区/缓存边界、周期分箱边界与多胞平移。",
        '- 另有280组、12种元素的优化路径检查，包含0、1、1.5、2、12倍起伏；原有224组WASM通用回归也通过。',
        '- 5组全笛卡尔有限差分验力通过。细步长最大误差2.705e-6 kcal/mol/Å；插值结点附近的有限差分收敛趋势与旧版一致。',
        f'- 16组2 ps成对NVE：最大坐标差 {nve_position:.3e} Å，最大能量差 {nve_energy:.3e} kcal/mol。',
        '- 真实Chrome 155 Worker完成1.5/2倍起伏的两段2 ps冒烟检查；准备门槛、模块加载、势表校验、积分与暂停均通过。完整10场比赛是在Node中执行同一个实际WASM。', '',
        '完整20 nm比赛：150 K、驱动0.6、场增长1、速度种子17329，自动准备；同旧版逐场配对。8场基础对照使用Langevin种子2718281；苯2倍再补两个种子。成绩按0.05 ps采样判定。', '',
        '| 分子 | 起伏 | 随机种子 | 旧成绩 / ps | 新成绩 / ps | 差 / ps | 最大COM分量差 / nm |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for r in sorted(comparison, key=lambda r: (r['car'], r['scale'], r['seed'])):
        lines.append(f"| {NAMES[r['car']]} | {r['scale']:g} | {r['seed']} | {r['old_finish_ps']:.2f} | {r['new_finish_ps']:.2f} | {r['finish_delta_ps']:+.2f} | {r['max_center_difference_nm']:.3e} |")
    lines += ['', '![完整比赛位移](race-trajectories.png)', '',
        f"最大配对成绩变化为 {max(abs(r['finish_delta_ps']) for r in comparison):.2f} ps；因此不能称为逐场成绩或逐帧轨迹完全一致。完整链路还包含数值松弛，准备阶段的微小差异和后续舍入可以在长程非线性运动中放大。前述固定相空间NVE对照用于分开检查求力实现本身。", '',
        f"对单次变化较大的苯2倍补充相同三种子后，旧成绩均值±样本SD为 {sensitivity['old']['mean_ps']:.2f}±{sensitivity['old']['sd_ps']:.2f} ps，新版为 {sensitivity['new']['mean_ps']:.2f}±{sensitivity['new']['sd_ps']:.2f} ps。三种子仅提供描述性对照，不能据此宣称通过统计等效性检验。", '',
        '全部10场通过准备门槛、有限数值、键长和持续离面检查，并由保存坐标独立重算COM。除苯2倍外，每条件仍只有一个共同种子；这是实现回归，不是种群层面的统计等效性证明。长时间随机非线性轨迹会放大浮点舍入，不能承诺任意设备逐帧相同。', '',
        '## 半分钟目标与交付', '',
        '当前Worker每100步维持约40 ms的播放节奏，100 ps名义播放下限约80秒；本轮未改变。后续应把计算批次与渲染节奏解耦，同时保留物理步长和轨迹采样。仅取消播放节流仍不足以让大分子达到30秒。', '',
        '再大幅降低方向相关表面项的成本，需要另行评估此前第四项的预计算/代理方案及误差验收。本轮没有启用该近似，也没有增加时间步。', '',
        '已更新`wasm/dist/`与`outputs/wasm-preparation/molecular-car-web-preview.zip`（约10.85 MB），55文件清单、源哈希、势表哈希和ZIP校验均通过。WASM增加13,299字节；首次非单位起伏额外生成约25 MB平均势表，另需候选表内存。尚未做手机与多浏览器性能矩阵。', '',
        '证据：`benchmark.json`、`numerical-validation.json`、`element-validation.json`、`browser-validation.json`、`summary.json`及`races/`。优化前源码与WASM在`reference/`；旧147场记录维持原哈希，未冒充新版重新运行结果。第一轮实验在`iteration-1/`，最终结论使用当前文件与源码哈希。', '']
    (OUT / 'REPORT.md').write_text('\n'.join(lines), encoding='utf-8')
    print(json.dumps({k: result[k] for k in ['status', 'static_cases', 'max_energy_error_kcal_mol', 'max_force_error_kcal_mol_A']}, indent=2))
    print(json.dumps(performance, indent=2))
    print(json.dumps(comparison, indent=2))


if __name__ == '__main__':
    main()
