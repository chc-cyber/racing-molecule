"""Finalize only after race audit, preparation comparison and real browser runs."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/default-candidate-acceptance'
read = lambda f: json.loads((OUT/f).read_text(encoding='utf-8'))
audit, browser, prep, force = [read(f) for f in ['audit.json','browser-v2-races.json','preparation.json','reused-force-validation.json']]
environment = read('environment.json')
assert audit['status'] == 'passed' and audit['race_stability_pass']
assert browser['passed'] and len(browser['cases']) == 3
assert prep['status'] == 'complete' and len(prep['rows']) == 4
assert force['status'] == 'passed'
for file,h in prep['sources'].items():
    assert hashlib.sha256(Path(file).read_bytes()).hexdigest() == h
for f in ['races.json','nve.json']:
    for file,h in read(f)['sources'].items():
        assert hashlib.sha256(Path(file).read_bytes()).hexdigest() == h
summary={'configuration':{'temperature_K':150,'strength':.6,'lattice_scale':1.5,'proxy':'v2','candidate_dt_fs':1,'reference_dt_fs':.5},
    'production_defaults_changed':False,'race_stability_pass':True,'race_count':24,
    'finished_races':audit['finished_races'],'nve_count':8,'browser_count':3,
    'dynamical_equivalence':'not established; n=3 thermostat seeds, one geometry and initial velocity draw per molecule',
    'precision_screen_1percent':audit['precision_screen_1percent'],
    'comparisons':audit['comparisons'],
    'max_preparation_component_difference_A':max(r['maximum_component_difference_A'] for r in prep['rows']),
    'coordinate_frames_audited':audit['coordinate_frames'],'max_com_error_nm':audit['max_com_error_nm'],
    'browser':[{**{k:r[k] for k in ['car','dt','scale','wall_seconds']},'finish_ps':r['step']*r['dt']/1000} for r in browser['cases']],
    'scope':'150 K, strength 0.6, 20 nm track, 400 ps race observation cap; no claim for arbitrary molecules, other temperatures or long-term trajectories'}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
browser_rows=''.join(f"<tr><td>{r['car']}</td><td>{r['dt']:g}</td><td>{r['finish_ps']:.2f}</td><td>{r['wall_seconds']:.2f}</td></tr>" for r in summary['browser'])
prep_rows=''.join(f"<tr><td>{r['car']}</td><td>{r['maximum_component_difference_A']:.6g}</td><td>{r['atom_rms_displacement_A']:.6g}</td></tr>" for r in prep['rows'])
section=f'''<h2>验收结果与适用范围</h2>
<p>24 场比赛均通过采样稳定性检查，其中 {audit['finished_races']} 场在 400 ps 内完赛。8 组逐步 NVE 均有限且通过结构检查。保守的 1% 初始动能精度筛选：0.5 fs 通过 {audit['precision_screen_1percent']['0.5']}/4，1 fs 通过 {audit['precision_screen_1percent']['1']}/4。游戏稳定性与精度筛选是不同结论。</p>
<p>支持将 1 fs 作为性能与玩法取舍的游戏候选，不宣称与 0.5 fs 动力学等价。平均成绩接近也不保证游走、滞留、成绩分布一致。本轮没有修改游戏默认参数，也没有修改势能或摩擦模型。</p>
<p>保存坐标共复核 {audit['coordinate_frames']} 帧；质心重算最大差异 {audit['max_com_error_nm']:.3g} nm；成对初始位置和速度最大差异 {audit['paired_initial_state_max_error']:.3g}。采样检查不能排除未保存时刻的瞬态异常；NVE 逐步检查仅覆盖 2 ps。</p>
<h2>真实 Chrome Worker 完整比赛</h2><p>测试设备：{environment['cpu']}；Node {environment['node']}；Chrome {browser['browser']}。生产 Worker 与分子文件，默认模型种子。包含调度消息开销；排除首次下载、准备和 Three.js 绘制。</p>
<table><tr><th>赛车</th><th>步长 fs</th><th>完赛 ps</th><th>实际秒</th></tr>{browser_rows}</table>
<h2>直接 v2 准备路径核对</h2><p>四种分子的两条准备路径均收敛并被接受。下表比较实际 v2 准备与成对试验所用周期势准备；非零差异仍可能随长时间混沌动力学放大，不表示逐轨迹等价。</p>
<table><tr><th>赛车</th><th>最大坐标分量差 Å</th><th>原子 RMS 位移 Å</th></tr>{prep_rows}</table>
<p><a href="summary.json">验收汇总</a> · <a href="preparation.json">准备路径原始数据</a> · <a href="manifest.json">复现文件校验清单</a></p>'''
page=OUT/'report.html'
body=page.read_text(encoding='utf-8')
assert '<h2>验收结果与适用范围</h2>' not in body, 'Regenerate the analysis report before finalizing again'
page.write_text(body.replace('</html>',section+'</html>'),encoding='utf-8')
files=[f for f in OUT.rglob('*') if f.is_file() and f.name!='manifest.json']
files += [ROOT/f for f in ['wasm/tools/acceptance_support.mjs','wasm/tools/acceptance_nve.mjs',
    'wasm/tools/acceptance_races.mjs','wasm/tools/acceptance_preparation.mjs',
    'scripts/analyze_candidate_acceptance.py','scripts/finalize_candidate_acceptance.py','scripts/verify_candidate_browser.cjs']]
files += [ROOT/'wasm/dist/browser/worker.mjs', ROOT/'frontend/browser/worker.mjs']
manifest={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
(OUT/'manifest.json').write_text(json.dumps({'status':'complete','files':manifest},indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False))
