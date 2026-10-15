"""Read-only independent audit of the saved 1.5x supplement; writes its own report only."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/crystal-contrast-guarded'
CARS = ('benzene-12', 'coronene-36', 'ethynyl-rigid', 'test-car-52',
        'fullerene-c60', 'pyridine-11', 'benzonitrile-13')
SEEDS = (2718281, 1234567, 31415927)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def identity(config):
    return (config['molecule'], config['growth'], config['seed'])


def main():
    target = OUT / '1p5-independent-verification.json'
    expected = {(car, growth, seed) for car in CARS for growth in (1., .25) for seed in SEEDS}
    report = {'status': 'running', 'expected_trials': 42, 'checked_trials': 0,
              'checks': [], 'input_sha256': {}, 'source_integrity': {'status': 'running'},
              'independent_verifier_sha256': digest(Path(__file__)),
              'sampling': {'integration_timestep_fs': .5, 'score_contact_sample_ps': .05,
                           'full_coordinate_sample_ps': .5, 'terminal_coordinates_saved': True,
                           'limitation': 'Recomputed evidence covers saved frames, not every integration step. Bond bounds at intervening score frames are scanner-recorded guards.'}}
    locked = {}
    try:
        native_path = OUT / '1p5-numerics.json'
        wasm_path = OUT / '1p5-wasm-numerics.json'
        native, wasm = read(native_path), read(wasm_path)
        assert native['status'] == wasm['status'] == 'passed'
        assert native['source_integrity']['status'] == wasm['source_integrity']['status'] == 'passed'
        assert native['source_integrity']['changed_files'] == wasm['source_integrity']['changed_files'] == []
        assert wasm['native_reference']['sha256'] == digest(native_path)
        assert wasm['native_reference']['source_sha256'] == native['source_sha256']
        assert wasm['native_reference']['source_integrity'] == native['source_integrity']
        for name, sha in native['source_sha256'].items():
            path = Path(native['resolved_source_paths'][name])
            assert digest(path) == sha, f'Stale native numerical dependency: {name}'
            locked[str(path)] = sha
        for name, sha in wasm['source_sha256'].items():
            path = ROOT / name
            assert digest(path) == sha, f'Stale WASM numerical dependency: {name}'
            locked[str(path)] = sha
        for path in (native_path, wasm_path, Path(__file__)):
            locked[str(path)] = digest(path)
        report['numerical_reference'] = {
            'native': {'path': str(native_path.relative_to(ROOT)), 'sha256': digest(native_path)},
            'wasm': {'path': str(wasm_path.relative_to(ROOT)), 'sha256': digest(wasm_path)},
            'source_reference_match': True,
            'native_current_dependencies_verified': len(native['source_sha256']),
            'wasm_current_dependencies_verified': len(wasm['source_sha256'])}

        progress = read(OUT / 'browser-confirm-report.json')
        reported = {identity(row['config']): row for row in progress['trials']
                    if row['config']['lattice_scale'] == 1.5}
        rows = []
        for path in sorted((OUT / 'browser-trials').glob('*.json')):
            row = read(path)
            if row['config']['lattice_scale'] == 1.5:
                rows.append((path, row))
        seen = set()
        for path, row in rows:
            config = row['config']
            key = identity(config)
            assert key in expected and key not in seen, key
            seen.add(key)
            assert config['phase'] == 'browser-confirm'
            assert config['temperature_K'] == 150 and config['strength'] == .6 and config['limit_ps'] == 400
            assert row['status'] == 'passed' and row['execution'] == 'actual-lammps-wasm'
            if key in reported:
                assert row == reported[key], f'Progress and saved trial differ: {row["key"]}'
            runtime = row['runtime']
            assert runtime['timestep_fs'] == .5 and runtime['langevin_seed'] == config['seed']
            assert runtime['velocity_seed'] == 17329
            assert runtime['surface']['lattice_scale'] == 1.5
            assert runtime['surface']['corrugation_control_version'] == 'six-phase-atomic-crystal-contrast-v3'
            prep = runtime['preparation']
            assert prep['accepted'] and prep['final_coordinate_force_verified']
            assert prep['global_force_tolerance_kcal_mol_A'] == 1e-4
            assert 0 <= prep['global_force_norm_kcal_mol_A'] <= 1e-4
            assert runtime['start_pose']['requested']['mode'] == 'auto'

            samples = row['samples']
            step = np.array([sample['step'] for sample in samples])
            time = np.array([sample['time_ps'] for sample in samples])
            center = np.array([sample['center_nm'] for sample in samples])
            height_min = np.array([sample['min_z_A'] for sample in samples])
            height_max = np.array([sample['max_z_A'] for sample in samples])
            energy = np.array([sample['energy_kcal_mol'] for sample in samples])
            assert np.array_equal(step, np.arange(len(samples)) * 100)
            assert np.allclose(time, step * .0005, rtol=0, atol=1e-12)
            assert all(np.isfinite(values).all() for values in (center, height_min, height_max, energy))
            assert np.array_equal(center[0], row['initial_center_nm'])
            assert abs(samples[0]['internal_temperature_K'] - 150) < 1e-10
            crossing = np.flatnonzero(center[:, 0] >= 20.)
            assert row['outcome'] == 'finish' and not row['right_censored']
            assert len(crossing) == 1 and crossing[0] == len(samples) - 1
            assert row['finish_time_ps'] == row['observed_time_ps'] == time[-1]
            assert np.allclose(row['finish_crossing_interval_ps'], time[-2:], rtol=0, atol=1e-12)
            assert row['final_x_nm'] == center[-1, 0]
            assert row['final_forward_nm'] == center[-1, 0] - center[0, 0]
            lateral = float(np.max(abs(center[:, 1] - center[0, 1])))
            assert abs(lateral - row['maximum_lateral_excursion_nm']) < 1e-12
            lift = maximum_lift = 0
            for value in height_min[1:]:
                lift = lift + 1 if value > 8 else 0
                maximum_lift = max(maximum_lift, lift)
            assert row['maximum_lift_duration_ps'] == maximum_lift * .05
            assert row['sustained_lift'] == (maximum_lift >= 100)
            assert row['competitive_valid'] == (maximum_lift < 100)
            assert row['competitive_valid']
            assert np.allclose(row['atom_height_range_A'], [height_min.min(), height_max.max()], rtol=0, atol=1e-12)

            coordinates_path = OUT / 'browser-trials' / row['coordinates_file']
            coordinates = json.loads(gzip.decompress(coordinates_path.read_bytes()))
            xyz = np.asarray(coordinates['positions']).reshape((-1, row['atoms'], 3))
            coordinate_steps = np.asarray(coordinates['steps'], dtype=int)
            expected_steps = sorted(set(range(0, int(step[-1]) + 1, 1000)) | {int(step[-1])})
            assert np.array_equal(coordinate_steps, expected_steps)
            mass = np.asarray(coordinates['masses'])
            assert xyz.shape[0] == len(coordinate_steps) and len(mass) == row['atoms']
            assert np.isfinite(xyz).all() and np.isfinite(mass).all() and (mass > 0).all()
            indices = coordinate_steps // 100
            coordinate_com = np.average(xyz, axis=1, weights=mass) / 10
            com_error = float(np.max(abs(coordinate_com - center[indices])))
            assert com_error < 1e-10
            assert np.allclose(xyz[:, :, 2].min(axis=1), height_min[indices], rtol=0, atol=1e-12)
            assert np.allclose(xyz[:, :, 2].max(axis=1), height_max[indices], rtol=0, atol=1e-12)
            bonds = np.asarray(coordinates['bonds'], dtype=int)
            assert bonds.shape[1] == 2 and bonds.min() >= 0 and bonds.max() < row['atoms']
            lengths = np.linalg.norm(xyz[:, bonds[:, 0]] - xyz[:, bonds[:, 1]], axis=2)
            assert np.all(lengths[0] > 0)
            ratios = lengths / lengths[0]
            assert .65 <= ratios.min() <= ratios.max() <= 1.45
            assert .65 <= row['bond_ratio_range'][0] <= ratios.min() + 1e-12
            assert ratios.max() <= row['bond_ratio_range'][1] + 1e-12 and row['bond_ratio_range'][1] <= 1.45

            for name, sha in row['source_sha256'].items():
                source_path = OUT / name if name.startswith('inputs/') else ROOT / name
                assert digest(source_path) == sha, f'Trial source changed: {name}'
                locked[str(source_path)] = sha
                if name.startswith('inputs/'):
                    assert native['source_sha256']['outputs/crystal-contrast-guarded/' + name] == sha
                elif name in wasm['source_sha256']:
                    assert wasm['source_sha256'][name] == sha
            for path_to_lock in (path, coordinates_path):
                sha = digest(path_to_lock)
                locked[str(path_to_lock)] = sha
                report['input_sha256'][str(path_to_lock.relative_to(ROOT)).replace('\\', '/')] = sha
            report['checks'].append({
                'key': row['key'], 'status': 'passed', 'molecule': config['molecule'],
                'growth': config['growth'], 'langevin_seed': config['seed'], 'velocity_seed': 17329,
                'integration_steps': int(step[-1]), 'finish_time_ps': float(time[-1]),
                'score_contact_frames': len(samples), 'coordinate_frames': len(xyz),
                'com_reconstruction_max_error_nm': com_error,
                'min_atom_z_A': float(height_min.min()), 'highest_whole_car_min_z_A': float(height_min.max()),
                'maximum_lift_duration_ps': maximum_lift * .05,
                'maximum_lateral_excursion_nm': lateral,
                'initial_COM_x_nm': float(center[0, 0]),
                'preparation_global_force_norm_kcal_mol_A': prep['global_force_norm_kcal_mol_A'],
                'saved_coordinate_bond_ratio_range': [float(ratios.min()), float(ratios.max())],
                'scanner_all_score_frame_bond_ratio_range': row['bond_ratio_range']})

        checks = report['checks']
        report['checked_trials'] = len(checks)
        report['missing_trials'] = [dict(zip(('molecule', 'growth', 'seed'), key)) for key in sorted(expected - seen)]
        report['progress_snapshot'] = {'status': progress['status'], 'planned': progress['planned'],
                                       'supplement_records': len(reported)}
        if checks:
            report['metrics'] = {
                'total_integration_steps': sum(row['integration_steps'] for row in checks),
                'total_simulated_time_ps': sum(row['finish_time_ps'] for row in checks),
                'minimum_atom_z_A': min(row['min_atom_z_A'] for row in checks),
                'highest_whole_car_min_z_A': max(row['highest_whole_car_min_z_A'] for row in checks),
                'maximum_lateral_excursion_nm': max(row['maximum_lateral_excursion_nm'] for row in checks),
                'initial_COM_x_range_nm': [min(row['initial_COM_x_nm'] for row in checks), max(row['initial_COM_x_nm'] for row in checks)],
                'maximum_preparation_global_force_norm_kcal_mol_A': max(row['preparation_global_force_norm_kcal_mol_A'] for row in checks),
                'maximum_sustained_lift_duration_ps': max(row['maximum_lift_duration_ps'] for row in checks),
                'saved_coordinate_bond_ratio_range': [min(row['saved_coordinate_bond_ratio_range'][0] for row in checks), max(row['saved_coordinate_bond_ratio_range'][1] for row in checks)],
                'scanner_all_score_frame_bond_ratio_range': [min(row['scanner_all_score_frame_bond_ratio_range'][0] for row in checks), max(row['scanner_all_score_frame_bond_ratio_range'][1] for row in checks)]}
        changed = [path for path, sha in locked.items() if digest(Path(path)) != sha]
        assert not changed, f'Audited files changed: {changed}'
        report['source_integrity'] = {'status': 'passed', 'files': len(locked), 'changed_files': changed,
                                      'check': 'Numerical references, sources, runtime, molecular inputs and saved trial records hashed before and after this independent audit.'}
        complete = len(checks) == 42 and set(reported) == expected and progress['status'] == 'complete'
        report['status'] = 'passed' if complete else 'incomplete'
    except Exception as error:
        report.update(status='failed', error=str(error))
        raise
    finally:
        temporary = target.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf8')
        temporary.replace(target)
    print(json.dumps({key: value for key, value in report.items() if key in ('status', 'expected_trials', 'checked_trials', 'metrics', 'missing_trials')}))


if __name__ == '__main__':
    main()
