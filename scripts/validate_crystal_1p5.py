"""Supplementary 1.5x crystal checks; never overwrite the original scan gates."""
from __future__ import annotations

import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '1'

import json
from pathlib import Path
import sys
import time

import numpy as np
import numba
from rdkit import rdBase

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
from validate_crystal_contrast import validation_sources as baseline_sources
from scan_crystal_contrast import CARS
from scan_structure_races import digest, write_json
from molecular_car.molecule import prepare_molecule
from molecular_car.periodic_au import create_periodic_surface, BONDED_NORMAL_POLICY
from molecular_car.lattice_corrugation import LatticeCorrugation
from molecular_car.external_fields import ExternalFieldEvaluator
from molecular_car.engine import Simulation


def validation_sources(out):
    paths = baseline_sources(out)
    paths['scripts/validate_crystal_1p5.py'] = Path(__file__).resolve()
    paths['outputs/crystal-contrast-guarded/numerics.json'] = out / 'numerics.json'
    return paths


def main():
    out = ROOT / 'outputs/crystal-contrast-guarded'
    target = out / '1p5-numerics.json'
    paths = validation_sources(out)
    hashes = {name: digest(path) for name, path in paths.items()}
    baseline = json.loads((out / 'numerics.json').read_text(encoding='utf8'))
    assert baseline['status'] == 'passed'
    assert baseline['source_integrity']['status'] == 'passed'
    for name, expected in baseline['source_sha256'].items():
        assert hashes[name] == expected, f'Baseline validation is stale: {name}'
    report = {
        'status': 'running', 'scope': '1.5x supplementary validation',
        'protocol': {'scale': 1.5, 'growths': [1., .25], 'strength': .6,
                     'nve_temperature_K': 150., 'nve_duration_ps': 2., 'timestep_fs': .5,
                     'static': 'All Cartesian total-energy finite differences on seven baseline fixed configurations',
                     'nve': 'Independent strict zero-drive preparation; frozen static drive; no thermostat'},
        'source_sha256': hashes,
        'resolved_source_paths': {name: str(path.resolve()) for name, path in paths.items()},
        'source_integrity': {'status': 'running', 'files': len(paths),
                             'check': 'Sources, original numerical reference, molecular inputs and native runtime hashed before and after'},
        'runtime_versions': {'python': sys.version.split()[0], 'rdkit': rdBase.rdkitVersion,
                             'numpy': np.__version__, 'numba': numba.__version__},
        'static': [], 'wasm_fixtures': [], 'nve': []}
    write_json(target, report)
    began = time.perf_counter()
    try:
        for name in CARS:
            path = out / 'inputs' / (name + '.mol')
            content = path.read_text(encoding='utf8')
            prepared = prepare_molecule(content, path.name, start_pose={'mode': 'auto'})
            base = create_periodic_surface(prepared.atoms, prepared.positions,
                    bonds=prepared.bonds, normal_policy=BONDED_NORMAL_POLICY)
            numpy_base = create_periodic_surface(prepared.atoms, prepared.positions,
                    bonds=prepared.bonds, normal_policy=BONDED_NORMAL_POLICY, accelerator='numpy')
            potential = LatticeCorrugation(base, 1.5)
            numpy_potential = LatticeCorrugation(numpy_base, 1.5)
            charges = np.array([atom['charge'] for atom in prepared.atoms])
            alphas = np.array([atom['polarizability'] for atom in prepared.atoms])
            fields = ExternalFieldEvaluator(charges, alphas, surface_enabled=False)
            prior = next(row for row in baseline['wasm_fixtures']
                         if row['filename'] == path.name and row['scale'] == 1.)
            positions = np.array(prior['positions']).reshape((-1, 3))
            for growth in (1., .25):
                def evaluate(x):
                    energy, force = potential.energy_and_forces(x)
                    flat = x.ravel().tolist()
                    energy += prepared.forcefield.CalcEnergy(flat)
                    force -= np.array(prepared.forcefield.CalcGrad(flat)).reshape(x.shape)
                    field_energy, field_force = fields.energy_and_forces(x, 1, .6, 6, 0, 7, 0., growth)
                    return float(energy + field_energy), force + field_force
                energy, force = evaluate(positions)
                step = 1e-5
                error = 0.
                for atom, axis in np.ndindex(positions.shape):
                    displacement = np.zeros_like(positions)
                    displacement[atom, axis] = step
                    finite_difference = -(evaluate(positions + displacement)[0] -
                                          evaluate(positions - displacement)[0]) / (2 * step)
                    error = max(error, abs(finite_difference - force[atom, axis]))
                surface_energy, surface_force = potential.energy_and_forces(positions)
                other_energy, other_force = numpy_potential.energy_and_forces(positions)
                periodic_energy, periodic_force = potential.energy_and_forces(
                    positions + [base.period[0] * 11, -base.period[1] * 13, 0])
                row = {'molecule': name, 'scale': 1.5, 'growth': growth,
                       'coordinates': int(positions.size), 'fd_step_A': step,
                       'max_total_force_fd_error': float(error),
                       'numpy_numba_energy_error': float(abs(other_energy - surface_energy)),
                       'numpy_numba_force_error': float(np.max(abs(other_force - surface_force))),
                       'periodic_energy_error': float(abs(periodic_energy - surface_energy)),
                       'periodic_force_error': float(np.max(abs(periodic_force - surface_force)))}
                assert error < 2e-5, row
                assert row['numpy_numba_force_error'] < 1e-8 and row['periodic_force_error'] < 1e-8, row
                report['static'].append(row)
                report['wasm_fixtures'].append({'filename': path.name, 'content': content,
                    'positions': positions.ravel().tolist(), 'scale': 1.5, 'growth': growth,
                    'strength': .6, 'energy': energy, 'forces': force.ravel().tolist()})
            write_json(target, report)
            print(json.dumps({'static_complete': name, 'max_force_error': max(
                row['max_total_force_fd_error'] for row in report['static'])}), flush=True)
        for name in ('benzene-12', 'test-car-52'):
            path = out / 'inputs' / (name + '.mol')
            for growth in (1., .25):
                simulation = Simulation(path.read_text(encoding='utf8'), path.name,
                    {'temperature': 150., 'lattice_scale': 1.5}, thermostat=False,
                    drive={'mode': 'gradient', 'strength': .6, 'growth': growth},
                    start_pose={'mode': 'auto'})
                try:
                    phase = {key: np.ctypeslib.as_array(simulation.lmp.gather_atoms(atom, 1, 3),
                                shape=(simulation.n * 3,)).copy().tolist()
                             for atom, key in (('x', 'positions'), ('v', 'velocities'))}
                    frames = [simulation.frame()]
                    for _ in range(40):
                        frames.append(simulation.advance(100))
                    energies = np.array([frame['energy'] for frame in frames])
                    row = {'molecule': name, 'scale': 1.5, 'growth': growth,
                        'timestep_fs': simulation.lmp.extract_global('dt'),
                        'steps': frames[-1]['step'], 'duration_ps': 2.,
                        'max_energy_deviation_kcal_mol': float(np.max(abs(energies - energies[0]))),
                        'final_energy_change_kcal_mol': float(energies[-1] - energies[0]),
                        'accepted': simulation.runtime['preparation']['accepted'],
                        'residual_force': simulation.runtime['preparation']['global_force_norm_kcal_mol_A'],
                        'phase_space': phase, 'frames': frames}
                    assert row['timestep_fs'] == .5 and row['steps'] == 4000, row
                    assert row['accepted'] and row['residual_force'] <= 1e-4, row
                    assert row['max_energy_deviation_kcal_mol'] < .5, row
                    report['nve'].append(row)
                    write_json(target, report)
                    print(json.dumps({key: value for key, value in row.items()
                                      if key not in ('phase_space', 'frames')}), flush=True)
                finally:
                    simulation.close()
        report['status'] = 'passed'
    except Exception as error:
        report.update(status='failed', error=str(error))
        raise
    finally:
        changed = [name for name, path in paths.items() if not path.is_file() or digest(path) != hashes[name]]
        report['source_integrity'].update(status='failed' if changed else 'passed', changed_files=changed)
        if changed:
            report.update(status='failed', error='Validation inputs changed: ' + ', '.join(changed))
        report['wall_seconds'] = time.perf_counter() - began
        write_json(target, report)
    assert report['status'] == 'passed', report.get('error')


if __name__ == '__main__':
    main()
