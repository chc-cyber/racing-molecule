"""Native v2 regression: real preparation, force parity, cost, short dynamics."""
from pathlib import Path
import os
import sys
import json
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from molecular_car.engine import Simulation, DEFAULT_SETTINGS, DEFAULT_DRIVE, TIMESTEP_FS
from molecular_car.fourier_atlas import FourierAtlasPotential
from molecular_car.lattice_corrugation import LatticeCorrugation
from molecular_car.relaxation import JointRelaxationError


def main():
    out = ROOT/'outputs/native-atlas-acceleration'
    out.mkdir(parents=True, exist_ok=True)
    inputs = [ROOT/'examples'/name for name in ('benzene.mol', 'tetraphenylbenzene.mol')]
    inputs += [ROOT/'outputs/structure-race-scan/inputs'/name for name in
               ('fullerene-c60.mol', 'folded-eicosane.mol', 'coronene-36.mol', 'ethynyl-rigid.mol')]
    inputs += [ROOT/'outputs/structure-race-polar/inputs'/name for name in
               ('pyridine-11.mol', 'acetone-10.mol', 'benzonitrile-13.mol')]
    report = dict(settings=DEFAULT_SETTINGS, drive=DEFAULT_DRIVE, timestep_fs=TIMESTEP_FS,
                  cost_scope='warmed surface force only; median of 5 calls, same relaxed coordinates',
                  molecules=[], paired_dynamics=[])

    def save():
        (out/'report.json').write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')

    for path in inputs:
        row = dict(molecule=path.name)
        sim = None
        started = time.perf_counter()
        try:
            sim = Simulation(path.read_text(), path.name)
            row.update(ready_seconds=time.perf_counter()-started, atoms=sim.n,
                       preparation=sim.runtime['preparation'], surface=sim.surface.metadata())
            reference = LatticeCorrugation(FourierAtlasPotential(
                sim.prepared.atoms, sim.prepared.positions, bonds=sim.prepared.bonds,
                normal_policy='bonded-saip-game-v1', accelerator='numpy'), 1.5)
            x = sim.coordinates()
            a = sim.surface.energy_and_forces(x); b = reference.energy_and_forces(x)
            row.update(energy_difference_kcal_mol=abs(a[0]-b[0]),
                       force_max_difference_kcal_mol_A=float(np.max(abs(a[1]-b[1]))))
            assert row['energy_difference_kcal_mol'] < 1e-8
            assert row['force_max_difference_kcal_mol_A'] < 1e-8
            for name, surface in [('compiled', sim.surface), ('python_reference', reference)]:
                samples = []
                for _ in range(5):
                    t = time.perf_counter(); surface.energy_and_forces(x)
                    samples.append((time.perf_counter()-t)*1000)
                row[name+'_force_ms'] = float(np.median(samples))
            row['surface_speedup'] = row['python_reference_force_ms']/row['compiled_force_ms']
            started = time.perf_counter()
            frame = sim.advance(1000)
            row.update(one_ps_seconds=time.perf_counter()-started, one_ps_frame=frame)
            assert np.isfinite(frame['positions']).all() and np.isfinite(frame['energy'])
        except JointRelaxationError as error:
            row.update(ready_seconds=time.perf_counter()-started, preparation=error.diagnostics)
        finally:
            if sim is not None:
                sim.close()
        report['molecules'].append(row); save()
        print(json.dumps({k:v for k,v in row.items() if k not in ('surface','one_ps_frame','preparation')},
                         ensure_ascii=False), flush=True)
        print('preparation accepted:', row['preparation']['accepted'], flush=True)

    # Identical prepared coordinates, velocities and random seed; replace only
    # the force evaluator. This checks both CH and hetero-Au Langevin paths.
    for path in (inputs[0], inputs[6]):
        trajectories = []
        for reference in (False, True):
            sim = Simulation(path.read_text(), path.name)
            try:
                initial = sim.frame()
                if reference:
                    sim.surface = LatticeCorrugation(FourierAtlasPotential(sim.prepared.atoms,
                        sim.prepared.positions, bonds=sim.prepared.bonds,
                        normal_policy='bonded-saip-game-v1', accelerator='numpy'), 1.5)
                trajectories.append([initial]+[sim.advance(10) for _ in range(10)])
            finally:
                sim.close()
        a, b = trajectories
        row = dict(molecule=path.name, duration_fs=100,
                   max_coordinate_difference_A=max(float(np.max(abs(np.array(x['positions'])-y['positions']))) for x,y in zip(a,b)),
                   max_energy_difference_kcal_mol=max(abs(x['energy']-y['energy']) for x,y in zip(a,b)))
        assert row['max_coordinate_difference_A'] < 1e-7
        assert row['max_energy_difference_kcal_mol'] < 1e-7
        report['paired_dynamics'].append(row); save(); print(json.dumps(row), flush=True)
    assert all(row['preparation']['accepted'] for row in report['molecules'])
    report['passed'] = True
    save()


if __name__ == '__main__':
    main()
