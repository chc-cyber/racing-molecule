"""Actual preparation and early-MD comparison of rigid start placement modes."""
from __future__ import annotations
import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[key] = "1"
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from molecular_car.engine import Simulation
from molecular_car.molecule import prepare_molecule
from molecular_car.start_pose import DEFAULT_START_POSE


def write(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def rotated_mol(content, axis, degrees):
    mol = Chem.MolFromMolBlock(content, removeHs=False)
    conf = mol.GetConformer();positions = conf.GetPositions();center = positions.mean(axis=0)
    radians = np.deg2rad(degrees);c,s=np.cos(radians),np.sin(radians)
    matrix = np.array([[c,-s,0],[s,c,0],[0,0,1]]) if axis == 'z' else np.array([[c,0,s],[0,1,0],[-s,0,c]])
    for index, position in enumerate((positions-center) @ matrix.T+center):
        conf.SetAtomPosition(index,tuple(float(v) for v in position))
    return Chem.MolToMolBlock(mol)


def body_rotation(initial, current, indices):
    a,b=initial[indices],current[indices];a=a-a.mean(axis=0);b=b-b.mean(axis=0)
    u,_,vt=np.linalg.svd(a.T@b)
    if np.linalg.det(vt.T@u.T)<0:vt[-1]*=-1
    matrix=vt.T@u.T
    return {"body_rotation_deg":float(np.degrees(np.arccos(np.clip((np.trace(matrix)-1)*.5,-1,1)))),
            "yaw_rotation_deg":float(np.degrees(np.arctan2(matrix[1,0],matrix[0,0]))),
            "tilt_rotation_deg":float(np.degrees(np.arccos(np.clip(matrix[2,2],-1,1))))}


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/start-pose-validation')
    parser.add_argument('--seeds',type=int,nargs='+',default=[98245165,2718281])
    parser.add_argument('--steps',type=int,default=10000)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    source_files=['molecular_car/start_pose.py','molecular_car/molecule.py','molecular_car/engine.py','molecular_car/relaxation.py','molecular_car/preparation_bfgs.py']
    sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in source_files}
    report={"passed":False,"source_sha256":sources,"cases":[],"failures":[],
            "qualification":"Finite paired-seed early trajectories; body rotation is a best fit of fixed heavy-atom identities and is affected by internal flexibility. No claim of universally fastest starting orientation."}
    original=(ROOT/'examples/tetraphenylbenzene.mol').read_text(encoding='utf8')
    fullerene=(ROOT/'tests/data/connected-fullerene.mol').read_text(encoding='utf8')
    inputs=[('test-car-original',original),('test-car-sideways',rotated_mol(original,'z',90)),
            ('test-car-upright',rotated_mol(original,'y',90)),('connected-fullerene',fullerene)]
    poses=[]
    for name,content in inputs:
        (args.output/(name+'.mol')).write_text(content,encoding='utf8')
        for config in ({'mode':'preserve'},DEFAULT_START_POSE,{'mode':'manual','yaw_deg':90,'pitch_deg':-30,'roll_deg':45},
                       {'mode':'auto','axis_atoms':[1,2]}):
            prepared=prepare_molecule(content,name+'.mol',start_pose=config)
            poses.append({'filename':name+'.mol','content':content,'config':config,
                          'positions':prepared.positions.ravel().tolist(),'report':prepared.start_pose})
    write(args.output/'pose-fixtures.json',poses)
    started=time.perf_counter()
    for name,content in inputs:
        for mode in ('preserve','auto'):
            for seed in args.seeds:
                row={'input':name,'mode':mode,'seed':seed,'surface':'periodic-au'};simulation=None
                try:
                    begin=time.perf_counter()
                    simulation=Simulation(content,name+'.mol',random_seed=seed,start_pose={'mode':mode})
                    row.update(prepare_wall_seconds=time.perf_counter()-begin,
                               relaxation_wall_seconds=simulation.runtime['preparation']['wall_seconds'],
                               conservative_force_norm=simulation.runtime['preparation']['global_force_norm_kcal_mol_A'],
                               pose=simulation.runtime['start_pose'],atom_count=simulation.n)
                    assert row['conservative_force_norm']<=1e-4 and simulation.lmp.get_thermo('step')==0
                    assert simulation.lmp.extract_global('dt')==.5
                    initial=simulation.coordinates();initial_center=simulation.center()
                    reference=[i for i,a in enumerate(simulation.prepared.atoms) if a['element']!='H']
                    endpoints=np.array([(b['a']-1,b['b']-1) for b in simulation.prepared.bonds])
                    lengths=np.linalg.norm(initial[endpoints[:,0]]-initial[endpoints[:,1]],axis=1)
                    frames=[];snapshots=[];run_started=time.perf_counter()
                    for step in range(0,args.steps+1,100):
                        if step:simulation.lmp.command('run 100 post no');simulation._check_callback()
                        x=simulation.coordinates();assert np.isfinite(x).all()
                        ratios=np.linalg.norm(x[endpoints[:,0]]-x[endpoints[:,1]],axis=1)/lengths
                        assert ratios.min()>.65 and ratios.max()<1.45
                        sample={'step':step,'time_ps':step*.0005,'forward_displacement_A':float(simulation.center()[0]-initial_center[0]),
                                'lateral_displacement_A':float(simulation.center()[1]-initial_center[1]),
                                'minimum_height_A':float(x[:,2].min()),**body_rotation(initial,x,reference)}
                        frames.append(sample);snapshots.append(x)
                    np.savez_compressed(args.output/f'{name}-{mode}-{seed}.npz',coordinates=np.array(snapshots),initial=initial)
                    row.update(passed=True,frames=frames,propagation_wall_seconds=time.perf_counter()-run_started,
                               maximum_body_rotation_first_2ps_deg=max(f['body_rotation_deg'] for f in frames if f['time_ps']<=2),
                               maximum_tilt_first_2ps_deg=max(f['tilt_rotation_deg'] for f in frames if f['time_ps']<=2),
                               final_forward_displacement_A=frames[-1]['forward_displacement_A'])
                except Exception as error:
                    row.update(passed=False,error=str(error));report['failures'].append({key:row[key] for key in ('input','mode','seed','error')})
                finally:
                    if simulation:simulation.close()
                report['cases'].append(row);write(args.output/'report.json',report)
                print(f"{name} {mode} seed {seed}: {'PASS' if row['passed'] else row['error']}",flush=True)
    report['source_unchanged']=all(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==value for name,value in sources.items())
    report.update(passed=not report['failures'] and report['source_unchanged'],wall_seconds=time.perf_counter()-started)
    write(args.output/'report.json',report)
    if not report['passed']:raise SystemExit('Start-pose verification failed; see report.json')


if __name__=='__main__':main()
