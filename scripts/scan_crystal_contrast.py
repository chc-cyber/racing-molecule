"""Seven-car actual LAMMPS scan, with fixed timestep, bath and neutral inputs."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import shutil
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from scan_structure_races import run_trial, digest, write_json, SOURCES
from molecular_car.periodic_au import parameter_file

CARS=('benzene-12','coronene-36','ethynyl-rigid','test-car-52','fullerene-c60','pyridine-11','benzonitrile-13')


def main():
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--jobs',type=Path)
    p.add_argument('--phase',default='screen')
    p.add_argument('--output',type=Path,default=ROOT/'outputs/crystal-contrast-guarded')
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--prepare-only',action='store_true')
    args=p.parse_args();out=args.output.resolve();(out/'inputs').mkdir(parents=True,exist_ok=True)
    fixtures=[]
    for name in CARS:
        old=ROOT/'scripts/data/crystal-race-inputs'/(name+'.mol')
        target=out/'inputs'/old.name
        shutil.copyfile(old,target)
        fixtures.append({'id':name,'file':target.name,'sha256':digest(target)})
    if args.prepare_only:
        write_json(out/'inputs/manifest.json',{'fixtures':fixtures})
        print(json.dumps({'prepared':len(fixtures),'output':str(out/'inputs')},ensure_ascii=False));return
    sources={f:digest(ROOT/f) for f in (*SOURCES,'molecular_car/lattice_corrugation.py',
        'molecular_car/hetero_au.py','molecular_car/hetero_au_kernel.py','molecular_car/data/hetero_au_v1.json',
        'scripts/scan_crystal_contrast.py')}
    sources['CHAu.ILP']=digest(parameter_file())
    sources.update({str((out/'inputs'/f['file']).relative_to(ROOT)):f['sha256'] for f in fixtures})
    jobs=json.loads(args.jobs.read_text(encoding='utf8')) if args.jobs else [
        {'phase':args.phase,'molecule':name,'temperature_K':150.,'strength':.6,'lattice_scale':scale,'growth':growth,
         'seed':98245165,'limit_ps':40.} for scale in (1.,3.,6.,12.) for growth in (1.,.25) for name in CARS]
    if not 1<=args.workers<=4:raise ValueError('Use 1..4 workers')
    for c in jobs:
        if (c['molecule'] not in CARS or not 0<c['temperature_K']<=1000 or not 0<=c['strength']<=1
            or not 0<=c['lattice_scale']<=12 or not .1<=c['growth']<=2 or not 0<c['limit_ps']<=1000):
            raise ValueError('Invalid scan configuration')
    write_json(out/(args.phase+'-jobs.json'),jobs)
    report={'status':'running','fixtures':fixtures,'source_sha256':sources,'planned':len(jobs),'trials':[],
            'protocol':{'timestep_fs':.5,'damping_fs':350.,'thermostat':'unchanged per-atom LAMMPS Langevin',
                        'track_nm':[20,20],'batch_steps':100,'sample_ps':.05,'full_coordinates_ps':.5,
                        'joint_force_gate':1e-4,'surface_guard':'all atoms >8 A for >=5 ps invalid',
                        'finish':'sampled absolute COM x >=200 A; never extrapolate unfinished times',
                        'scope':'seven cars excluding acetone and eicosane','phase':args.phase}}
    began=time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_trial,c,str(out),sources) for c in jobs]
        for f in as_completed(futures):
            row=f.result();report['trials'].append(row)
            write_json(out/(args.phase+'-report.json'),report)
            print(json.dumps({'completed':len(report['trials']),'planned':len(jobs),'key':row['key'],
                              'status':row['status'],'x_nm':row.get('final_x_nm'),'finish_ps':row.get('finish_time_ps'),
                              'lift':row.get('sustained_lift'),'error':row.get('error')},ensure_ascii=False),flush=True)
    report.update(status='complete',wall_seconds=time.perf_counter()-began)
    write_json(out/(args.phase+'-report.json'),report)


if __name__=='__main__':main()
