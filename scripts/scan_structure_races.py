"""Actual game races to separate molecular structure effects from seed noise.

Independent single-process LAMMPS instances; unchanged force model, start pose,
strict joint preparation, timestep, bath and 100-step game run cadence.
Unfinished races are right censored, never assigned a predicted finish time.
"""
from __future__ import annotations
import os
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[variable] = '1'
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from molecular_car.engine import Simulation, TIMESTEP_FS, FRAME_STEPS
from molecular_car.track import TRACK
from molecular_car.periodic_au import parameter_file

SOURCES = ('molecular_car/engine.py', 'molecular_car/molecule.py', 'molecular_car/start_pose.py',
           'molecular_car/relaxation.py', 'molecular_car/preparation_bfgs.py', 'molecular_car/track.py',
           'molecular_car/fields.py', 'molecular_car/fields_kernel.py', 'molecular_car/external_fields.py',
           'molecular_car/periodic_au.py', 'molecular_car/periodic_au_kernel.py',
           'molecular_car/element_models.py', 'molecular_car/data/au111_reference.json',
           'scripts/lammps_runtime.py', 'scripts/scan_structure_races.py')
FIXTURES = ('benzene-12', 'coronene-36', 'ethynyl-rigid', 'folded-eicosane', 'test-car-52', 'fullerene-c60')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, data):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf8')
    temp.replace(path)


def prepare_inputs(out):
    directory = out/'inputs'; directory.mkdir(parents=True, exist_ok=True)
    old = ROOT/'outputs/md-release-validation/inputs'
    old_manifest = json.loads((old/'manifest.json').read_text(encoding='utf8'))
    fixtures = []
    for name in FIXTURES:
        path = directory/(name+'.mol')
        if name == 'coronene-36':
            smiles = 'c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67'
            molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
            if rdMolDescriptors.CalcMolFormula(molecule) != 'C24H12' or sorted(len(r) for r in molecule.GetRingInfo().AtomRings()) != [6]*7:
                raise AssertionError('Coronene fixture graph failed formula/ring check')
            settings = AllChem.ETKDGv3(); settings.randomSeed=20261013;settings.numThreads=1
            if AllChem.EmbedMolecule(molecule, settings) != 0:
                raise ValueError('Could not generate coronene coordinates')
            if AllChem.UFFOptimizeMolecule(molecule, maxIters=3000) != 0:
                raise ValueError('Coronene fixture vacuum optimization did not converge')
            molecule.SetProp('_Name','Coronene C24H12 planar disk')
            path.write_text(Chem.MolToMolBlock(molecule),encoding='utf8')
            detail={'label':'Coronene disk','smiles':smiles,'source':'locally generated RDKit ETKDGv3 + UFF','embed_seed':20261013}
        else:
            reference = next(f for f in old_manifest['fixtures'] if f['id']==name)
            source = old/reference['file']
            if digest(source) != reference['file_sha256']:
                raise ValueError('Existing reference fixture changed: '+name)
            path.write_bytes(source.read_bytes())
            detail={'label':reference['label'],'source':str(source.relative_to(ROOT)),'reference_sha256':reference['file_sha256']}
        molecule=Chem.MolFromMolFile(str(path),removeHs=False)
        fixtures.append({'id':name,'file':path.name,'sha256':digest(path),'formula':rdMolDescriptors.CalcMolFormula(molecule),
                         'atoms':molecule.GetNumAtoms(), **detail})
    manifest={'rdkit_version':rdBase.rdkitVersion,'fixtures':fixtures}
    write_json(directory/'manifest.json',manifest)
    return manifest


def trial_key(c):
    extra = f"-L{c['lattice_scale']:g}-G{c['growth']:g}" if 'lattice_scale' in c else ''
    return f"{c['phase']}-{c['molecule']}-T{c['temperature_K']:g}-S{c['strength']:g}{extra}-seed{c['seed']}-{c['limit_ps']:g}ps"


def run_trial(c, folder_name, source_hashes):
    out=Path(folder_name); folder=out/'trials';folder.mkdir(exist_ok=True)
    key=trial_key(c);target=folder/(key+'.json');coords=folder/(key+'.npz')
    if target.exists() and coords.exists():
        old=json.loads(target.read_text(encoding='utf8'))
        if old.get('status')=='passed' and old['config']==c and old['source_sha256']==source_hashes:
            return old
    result={'key':key,'config':c,'source_sha256':source_hashes,'status':'running'}
    simulation=None;begin=time.perf_counter()
    try:
        path=out/'inputs'/(c['molecule']+'.mol')
        simulation=Simulation(path.read_text(encoding='utf8'),path.name,
            {'surface':'periodic-au','temperature':c['temperature_K'],'damping_fs':350.,'lattice_scale':c.get('lattice_scale',1.)},
            random_seed=c['seed'],drive={'mode':'gradient','strength':c['strength'],'growth':c.get('growth',1.)},start_pose={'mode':'auto'})
        preparation=simulation.runtime['preparation']
        assert preparation['accepted'] and preparation['global_force_norm_kcal_mol_A']<=1e-4
        assert simulation.lmp.get_thermo('step')==0 and simulation.lmp.extract_global('dt')==.5
        initial=simulation.coordinates();initial_center=simulation.center()
        bonds=np.array([(b['a']-1,b['b']-1) for b in simulation.prepared.bonds])
        reference_lengths=np.linalg.norm(initial[bonds[:,0]]-initial[bonds[:,1]],axis=1)
        result.update(preparation_wall_seconds=time.perf_counter()-begin, runtime=simulation.runtime,
                      atoms=simulation.n,mass_u=float(simulation.masses.sum()),alpha_sum=float(simulation.alphas.sum()),
                      initial_center_nm=(initial_center/10).tolist())
        samples=[];positions=[];position_steps=[]
        ratio_min,ratio_max,z_min,z_max=1.,1.,float('inf'),-float('inf')
        lift_ps=0.;max_lift_ps=0.;lower=None;outcome='time-limit'
        limit_steps=round(c['limit_ps']*1000/TIMESTEP_FS)
        propagate=time.perf_counter()
        for step in range(0,limit_steps+1,FRAME_STEPS):
            frame=simulation.advance(FRAME_STEPS) if step else simulation.frame()
            x=np.array(frame['positions']).reshape(-1,3)
            if not np.isfinite(frame['energy']) or not np.isfinite(x).all():
                raise AssertionError('Nonfinite trajectory')
            lengths=np.linalg.norm(x[bonds[:,0]]-x[bonds[:,1]],axis=1)/reference_lengths
            ratio_min=min(ratio_min,float(lengths.min()));ratio_max=max(ratio_max,float(lengths.max()))
            if ratio_min<.65 or ratio_max>1.45:
                raise AssertionError('Bond length diagnostic out of range')
            z_min=min(z_min,float(x[:,2].min()));z_max=max(z_max,float(x[:,2].max()))
            if step:
                lift_ps=lift_ps+FRAME_STEPS*TIMESTEP_FS/1000 if x[:,2].min()>8. else 0.
                max_lift_ps=max(max_lift_ps,lift_ps)
            entry={'step':step,'time_ps':step*TIMESTEP_FS/1000,'center_nm':(np.array(frame['center'])/10).tolist(),
                   'energy_kcal_mol':frame['energy'],'internal_temperature_K':frame['internal_temperature'],
                   'min_z_A':float(x[:,2].min()),'max_z_A':float(x[:,2].max())}
            samples.append(entry)
            finished=frame['center'][0]>=TRACK['finish_x']
            if step%(10*FRAME_STEPS)==0 or finished or step==limit_steps:
                positions.append(x);position_steps.append(step)
            if finished:
                outcome='finish';lower=(step-FRAME_STEPS)*TIMESTEP_FS/1000;break
        t=np.array([s['time_ps'] for s in samples]);centers=np.array([s['center_nm'] for s in samples]);relative=centers-centers[0]
        np.savez_compressed(coords,positions_A=np.stack(positions),steps=np.array(position_steps),masses_u=simulation.masses,bonds=bonds)
        result.update(status='passed',outcome=outcome,finish_time_ps=float(t[-1]) if outcome=='finish' else None,
                      finish_crossing_interval_ps=[lower,float(t[-1])] if outcome=='finish' else None,
                      right_censored=outcome!='finish',observed_time_ps=float(t[-1]),
                      final_x_nm=float(centers[-1,0]),final_forward_nm=float(relative[-1,0]),
                      maximum_lateral_excursion_nm=float(np.max(np.abs(relative[:,1]))),
                      bond_ratio_range=[ratio_min,ratio_max],atom_height_range_A=[z_min,z_max],
                      sustained_lift=max_lift_ps>=5.,maximum_lift_duration_ps=max_lift_ps,
                      competitive_valid=max_lift_ps<5.,samples=samples,coordinates_file=coords.name,
                      propagation_wall_seconds=time.perf_counter()-propagate)
    except BaseException as error:
        result.update(status='failed',error=str(error),traceback=traceback.format_exc())
    finally:
        if simulation is not None:simulation.close()
    result['total_wall_seconds']=time.perf_counter()-begin
    write_json(target,result)
    return result


def summarize(rows):
    groups={}
    for r in rows:
        c=r['config'];key=(c['temperature_K'],c['strength'])
        groups.setdefault(key,[]).append(r)
    summaries=[]
    for (temp,strength),conditions in sorted(groups.items()):
        structures=[]
        for name in FIXTURES:
            trials=[r for r in conditions if r['config']['molecule']==name]
            valid=[r for r in trials if r['status']=='passed' and r['competitive_valid']]
            finish=[r['finish_time_ps'] for r in valid if r['outcome']=='finish']
            # Capped time is a restricted observation, NOT a finish-time estimate.
            capped=[r['observed_time_ps'] for r in valid]
            structures.append({'molecule':name,'n':len(trials),'valid':len(valid),'finished':len(finish),
              'finish_mean_ps':float(np.mean(finish)) if finish else None,
              'finish_sd_ps':float(np.std(finish,ddof=1)) if len(finish)>1 else None,
              'finish_range_ps':[min(finish),max(finish)] if finish else None,
              'capped_mean_time_ps':float(np.mean(capped)) if capped else None,
              'mean_final_forward_nm':float(np.mean([r['final_forward_nm'] for r in valid])) if valid else None,
              'failed':len(trials)-sum(r['status']=='passed' for r in trials),
              'lifted':sum(r.get('sustained_lift',False) for r in trials),
              'seeds':[r['config']['seed'] for r in trials]})
        complete=all(g['n'] and g['finished']==g['n'] for g in structures)
        metrics={}
        if complete:
            means=np.array([g['finish_mean_ps'] for g in structures]);counts=np.array([g['n'] for g in structures])
            grand=float(np.average(means,weights=counts))
            between=float(np.sum(counts*(means-grand)**2))
            within=float(sum((g['n']-1)*g['finish_sd_ps']**2 for g in structures if g['n']>1))
            denominator=int(counts.sum()-len(counts))
            noise=float(np.sqrt(within/denominator)) if denominator else None
            metrics={'slow_fast_mean_ratio':float(means.max()/means.min()),'mean_gap_ps':float(np.ptp(means)),
                     'between_structure_mean_sd_ps':float(np.std(means,ddof=1)),
                     'pooled_within_structure_sd_ps':noise,'anova_eta_squared':between/(between+within) if between+within else None,
                     'gap_to_pooled_noise_ratio':float(np.ptp(means)/noise) if noise else None,
                     'grand_mean_ps':grand}
        summaries.append({'temperature_K':temp,'strength':strength,'planned':len(conditions),
                          'all_finished_valid':complete,'structures':structures,**metrics})
    return summaries


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/structure-race-scan')
    parser.add_argument('--jobs',type=Path)
    parser.add_argument('--phase',default='screen')
    parser.add_argument('--limit-ps',type=float,default=40.)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    if not 1<=args.workers<=6 or not 1<=args.limit_ps<=2000:
        parser.error('Use 1..6 workers and 1..2000 ps')
    if not (out/'inputs/manifest.json').exists():manifest=prepare_inputs(out)
    else:manifest=json.loads((out/'inputs/manifest.json').read_text(encoding='utf8'))
    if args.prepare_only:print(json.dumps(manifest,ensure_ascii=False,indent=2));return
    if args.jobs:
        tasks=json.loads(args.jobs.read_text(encoding='utf8'))
    else:
        tasks=[{'phase':args.phase,'molecule':name,'temperature_K':temp,'strength':strength,
                'seed':98245165,'limit_ps':args.limit_ps}
               for temp in (150.,300.,600.,900.) for strength in (.45,.60,.75) for name in FIXTURES]
    for c in tasks:
        if c['molecule'] not in FIXTURES or not 0<c['temperature_K']<=1000 or not 0<=c['strength']<=1 or not 1<=c['limit_ps']<=2000:
            raise ValueError('Invalid scan configuration')
    sources={name:digest(ROOT/name) for name in SOURCES}
    sources.update({str((out/'inputs'/f['file']).relative_to(ROOT)):digest(out/'inputs'/f['file']) for f in manifest['fixtures']})
    sources['CHAu.ILP']=digest(parameter_file())
    report={'status':'running','created_utc':datetime.now(timezone.utc).isoformat(),'source_sha256':sources,
            'planned':len(tasks),'fixtures':manifest['fixtures'],'phase':args.phase,'trials':[],
            'protocol':{'surface':'periodic-au','track_A':TRACK,'timestep_fs':TIMESTEP_FS,'damping_fs':350.,
                        'start_pose':{'mode':'auto'},'drive':'fixed gradient','warmup_ps':0.,'run_batch_steps':FRAME_STEPS,
                        'COM_sample_ps':FRAME_STEPS*TIMESTEP_FS/1000,'full_coordinates_sample_ps':.5,
                        'finish':'first sampled absolute COM x >=20 nm','censoring':'unfinished at explicit observation limit',
                        'lift_guard':'all atoms above 8 Å for >=5 ps; not a valid surface race'}}
    began=time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_trial,c,str(out),sources) for c in tasks]
        for future in as_completed(futures):
            row=future.result();report['trials'].append(row)
            report['summary']=summarize(report['trials'])
            write_json(out/(args.phase+'-report.json'),report)
            c=row['config'];print(f"{len(report['trials'])}/{len(tasks)} {c['molecule']} T={c['temperature_K']:g} S={c['strength']:g} seed={c['seed']} "
                + (f"{row['outcome']} {row['observed_time_ps']:.2f}ps x={row['final_x_nm']:.2f}nm" if row['status']=='passed' else row['error']),flush=True)
    unchanged=all(digest(ROOT/name)==expected for name,expected in sources.items() if name!='CHAu.ILP') and digest(parameter_file())==sources['CHAu.ILP']
    report.update(status='passed' if unchanged and all(r['status']=='passed' for r in report['trials']) else 'completed-with-errors',
                  source_unchanged=unchanged,elapsed_wall_seconds=time.perf_counter()-began)
    write_json(out/(args.phase+'-report.json'),report)
    fields=('key','molecule','temperature_K','strength','seed','status','outcome','finish_time_ps','observed_time_ps','final_x_nm','competitive_valid')
    with (out/(args.phase+'-trials.csv')).open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        for row in sorted(report['trials'],key=lambda r:r['key']):
            writer.writerow({k:row.get(k,row['config'].get(k,'')) for k in fields})
    print(json.dumps({'status':report['status'],'trials':len(tasks),'source_unchanged':unchanged,'wall_seconds':report['elapsed_wall_seconds']}))


if __name__=='__main__':main()
