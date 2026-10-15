"""Add neutral polar N/O controls using the existing game parameter extension."""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem,rdMolDescriptors
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'scripts'))
from scan_structure_races import run_trial,write_json,digest,SOURCES
from molecular_car.molecule import prepare_molecule
from molecular_car.periodic_au import parameter_file
from molecular_car.fields import GRADIENT_BASE_FIELD,GRADIENT_FIELD_SLOPE,GRADIENT_POLARIZATION_GAIN

OUT=ROOT/'outputs/structure-race-polar'
FIXTURES=('pyridine-11','benzonitrile-13','acetone-10')


def prepare():
    (OUT/'inputs').mkdir(parents=True,exist_ok=True);(OUT/'trials').mkdir(exist_ok=True)
    records=[]
    for name,label,smiles in (('pyridine-11','Pyridine','c1ccncc1'),('benzonitrile-13','Benzonitrile','N#Cc1ccccc1'),('acetone-10','Acetone','CC(=O)C')):
        molecule=Chem.AddHs(Chem.MolFromSmiles(smiles));settings=AllChem.ETKDGv3();settings.randomSeed=20261013;settings.numThreads=1
        if AllChem.EmbedMolecule(molecule,settings)!=0 or AllChem.UFFOptimizeMolecule(molecule,maxIters=3000)!=0:
            raise ValueError('Could not generate polar control '+name)
        molecule.SetProp('_Name',label+' neutral polar control')
        path=OUT/'inputs'/(name+'.mol');path.write_text(Chem.MolToMolBlock(molecule),encoding='utf8')
        p=prepare_molecule(path.read_text(encoding='utf8'),path.name,start_pose={'mode':'auto'})
        q=np.array([a['charge'] for a in p.atoms]);mass=np.array([a['mass'] for a in p.atoms]);center=np.average(p.positions,axis=0,weights=mass)
        dipole=np.sum(q[:,None]*(p.positions-center),axis=0)
        records.append({'id':name,'label':label,'smiles':smiles,'formula':rdMolDescriptors.CalcMolFormula(molecule),
                        'atoms':len(p.atoms),'file':path.name,'sha256':digest(path),'embed_seed':20261013,
                        'net_charge_e':float(q.sum()),'vacuum_model_dipole_D':float(np.linalg.norm(dipole)*4.8032047),
                        'dipole_qualification':'Fixed Gasteiger-charge estimate, not a measured or quantum-chemical dipole',
                        'Au_qualification':'C/H bonded SAIP adaptation; N/O use existing uncalibrated mixed van der Waals parameters'})
    manifest={'rdkit_version':rdBase.rdkitVersion,'fixtures':records};write_json(OUT/'inputs/manifest.json',manifest)
    return manifest


def annotate(row):
    if row['status']!='passed':return row
    c=row['config'];path=OUT/'inputs'/(c['molecule']+'.mol')
    p=prepare_molecule(path.read_text(encoding='utf8'),path.name,start_pose={'mode':'auto'})
    q=np.array([a['charge'] for a in p.atoms]);mass=np.array([a['mass'] for a in p.atoms]);alpha=np.array([a['polarizability'] for a in p.atoms])
    x=np.load(OUT/'trials'/row['coordinates_file'])['positions_A'][0];assert len(q)==len(x)
    dipole=np.sum(q[:,None]*(x-np.average(x,axis=0,weights=mass)),axis=0);norm=np.linalg.norm(dipole)
    g=GRADIENT_FIELD_SLOPE*c['strength'];base=GRADIENT_BASE_FIELD*c['strength']
    charged=float(np.sum(q*g*x[:,2]));induced=float(np.sum(GRADIENT_POLARIZATION_GAIN*alpha*g*(base+g*x[:,0])))
    row['initial_polar_response']={'net_charge_e':float(q.sum()),'model_dipole_D':float(norm*4.8032047),
        'model_dipole_vector_e_A':dipole.tolist(),'dipole_angle_to_z_deg':float(np.degrees(np.arccos(np.clip(dipole[2]/norm,-1,1)))) if norm>1e-12 else None,
        'charge_forward_force_kcal_mol_A':charged,'induced_forward_force_kcal_mol_A':induced,
        'direct_charge_to_induced_forward_force_ratio':charged/induced if induced else None,
        'qualification':'Initial forward-force decomposition only; partial charges can also affect torque/orientation and surface contact indirectly.'}
    write_json(OUT/'trials'/(row['key']+'.json'),row)
    return row


def main():
    parser=argparse.ArgumentParser(__doc__);parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    if not 1<=args.workers<=4:raise ValueError('Use 1..4 workers')
    OUT.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((OUT/'inputs/manifest.json').read_text(encoding='utf8')) if (OUT/'inputs/manifest.json').exists() else prepare()
    extra=('molecular_car/hetero_au.py','molecular_car/hetero_au_kernel.py','molecular_car/data/hetero_au_v1.json','scripts/scan_polar_races.py')
    sources={name:digest(ROOT/name) for name in (*SOURCES,*extra)}
    sources.update({str((OUT/'inputs'/r['file']).relative_to(ROOT)):digest(OUT/'inputs'/r['file']) for r in manifest['fixtures']})
    sources['CHAu.ILP']=digest(parameter_file())
    screen=[{'phase':'screen','molecule':name,'temperature_K':temp,'strength':strength,'seed':98245165,'limit_ps':40.}
            for temp in (150.,300.,600.,900.) for strength in (.45,.6,.75) for name in FIXTURES]
    confirm=[{'phase':'confirm','molecule':name,'temperature_K':temp,'strength':strength,'seed':seed,'limit_ps':400.}
             for temp,strength in ((150.,.6),(600.,.75),(600.,.6),(600.,.45))
             for seed in ((2718281,1234567,31415927,1618033,42424243) if temp==150 else (2718281,1234567,31415927)) for name in FIXTURES]
    report={'status':'running','source_sha256':sources,'fixtures':manifest['fixtures'],'planned':len(screen)+len(confirm),'trials':[],
            'protocol':'Same game protocol as scan_structure_races; polar controls retain standard masses, Gasteiger charges and existing N/O-Au approximations.'}
    begin=time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_trial,c,str(OUT),sources) for c in screen+confirm]
        for future in as_completed(futures):
            row=annotate(future.result());report['trials'].append(row);write_json(OUT/'report.json',report)
            c=row['config'];print(f"{len(report['trials'])}/{report['planned']} {c['phase']} {c['molecule']} T={c['temperature_K']:g} S={c['strength']:g} "
                 + (f"{row['outcome']} {row['observed_time_ps']:.2f}ps x={row['final_x_nm']:.2f}nm" if row['status']=='passed' else row['error']),flush=True)
    unchanged=all(digest(ROOT/name)==value for name,value in sources.items() if name!='CHAu.ILP') and digest(parameter_file())==sources['CHAu.ILP']
    report.update(status='passed' if unchanged and all(r['status']=='passed' for r in report['trials']) else 'completed-with-errors',
                  source_unchanged=unchanged,elapsed_wall_seconds=time.perf_counter()-begin)
    write_json(OUT/'report.json',report);print(json.dumps({'status':report['status'],'trials':report['planned'],'source_unchanged':unchanged}))


if __name__=='__main__':main()
