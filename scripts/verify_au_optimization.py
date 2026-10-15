"""Native OPT/safe-cell checks and short matched explicit-surface MD timings.

Only writes outputs/au-optimization-check.  Each timing case is a fresh
serial process.  No game server and no old 18-trial source are modified.
"""
from __future__ import annotations
import os
for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'

import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import sys
import time
import traceback

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
OUTPUT=ROOT/'outputs'/'au-optimization-check'
BIG=(57.69,59.96)
SMALL=(40.,39.9)


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def fingerprint():
    names=('molecular_car/explicit_au.py','scripts/verify_au_optimization.py',
           'molecular_car/surface_racing.py','molecular_car/au_surface.py',
           'examples/benzene.mol','outputs/au-surface-check/au111-relaxed.xyz')
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}


def native_eam_snapshot(style,coordinates,elements,cell):
    import numpy as np
    from molecular_car.explicit_au import _native_system,build_relaxed_supercell
    from scripts.lammps_runtime import open_lammps
    with open_lammps() as (lmp,runtime):
        installation=Path(runtime['lammps_library']).parents[1]
        pots={name:installation/'Potentials'/filename for name,filename in [('EAM','Au_u3.eam'),('SAIP','CHAu.ILP')]}
        slab=build_relaxed_supercell(cell)
        _native_system(lmp,slab,coordinates,elements,pots,style,mobile=True)
        lmp.command('run 0')
        n=len(slab.positions)+len(coordinates)
        forces=np.ctypeslib.as_array(lmp.gather_atoms('f',1,3),shape=(n*3,)).copy().reshape(n,3)
        return float(lmp.get_thermo('pe')),forces,runtime


def long_polyphenyl(rings=21):
    """Planar para-linked oligophenyl fixture: complete connected C/H graph."""
    import numpy as np
    r,ch,link=1.397,1.09,1.47
    step=2*r+link
    positions=[];elements=[]
    for ring in range(rings):
        for k in range(6):
            angle=k*np.pi/3
            direction=np.array([np.cos(angle),np.sin(angle),0.])
            center=np.array([ring*step,0.,3.34256139116])
            positions.append(center+r*direction);elements.append('C')
            bonded=(k==0 and ring<rings-1) or (k==3 and ring>0)
            if not bonded:
                positions.append(center+(r+ch)*direction);elements.append('H')
    return np.array(positions),elements


def atom_force_differences(positions,elements,cell,style):
    import numpy as np
    from molecular_car.explicit_au import FrozenAuOracle
    with FrozenAuOracle(positions,elements,cell=cell,style=style,allow_unapproved_opt=True) as oracle:
        result=oracle.evaluate(positions)
        numeric=np.zeros_like(positions)
        step=1e-5
        for atom in range(len(positions)):
            for axis in range(3):
                plus=positions.copy();minus=positions.copy()
                plus[atom,axis]+=step;minus[atom,axis]-=step
                numeric[atom,axis]=-(oracle.evaluate(plus)['energy_eV']-oracle.evaluate(minus)['energy_eV'])/(2*step)
        error=result['forces_eV_A']-numeric
        return {'style':style,'finite_difference_step_A':step,
            'energy_eV':result['energy_eV'],'analytic_forces_eV_A':result['forces_eV_A'].tolist(),
            'numeric_forces_eV_A':numeric.tolist(),'errors_eV_A':error.tolist(),
            'maximum_absolute_error_eV_A':float(np.max(np.abs(error))),
            'passed':bool(np.max(np.abs(error))<1e-6)}


def static_checks():
    import numpy as np
    from molecular_car.explicit_au import (FrozenAuOracle,choose_safe_cell,
        reference_geometry_metadata,build_relaxed_supercell)
    from molecular_car.surface_racing import benzene_input
    prepared,geometry,_=benzene_input()
    elements=[a['element'] for a in prepared.atoms]
    d=4.08/np.sqrt(2);h=np.sqrt(3)*d
    base=geometry+[5*d+.37,h/3+.23,3.34256139116+.15]
    big=choose_safe_cell(base,BIG);small=choose_safe_cell(base,SMALL)
    metadata=reference_geometry_metadata()
    report={'status':'running','source_sha256':fingerprint(),'reference':metadata,
            'cells':{'big':big.metadata(),'small':small.metadata()},'checks':[]}
    results={}
    for label,style,cell in [('base_big','base',big),('opt_big','opt',big),('base_small','base',small),('opt_small','opt',small)]:
        with FrozenAuOracle(base,elements,cell=cell,style=style,allow_unapproved_opt=True) as oracle:
            results[label]=oracle.evaluate(base,include_gold_forces=True)
    energy_error=abs(results['base_big']['energy_eV']-results['opt_big']['energy_eV'])
    force_error=float(np.max(np.abs(results['base_big']['all_forces_eV_A']-results['opt_big']['all_forces_eV_A'])))
    opt_match=energy_error<1e-8 and force_error<1e-8
    report['checks'].append({'name':'base_vs_opt_all_atoms','energy_error_eV':energy_error,'force_error_eV_A':force_error,'passed':opt_match,'required':False})
    report['force_finite_differences']=[]
    for name,positions in [('planar',base.copy()),('carbon_1_bent_z_0.12_A',base.copy())]:
        if name=='planar':positions[:,2]=positions[:,2].mean()
        else:positions[0,2]+=.12
        for style in ('base','opt'):
            check=atom_force_differences(positions,elements,big,style)
            check['geometry']=name;report['force_finite_differences'].append(check)
            if style=='base' and not check['passed']:raise AssertionError('Base native force failed full Cartesian finite differences')
    report['OPT_approved']=bool(opt_match and all(c['passed'] for c in report['force_finite_differences']))
    report['OPT_policy']='Disabled by default; any OPT timings are experimental, not an approved physical speedup.' if not report['OPT_approved'] else 'Validated'
    save(OUTPUT/'opt-validation.json',{'OPT_approved':report['OPT_approved'],'policy':report['OPT_policy'],'comparisons':report['checks'],'full_Cartesian_checks':report['force_finite_differences']})
    e_small=abs(results['base_big']['energy_eV']-results['base_small']['energy_eV'])
    f_small=float(np.max(np.abs(results['base_big']['forces_eV_A']-results['base_small']['forces_eV_A'])))
    if e_small>1e-7 or f_small>1e-7:raise AssertionError('6 nm and safe 4 nm same-registry cross interaction mismatch')
    report['checks'].append({'name':'same_registry_big_vs_small','energy_error_eV':e_small,'car_force_error_eV_A':f_small,'passed':True})
    e0,f0,runtime=native_eam_snapshot('base',base,elements,big)
    e1,f1,_=native_eam_snapshot('eam_opt',base,elements,big)
    e_error,f_error=abs(e1-e0),float(np.max(np.abs(f1-f0)))
    if e_error>1e-7 or f_error>1e-8:raise AssertionError('Base/OPT EAM+SAIP native all-atom mismatch')
    report['runtime']=runtime
    report['checks'].append({'name':'base_vs_EAMopt_with_base_SAIP_all_atoms','energy_error_eV':e_error,'force_error_eV_A':f_error,'passed':True})
    lengths=np.array(small.lengths_xy)
    poses={'seam_x':geometry+[lengths[0]-.2,h/3,3.34],
           'seam_y':geometry+[5*d,lengths[1]-.2,3.34],
           'corner':geometry+[lengths[0]-.2,lengths[1]-.2,3.34]}
    with FrozenAuOracle(base,elements,cell=small,style='base',auto_expand=False) as oracle:
        for name,pose in poses.items():
            reference=oracle.evaluate(pose)
            for ix,iy in [(1,0),(0,-1),(7,-5),(-11,9)]:
                translated=pose+[ix*lengths[0],iy*lengths[1],0]
                result=oracle.evaluate(translated)
                e=abs(reference['energy_eV']-result['energy_eV'])
                f=float(np.max(np.abs(reference['forces_eV_A']-result['forces_eV_A'])))
                identity=float(np.max(np.abs(result['unwrapped_positions']-translated)))
                if e>1e-8 or f>1e-8 or identity>1e-8:raise AssertionError(f'Periodic seam/translation failed {name}/{ix}/{iy}')
                report['checks'].append({'name':name,'integer_cells':[ix,iy],'energy_error_eV':e,'force_error_eV_A':f,'unwrapped_identity_error_A':identity,'car_atoms':len(pose),'passed':True})
    long_positions,long_elements=long_polyphenyl()
    rejected=False
    try:choose_safe_cell(long_positions,SMALL,auto_expand=False)
    except ValueError:rejected=True
    if not rejected:raise AssertionError('Unsafe complete long molecule was not rejected in strict mode')
    expanded=choose_safe_cell(long_positions,SMALL)
    with FrozenAuOracle(long_positions,long_elements,cell=expanded,style='base') as oracle:
        long_result=oracle.evaluate(long_positions)
        moved=long_positions+[3*expanded.lengths_xy[0],-2*expanded.lengths_xy[1],0]
        long_shifted=oracle.evaluate(moved)
        long_e=abs(long_result['energy_eV']-long_shifted['energy_eV'])
        long_f=float(np.max(np.abs(long_result['forces_eV_A']-long_shifted['forces_eV_A'])))
        if long_e>1e-7 or long_f>1e-7:raise AssertionError('Expanded long molecule periodic invariance failed')
    report['checks'].append({'name':'complete_long_molecule_auto_expand','atoms':len(long_positions),
        'span_xy_A':np.ptp(long_positions[:,:2],axis=0).tolist(),'cell':expanded.metadata(),
        'strict_mode_rejected':rejected,'energy_error_eV':long_e,'force_error_eV_A':long_f,'passed':True})
    # Geometry reconstruction is explicit and its approximation is measured.
    original=np.loadtxt(ROOT/'outputs/au-surface-check/au111-relaxed.xyz',skiprows=2,usecols=(1,2,3))
    canonical=build_relaxed_supercell(big).positions
    report['reference']['maximum_reconstruction_z_error_A']=float(np.max(np.abs(original[:,2]-canonical[:,2])))
    report['status']='passed'
    save(OUTPUT/'static-checks.json',report)
    print(f"PASS static: {len(report['checks'])} checks; frozen OPT f-error {force_error:.3g} eV/A; shrink f-error {f_small:.3g} eV/A",flush=True)
    return report


def performance_case(config):
    from molecular_car.explicit_au import ExplicitAuRace
    from scripts.verify_surface_racing import columns_and_row,actual_blas_threads
    simulation=None
    result={'status':'running','configuration':config,'source_sha256':fingerprint()}
    name=config['name'];path=OUTPUT/'timings'/name
    path.mkdir(parents=True,exist_ok=True)
    try:
        simulation=ExplicitAuRace(model=config['model'],style=config['style'],requested_xy=config['requested_xy'],allow_unapproved_opt=True)
        result.update(simulation.description())
        warmup=simulation.advance(1000)
        initial=simulation.sample();samples=[initial];timings=[]
        for batch in range(3):
            timing=simulation.advance(2000)
            timing['batch']=batch+1;timing['steps']=2000
            timing['seconds_per_step']=timing['wall_seconds']/2000
            timings.append(timing);samples.append(simulation.sample())
        trajectory=path/'trajectory.csv'
        with trajectory.open('w',newline='',encoding='utf-8') as handle:
            writer=csv.writer(handle);writer.writerow(columns_and_row(initial,initial)[0])
            for sample in samples:writer.writerow(columns_and_row(sample,initial)[1])
        result.update(status='passed',warmup_steps=1000,warmup_timing=warmup,
            measurement_steps=6000,measurement_duration_ps=3.0,batches=timings,
            median_wall_seconds_per_step=statistics.median(t['seconds_per_step'] for t in timings),
            initial=initial,final=samples[-1],production_monitor=simulation.monitor,
            trajectory_csv=str(trajectory),runtime_BLAS=actual_blas_threads(),
            timing_interpretation='Three contiguous 2000-step MD batches after 1000-step warmup; callback/monitor included; setup/sampling/CSV excluded')
        print(f"PASS timing {name}: {result['median_wall_seconds_per_step']*1e6:.2f} us/step; setup {simulation.preparation_wall_seconds:.3f}s",flush=True)
    except BaseException as error:
        result.update(status='failed',error=str(error),traceback=traceback.format_exc())
        traceback.print_exc()
    finally:
        if simulation is not None:simulation.close()
        save(path/'report.json',result)
    return 0 if result['status']=='passed' else 1


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--static-only',action='store_true')
    parser.add_argument('--worker',choices=('static','timing'))
    parser.add_argument('--configuration',type=Path)
    args=parser.parse_args()
    if args.worker=='static':static_checks();return 0
    if args.worker=='timing':return performance_case(json.loads(args.configuration.read_text(encoding='utf-8')))
    OUTPUT.mkdir(parents=True,exist_ok=True)
    process=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker','static'],cwd=ROOT)
    if process.returncode:return process.returncode
    report={'status':'running','source_sha256':fingerprint(),
        'static_checks':json.loads((OUTPUT/'static-checks.json').read_text(encoding='utf-8')),
        'timings':[],'protocol':{'temperature_K':600,'timestep_ps':.0005,'damping_ps':.35,
        'warmup_steps':1000,'batches':3,'steps_per_batch':2000,
        'cutoff_A':16,'skin_A':2,'layers':6,'bottom_frozen_layers':2,
        'limits':'Short matched-cost measurements, not new 18-trial race statistics; mobile small cell changes phonon modes'}}
    save(OUTPUT/'report.json',report)
    if not args.static_only:
        for model in ('frozen_au','mobile_au'):
            configurations=[('base','big',BIG),('base','small',SMALL)]
            if model=='frozen_au':configurations.insert(1,('opt','big',BIG))
            if model=='mobile_au':configurations += [('eam_opt','big',BIG)]
            for style,size,xy in configurations:
                config={'name':f'{model}-{style}-{size}','model':model,'style':style,'requested_xy':xy}
                path=OUTPUT/'timings'/config['name']/'configuration.json';save(path,config)
                print(f"START timing {config['name']}",flush=True)
                process=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker','timing','--configuration',str(path)],cwd=ROOT)
                result=json.loads((path.parent/'report.json').read_text(encoding='utf-8'))
                report['timings'].append(result)
                if process.returncode:report['status']='failed';save(OUTPUT/'report.json',report);return 1
                save(OUTPUT/'report.json',report)
    report['status']='passed';save(OUTPUT/'report.json',report)
    print(f"Report: {OUTPUT/'report.json'}",flush=True)
    return 0


if __name__=='__main__':raise SystemExit(main())
