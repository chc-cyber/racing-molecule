"""Independent Cartesian gradients, adhesion averages and fixed-step NVE gates."""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
import json
from pathlib import Path
import sys
import time
import numpy as np
import numba
from rdkit import Chem, rdBase
from rdkit.ForceField import rdForceField
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'scripts'))
from molecular_car.molecule import prepare_molecule
from molecular_car.periodic_au import create_periodic_surface, BONDED_NORMAL_POLICY, parameter_file
from molecular_car.lattice_corrugation import LatticeCorrugation
from molecular_car.external_fields import ExternalFieldEvaluator
from molecular_car.engine import Simulation
from molecular_car.relaxation import JointRelaxationError
from scan_crystal_contrast import CARS
from scan_structure_races import write_json, digest, SOURCES


def validation_sources(out):
    """Resolve the files actually used, including external UFF/SAIP inputs."""
    names=(*SOURCES, 'molecular_car/lattice_corrugation.py',
           'molecular_car/hetero_au.py', 'molecular_car/hetero_au_kernel.py',
           'molecular_car/data/hetero_au_v1.json',
           'scripts/scan_crystal_contrast.py', 'scripts/validate_crystal_contrast.py')
    paths={name:ROOT/name for name in names}
    if (ROOT/'.local/lammps.json').is_file():
        paths['.local/lammps.json']=ROOT/'.local/lammps.json'
    paths.update({str((out/'inputs'/(name+'.mol')).relative_to(ROOT)).replace('\\','/'):
                  out/'inputs'/(name+'.mol') for name in CARS})
    paths['CHAu.ILP']=parameter_file().resolve()
    # UFF parameters and evaluation are compiled into RDKit; record the actual
    # modules and shared libraries rather than only our Python caller.
    from rdkit.Chem import rdForceFieldHelpers
    for module in (rdForceField, rdForceFieldHelpers):
        path=Path(module.__file__).resolve()
        paths['runtime/'+path.name]=path
    libraries=Path(Chem.__file__).resolve().parents[2]/'rdkit.libs'
    for pattern in ('*ForceField*.dll', '*PartialCharges*.dll', '*GraphMol*.dll'):
        for path in sorted(libraries.glob(pattern)):
            paths['runtime/'+path.name]=path
    # parameter_file and the native loader select the same installation.
    installation=paths['CHAu.ILP'].parent.parent
    for relative in ('bin/liblammps.dll', 'Python/lammps/__init__.py', 'Python/lammps/core.py'):
        paths['runtime/lammps/'+relative]=installation/relative
    return paths


def main():
    out=ROOT/'outputs/crystal-contrast-guarded';report={'status':'running','static':[],'adhesion':[],'nve':[],'wasm_fixtures':[],'cache_node_fd_convergence':[]}
    source_paths=validation_sources(out)
    source_hashes={name:digest(path) for name,path in source_paths.items()}
    report.update(source_sha256=source_hashes,
        runtime_versions={'python':sys.version.split()[0], 'rdkit':rdBase.rdkitVersion,
                          'numpy':np.__version__, 'numba':numba.__version__},
        source_integrity={'status':'running', 'files':len(source_paths),
                          'check':'same resolved source/input/runtime files hashed before and after validation'})
    write_json(out/'numerics.json',report)
    rng=np.random.default_rng(20261014)
    for name in CARS:
        path=out/'inputs'/(name+'.mol');p=prepare_molecule(path.read_text(encoding='utf8'),path.name,start_pose={'mode':'auto'})
        base=create_periodic_surface(p.atoms,p.positions,bonds=p.bonds,normal_policy=BONDED_NORMAL_POLICY)
        numpy_base=create_periodic_surface(p.atoms,p.positions,bonds=p.bonds,normal_policy=BONDED_NORMAL_POLICY,accelerator='numpy')
        q=np.array([a['charge'] for a in p.atoms]);alpha=np.array([a['polarizability'] for a in p.atoms])
        fields=ExternalFieldEvaluator(q,alpha,surface_enabled=False)
        r=p.positions+rng.normal(0,.025,p.positions.shape)+[base.period[0]*17+.137,-base.period[1]*9+.281,.25]
        original_r=r.copy()
        for scale,growth,guard in ((1.,1.,False),(3.,.25,False),(12.,.1,False),(3.,.25,True)):
            r=original_r.copy()
            if guard:r[:,2]+=2.413-float(r[:,2].min())
            potential=LatticeCorrugation(base,scale);other=LatticeCorrugation(numpy_base,scale)
            def evaluate(x):
                e,f=potential.energy_and_forces(x)
                flat=x.ravel().tolist();e+=p.forcefield.CalcEnergy(flat);f-=np.array(p.forcefield.CalcGrad(flat)).reshape(x.shape)
                fe,ff=fields.energy_and_forces(x,1,.6,6,0,7,0.,growth)
                return float(e+fe),f+ff
            e,f=evaluate(r);epsilon=1e-5;error=0.
            for atom,axis in np.ndindex(r.shape):
                d=np.zeros_like(r);d[atom,axis]=epsilon
                numeric=-(evaluate(r+d)[0]-evaluate(r-d)[0])/(2*epsilon)
                error=max(error,abs(numeric-f[atom,axis]))
            se,sf=potential.energy_and_forces(r);ne,nf=other.energy_and_forces(r)
            periodic_e,periodic_f=potential.energy_and_forces(r+[base.period[0]*11,-base.period[1]*13,0])
            row={'molecule':name,'scale':scale,'growth':growth,'local_core_guard_case':guard,'coordinates':int(r.size),'max_total_force_fd_error':error,
                 'numpy_numba_energy_error':abs(ne-se),'numpy_numba_force_error':float(np.max(abs(nf-sf))),
                 'periodic_energy_error':abs(periodic_e-se),'periodic_force_error':float(np.max(abs(periodic_f-sf)))}
            assert error<2e-5,row
            assert row['numpy_numba_force_error']<1e-8 and row['periodic_force_error']<1e-8,row
            report['static'].append(row)
            report['wasm_fixtures'].append({'filename':path.name,'content':path.read_text(encoding='utf8'),'local_core_guard_case':guard,
                'positions':r.ravel().tolist(),'scale':scale,'growth':growth,'strength':.6,'energy':e,'forces':f.ravel().tolist()})
        if name=='benzene-12':
            r=original_r.copy();r[:,2]+=2.4-float(r[:,2].min())
            potential=LatticeCorrugation(base,3.);energy,force=potential.energy_and_forces(r)
            errors=[]
            for h in (1e-5,1e-6,1e-7):
                error=0.
                for atom,axis in np.ndindex(r.shape):
                    d=np.zeros_like(r);d[atom,axis]=h
                    numerical=-(potential.energy_and_forces(r+d)[0]-potential.energy_and_forces(r-d)[0])/(2*h)
                    error=max(error,abs(numerical-force[atom,axis]))
                errors.append({'h_A':h,'max_force_error':float(error)})
            assert errors[-1]['max_force_error']<2e-5 and errors[-1]['max_force_error']<errors[0]['max_force_error']/20
            report['cache_node_fd_convergence'].append({'molecule':name,'minimum_height_A':2.4,'steps':errors,
                'interpretation':'C1 Catmull-Rom knot: cross-knot centered FD is first-order; analytic derivative converges'})
        if name in ('benzene-12','fullerene-c60','pyridine-11','test-car-52'):
            for dz in (0.,.3,.8):
                sample=[];filtered=[];original_zforce=[];filtered_zforce=[]
                mean=LatticeCorrugation(base,0.);scaled=LatticeCorrugation(base,6.)
                for i in range(12):
                    for j in range(18):
                        x=p.positions+[i*base.period[0]/12,j*base.period[1]/18,dz]
                        a,af=base.energy_and_forces(x);b,bf=scaled.energy_and_forces(x);c,_=mean.energy_and_forces(x)
                        sample.append(a);filtered.append((b,c));original_zforce.append(af[:,2].sum());filtered_zforce.append(bf[:,2].sum())
                original=np.array(sample);altered=np.array(filtered)
                row={'molecule':name,'height_shift_A':dz,'cell_grid':[12,18],
                     'mean_energy_change':float(altered[:,0].mean()-original.mean()),
                     'mean_vertical_force_change':float(np.mean(filtered_zforce)-np.mean(original_zforce)),
                     'original_peak_to_peak_kcal_mol':float(np.ptp(original)),
                     'scale6_peak_to_peak_kcal_mol':float(np.ptp(altered[:,0])),
                     'filtered_residual_fraction':float(np.ptp(altered[:,1])/np.ptp(original))}
                assert abs(row['mean_energy_change'])<1e-8 and abs(row['mean_vertical_force_change'])<1e-8,row
                report['adhesion'].append(row)
        write_json(out/'numerics.json',report)
        print(json.dumps({'static_complete':name,'max_error':max(x['max_total_force_fd_error'] for x in report['static'])}),flush=True)
    for name in ('benzene-12','test-car-52','fullerene-c60','pyridine-11'):
        path=out/'inputs'/(name+'.mol')
        for scale,growth in ((1.,1.),(2.,.25),(3.,1.)):
            try:
                sim=Simulation(path.read_text(encoding='utf8'),path.name,{'temperature':150.,'lattice_scale':scale},
                    thermostat=False,drive={'mode':'gradient','strength':.6,'growth':growth},start_pose={'mode':'auto'})
            except JointRelaxationError as error:
                report['nve'].append({'molecule':name,'scale':scale,'growth':growth,'status':'preparation-rejected',
                                      'diagnostics':error.diagnostics,'error':str(error)})
                write_json(out/'numerics.json',report)
                print(json.dumps({'molecule':name,'scale':scale,'status':'preparation-rejected'}),flush=True)
                continue
            try:
                phase={k:np.ctypeslib.as_array(sim.lmp.gather_atoms(a,1,3),shape=(sim.n*3,)).copy().tolist()
                       for a,k in (('x','positions'),('v','velocities'))}
                initial=sim.frame()['energy'];energies=[initial];frames=[sim.frame()]
                for _ in range(40):
                    frame=sim.advance(100);energies.append(frame['energy']);frames.append(frame)
                error=max(abs(np.array(energies)-initial))
                row={'molecule':name,'scale':scale,'growth':growth,'timestep_fs':sim.lmp.extract_global('dt'),
                     'steps':frame['step'],'duration_ps':2.,'max_energy_deviation_kcal_mol':float(error),
                     'final_energy_change_kcal_mol':float(energies[-1]-initial),'accepted':sim.runtime['preparation']['accepted'],
                     'residual_force':sim.runtime['preparation']['global_force_norm_kcal_mol_A'],'phase_space':phase,'frames':frames}
                assert row['timestep_fs']==.5 and error<.5,row
                report['nve'].append(row);write_json(out/'numerics.json',report)
                print(json.dumps({k:v for k,v in row.items() if k not in ('phase_space','frames')}),flush=True)
            finally:sim.close()
    final_hashes={name:digest(path) if path.is_file() else None for name,path in source_paths.items()}
    changed=[name for name in source_hashes if final_hashes[name]!=source_hashes[name]]
    report['source_integrity'].update(status='failed' if changed else 'passed', changed_files=changed)
    if changed:
        report.update(status='failed', error='Validation source/input/runtime files changed during the run')
        write_json(out/'numerics.json',report)
        raise RuntimeError(report['error']+': '+', '.join(changed))
    report['status']='passed-with-preparation-rejections' if any(x.get('status')=='preparation-rejected' for x in report['nve']) else 'passed'
    write_json(out/'numerics.json',report)


if __name__=='__main__':main()
