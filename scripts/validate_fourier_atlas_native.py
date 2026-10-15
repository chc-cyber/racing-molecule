"""Independent full-surface comparisons, Cartesian FD, and atlas coverage."""
from pathlib import Path
import sys,json,os,time
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from molecular_car.molecule import prepare_molecule
from molecular_car.periodic_au import PeriodicAuPotential,BONDED_NORMAL_POLICY
from molecular_car.fourier_atlas import FourierAtlasPotential
from molecular_car.lattice_corrugation import LatticeCorrugation
out=ROOT/'outputs/c60-proxy-study'
fixtures=json.loads((ROOT/'outputs/crystal-contrast-guarded/1p5-numerics.json').read_text())['wasm_fixtures']
extra='--extra' in sys.argv
rows=[];fd=[];exports=[]
for car in (['benzene-12'] if extra else ['benzene-12','fullerene-c60','pyridine-11']):
    f=next(f for f in fixtures if f['filename']==car+'.mol' and f['growth']==1)
    mol=prepare_molecule(f['content'],f['filename'])
    args=dict(bonds=mol.bonds,normal_policy=BONDED_NORMAL_POLICY,parameter_path=ROOT/'.local/wasm/lammps/potentials/CHAu.ILP')
    atlas=FourierAtlasPotential(mol.atoms,mol.positions,**args,accelerator='numpy',fourier_cache_dir=out/'cache')
    exact=PeriodicAuPotential(mol.atoms,mol.positions,accelerator='numba',**args)
    x=np.asarray(f['positions']).reshape(-1,3)
    cases=[('start',x)]
    for z in [2.4,2.45,2.55,2.7,2.8,3.1,10.9,11.05,11.2]:
        a=x.copy();a[:,2]+=z-a[:,2].min();cases.append((f'height-{z}',a))
    for tilt in [np.arctan(1.),np.arctan(1.1),np.pi/2]:
        c,s=np.cos(tilt),np.sin(tilt);a=(x-x.mean(0))@np.array([[c,0,-s],[0,1,0],[s,0,c]])
        a[:,2]+=3.2-a[:,2].min();a[:,0]+=107*atlas.period[0];a[:,1]-=101*atlas.period[1]
        cases.append((f'tilt-{tilt}',a))
    if extra:
        cases=[]
        for label,normal in [('three-chart',[1.,1.,1.]),('side-y',[0.,1.,0.]),('negative-chart',[-1.,-1.,-1.])]:
            normal=np.array(normal);normal/=np.linalg.norm(normal)
            v=np.cross([0.,0.,1.],normal);c=normal[2]
            cross=np.array([[0.,-v[2],v[1]],[v[2],0.,-v[0]],[-v[1],v[0],0.]])
            rotation=np.eye(3)+cross+cross@cross/(1+c)
            a=(x-x.mean(0))@rotation.T;a[:,2]+=2.8-a[:,2].min()
            cases.append((label,a))
    for label,a in cases:
        for scale in [1.,2.]:
            p=LatticeCorrugation(atlas,scale);q=LatticeCorrugation(exact,scale)
            e,force=p.energy_and_forces(a);ee,ff=q.energy_and_forces(a)
            row=dict(car=car,label=label,scale=scale,energy_error=abs(e-ee),max_force_error=float(np.max(abs(force-ff))),
                rms_force_error=float(np.sqrt(np.mean((force-ff)**2))),counts=atlas.last_evaluation_counts.copy())
            rows.append(row)
            exports.append(dict(car=car,content=f['content'],positions=a.ravel().tolist(),scale=scale,energy=e,forces=force.ravel().tolist(),label=label))
            # Full Cartesian FD of benzene includes both chart overlap and height blend;
            # plus all 180 C60 components on its realistic contact configuration.
            if scale==2 and (extra or (car=='benzene-12' and label in ['height-2.55','tilt-0.8329812666744317','height-11.05']) or(car=='fullerene-c60' and label=='start')):
                eps=1e-5;numeric=np.zeros_like(a)
                for index in np.ndindex(a.shape):
                    plus=a.copy();minus=a.copy();plus[index]+=eps;minus[index]-=eps
                    numeric[index]=-(p.energy_and_forces(plus)[0]-p.energy_and_forces(minus)[0])/(2*eps)
                fd.append(dict(car=car,label=label,max_error=float(np.max(abs(numeric-force))),rms_error=float(np.sqrt(np.mean((numeric-force)**2)))))
                print('FD',fd[-1],flush=True)
        print(car,label,rows[-1]['max_force_error'],flush=True)
    del atlas
result=dict(cases=rows,finite_differences=fd,max_energy_error=max(r['energy_error'] for r in rows),max_force_error=max(r['max_force_error'] for r in rows),
    max_fd_error=max(r['max_error'] for r in fd))
(out/('native-extra-validation.json' if extra else 'native-validation.json')).write_text(json.dumps(result,indent=2)+'\n')
(out/('native-extra-fixtures.json' if extra else 'native-fixtures.json')).write_text(json.dumps(exports)+'\n')
print(json.dumps({k:v for k,v in result.items() if not isinstance(v,list)}),flush=True)
