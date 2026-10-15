"""Independent native six-query Fourier reference for the fused browser kernel."""
from pathlib import Path
import json,sys,os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='1'
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from molecular_car.molecule import prepare_molecule
from molecular_car.periodic_au import BONDED_NORMAL_POLICY
from molecular_car.fourier_au import FourierAuPotential
from molecular_car.lattice_corrugation import LatticeCorrugation

out=ROOT/'outputs/wasm-fourier-validation'
fixtures=json.loads((ROOT/'outputs/crystal-contrast-guarded/1p5-numerics.json').read_text())['wasm_fixtures']
rows=[]
for car in ['benzene-12','test-car-52','fullerene-c60','pyridine-11']:
    f=next(f for f in fixtures if f['filename']==car+'.mol' and f['growth']==1)
    mol=prepare_molecule(f['content'],f['filename'])
    p=FourierAuPotential(mol.atoms,mol.positions,bonds=mol.bonds,normal_policy=BONDED_NORMAL_POLICY,
        parameter_path=ROOT/'.local/wasm/lammps/potentials/CHAu.ILP',accelerator='numba',
        fourier_cache_dir=out/'cache')
    x=np.asarray(f['positions']).reshape(-1,3)
    cases=[('start',x)]
    for z in [2.8,2.95,3.1,7.85,8.]:
        a=x.copy();a[:,2]+=z-a[:,2].min();cases.append((f'height-{z}',a))
    for angle in [np.arctan(1.8),np.arctan(2.),np.pi/2]:
        c,s=np.cos(angle),np.sin(angle);a=(x-x.mean(0))@np.array([[c,0,-s],[0,1,0],[s,0,c]])
        a[:,2]+=3.4-a[:,2].min();a[:,0]-=103*p.period[0];a[:,1]+=97*p.period[1]
        cases.append((f'tilt-{angle}',a))
    for label,a in cases:
        for scale in [1.,1.5,2.]:
            e,force=LatticeCorrugation(p,scale).energy_and_forces(a)
            rows.append(dict(car=car,content=f['content'],label=label,scale=scale,positions=a.ravel().tolist(),energy=e,forces=force.ravel().tolist()))
    print(car,flush=True)
(out/'native-fixtures.json').write_text(json.dumps(rows)+'\n')
print('fixtures',len(rows),flush=True)
