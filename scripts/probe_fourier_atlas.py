"""Isolate angular interpolation / xy truncation before full atlas export."""
from pathlib import Path
import sys, json, os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']: os.environ[key]='1'
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from molecular_car.periodic_au import PeriodicAuPotential, _anisotropic, EV_TO_KCAL_MOL
from molecular_car.fourier_atlas import FourierAtlasPotential

p = FourierAtlasPotential.__new__(FourierAtlasPotential)
PeriodicAuPotential.__init__(p,['C'],np.array([[0.,0.,3.4]]),accelerator='numpy',
    parameter_path=ROOT/'.local/wasm/lammps/potentials/CHAu.ILP')
p.angular_nodes=np.linspace(-1.2,1.2,13);p.angular_step=.2;p.chart_bound=1.2;p.chart_blend=.2
p.fourier_nx,p.fourier_ny=8,18;p.modes=np.asarray([(m,n) for m in range(3) for n in range(-4,5) if (m+n)%2==0 and(m>0 or n>=0)])
p.wavevectors=p.modes*(2*np.pi/p.period);p.mode_multiplicity=np.where(np.any(p.modes!=0,axis=1),2.,1.)
p.height_step=.1
rng=np.random.default_rng(6014)
rows=[]
for z in [2.4,2.8,3.2,4.,7.,10.]:
    coeff=p._coefficients_at_height(z)
    p.atlas_coefficients=np.repeat(coeff[:,:,None],5,axis=2)
    p.fourier_z_start=z-.2
    for idx in range(100):
        n=rng.normal(size=3);n/=np.linalg.norm(n)
        if idx<3:n=np.eye(3)[idx]
        x=np.array([*rng.uniform(0,p.period),z])
        e,f,ng=p._proxy(0,x,n,0,0)
        ids,d,r=p._pairs(x[None],np.array([0]))
        ee,ff,gg=_anisotropic(d,r,np.repeat(n[None],len(r),axis=0),p.parameters['C','Au'])
        rows.append(dict(z=z,n=n.tolist(),energy_error=abs(e-ee.sum())*EV_TO_KCAL_MOL,
            lateral_force_error=float(np.max(np.abs(f[:2]-ff.sum(0)[:2]))*EV_TO_KCAL_MOL),
            normal_gradient_error=float(np.max(np.abs(ng-(gg.sum(0)-n*np.dot(n,gg.sum(0)))))*EV_TO_KCAL_MOL),energy=float(ee.sum()*EV_TO_KCAL_MOL)))
    print(z,max(r['energy_error'] for r in rows if r['z']==z),max(r['lateral_force_error'] for r in rows if r['z']==z),max(r['normal_gradient_error'] for r in rows if r['z']==z),flush=True)
out=ROOT/'outputs/c60-proxy-study';out.mkdir(exist_ok=True)
(out/'angular-probe.json').write_text(json.dumps(rows,indent=2)+'\n')

