"""Export optional v2 study assets; never overwrite deployed v1 assets."""
from pathlib import Path
import sys, json, gzip, hashlib, os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[key]='1'
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from molecular_car.fourier_atlas import FourierAtlasPotential, AXES

def main():
    out=ROOT/'outputs/c60-proxy-study';assets=out/'assets';assets.mkdir(parents=True,exist_ok=True)
    p=FourierAtlasPotential(['C'],np.array([[0.,0.,3.4]]),accelerator='numpy',
        parameter_path=ROOT/'.local/wasm/lammps/potentials/CHAu.ILP',fourier_cache_dir=out/'cache')
    modelpath=ROOT/'wasm/dist/browser/assets/model.json';model=json.loads(modelpath.read_text())
    assert p.parameter_sha256==model['surface']['parameter_sha256']
    raw=np.ascontiguousarray(p.atlas_coefficients[0],dtype='<c16').view('<f8').tobytes()
    packed=gzip.compress(raw,compresslevel=9,mtime=0);name='fourier-atlas-0.bin.gz';(assets/name).write_bytes(packed)
    descriptor=dict(path='assets/'+name,bytes=len(raw),download_bytes=len(packed),sha256=hashlib.sha256(raw).hexdigest())
    meta=dict(schema='fourier-saip-browser-v2',charts=3,chart_axes=AXES,
        parameter_sha256=p.parameter_sha256,base_model_sha256=hashlib.sha256(modelpath.read_bytes()).hexdigest(),
        nz=p.fourier_z_count,na=len(p.angular_nodes),modes=p.modes.tolist(),astep=p.angular_step,
        bound=p.chart_bound,blend=p.chart_blend,zstep=p.height_step,zstart=p.fourier_z_start,
        zmin=p.height_min,zmax=p.height_max,zblend=p.height_blend,arrays=[descriptor],
        storage='[species,chart,z,a,b,component,mode,re/im] float64 interleaved complex; lossless gzip')
    (assets/'fourier-atlas.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(meta),flush=True)
if __name__=='__main__':main()
