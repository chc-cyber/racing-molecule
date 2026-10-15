"""Developer-only lossless export; optional browser assets load on demand."""
from pathlib import Path
import sys, json, gzip, hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from molecular_car.molecule import prepare_molecule
from molecular_car.fourier_au import FourierAuPotential
from molecular_car.periodic_au import BONDED_NORMAL_POLICY

def main():
    assets=ROOT/'wasm/dist/browser/assets'
    model=json.loads((assets/'model.json').read_text())
    file=ROOT/'examples/benzene.mol'
    mol=prepare_molecule(file.read_text(),file.name)
    proxy=FourierAuPotential(mol.atoms,mol.positions,bonds=mol.bonds,normal_policy=BONDED_NORMAL_POLICY,
        parameter_path=ROOT/'.local/wasm/lammps/potentials/CHAu.ILP',accelerator='numpy',
        fourier_cache_dir=ROOT/'outputs/wasm-fourier-validation/cache')
    s=model['surface']
    assert proxy.parameter_sha256==s['parameter_sha256']
    assert np.allclose(proxy.period,s['period'],rtol=0,atol=1e-12)
    assert np.allclose(proxy.layer_heights,s['layer_heights_angstrom'],rtol=0,atol=1e-12)
    # Split species so the 75 MB native cache need not be copied simultaneously.
    descriptors=[]
    for code in range(2):
        raw=np.ascontiguousarray(proxy.coefficients[code],dtype='<c16').view('<f8').tobytes()
        packed=gzip.compress(raw,compresslevel=9,mtime=0)
        name=f'fourier-{code}.bin.gz';(assets/name).write_bytes(packed)
        descriptors.append(dict(path='assets/'+name,bytes=len(raw),download_bytes=len(packed),sha256=hashlib.sha256(raw).hexdigest()))
    info=dict(schema='fourier-saip-browser-v1',parameter_sha256=proxy.parameter_sha256,
        base_model_sha256=hashlib.sha256((assets/'model.json').read_bytes()).hexdigest(),
        nz=proxy.fourier_z_count,na=len(proxy.angular_nodes),modes=proxy.modes.tolist(),
        astep=proxy.angular_step,bound=proxy.chart_bound,blend=proxy.chart_blend,
        zstep=proxy.height_step,zstart=proxy.fourier_z_start,zmin=proxy.height_min,zmax=proxy.height_max,zblend=proxy.height_blend,
        arrays=descriptors,storage='float64 interleaved complex; lossless gzip')
    (assets/'fourier.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info))
if __name__=='__main__':main()
