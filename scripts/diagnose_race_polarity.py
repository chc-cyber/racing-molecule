"""Fixed-charge model dipoles and initial external-force components, not MD."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from molecular_car.molecule import prepare_molecule
from molecular_car.fields import GRADIENT_FIELD_SLOPE,GRADIENT_BASE_FIELD,GRADIENT_POLARIZATION_GAIN


def main():
    records=[]
    for folder in ('structure-race-scan','structure-race-polar'):
        out=ROOT/'outputs'/folder;manifest=json.loads((out/'inputs/manifest.json').read_text(encoding='utf8'))
        for m in manifest['fixtures']:
            path=out/'inputs'/m['file'];p=prepare_molecule(path.read_text(encoding='utf8'),path.name,start_pose={'mode':'auto'})
            charges=np.array([a['charge'] for a in p.atoms]);masses=np.array([a['mass'] for a in p.atoms]);alpha=np.array([a['polarizability'] for a in p.atoms])
            mu=np.sum(charges[:,None]*(p.positions-np.average(p.positions,axis=0,weights=masses)),axis=0)
            rows=[]
            for file in sorted((out/'trials').glob(f"confirm-{m['id']}-*.json")):
                r=json.loads(file.read_text(encoding='utf8'))
                if r['status']!='passed':continue
                x=np.load(out/'trials'/r['coordinates_file'])['positions_A'][0]
                relaxed=np.sum(charges[:,None]*(x-np.average(x,axis=0,weights=masses)),axis=0)
                c=r['config'];g=GRADIENT_FIELD_SLOPE*c['strength'];base=GRADIENT_BASE_FIELD*c['strength']
                qforce=float(np.sum(charges*g*x[:,2]));induced=float(np.sum(GRADIENT_POLARIZATION_GAIN*alpha*g*(base+g*x[:,0])))
                rows.append({'temperature_K':c['temperature_K'],'strength':c['strength'],'seed':c['seed'],
                    'initial_model_dipole_D':float(np.linalg.norm(relaxed)*4.8032047),
                    'charge_forward_force':qforce,'induced_forward_force':induced,
                    'charge_to_induced_forward_force_ratio':qforce/induced if induced else None})
            records.append({'molecule':m['id'],'formula':m['formula'],'atoms':len(p.atoms),'mass_u':float(masses.sum()),
                'alpha_sum':float(alpha.sum()),'alpha_per_mass':float(alpha.sum()/masses.sum()),
                'net_charge_e':float(charges.sum()),'vacuum_model_dipole_D':float(np.linalg.norm(mu)*4.8032047),
                'vacuum_model_dipole_vector_e_A':mu.tolist(),'initial_force_samples':rows})
    report={'qualification':'Fixed Gasteiger point charges, not measured dipoles. Force decomposition excludes surface/thermal forces; charges can affect orientation/contact indirectly.',
            'force_units':'kcal/mol/angstrom','molecules':records}
    target=ROOT/'outputs/structure-race-scan/polarity-diagnostics.json';target.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps([{k:r[k] for k in ('molecule','mass_u','vacuum_model_dipole_D','alpha_per_mass')} for r in records],indent=2))


if __name__=='__main__':main()
