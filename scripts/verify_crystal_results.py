"""Recompute scores/guards from saved frames and coordinates, without new MD."""
from __future__ import annotations
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/crystal-contrast-guarded'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(config):
    return tuple(config[k] for k in ('phase','molecule','temperature_K','strength','lattice_scale','growth','seed','limit_ps'))


def main():
    paths=sorted((OUT/'browser-trials').glob('*.json'))
    rows=[json.loads(f.read_text(encoding='utf8')) for f in paths]
    jobs=json.loads((OUT/'browser-confirm-jobs.json').read_text(encoding='utf8'))
    progress=json.loads((OUT/'browser-confirm-report.json').read_text(encoding='utf8'))
    expected={identity(c) for c in jobs}
    assert len(expected)==len(jobs) and len(rows)==len(jobs)
    assert {identity(r['config']) for r in rows}==expected
    assert progress['status']=='complete' and progress['planned']==len(jobs)
    assert {identity(r['config']) for r in progress['trials']}==expected
    by_key={identity(r['config']):r for r in progress['trials']}
    for row in rows:
        # Reuse metadata may exist only in the combined report; physical records must agree.
        combined=by_key[identity(row['config'])]
        assert {k:v for k,v in combined.items() if k!='reused_baseline'}=={k:v for k,v in row.items() if k!='reused_baseline'}
    baseline=json.loads((OUT/'baseline-compatibility.json').read_text(encoding='utf8'))
    assert baseline['status']=='passed' and baseline['checked_baseline_trials']==42
    for file,sha in baseline['source_sha256'].items():assert digest(ROOT/file)==sha,file
    for file,sha in baseline['archived_source_sha256'].items():assert digest(OUT/'development-v2/source'/file)==sha,file
    for file,sha in baseline['baseline_trials_sha256'].items():assert digest(OUT/file)==sha,file
    provenance={str(f.relative_to(OUT)):digest(f) for f in paths}
    for file in ('browser-confirm-jobs.json','browser-confirm-report.json','baseline-compatibility.json'):
        provenance[file]=digest(OUT/file)
    keys=set();checks=[]
    for row in rows:
        c=row['config'];key=(c['molecule'],c['temperature_K'],c['strength'],c['lattice_scale'],c['growth'],c['seed'])
        assert key not in keys,key
        keys.add(key)
        if row['status']!='passed':checks.append({'key':row['key'],'status':'rejected','competitive_valid':False,'reason':row.get('error')});continue
        if c['lattice_scale']!=1:assert row['runtime']['surface']['corrugation_control_version']=='six-phase-atomic-crystal-contrast-v3'
        assert row['runtime']['timestep_fs']==.5 and row['runtime']['velocity_seed']==17329
        assert row['runtime']['langevin_seed']==c['seed']
        prep=row['runtime']['preparation'];assert prep['accepted'] and prep['global_force_norm_kcal_mol_A']<=1e-4
        s=row['samples'];steps=np.array([x['step'] for x in s]);centers=np.array([x['center_nm'] for x in s])
        assert np.array_equal(steps,np.arange(len(s))*100)
        assert np.allclose([x['time_ps'] for x in s],steps*.0005,atol=1e-12,rtol=0)
        assert np.isfinite(centers).all() and np.isfinite([x['energy_kcal_mol'] for x in s]).all()
        crossings=np.flatnonzero(centers[:,0]>=20.)
        assert (row['outcome']=='finish')==bool(len(crossings))
        if len(crossings):assert crossings[0]==len(s)-1 and row['finish_time_ps']==s[-1]['time_ps']
        else:assert row['finish_time_ps'] is None and row['observed_time_ps']==c['limit_ps']
        lift=maximum=0
        for frame in s[1:]:lift=lift+1 if frame['min_z_A']>8 else 0;maximum=max(maximum,lift)
        assert row['competitive_valid']==(maximum<100)
        assert row['sustained_lift']==(maximum>=100)
        assert abs(row['maximum_lift_duration_ps']-maximum*.05)<1e-8
        coord=json.loads(gzip.decompress((OUT/'browser-trials'/row['coordinates_file']).read_bytes()))
        xyz=np.array(coord['positions']).reshape(-1,row['atoms'],3);mass=np.array(coord['masses'])
        assert np.isfinite(xyz).all() and np.isfinite(mass).all()
        com=np.average(xyz,axis=1,weights=mass)/10
        index=np.array(coord['steps'])//100
        assert np.max(abs(com-centers[index]))<1e-10
        assert np.allclose(xyz[:,:,2].min(axis=1),[s[i]['min_z_A'] for i in index],atol=1e-12,rtol=0)
        assert np.allclose(xyz[:,:,2].max(axis=1),[s[i]['max_z_A'] for i in index],atol=1e-12,rtol=0)
        bonds=np.array(coord['bonds'],dtype=int)
        lengths=np.linalg.norm(xyz[:,bonds[:,0]]-xyz[:,bonds[:,1]],axis=2)
        ratios=lengths/lengths[0]
        assert ratios.min()>=.65 and ratios.max()<=1.45
        assert row['bond_ratio_range'][0]<=ratios.min()+1e-12 and row['bond_ratio_range'][1]>=ratios.max()-1e-12
        provenance['browser-trials/'+row['coordinates_file']]=digest(OUT/'browser-trials'/row['coordinates_file'])
        for file,sha in row['source_sha256'].items():
            path=OUT/file if file.startswith('inputs/') else ROOT/file
            if digest(path)==sha:continue
            if c['lattice_scale']==1 and file in baseline['archived_source_sha256']:
                assert sha==baseline['archived_source_sha256'][file],file
                assert f"browser-trials/{row['key']}.json" in baseline['baseline_trials_sha256']
            elif file=='wasm/tools/scan_crystal_races.mjs':
                # The final scanner only fixes future exactly-5-ps guard rounding.
                assert digest(OUT/'completed-scan-source'/file)==sha,file
                assert maximum==0, 'Changed lift predicate requires explicit retrospective classification'
            else:raise AssertionError('Unverified changed source: '+file)
        checks.append({'key':row['key'],'status':'passed','competitive_valid':maximum<100,
            'outcome':row['outcome'],'steps':int(steps[-1]),'saved_coordinate_frames':len(xyz),
            'min_atom_z_A':float(min(x['min_z_A'] for x in s)),
            'highest_min_atom_z_A':float(max(x['min_z_A'] for x in s)),
            'maximum_lift_ps':maximum*.05})
    report={'status':'passed','planned_trials':len(jobs),'checked_trials':len(rows),
            'valid_surface_races':sum(x['status']=='passed' and x['competitive_valid'] for x in checks),
            'valid_finishes':sum(x['status']=='passed' and x['competitive_valid'] and x['outcome']=='finish' for x in checks),
            'right_censored_trials':sum(x['status']=='passed' and x['outcome']!='finish' for x in checks),
            'checks':checks,'saved_coordinate_frames_and_scores_recomputed':True,
            'sampling':{'md_timestep_fs':.5,'score_energy_contact_sample_ps':.05,'full_coordinate_sample_ps':.5,
                'limitation':'Checks cover saved frames, not every integration step; terminal coordinates also retained'},
            'input_sha256':provenance,'verification_script_sha256':digest(Path(__file__)),
            'baseline_compatibility':'42 specifically hashed older scale=1 trials checked against archived sources and paired old/current WASM',
            'scan_guard_fix':'Original completed scanner archived; current lift duration uses integer frame count. Legacy traces require zero lift for the rounding-only compatibility exception.'}
    (OUT/'verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('checks','input_sha256')}))


if __name__=='__main__':main()
