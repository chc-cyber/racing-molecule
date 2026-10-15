"""Combine the completed 42-race supplement with the preserved 105-race study."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/crystal-contrast-guarded'
sys.path.insert(0,str(ROOT/'scripts'))
from scan_structure_races import write_json
from verify_crystal_results import identity


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def verify_batch(report,jobs,expected_count,scale=None):
    assert report['status']=='complete' and report['planned']==expected_count
    assert len(jobs)==len(report['trials'])==expected_count
    planned={identity(c) for c in jobs}
    observed={identity(r['config']) for r in report['trials']}
    assert len(planned)==expected_count and planned==observed
    if scale is not None:assert all(c['lattice_scale']==scale for c in jobs)
    for row in report['trials']:
        saved=read(OUT/'browser-trials'/f"{row['key']}.json")
        assert {k:v for k,v in row.items() if k!='reused_baseline'}=={k:v for k,v in saved.items() if k!='reused_baseline'}


def main():
    old_dir=OUT/'pre-1p5-results'
    old=read(old_dir/'browser-confirm-report.json');old_jobs=read(old_dir/'browser-confirm-jobs.json')
    current=read(OUT/'browser-confirm-report.json')
    # Idempotent after combining: retain the separate raw supplement.
    if current.get('planned')==147:
        new=read(OUT/'1p5-browser-confirm-report.json');new_jobs=read(OUT/'1p5-browser-confirm-jobs.json')
    else:
        new=current;new_jobs=read(OUT/'browser-confirm-jobs.json')
    verify_batch(old,old_jobs,105);verify_batch(new,new_jobs,42,1.5)
    assert old['protocol']==new['protocol'],'The supplement must preserve the comparison protocol'
    trials=old['trials']+new['trials'];jobs=old_jobs+new_jobs
    assert len({identity(c) for c in jobs})==147
    if current.get('planned')!=147:
        write_json(OUT/'1p5-browser-confirm-report.json',new)
        write_json(OUT/'1p5-browser-confirm-jobs.json',new_jobs)
    combined=dict(new,planned=147,trials=trials,wall_seconds=old['wall_seconds']+new['wall_seconds'],
        aggregation={'original_trials':105,'supplement_trials':42,
            'source_note':'Each trial retains its original source hash; header sources describe the latest run.',
            'wall_seconds_note':'Sum of the two completed batch durations, not simultaneous CPU time.',
            'reports':{name:hashlib.sha256((OUT/name).read_bytes()).hexdigest() for name in (
                'pre-1p5-results/browser-confirm-report.json','1p5-browser-confirm-report.json')}})
    write_json(OUT/'browser-confirm-jobs.json',jobs)
    write_json(OUT/'browser-confirm-report.json',combined)
    print(json.dumps({'status':'complete','prior_trials':105,'added_trials':42,'combined_trials':147}))


if __name__=='__main__':main()
