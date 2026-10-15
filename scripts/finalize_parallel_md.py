"""Record the tested deliverable and its evidence without changing model assets."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/parallel-md-validation'

def read(path):
    return json.loads((ROOT / path).read_text(encoding='utf-8'))

def main():
    summary = read('outputs/parallel-md-validation/summary.json')
    assert summary['status'] == 'complete-experimental'
    assert sum(g['trials'] for g in summary['timestep']['race_groups']) == 68
    assert read('outputs/wasm-preparation/package-validation.json')['passed']
    assert read('outputs/wasm-preparation/wasm-validation.json')['passed']
    assert read('outputs/parallel-md-validation/browser-ui.json')['passed']
    assert read('outputs/parallel-md-validation/v1-regression.json')['status']=='passed'
    assert read('outputs/timestep-study/interface.json')['status']=='passed'
    integrity=read('outputs/timestep-study/integrity.json')
    assert integrity['status']=='passed' and integrity['verified']==68
    paths = set()
    for pattern in [
        'wasm/src/*', 'wasm/tools/*atlas*', 'wasm/tools/*timestep*',
        'wasm/tools/*parallel*', 'scripts/*parallel_md*',
        'scripts/verify_atlas_ui.cjs', 'wasm/tools/package_web.py',
        'wasm/tools/verify_package.py', 'wasm/tools/build.py',
        'scripts/*timestep*', 'frontend/browser/*.mjs',
        'molecular_car/fourier_atlas.py', 'frontend/app.js',
        'README.md', 'wasm/README.md',
        'outputs/parallel-md-validation/*.json',
        'outputs/parallel-md-validation/*.md',
        'outputs/parallel-md-validation/*.png',
        'outputs/timestep-study/*.json', 'outputs/timestep-study/*.md',
        'outputs/timestep-study/races*/*.json',
        'outputs/timestep-study/races*/*.json.coords.gz',
        'outputs/c60-proxy-study/*.json', 'outputs/c60-proxy-study/*.md',
        'outputs/c60-proxy-study/races/browser-confirm-report.json',
        'outputs/c60-proxy-study/races/browser-trials/*',
        'wasm/dist/site-manifest.json', 'wasm/dist/browser/runtime/*',
        'wasm/dist/browser/assets/*',
        'outputs/wasm-preparation/package-validation.json',
        'outputs/wasm-preparation/wasm-validation.json',
        'outputs/wasm-preparation/molecular-car-web-preview.zip',
    ]:
        paths.update(p for p in ROOT.glob(pattern) if p.is_file())
    target = OUT / 'completion-manifest.json'
    paths.discard(target)
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(paths)}
    result = {
        'status': 'complete-experimental',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'hash_scope': 'Final files; individual numerical reports additionally record tested input/build hashes and verify them at run completion.',
        'default_surface': 'periodic-au', 'default_timestep_fs': 0.5,
        'optional_browser_surface': 'fourier-au, v2 three-chart atlas',
        'timestep_1fs': 'diagnostic sensitivity only; not promoted',
        'publicly_deployed': False,
        'evidence_sha256': hashes,
    }
    target.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(json.dumps({'status': result['status'], 'files': len(hashes)}))

if __name__ == '__main__':
    main()
