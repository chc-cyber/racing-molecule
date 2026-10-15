"""Explicit approved layered matrix; split expensive races between two queues."""
from pathlib import Path
import json
import statistics
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/timestep-study'
seeds=[2718281,1234567,31415927]
jobs=[]
def add(car,surface,temp,seed):
    for dt in [.5,1]:jobs.append(dict(car=car,surface=surface,temperature_K=temp,lattice_scale=2,dt_fs=dt,seed=seed,strength=.6,limit_ps=400,sensitivity_only=True))
for car in ['fullerene-c60','folded-eicosane','benzene-12','pyridine-11']:
    for seed in seeds:add(car,'fourier-au',150,seed)
for car in ['benzene-12','pyridine-11']:
    for seed in seeds:add(car,'periodic-au',150,seed)
for car in ['benzene-12','folded-eicosane']:add(car,'fourier-au',600,seeds[0])
# Greedy estimated-cost allocation avoids putting every slow 0.5fs C60 on one queue.
queues=[[],[]];costs=[0.,0.]
est={'fullerene-c60':400.,'folded-eicosane':500.,'benzene-12':15.,'pyridine-11':15.}
# Use screening costs when available: a flexible car can cost more than C60
# over a full race because its physical finish time is longer.
if (OUT/'nve.json').exists():
    screening=json.loads((OUT/'nve.json').read_text(encoding='utf-8'))['rows']
    for car in est:
        rows=[r for r in screening if r['car']==car and r['dt_fs']==.5 and r['surface']=='fourier-au' and r.get('wall_seconds')]
        if rows:est[car]=statistics.mean(r['wall_seconds']*50 for r in rows)*({'fullerene-c60':1.05,'folded-eicosane':2,'benzene-12':1.3,'pyridine-11':1.3}[car])
def cost(j):return est[j['car']]*(.5/j['dt_fs'])*(5 if j['surface']=='periodic-au' and j['car'] in ['benzene-12','pyridine-11'] else 1)
for job in sorted(jobs,key=cost,reverse=True):
    idx=min(range(2),key=lambda i:costs[i]);queues[idx].append(job);costs[idx]+=cost(job)
OUT.mkdir(parents=True,exist_ok=True)
(OUT/'race-jobs-all.json').write_text(json.dumps(jobs,indent=2)+'\n')
for i,q in enumerate(queues):(OUT/f'race-jobs-{i+1}.json').write_text(json.dumps(q,indent=2)+'\n')
print(json.dumps({'total':len(jobs),'queues':[len(q) for q in queues],'cost_weights':costs}))
