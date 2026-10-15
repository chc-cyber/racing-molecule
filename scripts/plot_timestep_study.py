from pathlib import Path
import sys,os,json
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'.local/analysis-libs'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.local/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'outputs/timestep-study'
rows=json.loads((OUT/'nve.json').read_text(encoding='utf-8'))['rows']
peak_rows=json.loads((OUT/'sampling.json').read_text(encoding='utf-8'))['rows'] if (OUT/'sampling.json').exists() else []
colors={'benzene-12':'#2878B5','pyridine-11':'#F28E2B','fullerene-c60':'#59A14F','folded-eicosane':'#B44EAA'}
names={'benzene-12':'Benzene','pyridine-11':'Pyridine','fullerene-c60':'C60','folded-eicosane':'Folded eicosane'}
fig,axes=plt.subplots(1,2,figsize=(10,4.4),sharey=True,constrained_layout=True)
for ax,T in zip(axes,[150,600]):
 for car,color in colors.items():
  for surface,ls in [('periodic-au','-'),('fourier-au','--')]:
   for scale,alpha in [(1.5,.4),(2,1)]:
    values=sorted([r for r in rows if r['car']==car and r['surface']==surface and r['temperature_K']==T and r['lattice_scale']==scale and r.get('relative_to_initial_kinetic') is not None],key=lambda r:r['dt_fs'])
    if values:ax.plot([r['dt_fs'] for r in values],[100*r['relative_to_initial_kinetic'] for r in values],ls,marker='o',markersize=3,color=color,alpha=alpha,label=names[car] if scale==2 and surface=='periodic-au' else None)
 peaks=[r for r in peak_rows if r['temperature_K']==T and r['surface']=='periodic-au']
 if peaks:ax.scatter([r['dt_fs'] for r in peaks],[100*r['relative_every_step'] for r in peaks],marker='x',s=28,c='black',zorder=10,label='Every-step spot check (0.5 ps)')
 ax.axhline(1,color='#555',ls=':',lw=1);ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_xticks([.25,.5,1,2],['0.25','0.5','1','2']);ax.set_title(f'{T} K');ax.set_xlabel('Integration timestep (fs)');ax.grid(alpha=.15)
axes[0].set_ylabel('Maximum |E(t) - E(0)| / initial KE (%)')
axes[0].legend(frameon=False,fontsize=8)
fig.suptitle('2 ps paired-state NVE, sampled every 10 fs; dotted line = 1% budget')
fig.supxlabel('Solid: periodic reference; dashed: original Fourier; light/dark: corrugation 1.5 / 2.0',fontsize=8)
for ext in ['png','svg']:fig.savefig(OUT/f'nve-timestep.{ext}',dpi=180)
plt.close(fig)

trial_rows={}
for file in OUT.glob('races*.json'):
 report=json.loads(file.read_text(encoding='utf-8'))
 if isinstance(report,dict):
  for row in report.get('rows',[]):trial_rows[row['key']]=row
if trial_rows:
 import statistics
 fig,axes=plt.subplots(2,1,figsize=(11,7),constrained_layout=True,gridspec_kw={'height_ratios':[2,1]})
 for ax,T in zip(axes,[150,600]):
  keys=sorted({(r['config']['car'],r['config']['surface']) for r in trial_rows.values() if r['config']['temperature_K']==T})
  for dt,offset,color in [(.5,-.16,'#2878B5'),(1,.16,'#F28E2B')]:
   for index,key in enumerate(keys):
    group=[r for r in trial_rows.values() if (r['config']['car'],r['config']['surface'])==key and r['config']['temperature_K']==T and r['config']['dt_fs']==dt]
    valid=[r for r in group if r.get('competitive_valid',False) and r.get('outcome')=='finish']
    values=[r['finish_time_ps'] for r in valid]
    if values:
     mean=statistics.mean(values);sd=statistics.stdev(values) if len(values)>1 else 0
     ax.bar(index+offset,mean,.29,color=color,alpha=.65,label=f'{dt:g} fs' if index==0 else None,yerr=sd,capsize=3)
     ax.scatter([index+offset+(i-(len(values)-1)/2)*.025 for i in range(len(values))],values,s=15,c=color,zorder=5)
     ax.text(index+offset,max(mean+sd,max(values))+4,f'{len(values)}/{len(group)}',ha='center',fontsize=8)
    censored=[r for r in group if r.get('right_censored',False)]
    if censored:
     ax.scatter([index+offset]*len(censored),[r['observed_time_ps'] for r in censored],marker='^',s=40,c=color,zorder=6)
     ax.text(index+offset,400+8,f'{len(censored)} cens.',ha='center',fontsize=8)
    invalid=[r for r in group if r.get('outcome')=='finish' and not r.get('competitive_valid',False)]
    if invalid:ax.scatter([index+offset]*len(invalid),[r['finish_time_ps'] for r in invalid],marker='x',s=40,c=color,zorder=6)
  ax.set_xticks(range(len(keys)),[names[car]+'\n'+('Periodic' if surface=='periodic-au' else 'Original Fourier') for car,surface in keys],fontsize=9)
  ax.set_ylabel('Finish time (ps)');ax.set_title(f'{T} K, corrugation 2.0, drive 0.6');ax.grid(axis='y',alpha=.15);ax.legend(frameon=False)
 fig.suptitle('20 nm timestep sensitivity: bars = valid-finish mean +/- sample SD')
 fig.supxlabel('Dots: valid finishes; triangles: right-censored; labels: valid/total. Sparse seeds cannot establish equivalence.',fontsize=9)
 for ext in ['png','svg']:fig.savefig(OUT/f'race-timestep.{ext}',dpi=180)
 plt.close(fig)

# Equal physical-window behavior, so different finish times cannot bias exposure.
summary=json.loads((OUT/'summary.json').read_text(encoding='utf-8'))
if summary.get('paired_bootstrap',{}).get('rows'):
 fig,axes=plt.subplots(1,2,figsize=(10,4.5),constrained_layout=True)
 for ax,metric,label,factor in zip(axes,['first_50ps_lateral_rms_nm','first_50ps_low_progress_fraction'],['Lateral RMS over first 50 ps (nm)','Low-forward fraction over first 50 ps (%)'],[1,100]):
  for i,car in enumerate(['benzene-12','pyridine-11']):
   row=next((r for r in summary['paired_bootstrap']['rows'] if r['car']==car and r['metric']==metric),None)
   if row is None:continue
   for field,offset,color,dt in [('half_fs',-.16,'#2878B5',.5),('one_fs',.16,'#F28E2B',1)]:
    stat=row[field]
    ax.bar(i+offset,stat['mean']*factor,.29,yerr=(stat['sd'] or 0)*factor,capsize=4,color=color,alpha=.7,label=f'{dt:g} fs' if i==0 else None)
  ax.set_xticks([0,1],['Benzene','Pyridine']);ax.set_ylabel(label);ax.grid(axis='y',alpha=.15);ax.legend(frameon=False)
 fig.suptitle('Timestep behavior sensitivity: 150 K, corrugation 2, original Fourier')
 fig.supxlabel('10 seeds per timestep; bars = mean +/- sample SD. Fixed observation window = 50 ps.',fontsize=9)
 for ext in ['png','svg']:fig.savefig(OUT/f'behavior-timestep.{ext}',dpi=180)
 plt.close(fig)
