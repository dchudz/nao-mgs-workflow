"""Render standalone scientific figures from saved experiment tables."""
import argparse,csv,json,statistics
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('results');p.add_argument('output');a=p.parse_args();r=Path(a.results);out=Path(a.output);out.mkdir(exist_ok=True)
s=json.load(open(r/'summary.json'));rows=list(csv.DictReader(open(r/'sample-metrics.tsv'),delimiter='\t'));resources=json.load(open(r/'resource-summary.json')) if (r/'resource-summary.json').exists() else None
arms=['nt22','nt20','aa12','union'];labels=['NT22','NT20','AA12','NT24 + AA12'];colors=['#737373','#0072B2','#D55E00','#009E73']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
fig,ax=plt.subplots(1,3,figsize=(13,4.4),layout='constrained')
added=[s['arms'][k]['totals']['added_pairs'] for k in arms];lost=[s['arms'][k]['totals']['lost_pairs'] for k in arms]
ax[0].bar(range(4),added,color=colors);ax[0].bar(range(4),[-x for x in lost],color=colors,alpha=.45,hatch='//');ax[0].axhline(0,color='#555',lw=.7);ax[0].set_xticks(range(4),labels,rotation=25,ha='right');ax[0].set_title('A  Added and lost final read pairs',loc='left');ax[0].set_ylabel('Pairs relative to NT24')
for i,(gain,loss) in enumerate(zip(added,lost)):
 ax[0].annotate(f'+{gain:,}',(i,gain),xytext=(0,4),textcoords='offset points',ha='center',fontsize=9)
 if loss:ax[0].annotate(f'−{loss:,}',(i,-loss),xytext=(0,-12),textcoords='offset points',ha='center',fontsize=9)
ax[0].margins(y=.25)
by={arm:{x['sample']:int(x['added_pairs'])-int(x['lost_pairs']) for x in rows if x['arm']==arm} for arm in ['nt20','union']};names=sorted(by['nt20']);assert names==sorted(by['union'])
x=[by['nt20'][n] for n in names];y=[by['union'][n] for n in names];top=max(x+y+[1]);ax[1].scatter(x,y,s=18,alpha=.45,color='#009E73',edgecolors='none');ax[1].plot([0,top],[0,top],ls='--',lw=.8,color='#555');ax[1].set_xlabel('NT20 net added pairs');ax[1].set_ylabel('NT24 + AA12 net added pairs');ax[1].set_title('B  Paired sample comparison',loc='left');ax[1].set_xscale('symlog',linthresh=1);ax[1].set_yscale('symlog',linthresh=1);ax[1].grid(alpha=.15)
if resources:
 all_arms=['nt']+arms;rl=['NT24']+labels;scr=[resources['arms'][k]['screen_CPU_seconds']['median_per_sample'] for k in all_arms];post=[resources['arms'][k]['downstream_CPU_seconds_estimate']['median_per_sample'] for k in all_arms]
 ax[2].bar(range(5),scr,color='#56B4E9',label='Screen: GNU time');ax[2].bar(range(5),post,bottom=scr,color='#777777',label='Downstream: trace estimate');ax[2].set_xticks(range(5),rl,rotation=25,ha='right');ax[2].set_ylabel('Median CPU seconds per stage');ax[2].legend(frameon=False,fontsize=8);ax[2].set_title('C  Computational burden',loc='left')
else:ax[2].set_axis_off()
n=s['completed_samples'];fig.suptitle(f'Amino-acid screening with nucleotide downstream processing | {n:,} samples',fontsize=13)
fig.supxlabel('Five million input pairs per sample. Counts are pipeline assignments, not known-truth detections.\nCPU excludes extraction, indexing, union construction, orchestration and idle time.',fontsize=9)
for ext in ['png','svg','pdf']:fig.savefig(out/f'nucleotide-comparison.{ext}',dpi=180,bbox_inches='tight')
plt.close(fig)
if (r/'clustered-summary.json').exists():
 c=json.load(open(r/'clustered-summary.json'));means=[c['arms'][k]['mean_net_gain_per_sample'] for k in arms];intervals=[c['arms'][k]['location_cluster_bootstrap_95'] for k in arms]
 fig,ax=plt.subplots(figsize=(6,3.7),layout='constrained')
 for i,(m,ci,color) in enumerate(zip(means,intervals,colors)):
  ax.plot([ci[0],ci[1]],[i,i],color=color,lw=2);ax.scatter([m],[i],color=color,s=35,zorder=3)
 ax.axvline(0,color='#555',lw=.7);ax.set_yticks(range(4),labels);ax.invert_yaxis();ax.set_xlabel('Mean net added final read pairs per sample');ax.set_title(f'Paired changes in {n:,} samples');fig.supxlabel('Secondary 95% bootstrap intervals resample whole location labels.\nDescriptive for this cohort; locations can share cities.',fontsize=9)
 for ext in ['png','svg','pdf']:fig.savefig(out/f'clustered-effects.{ext}',dpi=180,bbox_inches='tight')
 plt.close(fig)
print(json.dumps({'samples':n,'figures':[str(p) for p in sorted(out.glob('*'))]}))
