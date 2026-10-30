import csv,json,sys,os
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR",str(Path(__file__).resolve().parent/"matplotlib-cache"))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter,ScalarFormatter
r=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=True)
rows=list(csv.DictReader(open(r/'screen-burden/sample-metrics.tsv'),delimiter='\t'));arms=['nt','nt22','nt20','aa12','union'];labels=['NT24','NT22','NT20','AA12','NT24 ∪ AA12'];by={a:{x['sample']:int(x['candidate_pairs']) for x in rows if x['arm']==a} for a in arms};n=len(by['nt']);assert all(set(v)==set(by['nt']) for v in by.values())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
fig,axes=plt.subplots(1,2,figsize=(11.7,4.7),layout='constrained')
data=[[v/5_000_000 for v in by[a].values()] for a in arms]
parts=axes[0].boxplot(data,tick_labels=labels,showfliers=True,patch_artist=True,flierprops={'markersize':2.5,'alpha':.35},medianprops={'color':'black'})
for patch,color in zip(parts['boxes'],['#72849a','#8aa7c2','#3977aa','#d6813c','#359d87']):patch.set_facecolor(color);patch.set_alpha(.65)
axes[0].set_yscale('log');axes[0].yaxis.set_major_formatter(FuncFormatter(lambda x,pos:f"{x*100:g}%"));axes[0].set_ylabel('Raw pairs admitted to downstream processing');axes[0].set_title('Typical burden and high-volume outliers');axes[0].grid(axis='y',alpha=.18)
ids=list(by['nt']);x=[by['nt20'][s] for s in ids];y=[by['aa12'][s] for s in ids];axes[1].scatter(x,y,s=17,alpha=.45,color='#3977aa',edgecolors='none');lim=[max(1,min(x+y)*.7),max(x+y)*1.3];axes[1].plot(lim,lim,color='#666',linestyle='--',linewidth=1,label='Same candidate count');axes[1].set(xscale='log',yscale='log',xlim=lim,ylim=lim,xlabel='NT20 candidate pairs',ylabel='AA12 candidate pairs',title='Paired comparison on identical raw reads');axes[1].legend(frameon=False,loc='lower right');axes[1].grid(alpha=.15)
s='SRR36995860'
if s in by['nt']:
 axes[1].scatter([by['nt20'][s]],[by['aa12'][s]],s=50,facecolors='none',edgecolors='#ad3d2d',linewidths=1.5);axes[1].annotate('Miami outlier',xy=(by['nt20'][s],by['aa12'][s]),xytext=(6,-17),textcoords='offset points',fontsize=9,color='#ad3d2d')
fig.suptitle(f'Candidate volume across {n:,} public samples · 5 million read pairs per sample',fontsize=14)
for suffix in ['png','svg','pdf']:fig.savefig(out/f'candidate-burden.{suffix}')
plt.close(fig)
p=r/'repeat-development/comparison.json';q=r/'repeat-development/screen-summary.json'
if p.exists() and q.exists():
 comparison={x['sample']:x for x in json.loads(p.read_text())['results']};screens=json.loads(q.read_text());names={'SRR39694948':'Riverside','SRR35946889':'Ottumwa','SRR33038139':'Boston North','SRR39694932':'Chicago','SRR39694979':'NYC Hospital B','SRR36995860':'Miami outlier'}
 fig,axes=plt.subplots(1,2,figsize=(11.7,4.8),layout='constrained');pos=list(range(len(screens)))
 for i,s in enumerate(screens):
  axes[0].plot([s['original_union_pairs'],s['union_pairs']],[i,i],color='#aaa',lw=2);axes[0].scatter(s['original_union_pairs'],i,color='#d6813c',s=35,label='Original AA12 union' if i==0 else None);axes[0].scatter(s['union_pairs'],i,color='#359d87',s=35,label='Repeat-filtered union' if i==0 else None)
  c=comparison[s['sample']];axes[1].plot([c['original_AA_union_pairs'],c['repeat_filtered_union_pairs']],[i,i],color='#aaa',lw=2);axes[1].scatter(c['original_AA_union_pairs'],i,color='#d6813c',s=60,facecolors='none');axes[1].scatter(c['repeat_filtered_union_pairs'],i,color='#359d87',s=20)
 for ax in axes:ax.set_yticks(pos,[names[s['sample']] for s in screens]);ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
 axes[0].set_xscale('log');axes[0].set_xlabel('Candidate pairs');axes[0].set_title('Candidate reduction');axes[0].legend(frameon=False,fontsize=9);axes[1].set_xlabel('Final assigned read pairs');axes[1].set_title('Final assignments: identical IDs and taxonomy')
 fig.suptitle('Exploratory repeat filter · five development samples and the triggering outlier',fontsize=13)
 for suffix in ['png','svg','pdf']:fig.savefig(out/f'repeat-filter-development.{suffix}')
