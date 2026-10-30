# Written by Codex, an AI, at David's request; David did not review it line by line.
"""Render frozen aggregate output counts; no sequence data or new analysis."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
P=Path(__file__).resolve().parent
rows=list(csv.DictReader((P/'funnel.tsv').open(),delimiter='\t'))
ranks=list(csv.DictReader((P/'taxonomic-resolution.tsv').open(),delimiter='\t'))
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','figure.facecolor':'white','savefig.facecolor':'white'})
INK='#223345';MUTED='#5D6B79';LINE='#AAB7C2';BLUE='#2866A2';ORANGE='#BD601A'

def render(depth):
    fig,ax=plt.subplots(figsize=(12,10.1));fig.subplots_adjust(0,0,1,1)
    ax.set(xlim=(0,12),ylim=(0,10.1));ax.axis('off')
    data={r['stage']:r for r in rows if int(r['depth_million_pairs'])==depth}
    n=data['input_pairs']['units']
    ax.text(.5,9.72,'More output, fewer species-resolved representatives',fontsize=22,weight='bold',color=INK,va='top')
    subtitle=f'{n} sequencing runs · {depth} million read pairs per run · same input for both pipelines'
    if depth==25: subtitle='10 of the same 50 runs · 25 million read pairs per run · same input for both pipelines'
    ax.text(.5,9.21,subtitle,fontsize=13,color=MUTED,va='top')
    def box(x,y,w,h,color,fill,label,value,labelsize=13):
        ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.015,rounding_size=.09',edgecolor=color,facecolor=fill,linewidth=1.05,zorder=2))
        ax.text(x,y+.17,label,ha='center',va='center',fontsize=labelsize,color=INK,linespacing=1.18,zorder=3)
        ax.text(x,y-.19,f'{value:,}',ha='center',va='center',fontsize=21,color=color,weight='bold',zorder=3)
    def arrow(a,b):
        ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=11,color=LINE,lw=1.3,zorder=1))
    for arm,title,x,color,fill in [('production','Production',3.05,BLUE,'#F1F6FB'),('protein_variant','Protein variant',8.95,ORANGE,'#FCF4EC')]:
        ax.text(x,8.7,title,ha='center',va='center',fontsize=19,weight='bold',color=color)
        stages=[('input_pairs','Input read pairs',8.03),('screened_pairs','Passed screening',6.88),('after_fastp_pairs','After FASTP trimming and filtering',5.73),('output_rows','Final output read pairs',4.22),('joint_exemplars','Selected by both duplicate checks',3.07)]
        for i,(key,label,y) in enumerate(stages):
            box(x,y,5.05,.81,color,fill,label,int(data[key][arm]))
            if i:
                arrow((x,stages[i-1][2]-.43),(x,y+.43))
        ax.text(x,4.99,'Alignment, filters and record conversion',ha='center',va='center',fontsize=11.5,color=MUTED,bbox=dict(facecolor='white',edgecolor='none',pad=3))
        r=next(r for r in ranks if int(r['depth_million_pairs'])==depth and r['arm']==arm)
        resolved=int(r['species_resolved']);unresolved=int(r['above_species_or_unresolved'])
        assert resolved+unresolved==int(data['joint_exemplars'][arm])
        for xx,label,value in [(x-1.3,'Assigned to\na species',resolved),(x+1.3,'Above species\nor unresolved',unresolved)]:
            ax.plot([x,x,xx],[2.64,2.15,2.15],color=LINE,lw=1.3,zorder=1)
            arrow((xx,2.15),(xx,1.83))
            box(xx,1.3,2.4,1.02,color,fill,label,value,12.5)
    ax.text(.5,.34,'The bottom boxes partition the representatives above. All counts are read pairs; box sizes are not to scale.',fontsize=11.5,color=MUTED)
    fig.savefig(P/f'output-flow-{depth}m.png',dpi=180)
    fig.savefig(P/f'output-flow-{depth}m.svg')
    plt.close(fig)
for d in (5,25):render(d)
