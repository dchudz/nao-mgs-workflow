# Written by Codex, an AI, at David's request; David did not review it line by line.
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
P=Path(__file__).parent
BLUE='#2866A2';ORANGE='#D87526';GREEN='#39816B';RED='#B54B57';INK='#263442';MUTED='#667684';GRID='#E6EBEF'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'text.color':INK,'axes.labelcolor':MUTED,'xtick.color':MUTED,'ytick.color':INK,'axes.spines.top':False,'axes.spines.right':False,'axes.spines.left':False,'axes.edgecolor':GRID,'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','svg.fonttype':'none'})
def save(fig,name):
 fig.savefig(P/(name+'.png'),dpi=200,bbox_inches='tight',pad_inches=.22)
 fig.savefig(P/(name+'.svg'),bbox_inches='tight',pad_inches=.22)
 plt.close(fig)
# Waterfalls preserve counts and show how gains and losses sum to the final count.
fig,axs=plt.subplots(1,2,figsize=(10.5,4.8));fig.subplots_adjust(left=.08,right=.98,top=.74,bottom=.22,wspace=.28)
fig.text(.02,.97,'More calls meet the threshold only in the protein variant',fontsize=16,weight='bold',va='top')
fig.text(.02,.89,'Run/species calls at the same three-exemplar threshold',fontsize=11,color=MUTED)
for ax,title,vals in zip(axs,['5M pairs × 50 runs','25M pairs × 10 runs'],[(661,391,85,-79,1058),(231,152,26,-29,380)]):
 start,gain,prom,loss,end=vals;levels=[start,start+gain,start+gain+prom,end]
 bottoms=[0,start,start+gain,end,0];heights=[start,gain,prom,abs(loss),end]
 colors=[BLUE,GREEN,GREEN,RED,ORANGE]
 for i,(bottom,height,col) in enumerate(zip(bottoms,heights,colors)):
  ax.bar(i,height,bottom=bottom,color=col,width=.62,zorder=3)
  label=f'{vals[i]:,}' if i in (0,4) else f'{vals[i]:+,}'.replace('-', '−')
  top=bottom+height;ax.text(i,top+end*.035,label,ha='center',fontsize=11,weight='bold',color=col)
 for i,level in enumerate(levels):ax.plot([i+.31,i+1-.31],[level,level],ls='--',lw=1,color='#AEBBC5')
 ax.set_title(f'{title}\nVariant has {end-start:,} more calls',loc='left',fontsize=12,pad=16)
 ax.set_xticks(range(5),['Production','Variant only\n(no support)','Variant only\n(low support)','Production\nonly','Variant'])
 ax.tick_params(axis='x',length=0,pad=9,labelsize=9)
 ax.yaxis.grid(True,color=GRID);ax.set_axisbelow(True);ax.set_ylim(0,max(start+gain+prom,end)*1.18)
 ax.tick_params(axis='y',length=0,labelsize=9)
axs[0].set_ylabel('Run/species calls')
fig.text(.02,.035,'Zero/low support refer to production. Panels use different scales. Calls are not confirmed detections.',fontsize=9.5,color=MUTED)
save(fig,'call-changes')
# Two non-complementary validation measures, plotted separately with exact values.
fig,axs=plt.subplots(1,2,figsize=(10.5,4.9));fig.subplots_adjust(left=.105,right=.97,top=.68,bottom=.25,wspace=.28)
fig.text(.02,.97,'DNA validation differs across call categories',fontsize=17,weight='bold',va='top')
handles=[Line2D([],[],marker='o',ls='',color=BLUE,label='Production · both pipelines',markersize=7),Line2D([],[],marker='o',ls='',color=ORANGE,label='Variant · both pipelines',markersize=7),Line2D([],[],marker='D',ls='',markerfacecolor='white',markeredgecolor=ORANGE,markeredgewidth=2,label='Variant only · zero support',markersize=7)]
fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.01,.91),ncol=3,frameon=False,fontsize=10,handletextpad=.4,columnspacing=1.2)
series=[[[.74,6.59,32.42],[.18,7.16,39.81]],[[49.26,53.22,.35],[49.04,51.11,.29]]]
for k,(ax,data) in enumerate(zip(axs,series)):
 for y,vals in zip([1,0],data):
  ax.plot(vals[:2],[y,y],color='#CED6DD',lw=2,zorder=1)
  for j,(val,c,marker) in enumerate(zip(vals,[BLUE,ORANGE,ORANGE],['o','o','D'])):
   ax.scatter(val,y,s=64,marker=marker,facecolor='white' if j==2 else c,edgecolor=c,lw=1.8,zorder=3)
   offset=13 if j!=1 else -22
   ax.annotate(f'{val:.2f}%',(val,y),xytext=(0,offset),textcoords='offset points',ha='center',fontsize=10,color=c,weight='bold')
 ax.set_yticks([1,0],['5M','25M']);ax.tick_params(axis='y',length=0,pad=12)
 ax.set_ylim(-.42,1.48);ax.set_xlim(-3,57 if k==0 else 66)
 ticks=[0,10,20,30,40,50] if k==0 else [0,20,40,60]
 ax.set_xticks(ticks,[f'{v}%' for v in ticks]);ax.xaxis.grid(True,color=GRID);ax.set_axisbelow(True)
 ax.set_title('No alignment' if k==0 else 'Compatible assignment',fontsize=12,loc='left',pad=15)
 ax.set_xlabel('Among all checked rows' if k==0 else 'Among rows with a numeric distance',fontsize=10,labelpad=12)
fig.text(.02,.09,'Different denominators: the panels are not complementary rates. Compatibility means zero taxonomic distance.',fontsize=9.3,color=MUTED)
fig.text(.02,.035,'Both pipelines: the call meets the threshold in both. Zero support: no production species-level support.',fontsize=9.3,color=MUTED)
save(fig,'validation')
print('Saved two PNG and SVG figures.')
