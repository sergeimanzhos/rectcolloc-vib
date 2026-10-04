import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ref=np.array([float(x) for x in open('reference.txt').read().split()])
def met(E):
    E=np.array(E); e=np.abs((E[1:]-E[0])-ref[1:100]); return np.median(e[:20]),e[:20].max(),np.median(e),e.max()
STD={2000:'std_N2000',5000:'std_N5000',10000:'std_N10000'}
GRD={2000:'p3_N2000',5000:'p3_long',10000:'p3_N10000'}
def series(tag,mmin=0):
    d=json.load(open(f'{tag}.json')); N=d['N']; R=[r for r in d['recs'] if r['M']/N>=mmin]
    return np.array([r['M']/N for r in R]), np.array([met(r['E']) for r in R])
plt.rcParams.update({'font.size':20,'axes.titlesize':20})
blue,orange='#2a78d6','#eb6834'
# ---- Figure 1: error vs N (largest M of each run) ----
Ns=np.array([2000,5000,10000]); fig,axs=plt.subplots(1,2,figsize=(14,6.9),sharey=True)
for col,(lab,(im,ix)) in enumerate([('levels 1–20',(0,1)),('levels 1–99',(2,3))]):
    ax=axs[col]
    for name,tags,c,mk in (('energies only',STD,blue,'o'),('energies + gradients',GRD,orange,'s')):
        v=np.array([series(tags[N])[1][-1] for N in Ns])
        ax.plot(Ns,v[:,im],'-',marker=mk,color=c,lw=2.8,ms=11,label=f'{name}, median')
        ax.plot(Ns,v[:,ix],'--',marker=mk,color=c,lw=2.2,ms=11,mfc='white',mew=2,label=f'{name}, max')
    for name,tags,ls in ((('energies only',STD,'-'),('energies + gradients',GRD,':')) if col==0 else ()):
        z=[abs(json.load(open(f'{tags[N]}.json'))['recs'][-1]['E'][0]-ref[0]) for N in Ns]
        ax.plot(Ns,z,ls,marker='o',color='0.72',mfc='0.72',lw=2.2,ms=10,label=f'{name}, |ZPE error|')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xticks(Ns); ax.set_xticklabels(['2000','5000','10000']); ax.minorticks_off()
    ax.set_xlabel('N'); ax.set_title(lab); ax.grid(alpha=0.3,which='both')
from matplotlib.ticker import FixedLocator, FixedFormatter
for a in axs:
    a.yaxis.set_major_locator(FixedLocator([0.05,0.1,0.2,0.5,1,2,5,10])); a.yaxis.set_major_formatter(FixedFormatter(['0.05','0.1','0.2','0.5','1','2','5','10'])); a.yaxis.set_minor_locator(matplotlib.ticker.NullLocator()); a.tick_params(labelleft=True)
axs[0].set_ylabel('error vs reference (cm$^{-1}$)')
h,l=axs[0].get_legend_handles_labels(); order=[0,1,4,2,3,5]; fig.legend([h[i] for i in order],[l[i] for i in order],loc='lower center',ncol=2,fontsize=16,frameon=False)
fig.tight_layout(rect=(0,0.24,1,1)); fig.savefig('cl2co_error_vs_N.png',dpi=200)
# ---- Figure 2: median error vs M/N ----
cols={2000:'#2a78d6',5000:'#eb6834',10000:'#1baf7a'}
fig,axs=plt.subplots(1,2,figsize=(14,6.6),sharey=True)
for col,(lab,im) in enumerate([('levels 1–20',0),('levels 1–99',2)]):
    ax=axs[col]
    for N in Ns:
        x,v=series(STD[N],mmin=2); ax.plot(x,v[:,im],'-o',color=cols[N],lw=2.8,ms=10,label=f'N = {N}, energies only')
        x,v=series(GRD[N]); ax.plot(x,v[:,im],'--s',color=cols[N],lw=2.5,ms=10,mfc='white',mew=2,label=f'N = {N}, energies + gradients')
    ax.set_yscale('log'); ax.set_xlabel('M / N'); ax.set_title(lab); ax.grid(alpha=0.3,which='both'); ax.set_xlim(0,6.3)
for a in axs:
    a.yaxis.set_major_locator(FixedLocator([0.05,0.1,0.2,0.5,1,2,5])); a.yaxis.set_major_formatter(FixedFormatter(['0.05','0.1','0.2','0.5','1','2','5'])); a.yaxis.set_minor_locator(matplotlib.ticker.NullLocator()); a.tick_params(labelleft=True)
axs[0].set_ylabel('median error vs reference (cm$^{-1}$)')
h,l=axs[0].get_legend_handles_labels(); order=[0,2,4,1,3,5]
fig.legend([h[i] for i in order],[l[i] for i in order],loc='lower center',ncol=2,fontsize=16,frameon=False)
fig.tight_layout(rect=(0,0.22,1,1)); fig.savefig('cl2co_error_vs_MN.png',dpi=200)
