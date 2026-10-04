import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':20})
cols={2000:'#2a78d6',5000:'#eb6834',10000:'#1baf7a'}
fig,ax=plt.subplots(figsize=(10,8))
for N in (2000,5000,10000):
    for mode,ls,mk in (('std','-','o'),('grad','--','s')):
        R=json.load(open(f'{mode}_N{N}.json'))['recs']; R=[r for r in R if r['M']>=N*0.5 or mode=='grad']
        if mode=='std': R=[r for r in R if r['M']>=2*N]
        x=[r['M']/N for r in R]; y=[r['med99'] for r in R]
        ax.plot(x,y,ls,marker=mk,color=cols[N],lw=2.5,ms=9,mfc='white' if mode=='grad' else cols[N],mew=2,
                label=f'N = {N}, '+('energies only' if mode=='std' else 'energies + gradients'))
ax.set_yscale('log'); ax.set_xlabel('collocation points / N'); ax.set_ylabel('median error, levels 1–99 (cm$^{-1}$)')
ax.set_xlim(0,6.3)
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
ax.yaxis.set_major_locator(FixedLocator([0.2,0.5,1,2,5])); ax.yaxis.set_major_formatter(FixedFormatter(['0.2','0.5','1','2','5'])); ax.yaxis.set_minor_locator(NullLocator())
ax.grid(alpha=0.3,which='both'); ax.legend(fontsize=15,frameon=False,ncol=2,loc='upper center',bbox_to_anchor=(0.5,-0.14))
fig.tight_layout(); fig.savefig('cl2co_convergence.png',dpi=200)
# spectrum (linear dipole, N=10000 energies-only)
d=json.load(open('ir_std_N10000.json')); rows=[r for r in d['rows'] if r['A']>0.05]
fig,ax=plt.subplots(figsize=(11,6))
E=np.array([r['E'] for r in rows]); A=np.array([r['A'] for r in rows])
ax.vlines(E,0,A,color='#2a78d6',lw=3,label='collocation, PES + linear DFT dipole')
h=d['harm']; ax.plot([x['omega'] for x in h],[x['A'] for x in h],'^',color='#eb6834',ms=11,label='double harmonic')
lab={7:r'$\nu_5$',67:r'$\nu_1$',68:r'$\nu_2\!+\!\nu_5\!+\!\nu_6$',3:r'$\nu_2$',48:r'$2\nu_5$',12:r'$\nu_2\!+\!\nu_6$'}
for r in rows:
    if r['n'] in lab:
        dx,ha=(-15,'right') if r['n'] in (7,67) else ((-25,'right') if r['n']==3 else (0,'center'))
        ax.text(r['E']+dx,r['A']+(12 if r['n'] not in (7,67) else -40),lab[r['n']],ha=ha,fontsize=17)
ax.set_xlabel(r'wavenumber (cm$^{-1}$)'); ax.set_ylabel('intensity (km/mol)'); ax.set_ylim(0,520); ax.set_xlim(200,2150)
ax.legend(fontsize=15,frameon=False,loc='upper center',bbox_to_anchor=(0.62,1.0)); ax.grid(alpha=0.3,axis='y')
fig.tight_layout(); fig.savefig('cl2co_IR_linear_dipole.png',dpi=200)
