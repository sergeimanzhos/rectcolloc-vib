"""Linear-dipole IR intensities from Cl2CO collocation wavefunctions (rho-weighted QMC on dense R2 grid)."""
import numpy as np, sys, json, itertools
from collo6 import *
from r2 import r2
import coords6 as c
from scipy.special import eval_hermite, factorial
from scipy import constants as K
H=H2CM; ea0=K.e*K.physical_constants['Bohr radius'][0]
def A_kmmol(nu,mu2): return 8*np.pi**3*K.N_A*nu*100*mu2*ea0**2/(3*K.h*K.c*4*np.pi*K.epsilon_0)/1000
tag=sys.argv[1]; N=int(sys.argv[2]); vmax=float(sys.argv[3]); NPTS=int(sys.argv[4]) if len(sys.argv)>4 else 400000
BASEW=np.array([0.857375,0.857375,0.95,0.9025,0.95,0.95])
cq=np.load(f'pts_20000000_1150_{N}.npz')['q']; bas=Basis(cq,BASEW)
E=np.load(f'{tag}_E.npy'); C=np.load(f'{tag}_C.npy')
Cr=np.zeros(C.shape); pairs=[]; used=set(); Ei=E.imag*H
for f in range(C.shape[1]):
    if f in used: continue
    cand=[j for j in range(C.shape[1]) if j!=f and j not in used and abs(E[j]-np.conj(E[f]))*H<1e-4]
    if abs(Ei[f])>1e-6 and cand:
        g=cand[0]; Cr[:,f]=C[:,f].real; Cr[:,g]=C[:,f].imag; used|={f,g}; pairs.append((f,g))
    else:
        ph=np.exp(-0.5j*np.angle((C[:,f]**2).sum())); Cr[:,f]=(C[:,f]*ph).real; used.add(f)
print('complex pairs (state indices):',pairs, ' max|Im E| %.2f cm-1'%np.abs(Ei).max())
import os
DIP=os.environ.get('DIP','dipder.npz'); dd=np.load(DIP)
POLY='mono' in dd.files                                       # polynomial dipole fit: mono (n,6) exponents, coef (n,3)
if POLY: MONO=dd['mono']; COEF=dd['coef']; dmu=np.array([COEF[[i for i,a in enumerate(MONO) if a.sum()==1 and a[k]==1][0]] for k in range(NQ)])
else: dmu=dd['dmu']                                           # (mode, xyz)
def dipole(q):
    if not POLY: return q@dmu
    F=np.prod(q[:,None,:]**MONO[None,:,:],axis=2); return F@COEF
b=np.load('box.npz'); lo,hi=b['lo'],b['hi']
# HO product functions for assignment: harmonic energy <= 2700 cm-1 above ZPE
qn=[t for t in itertools.product(*[range(int(2700/(w*H))+1) for w in c.omega]) if np.dot(t,c.omega)*H<=2700]
def ho1(n,x): return eval_hermite(n,x)*np.exp(-x**2/2)/np.sqrt(2.0**n*factorial(n)*np.sqrt(np.pi))
nf=C.shape[1]; Sq=np.zeros(nf); Sn=np.zeros(nf); S0f=np.zeros(nf); M0f=np.zeros((nf,3)); Mfull=np.zeros((3,nf,nf)); Ov=np.zeros((nf,nf)); HOov=np.zeros((nf,len(qn))); HOn=np.zeros(len(qn))
cnt=0; i=30_000_000
while cnt<NPTS:
    q=lo+(hi-lo)*r2(np.arange(i,i+100000),NQ); i+=100000; q=q[V(q)<=vmax/H]; cnt+=len(q)
    w=c.rho(q); P=bas(q)@Cr; mu=dipole(q)                       # dipole (constant removed later)
    Sq+=(P**2).sum(0); Sn+=(w[:,None]*P**2).sum(0); Ov+=P.T@(w[:,None]*P); Mfull+=np.einsum('nc,na,nb->cab',mu*w[:,None],P,P)
    hp=np.ones((len(q),len(qn)))
    for k in range(NQ):
        hk=np.array([ho1(n,q[:,k]) for n in range(max(t[k] for t in qn)+1)]).T
        hp*=hk[:,[t[k] for t in qn]]
    HOov+=P.T@hp; HOn+=(hp**2).sum(0)
nrm=np.sqrt(Sn); Ov/=np.outer(nrm,nrm); Mfull/=np.outer(nrm,nrm)[None]; HOov/=np.sqrt(Sq)[:,None]
O=Ov.copy(); np.fill_diagonal(O,0); print('before Loewdin: max|overlap| %.2e median %.2e'%(np.abs(O).max(),np.median(np.abs(O[np.triu_indices(nf,1)]))))
ev,U=np.linalg.eigh(Ov); T=U@np.diag(ev**-0.5)@U.T            # symmetric (Loewdin) orthonormalisation, rho metric
Mfull=np.einsum('ag,cab,bh->cgh',T,Mfull,T); HOov=T.T@HOov; Ov=T.T@Ov@T
M0f=(Mfull[:,0,:]-Mfull[:,0,0][:,None]*Ov[0][None,:]).T        # <0|mu - <0|mu|0>|f>
mu00=M0f[0].copy(); M0f-=mu00[None]*Ov[0][:,None]           # <0|mu-<0|mu|0>|f>
Pn=(P**2).sum(0)
Et=(E.real-E.real[0])*H; O=Ov.copy(); np.fill_diagonal(O,0)
print('points',cnt,' after Loewdin max|overlap| %.1e'%np.abs(O).max())
# normalise weights per state over dq: use unweighted norms
Sq=np.zeros(nf)
rows=[]
for f in range(1,nf):
    wf=HOov[f]**2/HOn                                       # weight of each normalised HO product function (dq metric)
    j=np.argmax(wf); lab=''.join(str(x) for x in qn[j])
    A=A_kmmol(Et[f],(M0f[f]**2).sum()); sym={0:'A1'}  # symmetry from dipole component
    comp='xyz'[int(np.argmax(np.abs(M0f[f])))] if A>1e-6 else '-'
    rows.append(dict(pair=any(f in p for p in pairs),Tdiag=float(T[f,f]),n=f,E=Et[f],label=lab,w=float(wf[j]),A=A,comp=comp,ov0=float(Ov[0,f])))
# double-harmonic reference
harm=[dict(mode=k,omega=c.omega[k]*H,A=A_kmmol(c.omega[k]*H,(dmu[k]**2).sum()/2)) for k in range(NQ)]
json.dump(dict(rows=rows,harm=harm,ZPE=E.real[0]*H),open(f'ir_{tag}_'+os.path.basename(DIP).replace('.npz','')+'.json','w'),indent=1)
print('double harmonic:',[(round(h['omega'],1),round(h['A'],2)) for h in harm])
for r in rows:
    if r['A']>0.01: print(f"{r['n']:3d} {r['E']:8.2f} {r['label']} w={r['w']:.2f} {r['comp']} A={r['A']:.3f} <0|f>={r['ov0']:.1e}")
