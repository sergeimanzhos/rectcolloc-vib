import numpy as np, json, sys, os, time
from collo6 import *
from r2 import r2
import coords6 as c
H=H2CM; ref=np.array([float(x) for x in open('reference.txt').read().split()])
BASEW=np.array([0.857375,0.857375,0.95,0.9025,0.95,0.95])
def box(E=6500/H):
    from scipy.optimize import brentq
    lo,hi=np.zeros(NQ),np.zeros(NQ)
    for k in range(NQ):
        for s in (-1,1):
            t=np.arange(0.05,25,0.05); v=V(np.outer(s*t,np.eye(NQ)[k]))-E; i=np.argmax(v>0)
            r=brentq(lambda tt: V(np.eye(NQ)[k][None]*s*tt)[0]-E,t[i-1],t[i]); (hi if s>0 else lo)[k]=s*r
    return lo,hi
if not os.path.exists('box.npz'):
    lo,hi=box(); np.savez('box.npz',lo=lo,hi=hi)
b=np.load('box.npz'); lo,hi=b['lo'],b['hi']
def points(start,n,vmax,grad=False):
    fn=f'pts_{start}_{vmax:.0f}_{n}.npz'
    if os.path.exists(fn): d=np.load(fn); return d['q'],d['V'],d['g']
    Q=[];Vs=[];i=start
    while sum(len(x) for x in Q)<n:
        q=lo+(hi-lo)*r2(np.arange(i,i+200000),NQ); i+=200000; v=V(q); w=v<=vmax/H
        assert v.min()>-1e-6, 'hole in PES'
        Q.append(q[w]); Vs.append(v[w])
    q=np.concatenate(Q)[:n]; v=np.concatenate(Vs)[:n]; g=gradV(q)
    np.savez(fn,q=q,V=v,g=g); return q,v,g
VG=None
def run(N,vmax,ec,Ms,grad,tag,lam=1.0,fac=1.0,fadep=None):
    out={'N':N,'vmax':vmax,'ec':ec,'grad':grad,'fadep':fadep,'recs':[]}
    qf,Vf,gf=points(1,max(Ms),vmax); qv,Vv,_=points(10_000_000,2000,vmax)
    if fadep is not None: lam=np.clip(1-Vf*H/vmax,0,1)**fadep
    cq,_,_=points(20_000_000,N,ec); bas=Basis(cq,BASEW*fac)
    S=H_=None; done=0; t0=time.time()
    for M in Ms:
        lamv=lam*np.ones(M-done) if np.isscalar(lam) else lam[done:M]
        if VG is not None: lamv=lamv*(Vf[done:M]<=VG/H)
        S,H_=accumulate(bas,qf[done:M],Vf[done:M],gf[done:M] if grad else None,lam=lamv,S=S,H=H_); done=M
        ta=time.time()-t0
        E,Cv=solve(S,H_,100,2250/H); r=residual(bas,qv,Vv,E,Cv)
        Er=E.real*H; tr=Er[1:]-Er[0]; d=np.abs(tr-ref[1:100])
        rec=dict(M=M,eqN=M*(1+NQ*grad)/N,zpe=Er[0],zpe_err=Er[0]-ref[0],med20=np.median(d[:19]),max20=d[:19].max(),
                 med99=np.median(d),max99=d.max(),res=float(np.median(r)),maxim=float(np.abs(E.imag).max()*H),t=time.time()-t0,E=Er.tolist())
        out['recs'].append(rec); json.dump(out,open(f'{tag}.json','w'))
        print(f"{tag} M={M} eq/N={rec['eqN']:.2f} ZPE err {rec['zpe_err']:+.3f} | 1-19 {rec['med20']:.3f}/{rec['max20']:.3f} | 1-99 {rec['med99']:.3f}/{rec['max99']:.3f} | res {rec['res']:.1f} | maxIm {rec['maxim']:.2f} | {rec['t']:.0f}s",flush=True)
        if M==Ms[-1]: np.save(f'{tag}_C.npy',Cv); np.save(f'{tag}_E.npy',E)
    return out
if __name__=='__main__':
    N=int(sys.argv[1]); vmax=float(sys.argv[2]); ec=1150; mode=sys.argv[3]
    if mode=='std': run(N,vmax,ec,[int(N*f) for f in (1,2,3,4,6)] if N<10000 else [N*f for f in (1,2,3,4)],False,f'std_N{N}')
    elif mode=='grad': run(N,vmax,ec,[int(N*f) for f in (0.25,0.5,0.75,1.0,1.5)] if N<10000 else [int(N*f) for f in (0.25,0.5,0.75,1.0)],True,f'grad_N{N}')
    elif mode=='lam':   # python3 run.py N vmax lam LAM [VG]
        lam=float(sys.argv[4]); VG=float(sys.argv[5]) if len(sys.argv)>5 else None
        run(N,vmax,ec,[int(N*f) for f in (0.5,1.0,1.5)],True,f'grad_N{N}_lam{lam}'+(f'_vg{VG:.0f}' if VG else ''),lam=lam)
    elif mode=='cfg':   # python3 run.py N vmax cfg p fac ec "M1,M2,..." tag
        p=float(sys.argv[4]); fac=float(sys.argv[5]); ec=float(sys.argv[6]); Ms=[int(x) for x in sys.argv[7].split(',')]
        run(N,vmax,ec,Ms,p>=0,sys.argv[8],fadep=p if p>=0 else None,fac=fac)
    elif mode=='fade':  # python3 run.py N vmax fade p   ->  lambda_i = (1 - V_i/Vmax)^p
        p=float(sys.argv[4]); Ms=[int(N*f) for f in ((0.5,1,1.5,2,3) if N<=2000 else (0.5,1.0,1.5))]
        run(N,vmax,ec,Ms,True,f'fade_N{N}_p{p:g}',fadep=p)
