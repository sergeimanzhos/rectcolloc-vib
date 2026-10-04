"""6D rectangular collocation for Cl2CO with optional potential-gradient rows.
Exact Cartesian 5-pt FD KEO through Eckart coords; generalized shift-invert eigensolver."""
import numpy as np, numexpr as ne, time, os
import scipy.linalg as sla
from scipy.sparse.linalg import eigs, LinearOperator
import coords6 as c
from pes import NNPES, rij
H2CM=c.H2CM; NQ=c.NQ; NC=3*c.NA
_P=NNPES(os.environ.get('PESFILE',os.path.join(os.path.dirname(__file__),'fitNN100NNs4.dat')))
def V(q):  # hartree, relative to PES minimum
    return (_P(rij(c.X_of_q(q)))-c.Vmin)/H2CM
def gradV(q,d=1e-3):  # dV/dq_k by 5-pt central FD in q (X(q) linear => exact direction)
    g=np.zeros((len(q),NQ))
    for k in range(NQ):
        f=[V(q+s*d*np.eye(NQ)[k]) for s in (-2,-1,1,2)]
        g[:,k]=(f[0]-8*f[1]+8*f[2]-f[3])/(12*d)
    return g
FDW=np.array([-1,16,-30,16,-1])/12.0; FDS=np.array([-2,-1,0,1,2]); D1W=np.array([1,-8,0,8,-1])/12.0
class Basis:
    def __init__(s,C,beta):
        s.beta=np.asarray(beta,float); s.Cs=np.asarray(C)/s.beta; s.c2=(s.Cs**2).sum(1)
    def __call__(s,q):
        qs=q/s.beta; G=qs@s.Cs.T; q2=(qs**2).sum(1)[:,None]; c2=s.c2[None,:]
        return ne.evaluate('exp(-0.5*(q2+c2-2*G))')
    def grad_k(s,q,k,B):
        qs=q/s.beta; return -(qs[:,k:k+1]-s.Cs[None,:,k])/s.beta[k]*B
def T_basis(bas,q,h=1e-3):
    n=len(q); X=c.X_of_q(q).reshape(n,NC)
    out=FDW[2]*(1/c.mm).sum()*bas(q)
    for i in range(NC):
        for sft,w in zip(FDS,FDW):
            if sft==0: continue
            Xd=X.copy(); Xd[:,i]+=sft*h; out+=(w/c.mm[i])*bas(c.q_of_X(Xd))
    return -0.5*out/h**2
def accumulate(bas,q,Vq,gq=None,lam=1.0,batch=None,S=None,H=None):
    N=len(bas.c2)
    if batch is None: batch=max(20,int(2e7/N/ (8 if gq is not None else 2)))
    if S is None: S=np.zeros((N,N)); H=np.zeros((N,N))
    for a in range(0,len(q),batch):
        qb=q[a:a+batch]; Vb=Vq[a:a+batch]
        B=bas(qb); A=T_basis(bas,qb)+Vb[:,None]*B
        S+=B.T@B; H+=B.T@A
        if gq is not None:
            gb=gq[a:a+batch]; lw=np.sqrt(np.broadcast_to(np.asarray(lam,float),(len(q),))[a:a+batch])[:,None]
            for k in range(NQ):
                Tchi=np.zeros_like(B)
                for sft,w in zip(FDS,D1W):
                    if sft==0: continue
                    qs=qb.copy(); qs[:,k]+=sft*1e-2; Tchi+=w*T_basis(bas,qs)
                Tchi/=1e-2
                chi=bas.grad_k(qb,k,B)*lw
                S+=chi.T@chi; H+=chi.T@((Tchi+Vb[:,None]*bas.grad_k(qb,k,B)+gb[:,k:k+1]*B)*lw)
                del Tchi
    return S,H
def solve(S,H,nev,sigma):
    """eigenpairs of H c = E S c nearest sigma (shift-invert; no inversion of S needed)."""
    lu=sla.lu_factor(H-sigma*S,overwrite_a=False,check_finite=False)
    op=LinearOperator(S.shape,matvec=lambda x: sla.lu_solve(lu,S@x,check_finite=False),dtype=float)
    nu,vec=eigs(op,k=nev,which='LM',tol=1e-12,maxiter=20000)
    E=sigma+1/nu; o=np.argsort(E.real); return E[o],vec[:,o]
def residual(bas,qv,Vv,E,C):
    B=bas(qv); A=T_basis(bas,qv)+Vv[:,None]*B
    Psi=B@C; R=A@C-Psi*E[None,:]
    return np.linalg.norm(R,axis=0)/np.linalg.norm(Psi,axis=0)*H2CM
