"""Eckart-frame dimensionless normal coordinates for Cl2CO (C,O,Cl1,Cl2), atomic units."""
import numpy as np
H2CM=219474.6313705
_d=np.load(__import__('os').path.join(__import__('os').path.dirname(__file__),'eq.npz'))
X0=_d['X0']; m=_d['m']; Vmin=float(_d['Vmin']); NA=len(m)
mm=np.repeat(m,3)
def _vib():
    Hm=_d['Hess']/np.sqrt(np.outer(mm,mm)); sq=np.sqrt(mm); V=[]
    for a in range(3):
        t=np.zeros((NA,3)); t[:,a]=1; V.append(t.ravel()*sq)
    for a in range(3):
        e=np.zeros(3); e[a]=1; V.append(np.cross(e,X0).ravel()*sq)
    Qv,_=np.linalg.qr(np.array(V).T); P=np.eye(3*NA)-Qv@Qv.T
    w2,L=np.linalg.eigh(P@Hm@P); idx=np.argsort(w2)[-(3*NA-6):]
    L=L[:,idx]; om=np.sqrt(w2[idx])
    for k in range(L.shape[1]):
        if L[np.argmax(np.abs(L[:,k])),k]<0: L[:,k]*=-1
    return om,L
omega,L=_vib(); NQ=len(omega)
Q=np.sqrt(omega)[:,None]*(L.T*np.sqrt(mm)[None,:])
A=(L/np.sqrt(mm)[:,None])/np.sqrt(omega)[None,:]
def X_of_q(q):
    q=np.atleast_2d(q); return X0[None]+(q@A.T).reshape(-1,NA,3)
def eckart_rotate(X):
    com=(m[None,:,None]*X).sum(1)/m.sum(); x=X-com[:,None,:]
    Hc=np.einsum('i,nia,ib->nab',m,x,X0); U,S,Vt=np.linalg.svd(Hc)
    d=np.sign(np.linalg.det(np.einsum('nba,ncb->nac',Vt,U)))
    D=np.ones((len(X),3)); D[:,2]=d
    R=np.einsum('nba,nb,ncb->nac',Vt,D,U)
    return np.einsum('nab,nib->nia',R,x)
def q_of_X(X,eckart=True):
    X=np.asarray(X).reshape(-1,NA,3); x=eckart_rotate(X) if eckart else X
    return (x-X0[None]).reshape(len(X),3*NA)@Q.T
def rho(q):
    q=np.atleast_2d(q); X=X_of_q(q); X=X-(m[None,:,None]*X).sum(1)[:,None]/m.sum()
    n=len(q); J=np.zeros((n,3*NA,3*NA))
    for a in range(3):
        t=np.zeros(3*NA); t[a::3]=1; J[:,:,a]=t
        e=np.zeros(3); e[a]=1; J[:,:,3+a]=np.cross(e[None,None,:],X).reshape(n,3*NA)
    J[:,:,6:]=A[None]
    return np.abs(np.linalg.det(J))
