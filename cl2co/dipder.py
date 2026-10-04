import numpy as np, coords6 as c, time, sys
BAS=sys.argv[1] if len(sys.argv)>1 else 'aug-cc-pvdz'; OUT=sys.argv[2] if len(sys.argv)>2 else 'dipder.npz'
from pyscf import gto, dft, lib
lib.num_threads(2)
def mu(X,dm0=None):
    mol=gto.M(atom=[(s,tuple(x)) for s,x in zip(['C','O','Cl','Cl'],X)],basis=BAS,unit='B',symmetry=False,verbose=0)
    mf=dft.RKS(mol); mf.xc='b3lyp'; mf.grids.level=3; mf.conv_tol=1e-10; mf.kernel(dm0=dm0)
    assert mf.converged; return mf.dip_moment(unit='AU',verbose=0), mf.make_rdm1()
t=time.time(); mu0,dm=mu(c.X0); print('mu0 (au)',mu0, '%.0fs'%(time.time()-t),flush=True)
d=0.02; D=[]
for k in range(c.NQ):
    mp,_=mu(c.X_of_q(d*np.eye(c.NQ)[k])[0],dm); mm,_=mu(c.X_of_q(-d*np.eye(c.NQ)[k])[0],dm)
    D.append((mp-mm)/(2*d)); print(k,np.round(D[-1],5),flush=True)
np.savez(OUT,basis=BAS,mu0=mu0,dmu=np.array(D),delta=d)
