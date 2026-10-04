"""Same stencil on the NN PES (validation of vp2.py against collocation on the same surface); dipole = B3LYP/aDZ linear."""
import numpy as np, itertools, coords6 as c
from collo6 import V, gradV
d=0.25; pts=[np.zeros(6)]
for i in range(6):
    for s in (-2,-1,1,2): v=np.zeros(6); v[i]=s*d; pts.append(v)
for i,j in itertools.combinations(range(6),2):
    for si,sj in itertools.product((-1,1),repeat=2): v=np.zeros(6); v[i]=si*d; v[j]=sj*d; pts.append(v)
for ss in itertools.product((-1,1),repeat=4): v=np.zeros(6); v[[5,2,4,1]]=np.array(ss)*d; pts.append(v)
q=np.array(pts); E=V(q); g=gradV(q,d=1e-3)
dd=np.load('dipder.npz'); mu=dd['mu0'][None]+q@dd['dmu']
X0=c.X0-(c.m[:,None]*c.X0).sum(0)/c.m.sum()
np.savez('vp_nn.npz',q=q,E=E,g=g,mu=mu,omega=c.omega,L=c.L,Dq=c.A,Xeq=X0,m=c.m,delta=d,basis='NN')
print('ok',len(q))
