import numpy as np
from pyscf import gto, dft, lib
lib.num_threads(1); HC=219474.6313705
d=np.load('ff_co2.npz'); X=d['X']; L=d['L']; w=d['omega']; m=d['m']; el=list(d['el'])
s=1/np.sqrt(np.repeat(m*1822.888486209,3)); N=3
def Hq(Xd):
    mol=gto.M(atom=[[e,x] for e,x in zip(el,Xd)],unit='Bohr',basis='aug-cc-pvdz',verbose=0)
    mf=dft.RKS(mol).set(xc='b3lyp',conv_tol=1e-12); mf.grids.level=5; mf.disp='d3bj'; mf.kernel()
    h=mf.Hessian().kernel().transpose(0,2,1,3).reshape(9,9)
    return L.T@(h*np.outer(s,s))@L
H0=d['Hq0']
for dQ in (1.0,):
    for k,(i,j) in ((0,(3,3)),(3,(0,0))):
        dd=(s*L[:,k]*dQ).reshape(N,3)
        D2=(Hq(X+dd)+Hq(X-dd)-2*H0)/dQ**2
        print(f'dQ={dQ} displace mode {k}, element ({i},{j}): phi = {D2[i,j]/np.sqrt(w[i]**2*w[k]**2)*HC:.2f} cm-1',flush=True)
for k,(i,j) in ((0,(3,3)),(3,(0,0))):
    D2=(d['Hqp'][k]+d['Hqm'][k]-2*H0)/2.0**2
    print(f'dQ=2.0 displace mode {k}, element ({i},{j}): phi = {D2[i,j]/np.sqrt(w[i]**2*w[k]**2)*HC:.2f} cm-1')
