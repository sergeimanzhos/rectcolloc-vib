import numpy as np
from pyscf import gto, dft, lib
lib.num_threads(1); HC=219474.6313705
d=np.load('ff_co2.npz'); X=d['X']; L=d['L']; w=d['omega']; m=d['m']; el=list(d['el'])
s=1/np.sqrt(np.repeat(m*1822.888486209,3))
Xb=X+(s*L[:,0]*1.0).reshape(3,3)          # bent geometry (dQ=1 along bend)
u3=(s*L[:,3]).reshape(3,3)                # Cartesian direction of Q3
LEV=8
def mk(Xd,disp):
    mol=gto.M(atom=[[e,x] for e,x in zip(el,Xd)],unit='Bohr',basis='aug-cc-pvdz',verbose=0)
    mf=dft.RKS(mol).set(xc='b3lyp',conv_tol=1e-12); mf.grids.level=LEV
    if disp: mf.disp='d3bj'
    mf.kernel(); return mf
for lev in (8,):
    disp=True
    mf=mk(Xb,disp); h=mf.Hessian().kernel().transpose(0,2,1,3).reshape(9,9)
    a=(s*L[:,3])@h@(s*L[:,3])
    g=lambda Xd: (lambda f: (f.nuc_grad_method().set(grid_response=True).kernel()))(mk(Xd,disp))
    hh=0.05
    fd=((g(Xb+hh*u3)-g(Xb-hh*u3)).ravel()@(s*L[:,3]))/(2*hh)
    print(f'disp={disp}: analytic Hq33={a:.6e}  FD-of-gradients Hq33={fd:.6e}  diff={(a-fd)/w[3]**2:.2e} (rel. to w3^2)',flush=True)
