import os
import numpy as np
from pyscf import gto, dft, lib, symm
from pyscf.tools import cubegen
lib.num_threads(1); BAS='aug-cc-pvdz'; EV=27.211386; OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','cubes')+'/'
def val_frac(mol, C):
    pmol=gto.M(atom=mol._atom, basis='minao', unit='Bohr', spin=mol.spin, verbose=0)
    pmol=pmol if not any('ghost' in a[0].lower() for a in mol._atom) else gto.M(
        atom=[a for a in mol._atom if 'ghost' not in a[0].lower()], basis='minao', unit='Bohr', verbose=0)
    S12=gto.intor_cross('int1e_ovlp', mol, pmol); S22=pmol.intor('int1e_ovlp')
    P=S12@np.linalg.solve(S22, S12.T)
    return np.einsum('pi,pq,qi->i', C, P, C)
d=np.load('ff_co2.npz')
mol=gto.M(atom=[[e,x] for e,x in zip(d['el'],d['X'])],basis=BAS,unit='Bohr',symmetry=True,verbose=0)
mf=dft.RKS(mol).set(xc='b3lyp',conv_tol=1e-11); mf.disp='d3bj'; mf.grids.level=5; mf.kernel()
lab=symm.label_orb_symm(mol,mol.irrep_name,mol.symm_orb,mf.mo_coeff); n=mol.nelectron//2
vf=val_frac(mol, mf.mo_coeff)
print('CO2 virtuals: idx, eV, irrep, valence fraction')
for i in range(n,n+12): print(i-n, f'{mf.mo_energy[i]*EV:7.3f}', lab[i], f'{vf[i]:.3f}')
pi=[i for i in range(n,n+30) if lab[i].startswith('E1u')]
best=max(pi, key=lambda i: vf[i]); pair=[i for i in pi if abs(mf.mo_energy[i]-mf.mo_energy[best])<1e-5]
print('valence pi* pair:', [p-n for p in pair], f'{mf.mo_energy[best]*EV:.3f} eV, valence frac {vf[best]:.3f}')
np.save('co2_pistar_idx.npy', np.array(pair))
import os
for f in ('CO2_LUMO_a.cube','CO2_LUMO_b.cube'):
    if os.path.exists(OUT+f): os.remove(OUT+f)
kw=dict(resolution=0.15, margin=6.0)
mol.symmetry=False
cubegen.orbital(mol,OUT+'CO2_pistar_x.cube',mf.mo_coeff[:,pair[0]],**kw)
cubegen.orbital(mol,OUT+'CO2_pistar_y.cube',mf.mo_coeff[:,pair[1]],**kw)
cubegen.orbital(mol,OUT+'CO2_KS_LUMO_diffuse_sigma.cube',mf.mo_coeff[:,n],**kw)
# LiCO2 SOMO and bent-CO2 LUMO valence character
dl=np.load('ff_lico2.npz')
mL=gto.M(atom=[[e,x] for e,x in zip(dl['el'],dl['X'])],basis=BAS,unit='Bohr',spin=1,verbose=0)
fL=dft.UKS(mL).set(xc='b3lyp',conv_tol=1e-11); fL.disp='d3bj'; fL.grids.level=5; fL.kernel()
na=mL.nelec[0]; print(f'LiCO2 SOMO {fL.mo_energy[0][na-1]*EV:.3f} eV, valence frac {val_frac(mL, fL.mo_coeff[0][:,na-1:na])[0]:.3f}')
el=list(dl['el']); X=dl['X']
mB=gto.M(atom=[[e,x] for e,x in zip(el[:3],X[:3])]+[['ghost-Li',X[3]]],
         basis={'C':BAS,'O':BAS,'ghost-Li':gto.basis.load(BAS,'Li')},unit='Bohr',verbose=0)
fB=dft.RKS(mB).set(xc='b3lyp',conv_tol=1e-11); fB.disp='d3bj'; fB.grids.level=5; fB.kernel(); nB=mB.nelectron//2
vB=val_frac(mB, fB.mo_coeff[:,nB:nB+4])
for k in range(4): print(f'bent CO2 virtual {k}: {fB.mo_energy[nB+k]*EV:7.3f} eV valence frac {vB[k]:.3f}')
S=mL.intor('int1e_ovlp'); somo=fL.mo_coeff[0][:,na-1]
for k in range(4): print(f'  |<SOMO|bent virt {k}>|^2 = {(somo@S@fB.mo_coeff[:,nB+k])**2:.3f}')
