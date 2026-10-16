import os
import numpy as np
from pyscf import gto, dft, lib
from pyscf.tools import cubegen
lib.num_threads(1); BAS='aug-cc-pvdz'; OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','cubes')+'/'
import os; os.makedirs(OUT, exist_ok=True)
def run(tag, bent=False):
    d=np.load(f'ff_{tag}.npz'); el=list(d['el']); X=d['X']
    if bent:
        mol=gto.M(atom=[[e,x] for e,x in zip(el[:3],X[:3])]+[['ghost-Li',X[3]]],
                  basis={'C':BAS,'O':BAS,'ghost-Li':gto.basis.load(BAS,'Li')},unit='Bohr',verbose=0); mf=dft.RKS(mol)
    else:
        spin=1 if tag=='lico2' else 0
        mol=gto.M(atom=[[e,x] for e,x in zip(el,X)],basis=BAS,unit='Bohr',spin=spin,verbose=0)
        mf=(dft.UKS if spin else dft.RKS)(mol)
    mf.xc='b3lyp'; mf.disp='d3bj'; mf.grids.level=5; mf.conv_tol=1e-11; mf.kernel()
    return mol,mf
kw=dict(resolution=0.15, margin=6.0)
mol,mf=run('co2'); n=mol.nelectron//2
cubegen.orbital(mol,OUT+'CO2_LUMO_a.cube',mf.mo_coeff[:,n],**kw)
cubegen.orbital(mol,OUT+'CO2_LUMO_b.cube',mf.mo_coeff[:,n+1],**kw)
cubegen.orbital(mol,OUT+'CO2_HOMO_a.cube',mf.mo_coeff[:,n-1],**kw)
print('CO2 LUMO/LUMO+1 eV:',mf.mo_energy[n:n+2]*27.2114, flush=True)
mol,mf=run('lico2'); na=mol.nelec[0]
cubegen.orbital(mol,OUT+'LiCO2_SOMO_alpha.cube',mf.mo_coeff[0][:,na-1],**kw)
dm=mf.make_rdm1(); cubegen.density(mol,OUT+'LiCO2_spin_density.cube',dm[0]-dm[1],**kw)
print('LiCO2 SOMO eV:',mf.mo_energy[0][na-1]*27.2114, flush=True)
mol,mf=run('lico2',bent=True); n=mol.nelectron//2
cubegen.orbital(mol,OUT+'bentCO2_LUMO_ghostLi.cube',mf.mo_coeff[:,n],**kw)
print('bent CO2 LUMO eV:',mf.mo_energy[n]*27.2114, flush=True)
