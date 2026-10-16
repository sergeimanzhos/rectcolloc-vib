"""In-plane contour maps of the LiCO2 HOMO (alpha SOMO) and of the valence LUMO (pi*) of free CO2, B3LYP-D3(BJ)/aug-cc-pVTZ.
LiCO2 at its aTZ equilibrium geometry; CO2 linear at its own aTZ equilibrium bond length.  usage: python orb_inplane.py out.png"""
import sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from pyscf import gto, dft, lib, symm
import coords as c
out = sys.argv[1]; BAS = 'aug-cc-pvtz'; EV = 27.211386; B = 0.529177210903; lib.num_threads(16)
def ks(mol):
    mf = (dft.UKS if mol.spin else dft.RKS)(mol).set(xc='b3lyp', conv_tol=1e-10, max_cycle=200); mf.grids.level = 5; mf.disp = 'd3bj'; mf.kernel(); assert mf.converged; return mf
def val_frac(mol, C):
    pm = gto.M(atom=mol._atom, basis='minao', unit='Bohr', spin=mol.spin, verbose=0)
    S12 = gto.intor_cross('int1e_ovlp', mol, pm); P = S12 @ np.linalg.solve(pm.intor('int1e_ovlp'), S12.T); return np.einsum('pi,pq,qi->i', C, P, C)
g = np.linspace(-4.2, 4.2, 281) / B; GX, GY = np.meshgrid(g, g, indexing='xy')          # plotting plane, bohr (4.2 A half-width)
def on_plane(mol, C, o, e1, e2):
    pts = o[None] + GX.ravel()[:, None] * e1[None] + GY.ravel()[:, None] * e2[None]; return (dft.numint.eval_ao(mol, pts) @ C).reshape(GX.shape + C.shape[1:])
# ---------------- LiCO2
X = c.X_of_q(np.load('cache_atz/qmin.npy')[None])[0]
mL = gto.M(atom=[[e, x] for e, x in zip(c.EL, X)], basis=BAS, unit='Bohr', spin=1, symmetry=True, verbose=0); fL = ks(mL); na = mL.nelec[0]
labL = symm.label_orb_symm(mL, mL.irrep_name, mL.symm_orb, fL.mo_coeff[0]); R = mL.atom_coords()
nrm = np.linalg.svd(R - R.mean(0))[2][2]; down = R[3] - R[0]; down /= np.linalg.norm(down); e1 = np.cross(down, nrm); oL = 0.5 * (R[0] + R[3])
somo = fL.mo_coeff[0][:, na - 1]; vL = on_plane(mL, somo, oL, e1, -down); AL = np.array([[(r - oL) @ e1, (r - oL) @ -down] for r in R]) * B
print('LiCO2: E %.8f, <S2> %.4f, HOMO (alpha SOMO) %s %.3f eV, valence fraction %.2f' % (fL.e_tot, fL.spin_square()[0], labL[na - 1], fL.mo_energy[0][na - 1] * EV, val_frac(mL, somo[:, None])[0]))
# ---------------- free CO2: optimise the bond length, then find the valence pi* pair
def co2(r): return gto.M(atom=[['C', (0, 0, 0)], ['O', (0, 0, r)], ['O', (0, 0, -r)]], basis=BAS, unit='Bohr', symmetry=True, verbose=0)
res = minimize_scalar(lambda r: ks(co2(r)).e_tot, bracket=(2.15, 2.20, 2.25), tol=1e-4); r0 = res.x
mC = co2(r0); fC = ks(mC); n = mC.nelectron // 2; labC = symm.label_orb_symm(mC, mC.irrep_name, mC.symm_orb, fC.mo_coeff); vf = val_frac(mC, fC.mo_coeff)
print('CO2: r(C-O) = %.4f A, E %.8f; KS LUMO: %s %.3f eV, valence fraction %.2f' % (r0 * B, fC.e_tot, labC[n], fC.mo_energy[n] * EV, vf[n]))
pi = [i for i in range(n, n + 60) if labC[i].startswith('E1u')]; best = max(pi, key=lambda i: vf[i]); pair = [i for i in pi if abs(fC.mo_energy[i] - fC.mo_energy[best]) < 1e-5]
print('CO2 valence pi* pair: virtual numbers %s, %.3f eV, valence fraction %.2f' % ([p - n + 1 for p in pair], fC.mo_energy[best] * EV, vf[best]))
print('CO2 E1u virtuals (eV, valence fraction):', [(round(float(fC.mo_energy[i] * EV), 2), round(float(vf[i]), 2)) for i in pi[:8]])
RC = mC.atom_coords(); ez = (RC[1] - RC[2]) / np.linalg.norm(RC[1] - RC[2]); ex = np.array([1.0, 0, 0]) if abs(ez[0]) < 0.9 else np.array([0, 1.0, 0]); ex -= ez * (ex @ ez); ex /= np.linalg.norm(ex)
V = on_plane(mC, fC.mo_coeff[:, pair], RC[0], ez, ex); M = np.einsum('xyi,xyj->ij', V, V); w, U = np.linalg.eigh(M); vC = V @ U[:, -1]     # the component lying in the plotting plane
AC = np.array([[(r - RC[0]) @ ez, (r - RC[0]) @ ex] for r in RC]) * B
# same phase convention in both panels: the lobe on carbon that points up (+y) is negative (orange)
def above_C(v, A): return v[np.abs(g * B - (A[0][1] + 0.9)).argmin(), np.abs(g * B - A[0][0]).argmin()]
vL = -vL * np.sign(above_C(vL, AL)); vC = -vC * np.sign(above_C(vC, AC))
# ---------------- figure
lev = np.array([0.02, 0.04, 0.08, 0.16, 0.32]); plt.rcParams.update({'font.size': 20}); fig, axs = plt.subplots(1, 2, figsize=(11.6, 6.6), sharey=True)
col = {'C': '#444444', 'O': '#d62728', 'Li': '#9467bd'}
def panel(ax, v, A, el, bonds, title):
    ax.contourf(GX * B, GY * B, v, levels=list(lev) + [10], colors=['#c6dbef', '#9ecae1', '#6baed6', '#3182bd', '#08519c']); ax.contourf(GX * B, GY * B, -v, levels=list(lev) + [10], colors=['#fdd0a2', '#fdae6b', '#fd8d3c', '#e6550d', '#a63603'])
    ax.contour(GX * B, GY * B, v, levels=lev, colors='#08306b', linewidths=0.8); ax.contour(GX * B, GY * B, v, levels=-lev[::-1], colors='#7f2704', linewidths=0.8, linestyles='-')
    for k, l, st in bonds: ax.plot(*zip(A[k], A[l]), color='k', lw=3 if st == '-' else 2, ls=st, zorder=3)
    for k, e in enumerate(el): ax.scatter(*A[k], s=330, color=col[e], edgecolor='k', zorder=4); ax.text(A[k][0], A[k][1], e, color='w', fontsize=13, fontweight='bold', ha='center', va='center', zorder=5)
    ax.set_aspect('equal'); ax.set_xlim(-2.9, 2.9); ax.set_ylim(-2.9, 3.6); ax.set_xlabel(r'x ($\mathrm{\AA}$)')
panel(axs[0], vL, AL, [mL.atom_symbol(k) for k in range(4)], [(0, 1, '-'), (0, 2, '-'), (3, 1, ':'), (3, 2, ':')], 'LiCO$_2$ HOMO (SOMO), %.2f eV' % (fL.mo_energy[0][na - 1] * EV))
panel(axs[1], vC, AC, [mC.atom_symbol(k) for k in range(3)], [(0, 1, '-'), (0, 2, '-')], r'free CO$_2$ LUMO ($\pi^*$), %.2f eV' % (fC.mo_energy[best] * EV))
axs[0].set_ylabel(r'y ($\mathrm{\AA}$)'); fig.tight_layout(pad=0.3, w_pad=0.6); fig.savefig(out, dpi=170, bbox_inches='tight', pad_inches=0.05); print('contour levels (a.u.): +-', lev.tolist()); print('saved', out)
