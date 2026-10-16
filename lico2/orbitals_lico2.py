"""HOMO and LUMO of LiCO2 at its equilibrium geometry (UB3LYP-D3(BJ)), as isosurfaces.  usage: python orbitals_lico2.py atz|adz out.png
HOMO = highest occupied spin-orbital (the alpha SOMO); LUMO = lowest unoccupied spin-orbital (alpha or beta, whichever is lower)."""
import sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure
from pyscf import gto, dft, lib, symm
import coords as c
tag, out = sys.argv[1], sys.argv[2]; BAS = {'atz': 'aug-cc-pvtz', 'adz': 'aug-cc-pvdz'}[tag]; EV = 27.211386; lib.num_threads(16)
X = c.X_of_q(np.load(f'cache_{tag}/qmin.npy')[None])[0]
mol = gto.M(atom=[[e, x] for e, x in zip(c.EL, X)], basis=BAS, unit='Bohr', spin=1, symmetry=True, verbose=0)
mf = dft.UKS(mol).set(xc='b3lyp', conv_tol=1e-10, max_cycle=200); mf.grids.level = 5; mf.disp = 'd3bj'; mf.kernel()
print('point group', mol.groupname, '| E = %.10f Eh (reference E0 %.10f) | converged %s | <S2> = %.4f' % (mf.e_tot, float(np.load(f'cache_{tag}/E0.npy')), mf.converged, mf.spin_square()[0]))
na, nb = mol.nelec
def val_frac(C):                                   # weight of the orbital in a minimal valence basis
    pm = gto.M(atom=mol._atom, basis='minao', unit='Bohr', spin=mol.spin, verbose=0)
    S12 = gto.intor_cross('int1e_ovlp', mol, pm); P = S12 @ np.linalg.solve(pm.intor('int1e_ovlp'), S12.T); return np.einsum('pi,pq,qi->i', C, P, C)
lab = [symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mf.mo_coeff[s]) for s in (0, 1)]
for s, n, name in ((0, na, 'alpha'), (1, nb, 'beta')):
    vf = val_frac(mf.mo_coeff[s][:, n - 2:n + 4])
    print(name, 'orbitals around the gap (index rel. to first virtual, eV, irrep, valence fraction):', [(i - 2, round(float(mf.mo_energy[s][n - 2 + i] * EV), 3), lab[s][n - 2 + i], round(float(vf[i]), 2)) for i in range(6)])
homo = (0, na - 1); lumo = (0, na) if mf.mo_energy[0][na] <= mf.mo_energy[1][nb] else (1, nb); SP = ('α', 'β')
orbs = [('HOMO (SOMO, %s)' % SP[homo[0]], homo), ('LUMO (%s)' % SP[lumo[0]], lumo)]
# frame for drawing: molecular plane = xy, C2 axis (C -> Li) pointing down
R = mol.atom_coords(); com = R.mean(0); nrm = np.linalg.svd(R - com)[2][2]; ax_ = R[3] - R[0]; ax_ /= np.linalg.norm(ax_); e1 = np.cross(ax_, nrm); F = np.array([e1, -ax_, nrm]); o = 0.5 * (R[0] + R[3])
def plot(ax, s, i, title, L, iso):
    g = np.linspace(-L, L, 91); gx, gy, gz = np.meshgrid(g, g, g, indexing='ij'); pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1) @ F + o
    val = (dft.numint.eval_ao(mol, pts) @ mf.mo_coeff[s][:, i]).reshape(gx.shape); sp = (g[1] - g[0],) * 3
    for sgn, col in ((1, '#1f77b4'), (-1, '#ff7f0e')):
        if (sgn * val).max() > iso:
            vs, fs, _, _ = measure.marching_cubes(sgn * val, iso, spacing=sp); ax.add_collection3d(Poly3DCollection(vs[fs] - L, alpha=0.42, facecolor=col, edgecolor='none'))
    A = (R - o) @ F.T; colr = {'C': '#444444', 'O': '#d62728', 'Li': '#9467bd'}; el = [mol.atom_symbol(k) for k in range(4)]
    for k, l, st in ((0, 1, '-'), (0, 2, '-'), (3, 1, ':'), (3, 2, ':')): ax.plot(*zip(A[k], A[l]), color='k', lw=3 if st == '-' else 2, ls=st, zorder=5)
    for k in range(4): ax.scatter(*A[k], s=260, color=colr[el[k]], edgecolor='k', depthshade=False, zorder=6); ax.text(A[k][0], A[k][1] + (0.9 if k != 3 else -1.2), A[k][2], el[k], fontsize=17, ha='center', zorder=7)
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(-L, L); ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(elev=28, azim=-62)
    ax.set_title(title, fontsize=18, y=0.97)
plt.rcParams.update({'font.size': 18}); fig = plt.figure(figsize=(14, 7.2))
for k, (name, (s, i)) in enumerate(orbs):
    vfi = float(val_frac(mf.mo_coeff[s][:, i:i + 1])[0]); diffuse = vfi < 0.5
    L, iso = (9.0, 0.012) if diffuse else (5.5, 0.04)
    plot(fig.add_subplot(1, 2, k + 1, projection='3d'), s, i, '%s: %s, %.2f eV\nisovalue %.3g a.u.' % (name, lab[s][i], mf.mo_energy[s][i] * EV, iso), L, iso)
    print('%s: spin %s, index %d, irrep %s, energy %.3f eV, valence fraction %.2f, box +-%.1f bohr, isovalue %.3g' % (name, 'ab'[s], i, lab[s][i], mf.mo_energy[s][i] * EV, vfi, L, iso))
fig.subplots_adjust(left=0, right=1, bottom=0, top=0.93, wspace=0); fig.savefig(out, dpi=170); print('saved', out)
