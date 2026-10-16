import os, numpy as np, json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from pyscf import gto, dft, lo, lib
lib.num_threads(1)
EV = 27.211386245988; BAS = 'aug-cc-pvdz'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures') + '/'

def scf(tag, ghost_li=False, bent_co2=False):
    d = np.load(f'ff_{tag}.npz'); el = list(d['el']); X = d['X']
    spin = 1 if tag == 'lico2' else 0
    if bent_co2:   # CO2 frozen at the LiCO2 geometry, Li as ghost (same AO basis)
        atoms = [[e, x] for e, x in zip(el[:3], X[:3])] + [['ghost-Li', X[3]]]
        basis = {'C': BAS, 'O': BAS, 'ghost-Li': gto.basis.load(BAS, 'Li')}
        mol = gto.M(atom=atoms, basis=basis, unit='Bohr', spin=0, verbose=0)
        mf = dft.RKS(mol)
    else:
        mol = gto.M(atom=[[e, x] for e, x in zip(el, X)], basis=BAS, unit='Bohr', spin=spin, verbose=0)
        mf = (dft.UKS if spin else dft.RKS)(mol)
    mf.xc = 'b3lyp'; mf.grids.level = 5; mf.conv_tol = 1e-11; mf.disp = 'd3bj'
    mf.kernel(dm0=None if bent_co2 else d['dm'])
    return mol, mf, X, el

def iao_charges(mol, mf):
    S = mol.intor('int1e_ovlp'); pmol = lo.iao.reference_mol(mol)
    mo = mf.mo_coeff if mf.mo_coeff.ndim == 3 else [mf.mo_coeff] * 1
    occs = mf.mo_occ if mf.mo_occ.ndim == 2 else [mf.mo_occ]
    pop = np.zeros(pmol.natm)
    for C, o in zip(mo, occs):
        Cocc = C[:, o > 0]; fac = o[o > 0]
        A = lo.iao.iao(mol, Cocc); A = lo.orth.vec_lowdin(A, S)
        Pm = A.T @ S @ (Cocc * fac) @ Cocc.T @ S @ A
        for i, (b0, b1) in enumerate(pmol.aoslice_by_atom()[:, 2:]): pop[i] += np.trace(Pm[b0:b1, b0:b1])
    return mol.atom_charges() - pop

res = {}
mol2, mf2, X2, el2 = scf('co2')
molL, mfL, XL, elL = scf('lico2')
molB, mfB, _, _ = scf('lico2', bent_co2=True)

# ---------- charges ----------
mulL = mfL.mulliken_pop(verbose=0)[1]
from pyscf.scf import uhf
spinpop = uhf.mulliken_spin_pop(molL, mfL.make_rdm1(), verbose=0)[1]
iaoL = iao_charges(molL, mfL); iao2 = iao_charges(mol2, mf2)
res['charges'] = {'LiCO2_mulliken': dict(zip(['C', 'O1', 'O2', 'Li'], np.round(mulL, 3).tolist())),
                  'LiCO2_IAO': dict(zip(['C', 'O1', 'O2', 'Li'], np.round(iaoL, 3).tolist())),
                  'LiCO2_mulliken_spin': dict(zip(['C', 'O1', 'O2', 'Li'], np.round(spinpop, 3).tolist())),
                  'CO2_IAO': dict(zip(['C', 'O1', 'O2'], np.round(iao2, 3).tolist()))}

# ---------- orbital energies and SOMO / LUMO character ----------
na = molL.nelec[0]
Ca, ea = mfL.mo_coeff[0], mfL.mo_energy[0]; Cb, eb = mfL.mo_coeff[1], mfL.mo_energy[1]
somo = Ca[:, na - 1]
SL = molL.intor('int1e_ovlp')
nB = molB.nelectron // 2
lumoB = mfB.mo_coeff[:, nB]                         # bent-CO2 LUMO (Li ghost basis = same AOs)
ovl = float((somo @ SL @ lumoB) ** 2)
# Li character of SOMO (Mulliken, per AO block)
li0, li1 = molL.aoslice_by_atom()[3, 2:]
SC = SL @ somo; li_frac = float(np.sum(somo[li0:li1] * SC[li0:li1]))
n2 = mol2.nelectron // 2
res['orbitals'] = {'CO2_HOMO_eV': float(mf2.mo_energy[n2 - 1] * EV),
                   'CO2_KS_LUMO_diffuse_sigma_eV': float(mf2.mo_energy[n2] * EV), 'CO2_pistar_eV': float(mf2.mo_energy[np.load('co2_pistar_idx.npy')[0]] * EV),
                   'bentCO2_LUMO_eV': float(mfB.mo_energy[nB] * EV),
                   'LiCO2_SOMO_eV': float(ea[na - 1] * EV),
                   'SOMO_overlap2_with_bentCO2_LUMO': round(ovl, 3),
                   'SOMO_Li_Mulliken_fraction': round(li_frac, 3)}

# ---------- PDOS ----------
def pdos(mol, C, e, occ, groups, grid, sig=0.25):
    S = mol.intor('int1e_ovlp'); W = C * (S @ C)
    sl = mol.aoslice_by_atom()[:, 2:]
    out = {}
    for name, atoms in groups.items():
        wt = sum(W[sl[a][0]:sl[a][1]].sum(0) for a in atoms)
        out[name] = (wt[None, :] * np.exp(-((grid[:, None] - e[None, :] * EV) / sig) ** 2 / 2)).sum(1) / (sig * np.sqrt(2 * np.pi))
    return out
g = np.linspace(-22, 6, 1400)
p2 = pdos(mol2, mf2.mo_coeff, mf2.mo_energy, mf2.mo_occ, {'C': [0], 'O': [1, 2]}, g)
pa = pdos(molL, Ca, ea, mfL.mo_occ[0], {'C': [0], 'O': [1, 2], 'Li': [3]}, g)
pb = pdos(molL, Cb, eb, mfL.mo_occ[1], {'C': [0], 'O': [1, 2], 'Li': [3]}, g)
fig, ax = plt.subplots(figsize=(10, 5.6))
co2tot = p2['C'] + p2['O']
ax.plot(g, co2tot, color='0.45', ls='--', lw=1.3, label='CO$_2$ total (per spin)')
ax.plot(g, -co2tot, color='0.45', ls='--', lw=1.3)
for sgn, p, lab in ((1, pa, r'$\alpha$'), (-1, pb, r'$\beta$')):
    ax.plot(g, sgn * (p['C'] + p['O']), color='#D85A30', lw=1.6, label=f'LiCO$_2$: CO$_2$ part' if sgn > 0 else None)
    ax.plot(g, sgn * p['Li'], color='#534AB7', lw=1.6, label='LiCO$_2$: Li part' if sgn > 0 else None)
ax.axhline(0, color='k', lw=0.6)
eS = ea[na - 1] * EV; eL2 = mf2.mo_energy[np.load('co2_pistar_idx.npy')[0]] * EV
ax.axvline(eS, color='#D85A30', lw=0.9, ls=':'); ax.axvline(eL2, color='0.45', lw=0.9, ls=':')
ax.annotate(f'LiCO$_2$ SOMO ({eS:.2f} eV)', (eS, ax.get_ylim()[1] * 0.85), xytext=(eS - 8.5, ax.get_ylim()[1] * 0.85),
            arrowprops=dict(arrowstyle='->', lw=0.7), fontsize=9)
ax.annotate(f'CO$_2$ $2\\pi_u^*$ ({eL2:.2f} eV)', (eL2, ax.get_ylim()[1] * 0.65), xytext=(eL2 + 0.6, ax.get_ylim()[1] * 0.65),
            fontsize=9)
ax.set_xlabel('Kohn–Sham orbital energy (eV, vacuum-referenced; both molecules neutral)')
ax.set_ylabel(r'PDOS (states/eV)   $\alpha$ up / $\beta$ down')
ax.set_title('B3LYP-D3(BJ)/aug-cc-pVDZ: Mulliken-projected DOS, Gaussian 0.25 eV')
ax.set_xlim(-22, 6); ax.legend(loc='lower left', fontsize=9, frameon=False)
fig.tight_layout(); fig.savefig(OUT + 'pdos_co2_lico2.png', dpi=160)

# ---------- orbital contour plots in the molecular plane ----------
def plane_frame(X, el):
    C = X[0]
    if len(X) == 4:
        e2 = X[3] - C; e2 /= np.linalg.norm(e2)            # C2 axis toward Li (vertical)
        e1 = X[1] - X[2]; e1 -= (e1 @ e2) * e2; e1 /= np.linalg.norm(e1)
    else:
        e1 = X[1] - X[2]; e1 /= np.linalg.norm(e1)
        t = np.array([1.0, 0, 0]) if abs(e1[0]) < 0.9 else np.array([0, 1.0, 0])
        e2 = t - (t @ e1) * e1; e2 /= np.linalg.norm(e2)
    return C, e1, e2
def grid_vals(mol, Cvec, C0, e1, e2, lim=4.0, ngr=220):
    u = np.linspace(-lim, lim, ngr) * 1.8897; v = np.linspace(-lim * 0.8, lim * 1.1, ngr) * 1.8897
    U_, V_ = np.meshgrid(u, v); pts = C0 + U_[..., None] * e1 + V_[..., None] * e2
    ao = mol.eval_gto('GTOval', pts.reshape(-1, 3))
    return U_ / 1.8897, V_ / 1.8897, (ao @ Cvec).reshape(U_.shape)
C0, e1, e2 = plane_frame(X2, el2)
# in-plane component of the degenerate CO2 LUMO pair (plane spanned by e1, e2)
pair = np.load('co2_pistar_idx.npy'); L1, L2 = mf2.mo_coeff[:, pair[0]], mf2.mo_coeff[:, pair[1]]
_, _, v1 = grid_vals(mol2, L1, C0, e1, e2, ngr=60); _, _, v2 = grid_vals(mol2, L2, C0, e1, e2, ngr=60)
M = np.array([[np.sum(v1 * v1), np.sum(v1 * v2)], [np.sum(v1 * v2), np.sum(v2 * v2)]])
c = np.linalg.eigh(M)[1][:, -1]; Lin = c[0] * L1 + c[1] * L2
panels = [('CO$_2$ valence $2\\pi_u^*$ (in-plane), empty', mol2, Lin, X2, el2, plane_frame(X2, el2)),
          ('LiCO$_2$ SOMO ($\\alpha$), 1 electron', molL, somo, XL, elL, plane_frame(XL, elL))]
fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
for ax, (title, mol, vec, X, el, fr) in zip(axs, panels):
    U_, V_, val = grid_vals(mol, vec, *fr)
    lv = np.array([0.01, 0.02, 0.04, 0.08, 0.16])
    ax.contour(U_, V_, val, levels=lv, colors='#D85A30', linewidths=0.9)
    ax.contour(U_, V_, val, levels=-lv[::-1], colors='#185FA5', linewidths=0.9, linestyles='solid')
    for e, x in zip(el, X):
        r = x - fr[0]; px, py = r @ fr[1] / 1.8897, r @ fr[2] / 1.8897
        ax.plot(px, py, 'o', ms={'C': 9, 'O': 10, 'Li': 9}[e],
                color={'C': '0.3', 'O': '#E24B4A', 'Li': '#7F77DD'}[e], mec='k', mew=0.5)
        ax.text(px + 0.18, py + 0.18, e, fontsize=10)
    ax.set_aspect('equal'); ax.set_title(title, fontsize=11); ax.set_xlabel('Å'); ax.set_ylabel('Å')
fig.suptitle('Molecular-plane contours (±0.01 … ±0.16 a.u.; orange +, blue −)', fontsize=10)
fig.tight_layout(); fig.savefig(OUT + 'orbitals_co2_lumo_vs_lico2_somo.png', dpi=160)

# ---------- 3D isosurfaces ----------
from skimage import measure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
def iso_panel(ax, mol, vec, X, el, title, iso=0.03):
    lo_, hi_ = X.min(0) - 6, X.max(0) + 6; ng = 70
    axes = [np.linspace(lo_[i], hi_[i], ng) for i in range(3)]
    G = np.stack(np.meshgrid(*axes, indexing='ij'), -1).reshape(-1, 3)
    val = (mol.eval_gto('GTOval', G) @ vec).reshape(ng, ng, ng)
    sp = [(a[1] - a[0]) for a in axes]
    for sgn, col in ((1, '#D85A30'), (-1, '#378ADD')):
        if (sgn * val).max() < iso: continue
        vs, fs, _, _ = measure.marching_cubes(sgn * val, iso, spacing=sp)
        vs = (vs + lo_) / 1.8897
        ax.add_collection3d(Poly3DCollection(vs[fs], alpha=0.35, facecolor=col, edgecolor='none'))
    Xa = X / 1.8897
    for e, x in zip(el, Xa):
        ax.scatter(*x, s={'C': 120, 'O': 140, 'Li': 120}[e], color={'C': '0.3', 'O': '#E24B4A', 'Li': '#7F77DD'}[e], edgecolor='k')
    for i, j in ([(0, 1), (0, 2)] + ([(1, 3), (2, 3)] if len(X) == 4 else [])):
        ax.plot(*np.array([Xa[i], Xa[j]]).T, color='0.4', lw=1.5 if j != 3 else 0.8, ls='-' if j != 3 else '--')
    c = Xa.mean(0); r = 2.6
    ax.set_xlim(c[0] - r, c[0] + r); ax.set_ylim(c[1] - r, c[1] + r); ax.set_zlim(c[2] - r, c[2] + r)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.set_title(title, fontsize=11)
fig = plt.figure(figsize=(11, 5.2))
iso_panel(fig.add_subplot(121, projection='3d'), mol2, Lin, X2, el2, 'CO$_2$ valence $2\\pi_u^*$ (in-plane)')
iso_panel(fig.add_subplot(122, projection='3d'), molL, somo, XL, elL, 'LiCO$_2$ SOMO (α)')
fig.suptitle('Isosurface ±0.03 a.u.', fontsize=10); fig.tight_layout()
fig.savefig(OUT + 'orbitals_3d.png', dpi=160)
json.dump(res, open('analysis.json', 'w'), indent=1); print(json.dumps(res, indent=1))
