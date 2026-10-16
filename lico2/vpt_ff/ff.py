"""(U)B3LYP-D3(BJ)/aug-cc-pVDZ: tight optimization, harmonic analysis, and a semi-quartic
force field from analytic Hessians at +/- displacements along each normal mode.
usage: python3 ff.py co2|lico2
Saves ff_<tag>.npz: el, m (amu), X (Bohr, principal axes), omega (Eh), L (mass-weighted Cartesian
normal vectors, 3N x nvib), Hq0/Hqp/Hqm (normal-coordinate Hessians), dQ (mass-weighted step).
"""
import os, sys, time, itertools, numpy as np
from pyscf import gto, dft, lib
from pyscf.geomopt.geometric_solver import optimize
lib.num_threads(int(os.environ.get('OMP_NUM_THREADS', '1')))
TAG = sys.argv[1]
BAS, XC, GRID = 'aug-cc-pvdz', 'b3lyp', 8
AMU = 1822.888486209; HC = 219474.6313705; ANG = 1 / 0.52917721092
MASS = {'Li': 7.016003437, 'C': 12.0, 'O': 15.994914620}
DQ = 2.0   # step in mass-weighted normal coordinate (Bohr amu_e^1/2), ~0.025 A amu^1/2
t0 = time.time()
def say(*a): print(f'[{time.time()-t0:6.0f}s]', *a, flush=True)

if TAG == 'co2':
    el = ['C', 'O', 'O']; X0 = np.array([[0, 0, 0], [0, 0, 1.167], [0, 0, -1.167]]); spin = 0
else:
    el = ['C', 'O', 'O', 'Li']
    X0 = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'start_geom', 'X_d3.npy')); spin = 1

def mkmf(X, unit='Angstrom'):
    mol = gto.M(atom=[[e, x] for e, x in zip(el, X)], unit=unit, basis=BAS, spin=spin,
                symmetry=False, verbose=0)
    mf = (dft.UKS if spin else dft.RKS)(mol).set(xc=XC, conv_tol=1e-12, max_cycle=200)
    mf.grids.level = GRID; mf.disp = 'd3bj'
    return mf

import os
CK = f'ckpt_{TAG}.npz'
if os.path.exists(CK):
    c = np.load(CK); X, H0, omega, L = c['X'], c['H0'], c['omega'], c['L']
    m = np.array([MASS[e] for e in el]); N = len(el)
    mw = np.repeat(m * AMU, 3); s = 1 / np.sqrt(mw); nvib = len(omega)
    dm_eq = c['dm']; E_eq = float(c['E'])
    def hess_cart(X, dm0=None):
        mf = mkmf(X, unit='Bohr'); mf.kernel(dm0=dm0); assert mf.converged
        h = mf.Hessian().kernel(); return h.transpose(0, 2, 1, 3).reshape(3 * N, 3 * N), mf
    def Hq(Hc): return L.T @ (Hc * np.outer(s, s)) @ L
    Hq0 = Hq(H0); Hqp, Hqm, done = c['Hqp'], c['Hqm'], list(c['done'])
    say('resumed from checkpoint; modes done:', done)
else:
    conv = {'convergence_energy': 1e-9, 'convergence_grms': 1e-6, 'convergence_gmax': 1.5e-6,
            'convergence_drms': 4e-6, 'convergence_dmax': 6e-6}
    g_ = mkmf(X0).nuc_grad_method(); g_.grid_response = True
    mol_eq = optimize(g_, maxsteps=100, **conv)
    X = mol_eq.atom_coords()                          # Bohr
    m = np.array([MASS[e] for e in el]); N = len(el)
    # principal axes
    X = X - (m[:, None] * X).sum(0) / m.sum()
    I = np.einsum('i,ia,ib->ab', m, X, X); I = np.eye(3) * np.trace(I) - I
    w_I, R = np.linalg.eigh(I); X = X @ R
    ops_ = []
    for sg in itertools.product((1, -1), repeat=3):
        for pm in itertools.permutations(range(N)):
            if any(el[p] != el[i] for i, p in enumerate(pm)): continue
            Y = (X * np.array(sg))[list(pm)]
            if np.abs(Y - X).max() < 1e-3: ops_.append(Y)
    X = np.mean(ops_, axis=0)
    say(f'symmetrized over {len(ops_)} operations')
    say('optimized; principal moments (amu Bohr^2):', np.round(w_I, 4))

    def hess_cart(X, dm0=None):
        mf = mkmf(X, unit='Bohr'); mf.kernel(dm0=dm0)
        assert mf.converged
        h = mf.Hessian().kernel()                     # (N,N,3,3)
        return h.transpose(0, 2, 1, 3).reshape(3 * N, 3 * N), mf

    H0, mf0 = hess_cart(X)
    dm_eq = mf0.make_rdm1()
    say('E_eq =', mf0.e_tot, ' <S2> =', mf0.spin_square()[0] if spin else 0)
    # mass-weighted Hessian, projected on translations/rotations
    mw = np.repeat(m * AMU, 3); s = 1 / np.sqrt(mw)
    Hm = H0 * np.outer(s, s)
    TR = []
    for a in range(3):
        t = np.zeros((N, 3)); t[:, a] = 1; TR.append((t * np.sqrt(mw.reshape(N, 3))).ravel())
        e = np.zeros(3); e[a] = 1; r = np.cross(e, X); TR.append((r * np.sqrt(mw.reshape(N, 3))).ravel())
    U, sv, _ = np.linalg.svd(np.array(TR).T, full_matrices=False); U = U[:, sv > 1e-6 * sv.max()]
    P = np.eye(3 * N) - U @ U.T
    ev, V = np.linalg.eigh(P @ Hm @ P)
    nvib = 3 * N - U.shape[1]
    idx = np.argsort(ev)[-nvib:]; ev, L = ev[idx], V[:, idx]
    omega = np.sqrt(ev)
    say(f'nvib={nvib} harmonic (cm-1):', np.round(omega * HC, 2))

    def Hq(Hc):  # Cartesian Hessian -> normal-coordinate Hessian (mass-weighted Q)
        return L.T @ (Hc * np.outer(s, s)) @ L
    Hq0 = Hq(H0)
    Hqp = np.zeros((nvib, nvib, nvib)); Hqm = np.zeros_like(Hqp)
    E_eq = mf0.e_tot; done = []
    np.savez(CK, X=X, H0=H0, omega=omega, L=L, dm=np.asarray(dm_eq), E=E_eq,
             Hqp=Hqp, Hqm=Hqm, done=np.array(done, int))
# symmetry operations (atom permutation + axis sign flips) that leave the equilibrium structure invariant
ops = []
for signs in itertools.product((1, -1), repeat=3):
    for perm in itertools.permutations(range(N)):
        if any(el[p] != el[i] for i, p in enumerate(perm)): continue
        Y = (X * np.array(signs))[list(perm)]
        if np.abs(Y - X).max() < 1e-4:
            R = np.zeros((3 * N, 3 * N))
            for i, p in enumerate(perm):          # new atom i = old atom p with signs
                R[3 * i:3 * i + 3, 3 * p:3 * p + 3] = np.diag(signs)
            ops.append(R)
say(f'{len(ops)} symmetry operations found')
Smw = np.diag(1 / s)   # mass-weighting commutes with ops (same element masses)
for k in range(nvib):
    if k in done: continue
    d = (s * L[:, k] * DQ).reshape(N, 3)
    Hp_c = hess_cart(X + d, dm_eq)[0]
    Hqp[k] = Hq(Hp_c)
    flip = [R for R in ops if np.abs(R @ L[:, k] + L[:, k]).max() < 1e-3]
    if flip:
        R = flip[0]; Hqm[k] = Hq(R @ Hp_c @ R.T); how = 'by symmetry'
    else:
        Hqm[k] = Hq(hess_cart(X - d, dm_eq)[0]); how = 'computed'
    done.append(k)
    np.savez(CK, X=X, H0=H0, omega=omega, L=L, dm=np.asarray(dm_eq), E=E_eq,
             Hqp=Hqp, Hqm=Hqm, done=np.array(done, int))
    say(f'mode {k} ({omega[k]*HC:.1f}) displaced Hessians done (minus side {how})')
# dipole at equilibrium (for reference) and Mulliken/IAO-ready density saved separately
np.savez(f'ff_{TAG}.npz', el=el, m=m, X=X, omega=omega, L=L, Hq0=Hq0, Hqp=Hqp, Hqm=Hqm,
         dQ=DQ, E=E_eq, dm=np.asarray(dm_eq))
say('ALL DONE')
