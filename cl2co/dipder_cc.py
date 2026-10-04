"""CCSD(T) finite-field dipole derivatives along Eckart normal coordinates.
dmu_a/dq_k = -d2E/(dq_k dF_a), mixed central differences. C2v: each mode couples to one
dipole component only, so 4 energies per mode (+2 at equilibrium for mu0).
usage: python3 dipder_cc.py basis out.npz"""
import numpy as np, coords6 as c, time, sys, os
from pyscf import gto, scf, cc, lib
BAS = sys.argv[1]; OUT = sys.argv[2]
lib.num_threads(int(os.environ.get('NT', 2)))
F = 1e-3; d = 0.02
SYM = ['C', 'O', 'Cl', 'Cl']; NFROZ = 1 + 1 + 5 + 5
log = open(OUT.replace('.npz', '.log'), 'a')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); log.write(s + '\n'); log.flush()

import json as _js
CACHE = OUT.replace('.npz', '_cache.json')
try: _C = _js.load(open(CACHE))
except Exception: _C = {}
def energy(X, fvec, dm0=None, test=False):
    key = None if test else ','.join('%.10f' % v for v in np.concatenate([np.ravel(X), fvec]))
    if key in _C: return tuple(_C[key]) + (None,)
    r = _energy(X, fvec, dm0, test)
    if not test:
        _C[key] = list(r[:3]); _js.dump(_C, open(CACHE + '.tmp', 'w')); os.replace(CACHE + '.tmp', CACHE)
    return r
def _energy(X, fvec, dm0=None, test=False):
    mol = gto.M(atom=[(s, tuple(x)) for s, x in zip(SYM, X)], basis=BAS, unit='B',
                symmetry=False, verbose=0, max_memory=6000)
    r = mol.intor('int1e_r')                       # <i|r|j>, electrons have charge -1
    h = mol.intor('int1e_kin') + mol.intor('int1e_nuc') + np.einsum('x,xij->ij', fvec, r)
    enuc_f = -np.dot(fvec, mol.atom_charges() @ mol.atom_coords())
    mf = scf.RHF(mol); mf.get_hcore = lambda *a: h
    mf.energy_nuc = lambda *a: mol.energy_nuc() + enuc_f
    mf.conv_tol = 1e-12; mf.conv_tol_grad = 1e-8; mf.kernel(dm0=dm0); assert mf.converged
    if test: return mf
    mc = cc.CCSD(mf, frozen=NFROZ); mc.conv_tol = 1e-11; mc.conv_tol_normt = 1e-8
    mc.kernel(); assert mc.converged
    et = mc.ccsd_t()
    return mf.e_tot, mc.e_tot, mc.e_tot + et, mf.make_rdm1()

# sign/convention check at HF level: finite-field HF dipole vs analytic
if not os.environ.get('SKIPTEST'):
    ez = np.array([0, 0, 1.])
    m0 = energy(c.X0, 0 * ez, test=True)
    ep = energy(c.X0, F * ez, test=True).e_tot; em = energy(c.X0, -F * ez, test=True).e_tot
    P('HF check: FF mu_z', -(ep - em) / (2 * F), ' analytic', m0.dip_moment(unit='AU', verbose=0)[2])

comp = {0: 2, 1: 1, 2: 2, 3: 0, 4: 1, 5: 2}       # C2v component per mode (from DFT dmu)
t0 = time.time()
E0 = {}
for s in (+1, -1):
    fv = np.zeros(3); fv[2] = s * F
    E0[s] = energy(c.X0, fv)
    P('eq F_z=%+d' % s, E0[s][:3], '%.0fs' % (time.time() - t0))
mu0 = np.zeros((3, 3)); mu0[:, 2] = [-(E0[1][i] - E0[-1][i]) / (2 * F) for i in range(3)]
P('mu0_z HF/CCSD/CCSD(T):', mu0[:, 2])
dmu = np.zeros((3, c.NQ, 3))                       # [method, mode, component]
for k in range(c.NQ):
    a = comp[k]; E = {}
    for sq in (+1, -1):
        X = c.X_of_q(sq * d * np.eye(c.NQ)[k])[0]; dm = None
        for sf in (+1, -1):
            fv = np.zeros(3); fv[a] = sf * F
            E[sq, sf] = energy(X, fv, dm); dm = E[sq, sf][3]
    for i in range(3):
        dmu[i, k, a] = -(E[1, 1][i] - E[1, -1][i] - E[-1, 1][i] + E[-1, -1][i]) / (4 * d * F)
    P('mode', k, 'comp', 'xyz'[a], 'HF/CCSD/CCSD(T):', dmu[:, k, a], '%.0fs' % (time.time() - t0))
    np.savez(OUT, basis=BAS, method='CCSD(T)', mu0=mu0[2], dmu=dmu[2], dmu_hf=dmu[0],
             dmu_ccsd=dmu[1], delta=d, field=F, done=k + 1)
P('done %.0fs' % (time.time() - t0))
