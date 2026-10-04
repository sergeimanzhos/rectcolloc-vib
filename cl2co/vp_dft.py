"""DFT data for numerical VPT2 of Cl2CO.  usage: python3 vp_dft.py basis tag
Stage 1: B3LYP geometry optimisation (BFGS on analytic gradients, starts at NN-PES minimum).
Stage 2: Hessian by central differences of analytic gradients -> normal modes (C2v-classified, ordered
         as the collocation labels: nu3, nu6, nu2, nu4, nu5, nu1).
Stage 3: E, gradient, dipole on a stencil in dimensionless normal coordinates:
         +-d, +-2d along each mode; +-d+-d for every pair; +-d+-d+-d+-d in modes (nu1,nu2,nu5,nu6) for f_1256.
Every DFT result is cached (json) so restarts resume."""
import numpy as np, json, os, sys, time, itertools
from pyscf import gto, dft, lib
import scipy.optimize as so
BAS, TAG = sys.argv[1], sys.argv[2]
lib.num_threads(int(os.environ.get('NT', 2)))
H2CM = 219474.6313705
SYM = ['C', 'O', 'Cl', 'Cl']; m = np.array([21874.66183451, 29156.94569791, 63744.31892996, 63744.31892996])
X0nn = np.load('eq.npz')['X0']
CACHE = f'vp_{TAG}_cache.json'
try: C = json.load(open(CACHE))
except Exception: C = {}
LOG = open(f'vp_{TAG}.log', 'a')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
_dm = [None]
def calc(X, key):
    if key in C: r = C[key]; return r['E'], np.array(r['g']), np.array(r['mu'])
    mol = gto.M(atom=[(s, tuple(x)) for s, x in zip(SYM, X)], basis=BAS, unit='B', verbose=0, max_memory=5000)
    mf = dft.RKS(mol); mf.xc = 'b3lyp'; mf.grids.level = 3; mf.conv_tol = 1e-11; mf.conv_tol_grad = 1e-7
    mf.kernel(dm0=_dm[0])
    if not mf.converged: mf.kernel()
    assert mf.converged, key
    _dm[0] = mf.make_rdm1()
    g = mf.nuc_grad_method(); g.grid_response = True; G = g.kernel()
    mu = mf.dip_moment(unit='AU', verbose=0)
    C[key] = dict(E=float(mf.e_tot), g=G.tolist(), mu=list(mu), X=X.tolist())
    json.dump(C, open(CACHE + '.tmp', 'w')); os.replace(CACHE + '.tmp', CACHE)
    return mf.e_tot, G, mu

t0 = time.time()
# ---------- stage 1: optimisation (C2v kept by symmetric start; 12 Cartesians) ----------
if 'Xeq' not in C:
    it = [0]
    def f(x):
        it[0] += 1; E, G, _ = calc(x.reshape(4, 3), 'opt%03d' % it[0]); return E, G.ravel()
    r = so.minimize(f, X0nn.ravel(), jac=True, method='BFGS', options=dict(gtol=2e-6, maxiter=200))
    E, G, _ = calc(r.x.reshape(4, 3), 'eq')
    P('opt: %d evals, E=%.10f, max|g|=%.2e  %.0fs' % (it[0], E, abs(G).max(), time.time() - t0))
    C['Xeq'] = r.x.reshape(4, 3).tolist(); json.dump(C, open(CACHE, 'w'))
Xeq = np.array(C['Xeq'])
# principal-axis frame, COM at origin (planar molecule stays in yz-plane, C2 along z)
Xeq = Xeq - (m[:, None] * Xeq).sum(0) / m.sum()
# ---------- stage 2: Hessian ----------
h = 5e-3; Hc = np.zeros((12, 12))
for i in range(12):
    gp = calc((Xeq.ravel() + h * np.eye(12)[i]).reshape(4, 3), 'h%02dp' % i)[1].ravel()
    gm = calc((Xeq.ravel() - h * np.eye(12)[i]).reshape(4, 3), 'h%02dm' % i)[1].ravel()
    Hc[:, i] = (gp - gm) / (2 * h)
Hc = (Hc + Hc.T) / 2
mm = np.repeat(m, 3); Hm = Hc / np.sqrt(np.outer(mm, mm))
# project translations / rotations
sq = np.sqrt(mm); Vtr = []
for a in range(3):
    v = np.zeros((4, 3)); v[:, a] = 1; Vtr.append((v.ravel() * sq))
for a in range(3):
    e = np.eye(3)[a]; Vtr.append((np.cross(e, Xeq)).ravel() * sq)
Vtr = np.linalg.qr(np.array(Vtr).T)[0]; Vtr = Vtr[:, np.linalg.norm(Vtr, axis=0) > 1e-8]
Pr = np.eye(12) - Vtr @ Vtr.T; w2, L = np.linalg.eigh(Pr @ Hm @ Pr)
idx = np.argsort(w2)[6:]; w2 = w2[idx]; L = L[:, idx]; om = np.sqrt(w2)
# C2v classification: in-plane yz; C2 axis z. Swap Cl3<->Cl4 and y->-y gives A1 (+) or B2 (-); out-of-plane x -> B1
def irrep(l):
    d = l.reshape(4, 3) / np.sqrt(m)[:, None]
    if np.abs(d[:, 0]).max() > 0.5 * np.abs(d).max(): return 'B1'
    s = d[[0, 1, 3, 2]].copy(); s[:, 1] *= -1
    return 'A1' if (s * d).sum() / (d * d).sum() > 0 else 'B2'
ir = [irrep(L[:, k]) for k in range(6)]
A1 = sorted([k for k in range(6) if ir[k] == 'A1'], key=lambda k: -om[k])   # nu1, nu2, nu3
B2 = sorted([k for k in range(6) if ir[k] == 'B2'], key=lambda k: -om[k])   # nu5, nu6
B1 = [k for k in range(6) if ir[k] == 'B1']                                # nu4
nu = {1: A1[0], 2: A1[1], 3: A1[2], 4: B1[0], 5: B2[0], 6: B2[1]}
order = [nu[3], nu[6], nu[2], nu[4], nu[5], nu[1]]          # collocation label order
L = L[:, order]; om = om[order]
for k in range(6):                                          # sign convention: largest component positive
    if L[np.argmax(np.abs(L[:, k])), k] < 0: L[:, k] *= -1
P('harmonic (cm-1) nu3,nu6,nu2,nu4,nu5,nu1:', np.round(om * H2CM, 2), [ir[i] for i in order], '%.0fs' % (time.time() - t0))
Dq = (L / sq[:, None]) / np.sqrt(om)[None, :]              # dX/dq (12 x 6), q dimensionless
# ---------- stage 3: stencil ----------
d = float(os.environ.get('DELTA', 0.25))
pts = [np.zeros(6)]
for i in range(6):
    for s in (-2, -1, 1, 2): v = np.zeros(6); v[i] = s * d; pts.append(v)
for i, j in itertools.combinations(range(6), 2):
    for si, sj in itertools.product((-1, 1), repeat=2): v = np.zeros(6); v[i] = si * d; v[j] = sj * d; pts.append(v)
quad = [5, 2, 4, 1]                                         # nu1, nu2, nu5, nu6 in label order
for ss in itertools.product((-1, 1), repeat=4):
    v = np.zeros(6); v[quad] = np.array(ss) * d; pts.append(v)
pts = np.array(pts); E = []; Gq = []; MU = []
for n, q in enumerate(pts):
    X = (Xeq.ravel() + Dq @ q).reshape(4, 3)
    e, g, mu = calc(X, 's%.3f_' % d + '_'.join('%+.0f' % (x / d) for x in q))
    E.append(e); Gq.append(g.ravel() @ Dq); MU.append(mu)
    if n % 10 == 0: P('stencil %d/%d  %.0fs' % (n + 1, len(pts), time.time() - t0))
np.savez(f'vp_{TAG}.npz', q=pts, E=np.array(E), g=np.array(Gq), mu=np.array(MU), omega=om, L=L, Dq=Dq,
         Xeq=Xeq, m=m, Hc=Hc, delta=d, basis=BAS)
P('done %.0fs' % (time.time() - t0))
