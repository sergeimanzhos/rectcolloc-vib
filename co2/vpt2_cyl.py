"""Numerical GVPT2 (and VCI on the same Hamiltonian) from ff_<tag>.npz.
Watson J=0 Hamiltonian in dimensionless normal coordinates q (atomic units):
  H = sum w_k/2 (p_k^2+q_k^2) + V3 + V4 + sum_a B_a pi_a^2 (+ U = -sum B/4, nonlinear only)
For a linear molecule the axial term is absent (J_z = pi_z).
GVPT2: states with <= 2 quanta; near-degenerate (|dE0| < THR) states up to 3 quanta form polyads
treated with a 2nd-order (Van Vleck) effective Hamiltonian; other states by plain 2nd-order RSPT.
usage: python3 vpt2.py co2|lico2
"""
import sys, itertools, json, numpy as np
HC = 219474.6313705; AMU = 1822.888486209
import os
TAG = sys.argv[1]; K = int(os.environ.get('KMAX', 7)); THR = 100 / HC
d = np.load(f'ff_{TAG}.npz')
w, L, m, X = d['omega'], d['L'], d['m'], d['X']
n, N, dQ = len(w), len(m), float(d['dQ'])
Hq0, Hqp, Hqm = d['Hq0'], d['Hqp'], d['Hqm']

# ---------- force constants in mass-weighted Q, then dimensionless q ----------
acc = {}
def add(key, val):
    key = tuple(sorted(key)); acc.setdefault(key, []).append(val)
for k in range(n):
    D1 = (Hqp[k] - Hqm[k]) / (2 * dQ); D2 = (Hqp[k] + Hqm[k] - 2 * Hq0) / dQ ** 2
    for i in range(n):
        for j in range(n):
            add((i, j, k), D1[i, j]); add((i, j, k, k), D2[i, j])
cub, qua = {}, {}
for key, v in acc.items():
    f = np.mean(v) / np.sqrt(np.prod(w[list(key)]))
    (cub if len(key) == 3 else qua)[key] = f
if TAG.endswith('cyl'):
    # cylindrically symmetric quartic force field for the linear molecule: modes 0, 1 = bend components, 2 = sym., 3 = asym. stretch.
    # Only terms allowed by symmetry are kept, bend constants are averaged over the two components and phi_aabb = phi_aaaa / 3,
    # so that V depends on the bend only through qa^2 + qb^2.
    kbb1 = 0.5 * (cub[(0, 0, 2)] + cub[(1, 1, 2)]); pb = 0.5 * (qua[(0, 0, 0, 0)] + qua[(1, 1, 1, 1)])
    p11b = 0.5 * (qua[(0, 0, 2, 2)] + qua[(1, 1, 2, 2)]); p33b = 0.5 * (qua[(0, 0, 3, 3)] + qua[(1, 1, 3, 3)])
    print('bend quartics before symmetrisation (cm-1): aaaa %.2f bbbb %.2f aabb %.2f (cylindrical value %.2f)' % (qua[(0, 0, 0, 0)] * HC, qua[(1, 1, 1, 1)] * HC, qua[(0, 0, 1, 1)] * HC, pb / 3 * HC))
    C = dict(w=[float(0.5 * (w[0] + w[1])), float(w[2]), float(w[3])], k111=float(cub[(2, 2, 2)]), k133=float(cub[(2, 3, 3)]), k1bb=float(kbb1), p1111=float(qua[(2, 2, 2, 2)]), p3333=float(qua[(3, 3, 3, 3)]),
             p1133=float(qua[(2, 2, 3, 3)]), p11bb=float(p11b), p33bb=float(p33b), pbbbb=float(pb))
    cub = {(0, 0, 2): kbb1, (1, 1, 2): kbb1, (2, 2, 2): C['k111'], (2, 3, 3): C['k133']}
    qua = {(0, 0, 0, 0): pb, (1, 1, 1, 1): pb, (0, 0, 1, 1): pb / 3, (0, 0, 2, 2): p11b, (1, 1, 2, 2): p11b, (0, 0, 3, 3): p33b, (1, 1, 3, 3): p33b,
           (2, 2, 2, 2): C['p1111'], (3, 3, 3, 3): C['p3333'], (2, 2, 3, 3): C['p1133']}
    w = np.array([C['w'][0], C['w'][0] * (1 + 1e-7), C['w'][1], C['w'][2]])      # tiny splitting: the perturbation sums need non-degenerate zeroth-order states
    json.dump(C, open('vff_cyl.json', 'w'), indent=1); print('cylindrical force field (cm-1):', {k: (np.round(np.array(v) * HC, 3).tolist()) for k, v in C.items()})
# harmonic consistency check: diagonal of Hq0 vs w^2
print('Hq0 off-diagonal max / w^2 max:', np.abs(Hq0 - np.diag(np.diag(Hq0))).max() / (w ** 2).max())
import math
def taylor(key, phi):
    c = 1.0
    for a in set(key): c *= math.factorial(key.count(a))
    return phi / c

# ---------- basis ----------
S = np.array([s for s in itertools.product(range(K + 1), repeat=n) if sum(s) <= K]); NS = len(S)
E0 = (S + 0.5) @ w
nb = K + 8
a_ = np.diag(np.sqrt(np.arange(1, nb)), 1)
q_ = (a_ + a_.T) / np.sqrt(2); dq_ = (a_ - a_.T) / np.sqrt(2)          # q and d/dq
qpow = [np.linalg.matrix_power(q_, p)[:K + 1, :K + 1] for p in range(5)]
dq1 = dq_[:K + 1, :K + 1]

def op(factors):
    """dense <I| prod_k M_k |J> for factors {mode: 1D matrix}."""
    M = np.ones((NS, NS))
    for k in range(n):
        A = factors.get(k)
        if A is None: M *= (S[:, k][:, None] == S[:, k][None, :])
        else: M *= A[S[:, k][:, None], S[:, k][None, :]]
    return M

def mono(key):
    f = {}
    for k in set(key): f[k] = qpow[key.count(k)]
    return op(f)

H1 = np.zeros((NS, NS)); H2 = np.zeros((NS, NS))
for key, phi in cub.items():
    if abs(phi) > 1e-12: H1 += taylor(key, phi) * mono(key)
for key, phi in qua.items():
    if abs(phi) > 1e-12: H2 += taylor(key, phi) * mono(key)

# ---------- Coriolis ----------
Lc = L.reshape(N, 3, n)
zeta = np.zeros((3, n, n))
for k in range(n):
    for l in range(n): zeta[:, k, l] = np.cross(Lc[:, :, k], Lc[:, :, l]).sum(0)
I = np.einsum('i,ia,ib->ab', m * AMU, X, X); I = np.eye(3) * np.trace(I) - I
Id = np.diag(I); assert np.abs(I - np.diag(Id)).max() < 1e-6 * Id.max(), 'not principal axes'
linear = Id.min() < 1e-6 * Id.max()
B = np.array([0.0 if x < 1e-6 * Id.max() else 1 / (2 * x) for x in Id])
U = 0.0 if linear else -B.sum() / 4
for ax in range(3):
    if B[ax] == 0: continue
    R = np.zeros((NS, NS))
    for k in range(n):
        for l in range(n):
            if k != l and abs(zeta[ax, k, l]) > 1e-10:
                R += zeta[ax, k, l] * np.sqrt(w[l] / w[k]) * op({k: qpow[1], l: dq1})
    H2 += -B[ax] * (R @ R)        # pi = -i R  ->  pi^2 = -R^2
H2 += U * np.eye(NS)
print('linear:', linear, ' B (cm-1):', np.round(B * HC, 4))

# ---------- GVPT2 ----------
idx = {tuple(s): i for i, s in enumerate(S)}
targets = [i for i in range(NS) if S[i].sum() <= 2]
cand = [i for i in range(NS) if S[i].sum() <= 3]
# polyads: connected components of near-degeneracy among candidates linked to targets
parent = {i: i for i in cand}
H1sq = H1 @ H1
def find(i):
    while parent[i] != i: parent[i] = parent[parent[i]]; i = parent[i]
    return i
for a in targets:
    for b in cand:
        if b != a and abs(E0[a] - E0[b]) < THR and (abs(H1[a, b]) > 1e-9 or abs(H2[a, b]) > 1e-9 or
                                                   abs(H1sq[a, b]) > 1e-12):
            parent[find(a)] = find(b)
groups = {}
for i in cand: groups.setdefault(find(i), []).append(i)
res = {}
def heff(P):
    P = list(P); Pset = set(P)
    Q = np.array([i for i in range(NS) if i not in Pset])
    He = np.zeros((len(P), len(P)))
    for x, a in enumerate(P):
        for y, b in enumerate(P):
            v = H1[a, b] + H2[a, b] + (E0[a] if a == b else 0)
            v += 0.5 * np.sum(H1[a, Q] * H1[Q, b] * (1 / (E0[a] - E0[Q]) + 1 / (E0[b] - E0[Q])))
            He[x, y] = v
    return He
for g in groups.values():
    if not any(t in g for t in targets): continue
    if len(g) == 1:
        a = g[0]; Q = np.array([i for i in range(NS) if i != a])
        res[a] = E0[a] + H2[a, a] + np.sum(H1[a, Q] ** 2 / (E0[a] - E0[Q]))
        continue
    He = heff(g); ev, V = np.linalg.eigh(0.5 * (He + He.T))
    for k in range(len(g)):
        dom = g[np.argmax(V[:, k] ** 2)]
        mix = {''.join(map(str, S[g[j]])): round(float(V[j, k] ** 2), 3) for j in range(len(g)) if V[j, k] ** 2 > 0.05}
        res.setdefault(('poly', dom, k), (ev[k], mix))
# collect
g0 = idx[tuple([0] * n)]
Ezp = res[g0] if g0 in res else [v[0] for kk, v in res.items() if isinstance(kk, tuple) and kk[1] == g0][0]
print(f'ZPE harmonic {E0[g0]*HC:.1f}  GVPT2 {Ezp*HC:.1f} cm-1')
rows = []
for kk, v in res.items():
    if isinstance(kk, tuple):
        st, (E, mix) = kk[1], v
    else:
        st, E, mix = kk, v, None
    if S[st].sum() == 0: continue
    rows.append((''.join(map(str, S[st])), (E0[st] - E0[g0]) * HC, (E - Ezp) * HC, mix))
rows.sort(key=lambda r: r[2])

# ---------- VCI on the same Hamiltonian (convergence check of PT) ----------
Hfull = np.diag(E0) + H1 + H2
evv, Vv = np.linalg.eigh(Hfull)
vci = []
for k in range(min(80, NS)):
    dom = np.argmax(Vv[:, k] ** 2)
    vci.append((''.join(map(str, S[dom])), (evv[k] - evv[0]) * HC, float(Vv[dom, k] ** 2)))
print(f'\nmodes (harmonic cm-1): {np.round(w*HC,2)}')
print(f'{"state":>8} {"harm":>9} {"GVPT2":>9}   polyad mixing')
for r in rows: print(f'{r[0]:>8} {r[1]:9.1f} {r[2]:9.1f}   {r[3] if r[3] else ""}')
print('\nVCI (same quartic Hamiltonian, K=%d, %d states), lowest levels:' % (K, NS))
for r in vci[:26]: print(f'{r[0]:>8} {r[1]:9.1f}  weight {r[2]:.2f}')
json.dump({'w': (w * HC).tolist(), 'gvpt2': [(r[0], r[1], r[2], r[3]) for r in rows],
           'vci': vci, 'zeta': zeta.tolist(), 'B': (B * HC).tolist()},
          open(f'vpt2_{TAG}.json', 'w'), indent=1)
