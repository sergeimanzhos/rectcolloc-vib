"""6-mode GVPT2 (numerical, sparse harmonic-oscillator basis) for Cl2CO.  usage: python3 vp2.py adz
H0 = sum w_k (n_k+1/2); H1 = linear + cubic; H2 = quartic + off-diagonal quadratic + Coriolis (Watson J=0); U = -sum B/4.
Non-resonant states: standard 2nd-order RSPT (energy) and wavefunction through 2nd order.
Resonant states (Martin test on cubic couplings, |W/dE| test on 2nd-order effective couplings, and the
nu1 / nu2+nu5+nu6 pair forced): grouped into polyads P; effective Hamiltonian by 2nd-order Van Vleck with
resonant terms removed from the perturbation sum, diagonalised; wavefunctions = eigenvector mix of
{|a> + c1_a + c2_a} with Q-space (non-resonant) denominators.  Intensities from the full polynomial dipole
(linear+quadratic+cubic terms; orders 0,1,2) and, for comparison, from its linear part only."""
import numpy as np, itertools, sys, json, scipy.sparse as sp
from vp_fit import fit, H2CM
TAG = sys.argv[1]; K = int(sys.argv[2]) if len(sys.argv) > 2 else 10
cV, cM, info, d = fit(f'vp_{TAG}.npz')
NQ = 6; w = np.array([2 * cV[tuple(2 * (np.arange(6) == k))] for k in range(NQ)])
print('fit:', {k: (round(v, 4) if isinstance(v, float) else v) for k, v in info.items()})
print('omega (fit, cm-1) nu3,nu6,nu2,nu4,nu5,nu1:', np.round(w * H2CM, 2), ' Hessian:', np.round(d['omega'] * H2CM, 2))
# ---- basis ----
S = np.array([s for s in itertools.product(range(K + 1), repeat=NQ) if sum(s) <= K]); NS = len(S)
base = (K + 1) ** np.arange(NQ); LUT = -np.ones((K + 1) ** NQ, int); LUT[S @ base] = np.arange(NS)
print('basis states:', NS)
nb = K + 1 + 8; a = np.diag(np.sqrt(np.arange(1, nb)), 1); qb = (a + a.T) / np.sqrt(2); Pb = (a.T - a) / np.sqrt(2)
def q1(k): return np.linalg.matrix_power(qb, k)[:K + 1, :K + 1]
P1 = Pb[:K + 1, :K + 1]
def op(factors):
    """sparse operator = product over (mode, 1D matrix) factors (distinct modes)."""
    src = np.arange(NS); tgt = S.copy(); val = np.ones(NS)
    for k, M in factors:
        offs = sorted({r - c for r, c in zip(*np.nonzero(np.abs(M) > 1e-14))})
        ns, nt, nv = [], [], []
        for o in offs:
            n2 = tgt[:, k] + o; ok = (n2 >= 0) & (n2 <= K)
            v = np.zeros(len(src)); v[ok] = M[n2[ok], tgt[ok, k]] * val[ok]; keep = np.abs(v) > 0
            t = tgt[keep].copy(); t[:, k] = n2[keep]; ns.append(src[keep]); nt.append(t); nv.append(v[keep])
        src = np.concatenate(ns); tgt = np.concatenate(nt); val = np.concatenate(nv)
        ok = tgt.sum(1) <= K; src, tgt, val = src[ok], tgt[ok], val[ok]
    return sp.csr_matrix((val, (LUT[tgt @ base], src)), shape=(NS, NS))
def mono(al): return op([(k, q1(al[k])) for k in range(NQ) if al[k]])
E0 = (S + 0.5) @ w
H1 = sp.csr_matrix((NS, NS)); H2 = sp.csr_matrix((NS, NS))
for al, c in cV.items():
    dg = sum(al)
    if dg in (1, 3): H1 = H1 + c * mono(al)
    elif dg == 2 and max(al) == 1: H2 = H2 + c * mono(al)
    elif dg == 4: H2 = H2 + c * mono(al)
# Coriolis (Watson, J=0)
m = d['m']; X = d['Xeq']; L = d['L'].reshape(4, 3, NQ) * 1.0
zeta = np.zeros((3, NQ, NQ))
for k in range(NQ):
    for l in range(NQ): zeta[:, k, l] = np.cross(L[:, :, k], L[:, :, l]).sum(0)
I = np.einsum('i,ia,ib->ab', m, X, X); I = np.eye(3) * np.trace(I) - I
assert np.abs(I - np.diag(np.diag(I))).max() < 1e-6 * np.abs(I).max(), 'not principal axes'
B = 1 / (2 * np.diag(I)); U = -B.sum() / 4
for ax in range(3):
    R = sp.csr_matrix((NS, NS))
    for k in range(NQ):
        for l in range(NQ):
            if k != l and abs(zeta[ax, k, l]) > 1e-10:
                R = R + zeta[ax, k, l] * np.sqrt(w[l] / w[k]) * op([(k, q1(1)), (l, P1)])
    H2 = H2 - B[ax] * (R @ R)
H1 = H1.tocsc(); H2 = H2.tocsc(); H1r = H1.tocsr()
print('B (cm-1):', np.round(B * H2CM, 4))
# ---- targets & resonance detection ----
lab = lambda s: ''.join(str(x) for x in s)
idx = lambda s: int(LUT[np.array(s) @ base])
ZPE0 = E0[0]
T = [i for i in range(NS) if S[i].sum() <= 5 and (E0[i] - ZPE0) * H2CM < 1900]
NU1, COMB = (0, 0, 0, 0, 0, 1), (0, 1, 1, 0, 1, 0)
def col(M, i): v = M[:, i].toarray().ravel(); return v
edges = set(); why = {}
for mI in T:
    h1 = col(H1, mI); dE = (E0[mI] - E0) * H2CM; dE[mI] = np.inf
    mart = np.where(np.abs(h1) > 0, (np.abs(h1) * H2CM) ** 4 / np.abs(dE) ** 3, 0)
    for n in np.nonzero((mart >= 1.0) & (np.abs(dE) < 200))[0]:
        e = tuple(sorted((mI, int(n)))); edges.add(e); why[e] = 'Fermi (Martin %.1f cm-1)' % mart[n]
# second-order effective couplings (excluding Fermi partners from the sums)
fermi_of = {}
for a_, b_ in edges: fermi_of.setdefault(a_, set()).add(b_); fermi_of.setdefault(b_, set()).add(a_)
def c1vec(i, excl):
    dE = E0[i] - E0; dE[list(excl | {i})] = np.inf
    small = np.abs(dE) < 1e-10; dE[small] = np.inf
    return col(H1, i) / dE
for mI in T:
    dE = (E0[mI] - E0) * H2CM
    near = np.nonzero((np.abs(dE) < 200) & (np.arange(NS) != mI))[0]
    if len(near) == 0: continue
    cm_ = c1vec(mI, fermi_of.get(mI, set()))
    Hc = H1r @ cm_
    for n in near:
        cn = c1vec(int(n), fermi_of.get(int(n), set()))
        W = col(H2, int(n))[mI] + 0.5 * (Hc[n] + (H1r @ cn)[mI])
        if abs(W) * H2CM >= 0.1 * abs(dE[n]) and abs(W) * H2CM > 0.5:
            e = tuple(sorted((mI, int(n))))
            if e not in edges: edges.add(e); why[e] = 'effective 2nd-order |W|=%.2f cm-1, dE=%.1f' % (abs(W) * H2CM, dE[n])
auto = None
# union-find
par = {}
def fnd(x):
    par.setdefault(x, x)
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for a_, b_ in edges: par[fnd(a_)] = fnd(b_)
groups = {}
for x in set(T) | set(par): groups.setdefault(fnd(x), set()).add(x)
polyads = [sorted(g) for g in groups.values()]
# ---- solve ----
def solve(Pset):
    Pl = list(Pset); Pm = np.zeros(NS, bool); Pm[Pl] = True
    c1s = []; c2s = []
    for a_ in Pl:
        dE = E0[a_] - E0; dE[Pm] = np.inf
        bad = (np.abs(dE) < 1e-9) & (np.abs(col(H1, a_)) > 1e-12)
        assert not bad.any(), 'exact degeneracy with coupling'
        dE[np.abs(dE) < 1e-9] = np.inf
        c1 = col(H1, a_) / dE
        c2 = (col(H2, a_) + H1r @ c1) / dE; c2[Pm] = 0; c2[a_] = -0.5 * c1 @ c1
        c1s.append(c1); c2s.append(c2)
    n = len(Pl); Heff = np.zeros((n, n))
    for i, a_ in enumerate(Pl):
        for j, b_ in enumerate(Pl):
            h = (E0[a_] if i == j else 0) + H1[a_, b_] + H2[a_, b_]
            h += 0.5 * (col(H1, b_) @ c1s[i] + col(H1, a_) @ c1s[j])     # sum over Q of H1_an H1_nb/(E_a-E_n) etc.
            Heff[i, j] = h
    Heff = (Heff + Heff.T) / 2; ev, Uv = np.linalg.eigh(Heff)
    out = []
    for k in range(n):
        c0 = np.zeros(NS); c0[Pl] = Uv[:, k]
        out.append(dict(E=ev[k] + U, c=(c0, sum(Uv[i, k] * c1s[i] for i in range(n)), sum(Uv[i, k] * c2s[i] for i in range(n))),
                        comp={lab(S[Pl[i]]): float(Uv[i, k] ** 2) for i in range(n) if Uv[i, k] ** 2 > 0.02},
                        dom=lab(S[Pl[int(np.argmax(Uv[:, k] ** 2))]]), P=[lab(S[x]) for x in Pl]))
    return out
res = []
for Pl in polyads: res += solve(Pl)
g0 = [r for r in res if r['dom'] == '000000' and len(r['P']) == 1][0]; ZPE = g0['E']
# dipole operators: order k = monomials of degree k+1
def dip(comp, deg):
    M = sp.csr_matrix((NS, NS))
    for al, v in cM.items():
        if sum(al) == deg and abs(v[comp]) > 1e-14: M = M + v[comp] * mono(al)
    return M
Mfull = [[dip(c, dg) for dg in (1, 2, 3)] for c in range(3)]
def trans(cA, cB, Ms):
    t = 0.0
    for a_, b_, k in itertools.product(range(3), repeat=3):
        if a_ + b_ + k <= 2 and Ms[k] is not None: t += cA[a_] @ (Ms[k] @ cB[b_])
    return t
from scipy import constants as Kc
ea0 = Kc.e * Kc.physical_constants['Bohr radius'][0]
def A_km(nu, mu2): return 8 * np.pi ** 3 * Kc.N_A * nu * 100 * mu2 * ea0 ** 2 / (3 * Kc.h * Kc.c * 4 * np.pi * Kc.epsilon_0) / 1000
rows = []
for r in sorted(res, key=lambda r: r['E']):
    if r is g0: continue
    nu = (r['E'] - ZPE) * H2CM
    if nu > 1750: continue
    mf = [trans(g0['c'], r['c'], Mfull[c]) for c in range(3)]
    ml = [trans(g0['c'], r['c'], [Mfull[c][0], None, None]) for c in range(3)]
    rows.append(dict(label=r['dom'], E=nu, A=A_km(nu, sum(x * x for x in mf)), A_lin=A_km(nu, sum(x * x for x in ml)),
                     comp=r['comp'], polyad=r['P'] if len(r['P']) > 1 else None))
harm = [dict(mode=k, omega=w[k] * H2CM, A=A_km(w[k] * H2CM, sum(cM[tuple(1 * (np.arange(6) == k))] ** 2) / 2)) for k in range(NQ)]
json.dump(dict(rows=rows, harm=harm, ZPE=ZPE * H2CM, polyads=[[lab(S[x]) for x in Pl] for Pl in polyads if len(Pl) > 1],
               why={f'{lab(S[a_])}-{lab(S[b_])}': v for (a_, b_), v in why.items()}, nu1_comb_auto=auto, info=info,
               omega=list(w * H2CM)), open(f'vpt2_{TAG}.json', 'w'), indent=1)
print('ZPE (VPT2) %.2f cm-1' % (ZPE * H2CM))
print('resonance edges:'); [print('  ', lab(S[a_]), lab(S[b_]), v) for (a_, b_), v in why.items()]
print('nu1/nu2+nu5+nu6 detected automatically:', auto)
for r in rows:
    if r['A'] > 0.05 or r['label'].count('0') >= 5:
        print('%s %9.2f  A=%9.3f  A_lin=%9.3f  %s' % (r['label'], r['E'], r['A'], r['A_lin'], r['comp'] if r['polyad'] else ''))
