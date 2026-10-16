"""IR intensities for the LiCO2 collocation states (skill section 8b).
  python intens.py <cache> <vec npz> <out tag> [vmax]
Routes: (a) quasi-Monte-Carlo over the DFT points (fit + val) with weights rho(q);
        (b) degree-6 polynomial dipole surface integrated on a dense R2 grid in the same region, weights rho(q).
Also: symmetry labels (C2v, from the sign of psi under the two reflections expressed in q), harmonic-oscillator assignment
with overlap weight, double-harmonic intensities for comparison. Output: results/intens_<tag>.csv / .json."""
import sys, os, json, math, itertools, numpy as np
import coords as c, collo as co
from r2 import r2
H = c.H2CM; NQ = c.NQ; NA = c.NA
cache, vecfn, tag = sys.argv[1], sys.argv[2], sys.argv[3]; VMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 5500.0
AU2D = 2.541746473; EA0 = 8.4783536255e-30; NAV = 6.02214076e23; HPL = 6.62607015e-34; CL = 2.99792458e8; EPS0 = 8.8541878128e-12
KCONST = 8 * np.pi ** 3 * NAV / (3 * HPL * CL * 4 * np.pi * EPS0) / 1000.0     # A[km/mol] = KCONST * nu[m^-1] * mu^2[C^2 m^2]
def A_kmmol(nu_cm, mu_au): return KCONST * nu_cm * 100.0 * (mu_au * EA0) ** 2
# ---------- volume element rho(q) = |det[dX/dt, dX/dphi, dX/dq]| in the Eckart frame
def rho(q):
    X = c.X_of_q(q); n = len(X); J = np.zeros((n, 3 * NA, 3 * NA))
    for a in range(3):
        t = np.zeros((NA, 3)); t[:, a] = 1; J[:, :, a] = t.ravel()[None]
        e = np.zeros(3); e[a] = 1; J[:, :, 3 + a] = np.cross(e[None, None, :], X).reshape(n, -1)
    J[:, :, 6:] = c.A[None]
    return np.abs(np.linalg.det(J))
# ---------- polynomial tools (degree D in 6 coordinates)
def expo(D): return np.array([e for d in range(D + 1) for e in itertools.product(range(d + 1), repeat=NQ) if sum(e) == d])
def Phi(Q, EX):
    out = np.ones((len(Q), len(EX))); D = EX.max()
    for k in range(NQ):
        pw = np.stack([Q[:, k] ** n for n in range(D + 1)], 1); out *= pw[:, EX[:, k]]
    return out
# ---------- data
def load(n):
    d = np.load(os.path.join(cache, f'{n}.npz')); return d['q'], d['E'], d['dip'], d['conv']
qF, EF, dF, cF = load('fit'); qV, EV, dV, cV = load('val')
inF = cF & (EF <= VMAX / H); inV = cV & (EV <= VMAX / H)
qa = np.concatenate([qF[inF], qV[inV]]); da = np.concatenate([dF[inF], dV[inV]])
print(f'{tag}: DFT points in region: fit {inF.sum()}, val {inV.sum()}', flush=True)
# ---------- symmetry operations in q (C2v): sigma_plane (molecular plane), sigma_perp (swaps the O atoms)
X0 = c.X0 - (c.m[:, None] * c.X0).sum(0) / c.m.sum()
nrm = np.linalg.svd(X0)[2][2]                                  # plane normal
ax = X0[3] - X0[0]; ax -= nrm * (ax @ nrm); ax /= np.linalg.norm(ax)   # C2 axis (C -> Li)
u = np.cross(nrm, ax)                                          # in-plane, perpendicular to the axis
def op_matrix(v, perm):
    R = np.eye(3) - 2 * np.outer(v, v); P = np.zeros((3 * NA, 3 * NA))
    for i, j in enumerate(perm): P[3 * j:3 * j + 3, 3 * i:3 * i + 3] = R
    return c.Q @ P @ c.A
Sp = op_matrix(nrm, [0, 1, 2, 3]); Sv = op_matrix(u, [0, 2, 1, 3])
sp, sv = np.round(np.diag(Sp)).astype(int), np.round(np.diag(Sv)).astype(int)
dev = max(np.abs(Sp - np.diag(sp)).max(), np.abs(Sv - np.diag(sv)).max()); print('symmetry operations in q: max deviation from diagonal +-1: %.1e' % dev)
if dev > 2e-2: print(np.round(Sp, 3)); print(np.round(Sv, 3)); raise SystemExit('modes are not symmetry-adapted')
IRR = {(1, 1): 'A1', (1, -1): 'B2', (-1, 1): 'B1', (-1, -1): 'A2'}
print('mode symmetries:', [IRR[(a, b)] for a, b in zip(sp, sv)], ' omega', np.round(c.omega * H, 1), flush=True)
# ---------- states
d = np.load(vecfn); bas = co.Basis(d['centres'], d['beta']); Cc = d['C']; E = d['E'].real * H; nst = Cc.shape[1]
# real wavefunctions; a complex-conjugate pair of eigenvalues is replaced by the real and imaginary parts of its vector (as cl2co/ir6.py)
Ecm = d['E'] * H; Cr = np.zeros(Cc.shape); used = set(); PAIRS = []
for f in range(nst):
    if f in used: continue
    cand = [j for j in range(nst) if j != f and j not in used and abs(Ecm[j] - np.conj(Ecm[f])) < 1e-4]
    if abs(Ecm[f].imag) > 1e-6 and cand:
        g = cand[0]; Cr[:, f] = Cc[:, f].real; Cr[:, g] = Cc[:, f].imag; used |= {f, g}; PAIRS.append((f, g))
    else:
        ph = np.exp(-0.5j * np.angle((Cc[:, f] ** 2).sum())); Cr[:, f] = (Cc[:, f] * ph).real; used.add(f)
PAIRSET = {x for p_ in PAIRS for x in p_}
print('complex pairs (state indices):', PAIRS, ' max|Im E| %.2f cm-1' % np.abs(Ecm.imag).max())
def psi(q, batch=1500): return np.concatenate([bas(q[a:a + batch]) @ Cr for a in range(0, len(q), batch)])
def moments(P, w, mu):
    """rho-weighted norms, overlaps with state 0, transition dipoles <0|mu-<mu>_0|f>"""
    Nf = (w[:, None] * P ** 2).sum(0); P = P / np.sqrt(Nf)[None]
    S0 = (w[:, None] * P[:, :1] * P).sum(0); Sfull = (w[:, None] * P).T @ P
    ev, U = np.linalg.eigh(Sfull); Po = P @ (U @ np.diag(ev ** -0.5) @ U.T)      # Loewdin orthonormalisation in the rho metric (as cl2co/ir6.py)
    mu0 = (w[:, None] * Po[:, :1] ** 2 * mu).sum(0)
    T = np.einsum('i,i,ia,if->fa', w, Po[:, 0], mu - mu0[None], Po)
    return T, S0, mu0, Sfull, P
# route (a): the DFT points themselves
Pa = psi(qa); wa = rho(qa); Ta, S0a, mu0a, Sfa, _ = moments(Pa, wa, da)
# route (b): polynomial dipole surface + dense grid
EX6 = expo(6); okF = cF & (EF <= 1.0 * VMAX / H) & np.isfinite(dF).all(1)
coef = np.linalg.lstsq(Phi(qF[okF], EX6), dF[okF], rcond=None)[0]
errV = Phi(qV[inV], EX6) @ coef - dV[inV]
print('dipole fit (degree 6, %d terms, %d points): RMS error on validation points (a.u.) %s ; dipole range (a.u.) %s' % (len(EX6), okF.sum(), np.round(np.sqrt((errV ** 2).mean(0)), 5), np.round(np.ptp(dV[inV], axis=0), 3)), flush=True)
EX4 = expo(4); okE = cF & (EF < 2.5 * 5500 / H); cV4 = np.linalg.lstsq(Phi(qF[okE], EX4), EF[okE], rcond=None)[0]   # region surrogate (same as the centre surrogate)
b = np.load(os.path.join(cache, 'box.npz')); lo, hi = b['lo'], b['hi']
NG = int(os.environ.get('NGRID', 300000)); qs = []; i0 = 40_000_000; nb = 0
while sum(len(x) for x in qs) < NG:
    q = lo + (hi - lo) * r2(np.arange(i0, i0 + 200000), NQ); i0 += 200000; nb += 200000
    vff = c.V_ff(q); keep = (Phi(q, EX4) @ cV4 <= VMAX / H) & (vff <= 2.5 * VMAX / H) & (vff > -1e-3); qs.append(q[keep])
qg = np.concatenate(qs); VOL = np.prod(hi - lo) * len(qg) / nb
print(f'dense grid: {len(qg)} points (acceptance {len(qg)/nb:.4f})', flush=True)
Pg = psi(qg); wg = rho(qg); mug = np.concatenate([Phi(qg[a:a + 20000], EX6) @ coef for a in range(0, len(qg), 20000)])
Tb, S0b, mu0b, Sfb, Pgn = moments(Pg, wg, mug)
iu = np.triu_indices(nst, 1); offa = np.abs(Sfa - np.eye(nst)); offb = np.abs(Sfb - np.eye(nst))
print('state overlaps with rho: route a median %.1e max %.2f ; route b median %.1e max %.2f' % (np.median(offa[iu]), offa.max(), np.median(offb[iu]), offb.max()))
top = np.argsort(offb[iu])[::-1][:6]; print('   largest overlaps (route b):', [(int(iu[0][t]), int(iu[1][t]), round(float(offb[iu][t]), 3)) for t in top])
Pd = Pg / np.sqrt((Pg ** 2).sum(0))[None]; offd = np.abs(Pd.T @ Pd - np.eye(nst)); print('   (plain dq, for comparison: median %.1e max %.2f)' % (np.median(offd[iu]), offd.max()))
print('permanent dipole <0|mu|0> (a.u.): route a', np.round(mu0a, 4), ' route b', np.round(mu0b, 4), ' | along C2 axis (D): %.3f' % (mu0b @ ax * AU2D), flush=True)
# ---------- symmetry of the states and HO assignment (dense grid)
def parity(s):
    Ps = psi(qg * s[None]); return np.einsum('i,if,if->f', wg, Pgn, Ps) / np.sqrt((wg[:, None] * Ps ** 2).sum(0))
par_p, par_v = parity(sp), parity(sv)
sym = [IRR[(1 if a > 0 else -1, 1 if b_ > 0 else -1)] for a, b_ in zip(par_p, par_v)]
def ho(n, x):
    h = np.polynomial.hermite.hermval(x, [0] * n + [1]); return h * np.exp(-x ** 2 / 2) / math.sqrt(2 ** n * math.factorial(n) * math.sqrt(math.pi))
ELIM = max(2400.0, (E[-1] - E[0]) + 800.0)                        # HO labels: harmonic energy up to 800 cm-1 above the top state
ns = [n for n in itertools.product(range(9), repeat=NQ) if sum(n) <= 8 and np.dot(n, c.omega) * H <= ELIM]
H1 = [[ho(n, qg[:, k]) for n in range(9)] for k in range(NQ)]
Pq = Pg / np.sqrt((Pg ** 2).sum(0) * VOL / len(qg))[None]            # dq-normalised states
ov = np.array([(np.prod([H1[k][n[k]] for k in range(NQ)], axis=0)[:, None] * Pq).sum(0) * VOL / len(qg) for n in ns]) ** 2
lab = [''.join(map(str, ns[i])) for i in ov.argmax(0)]; wgt = ov.max(0)
# ---------- double-harmonic reference
lin = [int(np.nonzero((EX6 == np.eye(NQ, dtype=int)[k]).all(1))[0][0]) for k in range(NQ)]
dmu = coef[lin]                                                    # d mu / d q_k at q = 0 (a.u.)
Aharm = [A_kmmol(c.omega[k] * H, np.linalg.norm(dmu[k]) / math.sqrt(2)) for k in range(NQ)]
# check of the constant against the double-harmonic formula 42.256 km/mol per (D/A)^2/amu:
# <0|Q|1>^2 = hbar/(2 omega) for a mass-weighted coordinate, so A = KCONST * nu * (dmu/dQ)^2 * hbar/(2 * 2 pi c nu) (nu in m^-1)
chk = KCONST * (3.33564095e-30 / 1e-10) ** 2 / 1.66053907e-27 * (HPL / (2 * np.pi)) / (4 * np.pi * CL)
print('constant check: %.3f km/mol per (D/A)^2/amu (expected 42.256)' % chk)
print('double-harmonic: ', [(round(c.omega[k] * H, 1), round(Aharm[k], 2)) for k in range(NQ)])
# ---------- table
axes = np.array([ax, u, nrm]); rows = []
for f in range(1, nst):
    nu = E[f] - E[0]; ta, tb = axes @ Ta[f], axes @ Tb[f]
    rows.append(dict(level=f, nu=round(float(nu), 2), sym=sym[f], assign=lab[f], weight=round(float(wgt[f]), 3), I_a=float(A_kmmol(nu, np.linalg.norm(ta))), I_b=float(A_kmmol(nu, np.linalg.norm(tb))),
                     I_b_axis=float(A_kmmol(nu, abs(tb[0]))), I_b_inplane=float(A_kmmol(nu, abs(tb[1]))), I_b_oop=float(A_kmmol(nu, abs(tb[2]))),
                     ovl0=float(S0b[f]), res=float(d['res'][f]), parity=(round(float(par_p[f]), 2), round(float(par_v[f]), 2))))
json.dump(dict(rows=rows, zpe=float(E[0]), harm=[dict(nu=float(c.omega[k] * H), I=float(Aharm[k]), sym=IRR[(sp[k], sv[k])]) for k in range(NQ)], mu0_D=float(mu0b @ ax * AU2D)), open(f'results/intens_{tag}.json', 'w'), indent=1)
json.dump(dict(rows=[dict(label=r['assign'], E=r['nu'], A=r['I_b'], sym=r['sym'], w=r['weight']) for r in rows], ZPE=float(E[0])), open(f'results/ir_colloc_{tag}.json', 'w'), indent=1)
with open(f'results/intens_{tag}.csv', 'w') as fo:
    fo.write('level,nu_cm-1,sym,assignment(v1..v6),HO_weight,I_points_km/mol,I_fit_km/mol,residual_cm-1,overlap_with_ground,parity_plane,parity_swap,flag\n')
    for r in rows:
        flag = ('complex_pair ' if r['level'] in PAIRSET else '') + ('mixed ' if r['weight'] < 0.4 else '') + ('routes_differ' if max(r['I_a'], r['I_b']) > 0.1 and abs(r['I_a'] - r['I_b']) > 0.3 * max(r['I_a'], r['I_b']) else '')
        fo.write(f"{r['level']},{r['nu']:.2f},{r['sym']},{r['assign']},{r['weight']:.2f},{r['I_a']:.4f},{r['I_b']:.4f},{r['res']:.1f},{r['ovl0']:.1e},{r['parity'][0]},{r['parity'][1]},{flag}\n")
print(open(f'results/intens_{tag}.csv').read())
