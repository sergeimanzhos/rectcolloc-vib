"""Benchmark of rectangular collocation with a GAUSSIAN basis for free CO2 on the cylindrically symmetric quartic force field
(vff_cyl.json, same potential as the VCI / GVPT2 of vpt2_cyl.py).  Energy rows only.  Exact Cartesian FD KEO.
Gaussians: products of 1D Gaussians in the two stretch coordinates x1, x3 (co2c.py) and, for the doubly degenerate bend, the
l = 0 projection of a 2D Gaussian in the bend plane:  exp(-(q^2 + c^2) / (2 b^2)) I0(q c / b^2),  q = sqrt(2 xb) the dimensionless
bend radius.  This is smooth at linearity (even in q).  J = 1 f block: the same functions times (a x b)_Z.
  python co2_ffbench.py holes | ref | grid N | series N ec b bb"""
import sys, json, itertools, numpy as np, numexpr as ne, scipy.linalg as sla
from scipy.special import ive
from multiprocessing import Pool
import co2c as c, co2ho as ho
H = c.H2CM; FF = json.load(open(f'{c.HERE}/vff_cyl.json')); _d = np.load(f'{c.HERE}/ff_co2.npz')
X0 = _d['X']; mA = _d['m'] * c.AMU; Lm = _d['L']; om = _d['omega']; sq = np.sqrt(np.repeat(mA, 3)); X0 = X0 - (mA[:, None] * X0).sum(0) / mA.sum()
def q_rect(X):
    """rectilinear dimensionless normal coordinates (q1, q3, qbend^2) of geometries X (n,3,3) in the Eckart frame of the linear reference"""
    x = X - (mA[None, :, None] * X).sum(1, keepdims=True) / mA.sum(); Hc = np.einsum('i,nia,ib->nab', mA, x, X0); U, S, Vt = np.linalg.svd(Hc)
    dd = np.sign(np.linalg.det(np.einsum('nba,ncb->nac', Vt, U))); D = np.ones((len(X), 3)); D[:, 2] = dd; R = np.einsum('nba,nb,ncb->nac', Vt, D, U)
    dq = ((np.einsum('nab,nib->nia', R, x) - X0[None]).reshape(len(X), 9) * sq[None]) @ Lm * np.sqrt(om)[None]; return dq[:, 2], dq[:, 3], dq[:, 0] ** 2 + dq[:, 1] ** 2
def V_ff(X):
    q1, q3, b2 = q_rect(X); w = FF['w']
    return (0.5 * w[0] * b2 + 0.5 * w[1] * q1 ** 2 + 0.5 * w[2] * q3 ** 2 + FF['k111'] * q1 ** 3 / 6 + FF['k133'] * q1 * q3 ** 2 / 2 + FF['k1bb'] * q1 * b2 / 2 + FF['p1111'] * q1 ** 4 / 24
            + FF['p3333'] * q3 ** 4 / 24 + FF['p1133'] * q1 ** 2 * q3 ** 2 / 4 + FF['p11bb'] * q1 ** 2 * b2 / 4 + FF['p33bb'] * q3 ** 2 * b2 / 4 + FF['pbbbb'] * b2 ** 2 / 24)
VBOX = 1.4 * 8000 / H; LO, HI = c.box(VBOX)
def points(start, n, vmax):
    out, i, tot = [], start, 0
    while tot < n:
        x = LO + (HI - LO) * c.r2seq(np.arange(i, i + 100000)); i += 100000; v = V_ff(c.geom(x)); k = (v <= vmax / H) & (v > -1e-6) & (c.V_harm(x) <= 4 * vmax / H); out.append((x[k], v[k])); tot += k.sum()
    return np.concatenate([o[0] for o in out])[:n], np.concatenate([o[1] for o in out])[:n]
class GL0:
    def __init__(s, C, beta, J=0):
        s.J = J; s.b = np.asarray(beta, float); s.c1 = C[:, 0] / s.b[0]; s.c3 = C[:, 1] / s.b[1]; s.cq = np.sqrt(2 * np.clip(C[:, 2], 0, None)) / s.b[2]; s.N = len(C)
    def __call__(s, X):
        x, axb = c.internal(X); u1 = (x[:, 0] / s.b[0])[:, None]; u3 = (x[:, 1] / s.b[1])[:, None]; uq = (np.sqrt(2 * np.clip(x[:, 2], 0, None)) / s.b[2])[:, None]
        c1, c3, cq = s.c1[None], s.c3[None], s.cq[None]; B = ne.evaluate('exp(-0.5*((u1-c1)**2+(u3-c3)**2+(uq-cq)**2))') * ive(0, uq * cq)
        return B if s.J == 0 else axb[:, None] * B
def solve(bas, xf, Ef, xv, Ev, nev, thr=1e-12, batch=4000):
    N = bas.N; S = np.zeros((N, N)); Hm = np.zeros((N, N))
    for a in range(0, len(xf), batch):
        X = c.geom(xf[a:a + batch]); B = bas(X); A = c.T_basis(bas, X) + Ef[a:a + batch, None] * B; S += B.T @ B; Hm += B.T @ A
    ev, U = np.linalg.eigh(S); k = ev > thr * ev.max(); Xo = U[:, k] / np.sqrt(ev[k]); E, W = sla.eig(Xo.T @ Hm @ Xo); o = np.argsort(E.real)[:nev]; E = E[o]; Cv = Xo @ W[:, o]
    Xv = c.geom(xv); Bv = bas(Xv); Av = c.T_basis(bas, Xv) + Ev[:, None] * Bv; P = Bv @ Cv; R = Av @ Cv - P * E[None]
    return (E.real * H), float(np.abs(E.imag).max() * H), np.linalg.norm(R, axis=0) / np.linalg.norm(P, axis=0) * H, int(k.sum())
_P = {}
def pts(vmax, M):
    if (vmax, M) not in _P: _P[(vmax, M)] = points(1, M, vmax) + points(5_000_000, 2000, min(vmax, 6000.0))
    return _P[(vmax, M)]
def centres(N, ec):
    return points(20_000_000, N, ec)[0]
def gauss(N, ec, beta, vmax=8000.0, MN=12, thr=1e-12):
    xf, Ef, xv, Ev = pts(vmax, MN * N); C = centres(N, ec)
    E0, im0, r0, k0 = solve(GL0(C, beta, 0), xf, Ef, xv, Ev, 8, thr); E1, im1, r1, k1 = solve(GL0(C, beta, 1), xf, Ef, xv, Ev, 3, thr)
    return dict(N=N, ec=ec, beta=list(beta), vmax=vmax, M=len(xf), thr=thr, zpe=float(E0[0]), nu2=float(E1[0] - E0[0] - c.BE * H), J0=(E0[1:] - E0[0]).tolist(), res0=r0.tolist(), resF=r1.tolist(), kept=(k0, k1), im=max(im0, im1),
                J=float(np.median(list(r0[:4]) + list(r1[:1]))))
def horef(ecut, vmax=8000.0, M=40000):
    xf, Ef, xv, Ev = pts(vmax, M); b0 = ho.HO(ecut, 0); b1 = ho.HO(ecut, 1); E0, _, r0, _ = solve(b0, xf, Ef, xv, Ev, 8, 1e-13); E1, _, r1, _ = solve(b1, xf, Ef, xv, Ev, 3, 1e-13)
    return dict(ecut=ecut, N=(b0.N, b1.N), zpe=float(E0[0]), nu2=float(E1[0] - E0[0] - c.BE * H), J0=(E0[1:] - E0[0]).tolist(), res0=r0.tolist(), resF=r1.tolist())
REF = None
def line(r):
    s = 'N=%-5d Ec=%-5g beta=%s Vmax=%-5g M=%-6d thr=%.0e kept=%s | ZPE %.3f nu2 %.3f dyad %.3f %.3f nu3 %.3f | res %s ; f %s | J %.2f' % (r['N'], r['ec'], np.round(r['beta'], 2).tolist(), r['vmax'], r['M'], r['thr'], r['kept'],
        r['zpe'], r['nu2'], r['J0'][0], r['J0'][1], r['J0'][2], ' '.join('%.2f' % x for x in r['res0'][:4]), ' '.join('%.2f' % x for x in r['resF'][:1]), r['J'])
    if REF: s += ' | err vs ref: ZPE %+.3f nu2 %+.3f dyad %+.3f %+.3f nu3 %+.3f' % (r['zpe'] - REF['zpe'], r['nu2'] - REF['nu2'], r['J0'][0] - REF['J0'][0], r['J0'][1] - REF['J0'][1], r['J0'][2] - REF['J0'][2])
    return s
def _g(a): return gauss(*a)
if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'holes':
        x = LO + (HI - LO) * c.r2seq(np.arange(1, 400001)); v = V_ff(c.geom(x)) * H; vh = c.V_harm(x) * H
        print('box lo', np.round(LO, 2), 'hi', np.round(HI, 2)); print('V_ff over the box: min %.0f cm-1; fraction below 8000: %.3f; fraction below 0: %.4f' % (v.min(), (v <= 8000).mean(), (v < 0).mean()))
        for lim in (8000, 6500, 5000): m = v <= lim; print('  points with V_ff <= %d: harmonic V up to %.0f (99th percentile %.0f); x ranges' % (lim, vh[m].max(), np.percentile(vh[m], 99)), np.round(x[m].min(0), 2), np.round(x[m].max(0), 2))
        xx = np.zeros((9, 3)); xx[:, 0] = np.linspace(LO[0], HI[0], 9); print('  cut along x1:', np.round(xx[:, 0], 1), np.round(V_ff(c.geom(xx)) * H))
        xx = np.zeros((9, 3)); xx[:, 2] = np.linspace(0, HI[2], 9); print('  cut along xb:', np.round(xx[:, 2], 1), np.round(V_ff(c.geom(xx)) * H))
    elif mode == 'ref':
        out = []
        for ecut in (13000, 17000, 21000, 25000): r = horef(ecut); out.append(r); print('HO reference Ecut=%d N=%s | ZPE %.3f nu2 %.3f dyad %.3f %.3f nu3 %.3f | next J=0: %s | res %s ; f %s' % (ecut, r['N'], r['zpe'], r['nu2'], r['J0'][0], r['J0'][1], r['J0'][2],
            ' '.join('%.2f' % x for x in r['J0'][3:6]), ' '.join('%.3f' % x for x in r['res0'][:4]), ' '.join('%.3f' % x for x in r['resF'][:1])), flush=True)
        json.dump(out, open(f'{c.HERE}/ffbench_ref.json', 'w'), indent=1)
    else:
        REF = json.load(open(f'{c.HERE}/ffbench_ref.json'))[-1]
        if mode == 'grid':
            N = int(sys.argv[2]); grid = [(N, ec, (b, b, bb)) for ec in (2500.0, 4000.0, 6000.0) for b in (0.6, 0.8, 1.0) for bb in (0.6, 0.8, 1.0, 1.3)]
            with Pool(12) as pool: R = pool.map(_g, grid, chunksize=1)
            for r in sorted(R, key=lambda r: r['J']): print(line(r), flush=True)
            json.dump(R, open(f'{c.HERE}/ffbench_grid_N{N}.json', 'w'), indent=1)
        elif mode == 'series':
            N, ec, b, bb = int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]); jobs = [(n, ec, (b, b, bb)) for n in (300, 600, 1200, 2400)] + [(N, ec, (b, b, bb), 8000.0, mn) for mn in (3, 6, 20)] + \
                [(N, ec, (b, b, bb), v) for v in (6500.0, 7250.0)] + [(N, ec, (b, b, bb), 8000.0, 12, t) for t in (1e-10, 1e-14)]
            with Pool(11) as pool: R = pool.map(_g, jobs, chunksize=1)
            for r in R: print(line(r), flush=True)
            json.dump(R, open(f'{c.HERE}/ffbench_series.json', 'w'), indent=1)
