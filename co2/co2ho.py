"""Free CO2, rectangular collocation with a pruned harmonic-oscillator product basis (skill section 4, "Harmonic-oscillator products")
in the coordinates of co2c.py:  phi(n1, n3, nb) = h_n1(x1 / s) h_n3(x3 / s) L_nb(2 xb / sb) exp(-xb / sb),
pruned by n1 w1 + n3 w3 + (2 nb + l) w2 <= Ecut, l = 0 for the J = 0 block and l = 1 for the J = 1 f block (which carries the
extra factor (a x b)_Z).  Exact Cartesian FD KEO; energy rows only; DFT energies at the points (pts_fit / pts_val).
  python co2ho.py scan"""
import sys, json, math, itertools, numpy as np, scipy.linalg as sla
import co2c as c
H = c.H2CM
def load(name, vmax):
    d = np.load(f'{c.HERE}/pts_{name}.npz'); o = np.argsort(d['idx']); m = d['conv'][o] & np.isfinite(d['E'][o]) & (d['E'][o] <= vmax / H); return d['x'][o][m], d['E'][o][m]
class HO:
    def __init__(s, ecut, J=0, sc=1.0, sb=1.0):
        s.J, s.sc, s.sb = J, sc, sb; e = ecut / H
        s.n = np.array([(n1, n3, nb) for n1 in range(40) for n3 in range(40) for nb in range(40) if n1 * c.W1 + n3 * c.W3 + (2 * nb + J) * c.W2 <= e]); s.N = len(s.n)
    def __call__(s, X):
        x, axb = c.internal(X); n1m, n3m, nbm = s.n.max(0)
        def herm(u, nmax):                                  # normalised Hermite functions by recurrence
            out = np.zeros((nmax + 1, len(u))); out[0] = np.pi ** -0.25 * np.exp(-u ** 2 / 2)
            if nmax > 0: out[1] = np.sqrt(2.0) * u * out[0]
            for k in range(2, nmax + 1): out[k] = np.sqrt(2.0 / k) * u * out[k - 1] - np.sqrt((k - 1) / k) * out[k - 2]
            return out
        def lag(t, nmax):                                   # Laguerre L_k(t) exp(-t / 2) by recurrence
            out = np.zeros((nmax + 1, len(t))); out[0] = np.exp(-t / 2)
            if nmax > 0: out[1] = (1 - t) * out[0]
            for k in range(1, nmax): out[k + 1] = ((2 * k + 1 - t) * out[k] - k * out[k - 1]) / (k + 1)
            return out
        h1 = herm(x[:, 0] / s.sc, n1m); h3 = herm(x[:, 1] / s.sc, n3m); lb = lag(2 * x[:, 2] / s.sb, nbm)
        B = (h1[s.n[:, 0]] * h3[s.n[:, 1]] * lb[s.n[:, 2]]).T
        return B if s.J == 0 else axb[:, None] * B
def solve(J, ecut, vmax, M, nev, sc=1.0, sb=1.0, thr=1e-13, vval=6000):
    xf, Ef = load('fit', vmax); xf, Ef = xf[:M], Ef[:M]; xv, Ev = load('val', vval); bas = HO(ecut, J, sc, sb)
    X = c.geom(xf); B = bas(X); A = c.T_basis(bas, X) + Ef[:, None] * B; S = B.T @ B; Hm = B.T @ A
    ev, U = np.linalg.eigh(S); k = ev > thr * ev.max(); Xo = U[:, k] / np.sqrt(ev[k]); E, W = sla.eig(Xo.T @ Hm @ Xo); o = np.argsort(E.real)[:nev]; E = E[o]; Cv = Xo @ W[:, o]
    Xv = c.geom(xv); Bv = bas(Xv); Av = c.T_basis(bas, Xv) + Ev[:, None] * Bv; P = Bv @ Cv; R = Av @ Cv - P * E[None]
    return dict(E=(E.real * H).tolist(), im=float(np.abs(E.imag).max() * H), res=(np.linalg.norm(R, axis=0) / np.linalg.norm(P, axis=0) * H).tolist(), N=bas.N, kept=int(k.sum()), M=len(xf), cond=float(ev.max() / ev[k].min()))
def both(ecut, vmax, M, sc=1.0, sb=1.0, tag=''):
    a = solve(0, ecut, vmax, M, 8, sc, sb); f = solve(1, ecut, vmax, M, 3, sc, sb); E0 = a['E'][0]
    r = dict(tag=tag, ecut=ecut, vmax=vmax, M=a['M'], sc=sc, sb=sb, N=(a['N'], f['N']), kept=(a['kept'], f['kept']), zpe=E0, J0=[e - E0 for e in a['E']], F=[e - E0 - c.BE * H for e in f['E']], res0=a['res'], resF=f['res'], im=max(a['im'], f['im']))
    print('%-16s Ecut=%-6g Vmax=%-5g M=%-5d s=%.2f sb=%.2f N=%s kept=%s | ZPE %.3f | nu2 %.3f | J=0: %s | res J0 %s ; f %s | maxIm %.1e' % (tag, ecut, vmax, r['M'], sc, sb, r['N'], r['kept'], E0, r['F'][0],
          ' '.join('%.3f' % x for x in r['J0'][1:6]), ' '.join('%.2f' % x for x in r['res0'][:5]), ' '.join('%.2f' % x for x in r['resF'][:2]), r['im']), flush=True)
    return r
if __name__ == '__main__':
    print('re = %.5f A; harmonic: bend %.2f, sym. stretch %.2f, asym. stretch %.2f; Be = %.5f cm-1' % (c.RE * c.BOHR, c.W2 * H, c.W1 * H, c.W3 * H, c.BE * H))
    MALL = len(load('fit', 8000)[0]); out = []
    for ecut in (5000, 7000, 9000, 11000, 13000): out.append(both(ecut, 8000, MALL, tag='Ecut'))
    for m in (1000, 2000, 4000, 6000): out.append(both(9000, 8000, m, tag='M'))
    for v in (6500, 7250): out.append(both(9000, v, MALL, tag='Vmax'))
    for s_ in (0.9, 1.1): out.append(both(9000, 8000, MALL, sc=s_, tag='stretch scale'))
    for s_ in (0.9, 1.1): out.append(both(9000, 8000, MALL, sb=s_, tag='bend scale'))
    json.dump(dict(runs=out, harm=dict(bend=c.W2 * H, sym=c.W1 * H, asym=c.W3 * H), Be=c.BE * H, re=c.RE * c.BOHR), open(f'{c.HERE}/co2_ho.json', 'w'), indent=1); print('saved co2_ho.json')
