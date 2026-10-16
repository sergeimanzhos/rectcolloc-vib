"""Rectangular collocation for free CO2 on B3LYP-D3(BJ)/aug-cc-pVTZ points (energy rows only), exact Cartesian FD KEO.
J = 0 block (l = 0 states: ground state, Fermi dyad nu1 / 2 nu2, nu3, ...) and J = 1 f-parity block (lowest state = nu2, l = 1;
band origin = E - E0 - B, since the rotational energy is B [J(J+1) - l^2] = B).
  python co2_colloc.py check                      mandatory KEO checks
  python co2_colloc.py run N ec b1 b3 bb vmax M   one calculation (prints levels)
  python co2_colloc.py scan                       parameter scans + convergence, writes co2_colloc.json"""
import sys, json, itertools, numpy as np, scipy.linalg as sla
import co2c as c
H = c.H2CM
def load(name, vmax):
    d = np.load(f'{c.HERE}/pts_{name}.npz'); o = np.argsort(d['idx']); m = d['conv'][o] & np.isfinite(d['E'][o]) & (d['E'][o] <= vmax / H); return d['x'][o][m], d['E'][o][m]
EX = np.array([e for dg in range(7) for e in itertools.product(range(dg + 1), repeat=3) if sum(e) == dg])
def Phi(x): return np.prod(x[:, None, :] ** EX[None], axis=2)
_S = {}
def surrogate():
    if 'f' not in _S:
        d = np.load(f'{c.HERE}/pts_fit.npz'); ok = d['conv'] & np.isfinite(d['E']) & (d['E'] < 2.5 * 8000 / H); coef = np.linalg.lstsq(Phi(d['x'][ok]), d['E'][ok], rcond=None)[0]
        _S['f'] = lambda q: np.where(c.V_harm(q) <= 3 * 8000 / H, Phi(q) @ coef, np.inf); _S['rms'] = float(np.sqrt(np.mean((Phi(d['x'][ok]) @ coef - d['E'][ok]) ** 2)) * H)
    return _S['f']
_C = {}
def centres(N, ec):
    if (N, ec) not in _C:
        V = surrogate(); lo, hi = c.box(1.4 * 8000 / H); out, i, n = [], 20_000_000, 0
        while n < N:
            q = lo + (hi - lo) * c.r2seq(np.arange(i, i + 100000)); i += 100000; out.append(q[V(q) <= ec / H]); n += len(out[-1])
        _C[(N, ec)] = np.concatenate(out)[:N]
    return _C[(N, ec)]
def solve(J, N, ec, beta, vmax, M, nev, thr=1e-12, vval=None):
    xf, Ef = load('fit', vmax); xf, Ef = xf[:M], Ef[:M]; xv, Ev = load('val', vval or vmax)
    bas = c.Basis(centres(N, ec), beta, J); X = c.geom(xf); B = bas(X); A = c.T_basis(bas, X) + Ef[:, None] * B
    S = B.T @ B; Hm = B.T @ A; ev, U = np.linalg.eigh(S); k = ev > thr * ev.max(); Xo = U[:, k] / np.sqrt(ev[k]); E, W = sla.eig(Xo.T @ Hm @ Xo); o = np.argsort(E.real)[:nev]; E = E[o]; Cv = Xo @ W[:, o]
    Xv = c.geom(xv); Bv = bas(Xv); Av = c.T_basis(bas, Xv) + Ev[:, None] * Bv; P = Bv @ Cv; R = Av @ Cv - P * E[None]
    return dict(E=(E.real * H).tolist(), im=float(np.abs(E.imag).max() * H), res=(np.linalg.norm(R, axis=0) / np.linalg.norm(P, axis=0) * H).tolist(), kept=int(k.sum()), cond=float(ev.max() / max(ev[k].min(), 1e-300)), M=len(xf), nval=len(xv))
def both(N, ec, beta, vmax, M, vval=None):
    a = solve(0, N, ec, beta, vmax, M, 10, vval=vval); f = solve(1, N, ec, beta, vmax, M, 4, vval=vval); E0 = a['E'][0]
    return dict(N=N, ec=ec, beta=list(beta), vmax=vmax, M=a['M'], zpe=E0, J0=[e - E0 for e in a['E']], F=[e - E0 - c.BE * H for e in f['E']], res0=a['res'], resF=f['res'], kept=(a['kept'], f['kept']), im=max(a['im'], f['im']),
                J=float(np.median(a['res'][:6] + f['res'][:2])))
def line(r): return 'N=%-5d Ec=%-5g beta=%s Vmax=%-5g M=%-5d | ZPE %.2f | nu2 %.2f | J=0: %s | res med %.2f (J0 %s; f %s) | kept %s | maxIm %.1e' % (
    r['N'], r['ec'], np.round(r['beta'], 2).tolist(), r['vmax'], r['M'], r['zpe'], r['F'][0], ' '.join('%.2f' % x for x in r['J0'][1:7]), r['J'], ' '.join('%.1f' % x for x in r['res0'][:6]), ' '.join('%.1f' % x for x in r['resF'][:2]), r['kept'], r['im'])
if __name__ == '__main__':
    mode = sys.argv[1]
    print('re = %.5f A; harmonic: bend %.2f, sym. stretch %.2f, asym. stretch %.2f cm-1; Be = %.5f cm-1; kappa = %.2f' % (c.RE * c.BOHR, c.W2 * H, c.W1 * H, c.W3 * H, c.BE * H, c.KAP))
    if mode == 'check':
        rng = np.random.default_rng(1); x = rng.normal(size=(40, 3)) * [1.0, 1.0, 0.0] + [0, 0, 1]; x[:, 2] = rng.uniform(0.1, 6, 40); X = c.geom(x)
        Rr = np.linalg.qr(rng.normal(size=(3, 3)))[0]; Xr = X @ Rr.T + rng.normal(size=3); print('invariance of x under rotation + translation: %.1e' % np.abs(c.internal(Xr)[0] - c.internal(X)[0]).max(), '; round trip x -> X -> x: %.1e' % np.abs(c.internal(X)[0] - x).max())
        class F:                                   # function of one bond length: T f = -(f'' + 2 f'/r) / (2 mu)
            def __call__(s, X): r = np.linalg.norm(X[:, 1] - X[:, 0], axis=1); return np.exp(-(r - c.RE) ** 2)[:, None]
        r = np.linalg.norm(X[:, 1] - X[:, 0], axis=1); g = np.exp(-(r - c.RE) ** 2); mu = c.mC * c.mO / (c.mC + c.mO); ana = -((4 * (r - c.RE) ** 2 - 2) * g + 2 * (-2 * (r - c.RE) * g) / r) / (2 * mu)
        num = c.T_basis(F(), X)[:, 0]; print('exact-KEO check on f(r_CO): max relative error %.1e' % (np.abs(num - ana).max() / np.abs(ana).max()))
        bas = c.Basis(x[:5], [0.9, 0.9, 0.9]); t1 = c.T_basis(bas, X, 1e-3); t2 = c.T_basis(bas, X, 2e-3); print('h-sensitivity of T phi (h vs 2h), J = 0 basis: %.1e relative' % (np.abs(t1 - t2).max() / np.abs(t1).max()))
        bas = c.Basis(x[:5], [0.9, 0.9, 0.9], 1); t1 = c.T_basis(bas, X, 1e-3); t2 = c.T_basis(bas, X, 2e-3); print('h-sensitivity, J = 1 f basis: %.1e relative' % (np.abs(t1 - t2).max() / np.abs(t1).max()))
        xl = x.copy(); xl[:, 2] = 1e-6; Xl = c.geom(xl); bas = c.Basis(x[:5], [0.9, 0.9, 0.9]); t1 = c.T_basis(bas, Xl, 1e-3); t2 = c.T_basis(bas, Xl, 2e-3); print('h-sensitivity at linear geometries, J = 0 basis: %.1e relative' % (np.abs(t1 - t2).max() / np.abs(t1).max()))
    elif mode == 'run':
        N, ec, b1, b3, bb, vmax, M = int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), int(sys.argv[8]); print(line(both(N, ec, [b1, b3, bb], vmax, M)))
    elif mode == 'scan':
        xf, Ef = load('fit', 8000); xv, _ = load('val', 8000); surrogate(); print('fit points V <= 8000: %d; validation: %d; surrogate rms %.1f cm-1' % (len(xf), len(xv), _S['rms']))
        print('in-region points near the box faces (fraction within 3%% of lower / upper face):', np.round(((xf - c.box(1.4 * 8000 / H)[0]) / (c.box(1.4 * 8000 / H)[1] - c.box(1.4 * 8000 / H)[0]) < 0.03).mean(0), 3), np.round(((xf - c.box(1.4 * 8000 / H)[0]) / (c.box(1.4 * 8000 / H)[1] - c.box(1.4 * 8000 / H)[0]) > 0.97).mean(0), 3))
        MALL = len(xf); out = []; P = dict(N=800, ec=3500.0, beta=[0.9, 0.9, 0.9], vmax=8000.0)
        def run(tag, **kw):
            q = dict(P, **kw); r = both(q['N'], q['ec'], q['beta'], q['vmax'], kw.get('M', MALL), vval=6000); r['tag'] = tag; out.append(r); print('%-22s' % tag, line(r), flush=True); return r
        best = run('start')
        def greedy(name, values, setter):          # keep a change only if it lowers the median residual by more than 2 %
            global P
            cur = best
            for v in values:
                kw = setter(v); r = run(f'{name}={v:g}', **kw)
                if r['J'] < cur['J'] * 0.98: cur = r; P = dict(P, **kw)
            return cur
        best = greedy('Ec', (2000.0, 2750.0, 4500.0, 6000.0), lambda v: dict(ec=v))
        b0 = list(P['beta']); best = greedy('beta_all', (0.7, 0.8, 1.0, 1.1, 1.25), lambda v: dict(beta=[v, v, P['beta'][2]]))
        best = greedy('beta_bend', (0.5, 0.7, 1.1, 1.4, 1.8), lambda v: dict(beta=[P['beta'][0], P['beta'][1], v]))
        best = greedy('Ec2', (0.75 * P['ec'], 1.3 * P['ec']), lambda v: dict(ec=v))
        print('--- tuned parameters:', P)
        conv = [run(f'M={m}', M=m) for m in (2 * P['N'], 4 * P['N'], 6 * P['N'], MALL) if m <= MALL]
        nser = [run(f'N={n}', N=n) for n in (400, 600, 1200, 1600)]
        vser = [run(f'Vmax={v}', vmax=v) for v in (6500.0, 7250.0)]
        final = run('final'); json.dump(dict(params=P, final=final, runs=out, harm=dict(bend=c.W2 * H, sym=c.W1 * H, asym=c.W3 * H), Be=c.BE * H, re=c.RE * c.BOHR), open(f'{c.HERE}/co2_colloc.json', 'w'), indent=1)
        print('saved co2_colloc.json')
