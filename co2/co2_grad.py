"""Free CO2: rectangular collocation with potential-gradient (derivative) rows (skill section 5b), Gaussian basis of co2_ffbench.py.
For each point and each of the three coordinates u_k = (x1, x3, q), q = sqrt(2 xb) the bend radius, a CONSTANT Cartesian direction
a_k is used (mass-weighted right inverse of du/dX at that point, so a_k . grad u_j = delta_jk).  D_k = a_k . grad_X then commutes
exactly with the Cartesian Laplacian, and the extra rows are
    B_k = D_k phi_j ,   A_k = T(D_k phi_j) + V D_k phi_j + (D_k V) phi_j ,
with D_k phi and D_k (T phi) taken by 5-point finite differences along the straight line X + eps a_k (eps = 1e-2).
  python co2_grad.py check                 commutation / direction checks
  python co2_grad.py ff                    benchmark on the cylindrical quartic force field, energies vs energies + gradients
  python co2_grad.py dftpts <n> <nproc>    B3LYP-D3(BJ)/aug-cc-pVTZ gradients at the first n fit points
  python co2_grad.py dft                   collocation on the DFT points, energies vs energies + gradients"""
import os, sys, json, time, numpy as np, scipy.linalg as sla
from multiprocessing import Pool
import co2c as c, co2ho as ho, co2_ffbench as fb
H = c.H2CM; D1W = np.array([1, -8, 0, 8, -1]) / 12.0; EPS = 1e-2
def directions(X, h=1e-4):
    """a (n,3,9): constant Cartesian directions with a_k . grad u_j = delta_jk for u = (x1, x3, q)"""
    n = len(X); Xf = X.reshape(n, 9); J = np.zeros((n, 3, 9))
    for i in range(9):
        Xp = Xf.copy(); Xp[:, i] += h; Xm = Xf.copy(); Xm[:, i] -= h; J[:, :, i] = (c.internal(Xp.reshape(n, 3, 3))[0] - c.internal(Xm.reshape(n, 3, 3))[0]) / (2 * h)
    q = np.sqrt(2 * np.clip(c.internal(X)[0][:, 2], 1e-300, None)); J[:, 2, :] /= q[:, None]            # grad q = grad xb / q
    W = 1 / c.MM; JW = J * W[None, None, :]; G = np.einsum('nki,nli->nkl', JW, J); return np.einsum('nki,nkl->nli', JW, np.linalg.inv(G)), q
def dline(f, X, a):
    """d/d eps f(X + eps a) by 5-point central differences (f returns an array with leading dimension n)"""
    n = len(X); Xf = X.reshape(n, 9); out = 0
    for s, w in zip((-2, -1, 1, 2), (D1W[0], D1W[1], D1W[3], D1W[4])): out = out + w * f((Xf + s * EPS * a).reshape(n, 3, 3))
    return out / EPS
def accumulate(bas, x, V, GV=None, lam=1.0, batch=1500):
    """S = sum B^T B, H = sum B^T A over energy rows and (if GV is given, shape (n,3): dV/d(x1, x3, q)) gradient rows weighted by lam (scalar or one value per point)"""
    N = bas.N; S = np.zeros((N, N)); Hm = np.zeros((N, N))
    for b0 in range(0, len(x), batch):
        X = c.geom(x[b0:b0 + batch]); v = V[b0:b0 + batch]; B = bas(X); A = c.T_basis(bas, X) + v[:, None] * B; S += B.T @ B; Hm += B.T @ A
        if GV is not None:
            a, _ = directions(X)
            for k in range(3):
                Bk = dline(bas, X, a[:, k]); Ak = dline(lambda Y: c.T_basis(bas, Y), X, a[:, k]) + v[:, None] * Bk + GV[b0:b0 + batch, k][:, None] * B; Bw = Bk * (lam[b0:b0 + batch, None] if np.ndim(lam) else lam); S += Bw.T @ Bk; Hm += Bw.T @ Ak
    return S, Hm
def finish(bas, S, Hm, xv, Ev, nev, thr=1e-12):
    ev, U = np.linalg.eigh(S); k = ev > thr * ev.max(); Xo = U[:, k] / np.sqrt(ev[k]); E, W = sla.eig(Xo.T @ Hm @ Xo); o = np.argsort(E.real)[:nev]; E = E[o]; Cv = Xo @ W[:, o]
    Xv = c.geom(xv); Bv = bas(Xv); Av = c.T_basis(bas, Xv) + Ev[:, None] * Bv; P = Bv @ Cv; R = Av @ Cv - P * E[None]
    return E.real * H, float(np.abs(E.imag).max() * H), np.linalg.norm(R, axis=0) / np.linalg.norm(P, axis=0) * H, int(k.sum())
def levels(C, beta, x, V, GV, xv, Ev, lam=1.0, thr=1e-12):
    out = {}
    for J, nev in ((0, 8), (1, 3)):
        bas = fb.GL0(C, beta, J); S, Hm = accumulate(bas, x, V, GV, lam); out[J] = finish(bas, S, Hm, xv, Ev, nev, thr)
    E0, E1 = out[0][0], out[1][0]
    return dict(zpe=float(E0[0]), nu2=float(E1[0] - E0[0] - c.BE * H), J0=(E0[1:] - E0[0]).tolist(), res0=out[0][2].tolist(), resF=out[1][2].tolist(), kept=(out[0][3], out[1][3]), im=max(out[0][1], out[1][1]))
def gradV_ff(x):
    X = c.geom(x); a, _ = directions(X); return np.stack([dline(fb.V_ff, X, a[:, k]) for k in range(3)], 1)
def show(tag, r, ref=None):
    s = '%-34s kept=%s | ZPE %.3f nu2 %.3f dyad %.3f %.3f nu3 %.3f | res %s ; f %s | maxIm %.0e' % (tag, r['kept'], r['zpe'], r['nu2'], r['J0'][0], r['J0'][1], r['J0'][2], ' '.join('%.2f' % v for v in r['res0'][:4]), ' '.join('%.2f' % v for v in r['resF'][:1]), r['im'])
    if ref: e = [r['zpe'] - ref['zpe'], r['nu2'] - ref['nu2'], r['J0'][0] - ref['J0'][0], r['J0'][1] - ref['J0'][1], r['J0'][2] - ref['J0'][2]]; s += ' | err %s | max|err| %.3f' % (' '.join('%+.3f' % v for v in e), max(abs(v) for v in e[1:]))
    print(s, flush=True)
BETA = (0.9, 0.9, 0.8); EC = 4000.0
def job_ff(a):
    N, M, grad, lam = a; x, V = fb.points(1, M, 8000.0); xv, Ev = fb.points(5_000_000, 2000, 6000.0); C = fb.centres(N, EC)
    r = levels(C, BETA, x, V, gradV_ff(x) if grad else None, xv, Ev, lam); r.update(N=N, M=M, grad=grad, lam=lam); return r
def job_dft(a):
    N, M, grad, lam = a; d = np.load(f'{c.HERE}/pts_grad.npz'); x, V, GV = d['x'][:M], d['E'][:M], d['gu'][:M]; xv, Ev = ho.load('val', 6000); C = fb.centres(N, EC)
    r = levels(C, BETA, x, V, GV if grad else None, xv, Ev, lam); r.update(N=N, M=M, grad=grad, lam=lam); return r
def dft_one(a):
    i, x, Eref = a
    from pyscf import gto, dft, lib
    lib.num_threads(1); X = c.geom(x[None])[0]; mol = gto.M(atom=[[e, y] for e, y in zip(('C', 'O', 'O'), X)], basis='aug-cc-pvtz', unit='Bohr', symmetry=False, verbose=0)
    mf = dft.RKS(mol).set(xc='b3lyp', conv_tol=1e-10, max_cycle=150); mf.grids.level = 5; mf.disp = 'd3bj'; mf.kernel(dm0=np.load(f'{c.HERE}/dm0.npy')); g = mf.nuc_grad_method(); g.grid_response = True
    return i, mf.e_tot, g.kernel().ravel(), bool(mf.converged)
if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'check':
        rng = np.random.default_rng(2); x = np.column_stack([rng.normal(size=30), rng.normal(size=30), rng.uniform(0.05, 6, 30)]); X = c.geom(x); a, q = directions(X)
        u = lambda Y: np.column_stack([c.internal(Y)[0][:, 0], c.internal(Y)[0][:, 1], np.sqrt(2 * c.internal(Y)[0][:, 2])]); d = np.stack([dline(u, X, a[:, k]) for k in range(3)], 2)
        print('a_k . grad u_j = delta_jk: max deviation %.1e' % np.abs(d - np.eye(3)[None]).max())
        class F:
            N = 1
            def __call__(s, Y): r = np.linalg.norm(Y[:, 1] - Y[:, 0], axis=1); return np.exp(-(r - c.RE) ** 2)[:, None]
        mu = c.mC * c.mO / (c.mC + c.mO); g = lambda r: -((4 * (r - c.RE) ** 2 - 2) * np.exp(-(r - c.RE) ** 2) + 2 * (-2 * (r - c.RE) * np.exp(-(r - c.RE) ** 2)) / r) / (2 * mu)
        r = np.linalg.norm(X[:, 1] - X[:, 0], axis=1); hh = 1e-5; dg = (g(r + hh) - g(r - hh)) / (2 * hh); bond = lambda Y: np.linalg.norm(Y[:, 1] - Y[:, 0], axis=1)
        for k in range(3):
            num = dline(lambda Y: c.T_basis(F(), Y), X, a[:, k])[:, 0]; ana = dg * dline(bond, X, a[:, k]); print('D_%d (T f) for f(r_CO): max relative error %.1e' % (k, np.abs(num - ana).max() / np.abs(ana).max()))
        bas = fb.GL0(fb.centres(5, EC), BETA, 0); t1 = dline(lambda Y: c.T_basis(bas, Y), X, a[:, 2]); t2 = c.T_basis(type('D', (), {'__call__': lambda s, Y: dline(bas, Y, a[:, 2])})(), X)
        print('commutation T D_q phi = D_q T phi (Gaussian basis): max relative difference %.1e' % (np.abs(t1 - t2).max() / np.abs(t1).max()))
    elif mode == 'ff':
        REF = json.load(open(f'{c.HERE}/ffbench_ref.json'))[-1]; N = 600
        jobs = [(N, m, False, 1.0) for m in (300, 600, 1200, 2400, 7200)] + [(N, m, True, 1.0) for m in (150, 300, 600, 1200, 2400)] + [(N, 600, True, lam) for lam in (0.25, 4.0)] + [(300, m, True, 1.0) for m in (150, 300, 600)] + [(300, m, False, 1.0) for m in (600, 3600)]
        with Pool(12) as pool: R = pool.map(job_ff, jobs, chunksize=1)
        for r in R: show('N=%d M=%d %s lam=%g (eq/N %.1f)' % (r['N'], r['M'], 'E+grad' if r['grad'] else 'E only', r['lam'], r['M'] * (4 if r['grad'] else 1) / r['N']), r, REF)
        json.dump(R, open(f'{c.HERE}/co2_grad_ff.json', 'w'), indent=1)
    elif mode == 'dftpts':
        n, NP = int(sys.argv[2]), int(sys.argv[3]); x, E = ho.load('fit', 8000); x, E = x[:n], E[:n]; E0 = float(np.load(f'{c.HERE}/E0.npy')); t0 = time.time()
        with Pool(NP) as pool: R = pool.map(dft_one, [(i, x[i], E[i]) for i in range(n)], chunksize=4)
        R.sort(key=lambda t: t[0]); Eg = np.array([t[1] for t in R]) - E0; g = np.array([t[2] for t in R]); a, q = directions(c.geom(x)); gu = np.einsum('ni,nki->nk', g, a)
        print('%d gradients in %.0f s; converged %d; max |E - cached E| = %.1e cm-1' % (n, time.time() - t0, sum(t[3] for t in R), np.abs(Eg - E).max() * H))
        np.savez(f'{c.HERE}/pts_grad.npz', x=x, E=E, g=g, gu=gu)
    elif mode == 'dft':
        REF = dict(zpe=2547.696, nu2=669.679, J0=[1293.594, 1401.139, 2354.930])         # converged HO-basis collocation on the same DFT surface (co2_ho2.json, Ecut 25000)
        nmax = len(np.load(f'{c.HERE}/pts_grad.npz')['E']); N = 600
        jobs = [(N, m, False, 1.0) for m in (300, 600, 1200, 2400) if m <= nmax] + [(N, m, True, 1.0) for m in (150, 300, 600, 1200, 2400) if m <= nmax] + [(300, m, True, 1.0) for m in (150, 300, 600)] + [(N, 600, True, lam) for lam in (0.25, 4.0)]
        with Pool(12) as pool: R = pool.map(job_dft, jobs, chunksize=1)
        for r in R: show('N=%d M=%d %s lam=%g (eq/N %.1f)' % (r['N'], r['M'], 'E+grad' if r['grad'] else 'E only', r['lam'], r['M'] * (4 if r['grad'] else 1) / r['N']), r, REF)
        json.dump(R, open(f'{c.HERE}/co2_grad_dft.json', 'w'), indent=1)
