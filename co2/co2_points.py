"""B3LYP-D3(BJ)/aug-cc-pVTZ energies of free CO2 at R2 quasi-random points in (x1, x3, xb) for collocation (skill section 1).
  python co2_points.py fit|val <n in region> <vmax cm-1> <nproc>
Cache: pts_<set>.npz (idx, x, E relative to the energy at the reference geometry, conv); resumable. Candidates are screened with
the harmonic potential at first (generously) and, once 600 energies exist, with a degree-6 polynomial fit at 1.3 Vmax."""
import os, sys, time, itertools, numpy as np
import co2c as c
SET, NIN, VMAX, NP = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]); H = c.H2CM; vmax = VMAX / H
START = {'fit': 1, 'val': 5_000_000}; FN = os.path.join(c.HERE, f'pts_{SET}.npz'); VBOX = 1.4 * 8000 / H; LO, HI = c.box(VBOX)
BAS = 'aug-cc-pvtz'
def mf_at(X):
    from pyscf import gto, dft
    mol = gto.M(atom=[[e, x] for e, x in zip(('C', 'O', 'O'), X)], basis=BAS, unit='Bohr', symmetry=False, verbose=0)
    mf = dft.RKS(mol).set(xc='b3lyp', conv_tol=1e-10, max_cycle=150); mf.grids.level = 5; mf.disp = 'd3bj'; return mf
def init():
    global DM0
    from pyscf import lib; lib.num_threads(1); DM0 = np.load(os.path.join(c.HERE, 'dm0.npy'))
def job(a):
    idx, x = a
    try:
        mf = mf_at(c.geom(x[None])[0]); mf.kernel(dm0=DM0)
        if not mf.converged: mf = mf_at(c.geom(x[None])[0]); mf.level_shift = 0.2; mf.max_cycle = 300; mf.kernel(dm0=DM0)
        return idx, x, mf.e_tot, bool(mf.converged)
    except Exception as ex: print('job', idx, 'failed', repr(ex)[:100], flush=True); return idx, x, np.nan, False
EX = np.array([e for d in range(7) for e in itertools.product(range(d + 1), repeat=3) if sum(e) == d])
def Phi(x): return np.prod(x[:, None, :] ** EX[None], axis=2)
if __name__ == '__main__':
    from pyscf import lib; lib.num_threads(min(NP, 32))
    e0f = os.path.join(c.HERE, 'E0.npy')
    if not os.path.exists(e0f):
        mf = mf_at(c.geom(np.zeros((1, 3)))[0]); mf.kernel(); assert mf.converged; np.save(e0f, mf.e_tot); np.save(os.path.join(c.HERE, 'dm0.npy'), mf.make_rdm1()); print('E0 = %.10f' % mf.e_tot, flush=True)
    E0 = float(np.load(e0f)); t0 = time.time()
    if os.path.exists(FN): d = np.load(FN); idx, xs, Es, cv, nxt = list(d['idx']), list(d['x']), list(d['E']), list(d['conv']), int(d['next'])
    else: idx, xs, Es, cv, nxt = [], [], [], [], START[SET]
    other = os.path.join(c.HERE, 'pts_fit.npz'); extra = np.load(other) if (SET == 'val' and os.path.exists(other)) else None
    from multiprocessing import Pool
    nin = lambda: int(sum(1 for e, k in zip(Es, cv) if k and e <= vmax))
    with Pool(NP, init) as pool:
        while nin() < NIN:
            X_ = np.array(xs + (list(extra['x']) if extra is not None else [])); E_ = np.array(Es + (list(extra['E']) if extra is not None else [])); ok = np.isfinite(E_) & (E_ < 2.5 * vmax) if len(E_) else np.array([], bool)
            if ok.sum() >= 600: coef = np.linalg.lstsq(Phi(X_[ok]), E_[ok], rcond=None)[0]; surr = lambda q: Phi(q) @ coef; lim = 1.3 * vmax; rms = np.sqrt(np.mean((Phi(X_[ok]) @ coef - E_[ok]) ** 2)) * H
            else: surr = c.V_harm; lim = 2.2 * vmax; rms = -1
            ids = np.arange(nxt, nxt + 60 * NP); q = LO + (HI - LO) * c.r2seq(ids); keep = np.nonzero(surr(q) <= lim)[0][:8 * NP]
            if len(keep) == 0: nxt = int(ids[-1]) + 1; continue
            nxt = int(ids[keep[-1]]) + 1
            for i, x, e, k in pool.imap_unordered(job, [(int(ids[j]), q[j]) for j in keep], chunksize=1): idx.append(i); xs.append(x); Es.append(e - E0); cv.append(k)
            np.savez(FN + '.tmp.npz', idx=np.array(idx), x=np.array(xs), E=np.array(Es), conv=np.array(cv), next=nxt); os.replace(FN + '.tmp.npz', FN)
            print('[%6.0fs] %s: computed %d, in region %d/%d, unconverged %d, surrogate rms %.0f cm-1' % (time.time() - t0, SET, len(Es), nin(), NIN, int(np.sum(~np.array(cv))), rms), flush=True)
    print('done', SET, flush=True)
