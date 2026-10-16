"""DFT points for LiCO2 rectangular collocation (workstation production).

UB3LYP-D3(BJ)/aug-cc-pVDZ single points at R2 quasi-random points in Eckart-frame dimensionless
normal coordinates. Energy first; gradient (with grid response) and dipole only if V <= Vmax.
Every result is cached by its R2 index in cache/<set>.npz, so the run can be stopped and resumed,
and changing N, M, Ec, widths or Vmax (downward) reuses the same data.

Pre-screening (skill section 1): candidates with surrogate V > SCREEN*Vmax are skipped. The surrogate
is the quartic force field until FIT_AFTER energies exist, then a degree-4 polynomial fitted to all
computed energies. A random AUDIT fraction of rejected candidates is computed anyway and any in-region
point among them is reported (it should be none).

usage:
  python dft_points.py check-grid  --nproc 8                  # grid 3/1e-9 vs grid 4/1e-11 on 6 points
  python dft_points.py fit  --n 3000 --vmax 5000 --nproc 12 --threads 4
  python dft_points.py val  --n 1000 --vmax 5000 --nproc 12 --threads 4
Sets: fit = R2 indices from 1, val = from 10,000,000 (never used for fitting).
Cache arrays: idx, q, E (hartree, relative to E0 at identical settings), gq (dV/dq, NaN if no gradient),
dip (a.u., Eckart frame), conv, s2, cand (number of R2 candidates consumed).
"""
import os, sys, time, argparse, itertools, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import coords as c
from r2 import r2
H = c.H2CM
p = argparse.ArgumentParser()
p.add_argument('mode', choices=['fit', 'val', 'check-grid'])
p.add_argument('--n', type=int, default=1000, help='number of in-region points wanted (V <= vmax)')
p.add_argument('--vmax', type=float, default=5000.0)
p.add_argument('--box', type=float, default=1.15, help='box energy = box * vmax (turning points)')
p.add_argument('--nproc', type=int, default=1); p.add_argument('--threads', type=int, default=1)
p.add_argument('--grid', type=int, default=3); p.add_argument('--conv', type=float, default=1e-9)
p.add_argument('--screen', type=float, default=1.35); p.add_argument('--audit', type=float, default=0.02)
p.add_argument('--fit_after', type=int, default=400)
p.add_argument('--batch', type=int, default=0, help='DFT jobs per batch (default 50*nproc)')
args = p.parse_args()
BAS = 'aug-cc-pvdz'
START = {'fit': 1, 'val': 10_000_000}

def worker_init(threads):
    from pyscf import lib; lib.num_threads(threads); load_dm0()

def dft(job):
    """job = (idx, q, vmax_hartree, grid, conv, want_grad_if_inside). Returns dict."""
    idx, q, vmax, grid, conv, E0 = job
    from pyscf import gto, dft as D
    X = c.X_of_q(q[None])[0]
    mol = gto.M(atom=[[e, x] for e, x in zip(c.EL, X)], basis=BAS, unit='Bohr', spin=1,
                symmetry=False, verbose=0)
    mf = D.UKS(mol).set(xc='b3lyp', conv_tol=conv, max_cycle=150); mf.grids.level = grid; mf.disp = 'd3bj'
    mf.kernel(dm0=DM0)
    if not mf.converged:                       # one retry with level shift + fresh guess
        mf = D.UKS(mol).set(xc='b3lyp', conv_tol=conv, max_cycle=300, level_shift=0.2)
        mf.grids.level = grid; mf.disp = 'd3bj'; mf.kernel(dm0=DM0)
    E = mf.e_tot - (E0 if E0 is not None else 0.0)
    out = dict(idx=idx, q=q, E=E, conv=bool(mf.converged), s2=float(mf.spin_square()[0]),
               gq=np.full(c.NQ, np.nan), dip=np.full(3, np.nan))
    if E0 is not None and E <= vmax and mf.converged:
        g = mf.nuc_grad_method(); g.grid_response = True
        gc = g.kernel().ravel()
        out['gq'] = gc @ c.A                   # dV/dq_k = sum_i g_i dX_i/dq_k
        out['dip'] = mf.dip_moment(unit='AU', verbose=0)
    return out

DM0 = None
def load_dm0():
    global DM0
    DM0 = np.load(os.path.join(HERE, '..', 'vpt', 'ff_lico2.npz'))['dm']

def E0_ref(grid, conv):
    fn = os.path.join(HERE, 'cache', f'E0_g{grid}_c{conv:g}.npy')
    if os.path.exists(fn): return float(np.load(fn))
    r = dft((0, np.zeros(c.NQ), 0.0, grid, conv, None)); assert r['conv']
    np.save(fn, r['E']); return r['E']

def box(vmax_h):
    """1D turning points of the quartic-FF surrogate at the box energy, both signs, per coordinate."""
    lo, hi = np.zeros(c.NQ), np.zeros(c.NQ)
    for k in range(c.NQ):
        harm = np.sqrt(2 * vmax_h / c.omega[k])
        for s in (-1, 1):
            t = np.linspace(0.01, 3 * harm, 3000); qq = np.zeros((len(t), c.NQ)); qq[:, k] = s * t
            v = c.V_ff(qq); above = np.nonzero(v > vmax_h)[0]
            r = t[above[0]] if len(above) else harm * 1.5     # surrogate turned over: fall back
            if not len(above) or v[:above[0]].min() < -1e-6:
                print(f'  box: coord {k} sign {s:+d}: surrogate not monotone, using {r:.2f}', flush=True)
            (hi if s > 0 else lo)[k] = s * r
    return lo, hi

def polyfit(q, E, deg=4):
    exps = [e for d in range(deg + 1) for e in itertools.product(range(d + 1), repeat=c.NQ) if sum(e) == d]
    exps = np.array(exps)
    Phi = lambda Q: np.prod(Q[:, None, :] ** exps[None], axis=2)
    coef = np.linalg.lstsq(Phi(q), E, rcond=None)[0]
    rms = np.sqrt(np.mean((Phi(q) @ coef - E) ** 2)) * H
    return (lambda Q: Phi(np.atleast_2d(Q)) @ coef), rms

def main():
    os.makedirs(os.path.join(HERE, 'cache'), exist_ok=True)
    load_dm0()
    from multiprocessing import Pool
    if args.mode == 'check-grid':
        rng = np.random.default_rng(1); qs = rng.normal(size=(6, c.NQ)) * 0.8
        res = {}
        for grid, conv in ((3, 1e-9), (4, 1e-11)):
            E0 = E0_ref(grid, conv)
            with Pool(args.nproc, worker_init, (args.threads,)) as pool:
                rr = pool.map(dft, [(i, q, 1.0, grid, conv, E0) for i, q in enumerate(qs)])
            res[grid] = np.array([[r['E'] * H, *r['gq']] for r in rr])
        d = res[3] - res[4]
        print('grid 3 vs 4: |dE| max %.3f cm-1, |dgrad| max %.2e Eh/q' % (np.abs(d[:, 0]).max(), np.abs(d[:, 1:]).max() / H))
        return
    vmax_h = args.vmax / H
    fn = os.path.join(HERE, 'cache', f'{args.mode}.npz')
    if os.path.exists(fn):
        d = dict(np.load(fn)); recs = [dict(idx=i, q=q, E=e, gq=g, dip=dp, conv=cv, s2=s) for i, q, e, g, dp, cv, s in
                                       zip(d['idx'], d['q'], d['E'], d['gq'], d['dip'], d['conv'], d['s2'])]
        nxt = int(d['next'])
    else:
        recs, nxt = [], START[args.mode]
    lo, hi = box(args.box * vmax_h)
    np.savez(os.path.join(HERE, 'cache', 'box.npz'), lo=lo, hi=hi, vmax=args.vmax, boxfac=args.box)
    E0 = E0_ref(args.grid, args.conv)
    print(f'{args.mode}: {len(recs)} cached, next R2 index {nxt}; box lo {np.round(lo,2)} hi {np.round(hi,2)}', flush=True)
    rng = np.random.default_rng(nxt); surr = c.V_ff; t0 = time.time(); audit_bad = 0
    def save():
        np.savez(fn, idx=np.array([r['idx'] for r in recs]), q=np.array([r['q'] for r in recs]),
                 E=np.array([r['E'] for r in recs]), gq=np.array([r['gq'] for r in recs]),
                 dip=np.array([r['dip'] for r in recs]), conv=np.array([r['conv'] for r in recs]),
                 s2=np.array([r['s2'] for r in recs]), next=nxt)
    with Pool(args.nproc, worker_init, (args.threads,)) as pool:
        while sum(1 for r in recs if r['E'] <= vmax_h and r['conv']) < args.n:
            if len(recs) >= args.fit_after and surr is c.V_ff:
                ok = [r for r in recs if r['conv'] and r['E'] < 3 * vmax_h]
                f, rms = polyfit(np.array([r['q'] for r in ok]), np.array([r['E'] for r in ok]))
                print(f'  surrogate -> degree-4 fit to {len(ok)} energies, rms {rms:.1f} cm-1', flush=True)
                surr = f
            batch = args.batch or 50 * args.nproc
            ids = np.arange(nxt, nxt + 20 * batch); qs = lo + (hi - lo) * r2(ids, c.NQ)
            vs = surr(qs); keep = vs <= args.screen * vmax_h
            audit = (~keep) & (rng.random(len(ids)) < args.audit)
            sel = np.nonzero(keep | audit)[0][:batch]
            if len(sel) == 0: nxt = int(ids[-1]) + 1; continue
            nxt = int(ids[sel[-1]]) + 1
            jobs = [(int(ids[i]), qs[i], vmax_h, args.grid, args.conv, E0) for i in sel]
            for r in pool.imap_unordered(dft, jobs, chunksize=1):
                recs.append(r)
                i = int(np.nonzero(ids == r['idx'])[0][0])
                if audit[i] and r['E'] <= vmax_h: audit_bad += 1
            save()
            nin = sum(1 for r in recs if r['E'] <= vmax_h and r['conv'])
            nbad = sum(1 for r in recs if not r['conv']); s2max = max(r['s2'] for r in recs)
            print(f'  [{time.time()-t0:7.0f}s] computed {len(recs)}, in region {nin}/{args.n}, unconverged {nbad}, '
                  f'max <S2> {s2max:.4f}, audit misses {audit_bad}', flush=True)
    save(); print('done', flush=True)

if __name__ == '__main__':
    main()
