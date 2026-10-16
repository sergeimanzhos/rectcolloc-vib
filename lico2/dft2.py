"""DFT points for LiCO2 rectangular collocation, any basis set (extends dft_points.py; workstation production).

UB3LYP-D3(BJ)/<basis> at R2 quasi-random points in the Eckart-frame dimensionless normal coordinates of
coords.py (aDZ grid-8 reference geometry and modes; used as the coordinate system for every level of theory).
Cache: cache_<tag>/{fit,val}.npz (same arrays as dft_points.py), E relative to cache_<tag>/E0.npy.

  python dft2.py ref --tag adz --basis aug-cc-pvdz            # minimum (Newton in q), E0, dm0, qmin
  python dft2.py fit --tag adz --basis aug-cc-pvdz --n 8000 --vmax 5500 --nproc 48
  python dft2.py val --tag adz --basis aug-cc-pvdz --n 1300 --vmax 5500 --nproc 48
  python dft2.py fit --tag atz --basis aug-cc-pvtz --src adz --vmax 5000 --margin 1.06 [--limit 3000]
--src: take the candidate points from another level's cache (those with E_src <= margin*vmax); no screening needed.
"""
import os, sys, time, argparse, itertools, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import coords as c
from r2 import r2
H = c.H2CM
p = argparse.ArgumentParser()
p.add_argument('mode', choices=['fit', 'val', 'ref'])
p.add_argument('--tag', required=True); p.add_argument('--basis', required=True)
p.add_argument('--n', type=int, default=1000); p.add_argument('--vmax', type=float, default=5000.0)
p.add_argument('--box', type=float, default=1.15)
p.add_argument('--nproc', type=int, default=48)
p.add_argument('--grid', type=int, default=4); p.add_argument('--conv', type=float, default=1e-9)
p.add_argument('--screen', type=float, default=1.3); p.add_argument('--audit', type=float, default=0.03)
p.add_argument('--fit_after', type=int, default=400)
p.add_argument('--src', default=None); p.add_argument('--margin', type=float, default=1.06)
p.add_argument('--limit', type=int, default=0); p.add_argument('--srcmax', type=float, default=0, help='cm-1; default margin*vmax')
args = p.parse_args()
CD = os.path.join(HERE, f'cache_{args.tag}')
START = {'fit': 1, 'val': 10_000_000}
DM0 = None

def worker_init():
    global DM0
    from pyscf import lib; lib.num_threads(1)
    DM0 = np.load(os.path.join(CD, 'dm0.npy'))

def make_mf(q, level_shift=0.0, max_cycle=150):
    from pyscf import gto, dft as D
    X = c.X_of_q(np.asarray(q)[None])[0]
    mol = gto.M(atom=[[e, x] for e, x in zip(c.EL, X)], basis=args.basis, unit='Bohr', spin=1, symmetry=False, verbose=0)
    mf = D.UKS(mol).set(xc='b3lyp', conv_tol=args.conv, max_cycle=max_cycle, level_shift=level_shift)
    mf.grids.level = args.grid; mf.disp = 'd3bj'
    return mf

def grad_q(mf):
    g = mf.nuc_grad_method(); g.grid_response = True
    return g.kernel().ravel() @ c.A

def dft(job):
    idx, q, vmax, E0 = job
    out = dict(idx=idx, q=q, E=np.nan, conv=False, s2=np.nan, gq=np.full(c.NQ, np.nan), dip=np.full(3, np.nan))
    try:
        mf = make_mf(q); mf.kernel(dm0=DM0)
        if not mf.converged:
            mf = make_mf(q, 0.2, 300); mf.kernel(dm0=DM0)
        out.update(E=mf.e_tot - E0, conv=bool(mf.converged), s2=float(mf.spin_square()[0]))
        if out['E'] <= vmax and mf.converged:
            out['gq'] = grad_q(mf); out['dip'] = mf.dip_moment(unit='AU', verbose=0)
    except Exception as ex:
        print('  job', idx, 'failed:', repr(ex)[:200], flush=True)
    return out

def box(vmax_h):
    lo, hi = np.zeros(c.NQ), np.zeros(c.NQ)
    for k in range(c.NQ):
        harm = np.sqrt(2 * vmax_h / c.omega[k])
        for s in (-1, 1):
            t = np.linspace(0.01, 3 * harm, 3000); qq = np.zeros((len(t), c.NQ)); qq[:, k] = s * t
            v = c.V_ff(qq); above = np.nonzero(v > vmax_h)[0]
            r = t[above[0]] if len(above) else harm * 1.5
            (hi if s > 0 else lo)[k] = s * r
    return lo, hi

EXPS = np.array([e for d in range(5) for e in itertools.product(range(d + 1), repeat=c.NQ) if sum(e) == d])
def Phi(Q):
    Q = np.atleast_2d(Q); out = np.ones((len(Q), len(EXPS)))
    for k in range(c.NQ):
        pw = np.stack([Q[:, k] ** n for n in range(5)], 1); out *= pw[:, EXPS[:, k]]
    return out
def polyfit(q, E):
    coef = np.linalg.lstsq(Phi(q), E, rcond=None)[0]
    rms = np.sqrt(np.mean((Phi(q) @ coef - E) ** 2)) * H
    return (lambda Q: Phi(Q) @ coef), rms

def load(fn):
    d = dict(np.load(fn))
    recs = [dict(idx=int(i), q=q, E=e, gq=g, dip=dp, conv=bool(cv), s2=s) for i, q, e, g, dp, cv, s in
            zip(d['idx'], d['q'], d['E'], d['gq'], d['dip'], d['conv'], d['s2'])]
    return recs, int(d['next'])

def save(fn, recs, nxt):
    tmp = fn + '.tmp.npz'
    np.savez(tmp, idx=np.array([r['idx'] for r in recs]), q=np.array([r['q'] for r in recs]),
             E=np.array([r['E'] for r in recs]), gq=np.array([r['gq'] for r in recs]),
             dip=np.array([r['dip'] for r in recs]), conv=np.array([r['conv'] for r in recs]),
             s2=np.array([r['s2'] for r in recs]), next=nxt)
    os.replace(tmp, fn)

def ref():
    """Minimum at this level by Newton steps in q (diagonal aDZ Hessian), E0, dm0."""
    from pyscf import lib; lib.num_threads(args.nproc)
    q = np.zeros(c.NQ); dm = None
    if os.path.exists(os.path.join(CD, 'qmin.npy')): q = np.load(os.path.join(CD, 'qmin.npy'))
    for it in range(15):
        mf = make_mf(q); mf.kernel(dm0=dm); dm = mf.make_rdm1(); g = grad_q(mf)
        print(f'  ref it {it}: E {mf.e_tot:.10f} |g| {np.abs(g).max():.2e} q {np.round(q, 5)}', flush=True)
        if np.abs(g).max() < 1e-6: break
        q = q - g / c.omega
    np.save(os.path.join(CD, 'qmin.npy'), q); np.save(os.path.join(CD, 'E0.npy'), mf.e_tot); np.save(os.path.join(CD, 'dm0.npy'), dm)
    np.save(os.path.join(CD, 'dip0.npy'), mf.dip_moment(unit='AU', verbose=0))
    print(f'ref {args.tag}: E0 {mf.e_tot:.10f} qmin {np.round(q, 5)} <S2> {mf.spin_square()[0]:.4f}', flush=True)

def main():
    os.makedirs(CD, exist_ok=True)
    if args.mode == 'ref': return ref()
    from multiprocessing import Pool
    vmax_h = args.vmax / H; fn = os.path.join(CD, f'{args.mode}.npz')
    recs, nxt = load(fn) if os.path.exists(fn) else ([], START[args.mode])
    E0 = float(np.load(os.path.join(CD, 'E0.npy'))); t0 = time.time()
    def status():
        nin = sum(1 for r in recs if r['conv'] and r['E'] <= vmax_h); nbad = sum(1 for r in recs if not r['conv'])
        s2 = [r['s2'] for r in recs if r['conv']]
        return f'[{time.time()-t0:7.0f}s] computed {len(recs)}, in region {nin}, unconverged {nbad}, max <S2> {max(s2) if s2 else 0:.4f}'
    if args.src:                                   # candidates from another level's cache
        s, _ = load(os.path.join(HERE, f'cache_{args.src}', f'{args.mode}.npz'))
        have = {r['idx'] for r in recs}
        smax = (args.srcmax or args.margin * args.vmax) / H
        todo = sorted([r for r in s if r['conv'] and r['E'] <= smax], key=lambda r: r['idx'])
        if args.limit: todo = todo[:args.limit]
        todo = [r for r in todo if r['idx'] not in have]
        print(f'{args.tag} {args.mode}: {len(recs)} cached, {len(todo)} to do (from {args.src})', flush=True)
        with Pool(args.nproc, worker_init) as pool:
            for a in range(0, len(todo), 5 * args.nproc):
                jobs = [(r['idx'], r['q'], vmax_h, E0) for r in todo[a:a + 5 * args.nproc]]
                for r in pool.imap_unordered(dft, jobs, chunksize=1): recs.append(r)
                save(fn, recs, nxt); print('  ' + status(), flush=True)
        print('done', flush=True); return
    bf = os.path.join(CD, 'box.npz')
    if os.path.exists(bf): b = np.load(bf); lo, hi = b['lo'], b['hi']
    else:
        lo, hi = box(args.box * vmax_h); np.savez(bf, lo=lo, hi=hi, vmax=args.vmax, boxfac=args.box)
    extra = load(os.path.join(CD, 'fit.npz'))[0] if args.mode == 'val' and os.path.exists(os.path.join(CD, 'fit.npz')) else []
    print(f'{args.tag} {args.mode}: {len(recs)} cached, next R2 index {nxt}; box lo {np.round(lo,2)} hi {np.round(hi,2)}', flush=True)
    rng = np.random.default_rng(nxt); audit_n = audit_bad = 0
    with Pool(args.nproc, worker_init) as pool:
        while sum(1 for r in recs if r['conv'] and r['E'] <= vmax_h) < args.n:
            ok = [r for r in recs + extra if r['conv'] and r['E'] < 2.5 * vmax_h]
            if len(ok) >= args.fit_after:
                surr, rms = polyfit(np.array([r['q'] for r in ok]), np.array([r['E'] for r in ok]))
            else: surr, rms = c.V_ff, -1
            batch = 10 * args.nproc
            ids = np.arange(nxt, nxt + 40 * batch); qs = lo + (hi - lo) * r2(ids, c.NQ)
            keep = surr(qs) <= args.screen * vmax_h
            audit = (~keep) & (rng.random(len(ids)) < args.audit)
            sel = np.nonzero(keep | audit)[0][:batch]
            if len(sel) == 0: nxt = int(ids[-1]) + 1; continue
            nxt = int(ids[sel[-1]]) + 1; isaud = {int(ids[i]) for i in sel if audit[i]}
            jobs = [(int(ids[i]), qs[i], vmax_h, E0) for i in sel]
            for r in pool.imap_unordered(dft, jobs, chunksize=1):
                recs.append(r)
                if r['idx'] in isaud:
                    audit_n += 1; audit_bad += int(r['conv'] and r['E'] <= vmax_h)
            save(fn, recs, nxt)
            print(f'  {status()}/{args.n}, surrogate rms {rms:.0f}, audit misses {audit_bad}/{audit_n}', flush=True)
    print('done', flush=True)

if __name__ == '__main__':
    main()
