"""Rectangular-collocation driver for LiCO2 (J=0 levels), following the rectangular-collocation-vib skill.

  python run_rc.py --pes dft  --N 3000 --ec 900 --beta 0.9 --vmax 5000 --p 3 --Ms 750,1500,2250,3000 --nev 30 --workers 12 --threads 4
  python run_rc.py --pes ff   --N 1500 ... --ref ../vpt/vpt2_lico2_big.json   # validation on the quartic FF (compare with VCI)
  python run_rc.py --pes harm --N 1500 ...                                    # harmonic sanity test (exact levels known)

--pes dft  reads cache/fit.npz and cache/val.npz written by dft_points.py (points with V <= vmax only).
--pes ff/harm  generates R2 points on the fly with the surrogate (fast; no DFT).
Basis centres: R2 points (indices from 20,000,000) with surrogate V <= Ec; the surrogate is the degree-4
polynomial fitted to the DFT cache when --pes dft, otherwise the quartic FF / harmonic potential.
Gradient rows: on when --p >= 0, weights lam_i = (1 - V_i/Vmax)^p  (p = 3 recommended).
--beta: one value (isotropic) or six comma-separated values (per normal coordinate).
Outputs: results/levels_<tag>_M<M>.csv, results/<tag>.json, results/SH_<tag>_M<M>.npz (checkpoint).
"""
import os, sys, json, time, argparse, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ap = argparse.ArgumentParser()
ap.add_argument('--pes', choices=['dft', 'ff', 'harm'], required=True)
ap.add_argument('--N', type=int, required=True); ap.add_argument('--ec', type=float, default=900)
ap.add_argument('--beta', default='0.9'); ap.add_argument('--vmax', type=float, default=5000)
ap.add_argument('--p', type=float, default=3.0); ap.add_argument('--Ms', required=True)
ap.add_argument('--nev', type=int, default=30); ap.add_argument('--workers', type=int, default=1)
ap.add_argument('--threads', type=int, default=1); ap.add_argument('--ref', default=None)
ap.add_argument('--tag', default=None); ap.add_argument('--cache', default='cache')
ap.add_argument('--nvalraw', type=int, default=0, help='use only the first n computed validation records (fixed set while the cache grows)')
ap.add_argument('--batch', type=int, default=1500); ap.add_argument('--vval', type=float, default=0, help='validation region energy (default vmax)')
a = ap.parse_args()
os.environ.setdefault('OMP_NUM_THREADS', str(a.threads))
import coords as c, collo as co
from r2 import r2
H = c.H2CM; vmax_h = a.vmax / H
beta = np.array([float(x) for x in a.beta.split(',')]); beta = beta if len(beta) == c.NQ else np.full(c.NQ, beta[0])
tag = a.tag or f"{a.pes}_N{a.N}_ec{a.ec:g}_b{a.beta.replace(',','-')}_v{a.vmax:g}_p{a.p:g}"
os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
if a.pes == 'harm': co.ECKART = False          # rectilinear KEO: exact harmonic levels known

def surrogate_points(start, n, vlim, V):
    lo, hi = (lambda b: (b['lo'], b['hi']))(np.load(os.path.join(HERE, a.cache, 'box.npz'))) \
        if a.pes == 'dft' else (-np.sqrt(2 * vlim * 1.15 / c.omega), np.sqrt(2 * vlim * 1.15 / c.omega))
    out, i = [], start
    while sum(len(x) for x in out) < n:
        q = lo + (hi - lo) * r2(np.arange(i, i + 200000), c.NQ); i += 200000; out.append(q[V(q) <= vlim])
    return np.concatenate(out)[:n]

# ---------------- points ----------------
if a.pes == 'dft':
    def load(name):
        d = np.load(os.path.join(HERE, a.cache, f'{name}.npz'))
        o = np.argsort(d['idx']); m = (d['E'][o] <= vmax_h) & d['conv'][o] & np.isfinite(d['gq'][o]).all(1)
        return d['q'][o][m], d['E'][o][m], d['gq'][o][m]
    qf, Vf, gf = load('fit')
    def loadval():
        d = np.load(os.path.join(HERE, a.cache, 'val.npz')); o = np.argsort(d['idx'])[:a.nvalraw or None]
        m = (d['E'][o] <= vmax_h) & d['conv'][o]; return d['q'][o][m], d['E'][o][m]
    qv, Vv = loadval()
    from surr import polyfit
    d = np.load(os.path.join(HERE, a.cache, 'fit.npz')); ok = d['conv'] & (d['E'] < 2.5 * 5500 / H)
    Vs, rms = polyfit(d['q'][ok], d['E'][ok]); print(f'centre surrogate: degree-4 fit, rms {rms:.1f} cm-1')
else:
    V = c.V_ff if a.pes == 'ff' else c.V_harm
    Mmax = max(int(x) for x in a.Ms.split(','))
    qf = surrogate_points(1, Mmax, vmax_h, V); Vf = V(qf); qv = surrogate_points(10_000_000, 1000, vmax_h, V); Vv = V(qv)
    eps = 1e-4; gf = np.array([(V(qf + eps * e) - V(qf - eps * e)) / (2 * eps) for e in np.eye(c.NQ)]).T
    Vs = V
# centres: surrogate V <= Ec; for the DFT polynomial fit also require the quartic FF to be low (guards against
# spurious minima of the fitted polynomial far from the data)
Vc = (lambda q: np.where(c.V_ff(q) <= 2.5 * a.ec / H, Vs(q), np.inf)) if a.pes == 'dft' else Vs
# centres depend only on (PES/cache, Ec, N): cache them, since rejection sampling in the DFT box takes many minutes
cfn = os.path.join(HERE, 'results', f"centres_{a.pes}_{os.path.basename(a.cache) if a.pes == 'dft' else 'x'}_ec{a.ec:g}_N{a.N}.npy")
if os.path.exists(cfn): C = np.load(cfn)
else:
    C = surrogate_points(20_000_000, a.N, a.ec / H, Vc); np.save(cfn + f'.{os.getpid()}.npy', C); os.replace(cfn + f'.{os.getpid()}.npy', cfn)
bas = co.Basis(C, beta)
if a.vval: _m = Vv <= a.vval / H; qv, Vv = qv[_m], Vv[_m]
lam = np.clip(1 - Vf / vmax_h, 0, 1) ** a.p if a.p >= 0 else None
print(f'{tag}: {len(qf)} fit points available, {len(qv)} validation points, N={a.N}, beta={np.round(beta,3)}', flush=True)

# ---------------- reference levels (optional) ----------------
ref = None
if a.pes == 'harm':
    import itertools
    n_ = np.array(list(itertools.product(range(12), repeat=c.NQ))); e_ = np.sort((n_ + 0.5) @ c.omega)
    ref = e_[:a.nev] * H; ref = ref - ref[0]
elif a.ref:
    j = json.load(open(a.ref)); ref = np.array([0.0] + [e for _, e, _ in j['vci']][1:a.nev])

# ---------------- accumulation (parallel over chunks) ----------------
def chunk(args_):
    lo_, hi_ = args_
    import collo as cc
    if a.pes == 'harm': cc.ECKART = False
    s_, h_ = cc.accumulate(bas, qf[lo_:hi_], Vf[lo_:hi_], gf[lo_:hi_] if lam is not None else None,
                           lam=lam[lo_:hi_] if lam is not None else 1.0, batch=a.batch)
    # hand the N x N matrices back through a file: pickling them through the Pool pipe takes hours at N >= 5000 here
    fn = os.path.join(HERE, 'results', f'tmp_{tag}_{lo_}_{hi_}.npy'); np.save(fn, np.stack([s_, h_])); return fn
S = Hm = None; done = 0; t0 = time.time(); recs = []
from multiprocessing import Pool
for M in [int(x) for x in a.Ms.split(',')]:
    M = min(M, len(qf))
    if M <= done: continue
    edges = np.linspace(done, M, max(2, min(a.workers, M - done) + 1)).astype(int)
    parts = list(zip(edges[:-1], edges[1:]))
    if a.workers > 1:
        with Pool(a.workers) as pool: res = pool.map(chunk, parts)
    else:
        res = [chunk(x) for x in parts]
    for fn in res:
        s_, h_ = np.load(fn); os.remove(fn)
        S = s_ if S is None else S + s_; Hm = h_ if Hm is None else Hm + h_
    done = M; ta = time.time() - t0
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=max(a.workers, a.threads)):      # eigen stage: BLAS threads in the parent only
        ev = np.linalg.eigvalsh(S) if a.N <= 12000 else np.array([np.nan, np.nan])
        E, Cv = co.solve(S, Hm, a.nev, 0.9 * 0.5 * c.omega.sum())
    r = co.residual(bas, qv, Vv, E, Cv); Er = E.real * H; tr = Er - Er[0]
    np.savez(os.path.join(HERE, 'results', f'vec_{tag}_M{M}.npz'), E=E, C=Cv, centres=C, beta=beta, res=r)
    rec = dict(M=M, eq_per_N=M * (1 + (c.NQ if lam is not None else 0)) / a.N, zpe=Er[0], trans=tr.tolist(),
               res_med20=float(np.median(r[:20])), res_med=float(np.median(r)), res_max=float(r.max()), res=r.tolist(), maxIm=float(np.abs(E.imag).max() * H),
               condS=float(ev.max() / max(ev.min(), 1e-300)), t=time.time() - t0)
    if recs:
        dd = np.abs(tr - np.array(recs[-1]['trans'])); rec.update(chg_med=float(np.median(dd[1:])), chg_max=float(dd[1:].max()))
    if ref is not None:
        dr = np.abs(tr[1:len(ref)] - ref[1:]); rec.update(ref_med=float(np.median(dr)), ref_max=float(dr.max()))
    recs.append(rec); json.dump(dict(args=vars(a), recs=recs), open(os.path.join(HERE, 'results', f'{tag}.json'), 'w'), indent=1)
    with open(os.path.join(HERE, 'results', f'levels_{tag}_M{M}.csv'), 'w') as f:
        f.write('level,E_minus_ZPE_cm-1,abs_E_cm-1,Im_cm-1,residual_cm-1\n')
        for i in range(a.nev): f.write(f'{i},{tr[i]:.4f},{Er[i]:.4f},{E[i].imag*H:.2e},{r[i]:.3f}\n')
    s = f"M={M} eq/N={rec['eq_per_N']:.2f} ZPE={rec['zpe']:.2f} | res med20 {rec['res_med20']:.2f} all {rec['res_med']:.2f} | maxIm {rec['maxIm']:.2e} | cond {rec['condS']:.1e}"
    if 'chg_med' in rec: s += f" | change med {rec['chg_med']:.3f} max {rec['chg_max']:.3f}"
    if 'ref_med' in rec: s += f" | vs ref med {rec['ref_med']:.3f} max {rec['ref_max']:.3f}"
    print(s + f" | {rec['t']:.0f}s", flush=True)
