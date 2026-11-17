"""Max-residual parameter tuning for LiCO2 collocation on the UB3LYP-D3(BJ)/aug-cc-pVDZ DFT points, at N = 10000 (14 Nov 2026).
Ec = 1200 fixed (user decision after the N = 10000 Ec scan, scan_ec_N10000.py).
J = MAXIMUM over the 100 computed levels of the unweighted validation residual r_i = ||(H - E_i) psi_i|| / ||psi_i|| (energy rows,
cm-1) on a FIXED validation set: the first 3920 validation records of data/cache_adz with V <= 5000 cm-1.
Fixed setting of every trial: run_rc.py --pes dft --cache data/cache_adz --N 10000 --ec 1200 --Ms 12000 (first 12000 fit points by
index with V <= Vmax, energy + gradient rows, weights (1 - V/Vmax)^p) --nev 100 --nvalraw 3920 --vval 5000.
Start: Ec = 1200, Vmax = 5500, p = 0, beta = 1.00 x 6 (the N = 10000 Ec-scan trial is reused from the cache).
Stage 1b: one parameter at a time, accepting a value if it lowers J: Vmax in {5000, 5500}; then p in {0, 1, 3}.
Stage 1c: isotropic beta0 scanned around 1.00 in steps of 0.05 (0.95 and 1.05, then further while J decreases).
Stage 2: coordinate descent over k = 1..6 (nu1 ... nu6), beta_k x 0.95 / x 1.05, accept the better if it lowers J, sweeps until a sweep
         accepts nothing (at most 5).  Stage 3: the same with x 0.99 / x 1.01 (at most 20 sweeps).  (As for Cl2CO.)
One trial at a time (1 process, 2 threads; two at a time ran out of memory at N = 10000).
Trials are cached in results/<tag>.json; the script is restartable.  It exits (code 3) instead of starting a new trial when fewer than
45 min remain of the 7000 s run window, so that no trial is cut off; it is then relaunched."""
import numpy as np, json, os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, 'tune_li_max10k.log')
N, M, EC, NEV, NVALRAW, VVAL, CACHE = 10000, 12000, 1200.0, 100, 3920, 5000, 'data/cache_adz'
T0, WINDOW, NEED = time.time(), 7000, 2700
def log(s): print(s, flush=True); open(LOG, 'a').write(s + '\n')
def tag(P): return f"tmax_N{N}_M{M}_ec{EC:g}_v{P['vmax']:g}_p{P['p']:g}_b" + '-'.join(f"{b:.6f}" for b in P['beta'])
def trial(P):
    t = tag(P); fn = os.path.join(HERE, 'results', t + '.json')
    if not os.path.exists(fn):
        if time.time() - T0 + NEED > WINDOW: log('--- run window nearly used up; exiting before the next trial (relaunch to continue)'); sys.exit(3)
        env = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2')
        cmd = ['python3', 'run_rc.py', '--pes', 'dft', '--cache', CACHE, '--N', str(N), '--ec', f"{EC:g}",
               '--beta', ','.join(repr(float(b)) for b in P['beta']), '--vmax', f"{P['vmax']:g}", '--p', f"{P['p']:g}", '--Ms', str(M),
               '--nev', str(NEV), '--workers', '1', '--threads', '2', '--vval', str(VVAL), '--nvalraw', str(NVALRAW), '--batch', '1500', '--tag', t]
        subprocess.run(cmd, cwd=HERE, env=env, check=True, stdout=subprocess.DEVNULL)
        v = os.path.join(HERE, 'results', f'vec_{t}_M{M}.npz')
        if os.path.exists(v): os.remove(v)
    return json.load(open(fn))['recs'][-1]
def J(P):
    r = trial(P); res = np.array(r['res']); tr = np.array(r['trans'])
    log(f"vmax={P['vmax']:g} p={P['p']:g} beta={[round(b, 6) for b in P['beta']]} J={res.max():.4f} (level #{res.argmax()}) | 2nd {np.sort(res)[-2]:.3f} | "
        f"med {np.median(res):.3f} | med20 {np.median(res[:20]):.3f} | maxIm {r['maxIm']:.2f} | ZPE {r['zpe']:.2f} | top {tr[-1]:.1f} | t {r['t']:.0f}s")
    return float(res.max())
def iso(P, b0): return dict(P, beta=[b0] * 6)
log(f'=== tune LiCO2 aDZ: N={N}, M={M} (energy + gradient rows), Ec={EC:g} fixed, nev={NEV}, validation = first {NVALRAW} records with V <= {VVAL}; J = MAXIMUM residual')
P = dict(vmax=5500.0, p=0.0, beta=[1.0] * 6); Jc = J(P)
# ---- stage 1b
for key, vals in (('vmax', [5000.0, 5500.0]), ('p', [0.0, 1.0, 3.0])):
    cand = [dict(P, **{key: v}) for v in vals if v != P[key]]; js = [J(c) for c in cand]
    log(f"stage 1b {key}: current {P[key]:g} J={Jc:.4f}; " + ', '.join(f"{c[key]:g}: {j:.4f}" for c, j in zip(cand, js)))
    i = int(np.argmin(js))
    if js[i] < Jc: P, Jc = cand[i], js[i]; log(f'  accepted {key} = {P[key]:g}, J = {Jc:.4f}')
# ---- stage 1c
b0 = 1.0; Js = {b0: Jc}
for nb in (0.95, 1.05): Js[nb] = J(iso(P, nb))
while True:
    b = min(Js, key=Js.get)
    if b == min(Js): nb = round(b - 0.05, 2)
    elif b == max(Js): nb = round(b + 0.05, 2)
    else: break
    Js[nb] = J(iso(P, nb))
    if Js[nb] >= Js[b]: break
b0 = min(Js, key=Js.get); Jc = Js[b0]; P = iso(P, b0)
log(f'stage 1c: J(beta0) = { {k: round(v, 4) for k, v in sorted(Js.items())} } -> beta0 = {b0}, J = {Jc:.4f}')
log(f"stage 1 result: ec={EC:g} vmax={P['vmax']:g} p={P['p']:g} beta0={b0} J={Jc:.4f}")
# ---- stages 2 and 3
def descend(P, Jc, step, maxsw, name):
    for sweep in range(1, maxsw + 1):
        ch = False
        for k in range(6):
            cand = []
            for f in (1 - step, 1 + step):
                b = list(P['beta']); b[k] *= f; cand.append(dict(P, beta=b))
            js = [J(c) for c in cand]; i = int(np.argmin(js))
            if js[i] < Jc:
                P, Jc, ch = cand[i], js[i], True; log(f"  {name} sweep {sweep}: accepted beta[{k + 1}] x{[1 - step, 1 + step][i]} -> {P['beta'][k]:.6f}, J = {Jc:.4f}")
            else:
                log(f"  {name} sweep {sweep}: beta[{k + 1}] unchanged (J x{1 - step} = {js[0]:.4f}, x{1 + step} = {js[1]:.4f}, current {Jc:.4f})")
        if not ch: break
    return P, Jc, sweep
P, Jc, s2 = descend(P, Jc, 0.05, 5, 'stage 2 (5%)')
log(f"stage 2 result: beta = {[round(b, 6) for b in P['beta']]}, J = {Jc:.4f}")
P, Jc, s3 = descend(P, Jc, 0.01, 20, 'stage 3 (1%)')
json.dump(dict(P, ec=EC, J=Jc, beta0=b0, sweeps2=s2, sweeps3=s3), open(os.path.join(HERE, 'results', 'best_li_adz_max10k.json'), 'w'))
log(f"=== result: ec={EC:g} vmax={P['vmax']:g} p={P['p']:g} beta={[round(b, 6) for b in P['beta']]} J={Jc:.4f}")
