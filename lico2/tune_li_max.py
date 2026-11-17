"""Max-residual parameter tuning for LiCO2 collocation on the UB3LYP-D3(BJ)/aug-cc-pVDZ DFT points (14 Nov 2026).
J = MAXIMUM over the 100 computed levels of the unweighted validation residual r_i = ||(H - E_i) psi_i|| / ||psi_i|| (energy rows,
cm-1) on a FIXED validation set: the first 3920 validation records of data/cache_adz with V <= 5000 cm-1.
Fixed setting of every trial: run_rc.py --pes dft --cache data/cache_adz --N 5000 --Ms 6000 (first 6000 fit points by index with
V <= Vmax, energy + gradient rows, weights (1 - V/Vmax)^p) --nev 100 --nvalraw 3920 --vval 5000.
Stage 1a: isotropic beta0 in {0.80, 0.85, 0.90, 0.95, 1.00} (extended by 0.05 while the minimum is at an end), Ec = 900, Vmax = 5500, p = 0
          (the values of the earlier LiCO2 runs).
Stage 1b: at that beta0, one parameter at a time, accepting a value if it lowers J:
          Ec (V'max, centre cutoff) in {700, 900, 1200, 1500}; Vmax in {5000, 5500}; p in {0, 1, 3}.
          (Vmax > 5500 is not scanned: the box of this point set was built for 5500 cm-1.)
Stage 1c: if stage 1b changed anything, beta0 is scanned again around its value (steps of 0.05, while J decreases).
Stage 2: coordinate descent over k = 1..6 (nu1 ... nu6), beta_k x 0.95 / x 1.05, accept the better if it lowers J, sweeps until a sweep
         accepts nothing (at most 5).  Stage 3: the same with x 0.99 / x 1.01 (at most 20 sweeps).  (As for Cl2CO.)
Two trials run at a time (each 1 process, 1 thread). Trials are cached in results/<tag>.json; the script can be restarted."""
import numpy as np, json, os, subprocess
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, 'tune_li_max.log')
N, M, NEV, NVALRAW, VVAL, CACHE = 5000, 6000, 100, 3920, 5000, 'data/cache_adz'
def log(s): print(s, flush=True); open(LOG, 'a').write(s + '\n')
def tag(P): return f"tmax_N{N}_M{M}_ec{P['ec']:g}_v{P['vmax']:g}_p{P['p']:g}_b" + '-'.join(f"{b:.6f}" for b in P['beta'])
def trial(P):
    t = tag(P); fn = os.path.join(HERE, 'results', t + '.json')
    if not os.path.exists(fn):
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        cmd = ['python3', 'run_rc.py', '--pes', 'dft', '--cache', CACHE, '--N', str(N), '--ec', f"{P['ec']:g}",
               '--beta', ','.join(repr(float(b)) for b in P['beta']), '--vmax', f"{P['vmax']:g}", '--p', f"{P['p']:g}", '--Ms', str(M),
               '--nev', str(NEV), '--workers', '1', '--threads', '1', '--vval', str(VVAL), '--nvalraw', str(NVALRAW), '--batch', '1500', '--tag', t]
        subprocess.run(cmd, cwd=HERE, env=env, check=True, stdout=subprocess.DEVNULL)
        v = os.path.join(HERE, 'results', f'vec_{t}_M{M}.npz')
        if os.path.exists(v): os.remove(v)
    return json.load(open(fn))['recs'][-1]
def J(P):
    r = trial(P); res = np.array(r['res'])
    log(f"ec={P['ec']:g} vmax={P['vmax']:g} p={P['p']:g} beta={[round(b, 6) for b in P['beta']]} J={res.max():.4f} | res med {np.median(res):.3f} "
        f"med20 {np.median(res[:20]):.3f} | maxIm {r['maxIm']:.2f} | ZPE {r['zpe']:.2f} | t {r['t']:.0f}s")
    return float(res.max())
def Jmany(Ps):
    with ThreadPoolExecutor(2) as ex: return list(ex.map(J, Ps))
def iso(P, b0): return dict(P, beta=[b0] * 6)
log(f'=== tune LiCO2 aDZ: N={N}, M={M} (energy + gradient rows), nev={NEV}, validation = first {NVALRAW} records with V <= {VVAL}; J = MAXIMUM residual')
P0 = dict(ec=900.0, vmax=5500.0, p=0.0, beta=[1.0] * 6)
# ---- stage 1a
grid = [0.80, 0.85, 0.90, 0.95, 1.00]; Js = dict(zip(grid, Jmany([iso(P0, b) for b in grid])))
def walk(P, Js):
    while True:
        b = min(Js, key=Js.get)
        if b == min(Js): nb = round(b - 0.05, 2)
        elif b == max(Js): nb = round(b + 0.05, 2)
        else: return Js
        Js[nb] = J(iso(P, nb))
        if Js[nb] >= Js[b]: return Js
Js = walk(P0, Js); b0 = min(Js, key=Js.get); P = iso(P0, b0); Jc = Js[b0]
log(f'stage 1a: J(beta0) = { {k: round(v, 4) for k, v in sorted(Js.items())} } -> beta0 = {b0}, J = {Jc:.4f}')
# ---- stage 1b
changed = False
for key, vals in (('ec', [700.0, 900.0, 1200.0, 1500.0]), ('vmax', [5000.0, 5500.0]), ('p', [0.0, 1.0, 3.0])):
    cand = [dict(P, **{key: v}) for v in vals if v != P[key]]
    js = Jmany(cand)
    log(f"stage 1b {key}: current {P[key]:g} J={Jc:.4f}; " + ', '.join(f"{c[key]:g}: {j:.4f}" for c, j in zip(cand, js)))
    i = int(np.argmin(js))
    if js[i] < Jc: P, Jc, changed = cand[i], js[i], True; log(f'  accepted {key} = {P[key]:g}, J = {Jc:.4f}')
# ---- stage 1c
if changed:
    Js = {b0: Jc}
    for nb in (round(b0 - 0.05, 2), round(b0 + 0.05, 2)): Js[nb] = J(iso(P, nb))
    Js = walk(P, Js); b0 = min(Js, key=Js.get); Jc = Js[b0]; P = iso(P, b0)
    log(f'stage 1c: J(beta0) = { {k: round(v, 4) for k, v in sorted(Js.items())} } -> beta0 = {b0}, J = {Jc:.4f}')
log(f"stage 1 result: ec={P['ec']:g} vmax={P['vmax']:g} p={P['p']:g} beta0={b0} J={Jc:.4f}")
# ---- stages 2 and 3
def descend(P, Jc, step, maxsw, name):
    for sweep in range(1, maxsw + 1):
        ch = False
        for k in range(6):
            cand = []
            for f in (1 - step, 1 + step):
                b = list(P['beta']); b[k] *= f; cand.append(dict(P, beta=b))
            js = Jmany(cand); i = int(np.argmin(js))
            if js[i] < Jc:
                P, Jc, ch = cand[i], js[i], True; log(f"  {name} sweep {sweep}: accepted beta[{k + 1}] x{[1 - step, 1 + step][i]} -> {P['beta'][k]:.6f}, J = {Jc:.4f}")
            else:
                log(f"  {name} sweep {sweep}: beta[{k + 1}] unchanged (J x{1 - step} = {js[0]:.4f}, x{1 + step} = {js[1]:.4f}, current {Jc:.4f})")
        if not ch: break
    return P, Jc, sweep
P, Jc, s2 = descend(P, Jc, 0.05, 5, 'stage 2 (5%)')
log(f"stage 2 result: beta = {[round(b, 6) for b in P['beta']]}, J = {Jc:.4f}")
P, Jc, s3 = descend(P, Jc, 0.01, 20, 'stage 3 (1%)')
json.dump(dict(P, J=Jc, beta0=b0, sweeps2=s2, sweeps3=s3), open(os.path.join(HERE, 'results', 'best_li_adz_max.json'), 'w'))
log(f"=== result: ec={P['ec']:g} vmax={P['vmax']:g} p={P['p']:g} beta={[round(b, 6) for b in P['beta']]} J={Jc:.4f}")
