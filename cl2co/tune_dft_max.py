"""Max-residual width tuning for Cl2CO collocation directly from B3LYP/aug-cc-pVDZ DFT data (option 4 of 13 Nov 2026).
Same procedure as the PES tuning (PES max-residual tuning; widths in results/best_widths_pes_max.json): J = MAXIMUM over the 100 computed levels of the
UNWEIGHTED validation residual r_i = ||(H - E_i) psi_i|| / ||psi_i|| (energy rows, cm-1) on the DFT validation points.
Stage 1: isotropic beta0 in {0.80, 0.85, 0.90, 0.95, 1.00} (extended by 0.05 if the minimum is at an end).
Stage 2: coordinate descent over k = 1..6 (nu3, nu6, nu2, nu4, nu5, nu1), alpha_k x 0.95 / x 1.05, accept if J decreases, sweeps until
         a sweep accepts nothing (at most 5).  Stage 3: the same with x 0.99 / x 1.01 (at most 20 sweeps).
Differences from the PES tuning (the PES setting N = 5000, M = 15000, energy rows, 2000 validation points, is not possible with the
DFT data): each trial is the production DFT run itself, run_dft.py adz 5000 5000 3 "5000" (N = 5000 centres, M = 5000 DFT collocation
points with gradient rows, weights (1 - V/Vmax)^3, Vmax = 5000 cm-1, Ec = 1150 cm-1, nev = 100, sigma = 2250 cm-1), validation on the
1000 DFT validation points (R6 index from 10,000,000). No reference levels exist or are used.
Each trial runs as a separate process with OMP_NUM_THREADS=1; the trials of a step (and the stage-1 grid) run two at a time.
Trials are cached in runs/<tag>.json, so the script can be restarted and re-reads finished trials."""
import numpy as np, json, os, subprocess
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); os.makedirs(os.path.join(HERE, 'runs'), exist_ok=True)
TAG = os.environ.get('TAG', 'adz')
LOG = os.path.join(HERE, f'tune_dft_{TAG}_max.log')
def log(s): print(s, flush=True); open(LOG, 'a').write(s + '\n')
def trial(beta):
    tag = f'runs/tune_{TAG}_' + '-'.join(f'{b:.6f}' for b in beta)
    fn = os.path.join(HERE, tag + '.json')
    if not os.path.exists(fn):
        env = dict(os.environ, WIDTHS=','.join(repr(float(b)) for b in beta), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        subprocess.run(['python3', 'run_dft.py', TAG, '5000', '5000', '3', '5000', tag], cwd=HERE, env=env, check=True,
                       stdout=subprocess.DEVNULL)
        for x in ('_C.npy', '_E.npy'):
            if os.path.exists(os.path.join(HERE, tag + x)): os.remove(os.path.join(HERE, tag + x))
    return json.load(open(fn))['recs'][-1]
def J(beta):
    r = trial(beta); res = np.array(r['res'])
    log(f"beta={[round(b, 6) for b in beta]} J={res.max():.4f} | res med {np.median(res):.3f} med20 {np.median(res[:20]):.3f} | maxIm {r['maxim']:.2f} | ZPE {r['zpe']:.2f} | t {r['t']:.0f}s")
    return float(res.max())
def Jmany(betas):
    with ThreadPoolExecutor(2) as ex: return list(ex.map(J, betas))
log(f'=== tune DFT {TAG}: N=5000, M=5000 DFT points with gradient rows (p=3), Vmax=5000, Ec=1150, validation = DFT validation points; J = MAXIMUM UNWEIGHTED validation residual over 100 levels')
grid = [0.80, 0.85, 0.90, 0.95, 1.00]; Js = dict(zip(grid, Jmany([[b] * 6 for b in grid])))
while True:
    b = min(Js, key=Js.get)
    if b == min(Js): nb = round(b - 0.05, 2)
    elif b == max(Js): nb = round(b + 0.05, 2)
    else: break
    Js[nb] = J([nb] * 6)
    if Js[nb] >= Js[b]: break
b0 = min(Js, key=Js.get); Jc = Js[b0]; alpha = [1.0] * 6
log(f'stage 1: J(beta0) = { {k: round(v, 4) for k, v in sorted(Js.items())} } -> beta0 = {b0}')
def descend(alpha, Jc, step, maxsw, name):
    for sweep in range(1, maxsw + 1):
        changed = False
        for k in range(6):
            A2 = []
            for f in (1 - step, 1 + step):
                a = list(alpha); a[k] *= f; A2.append((a, f))
            js = Jmany([[b0 * x for x in a] for a, f in A2])
            trials = [(j, a, f) for j, (a, f) in zip(js, A2)]
            jt, a, f = min(trials, key=lambda t: t[0])
            if jt < Jc:
                alpha, Jc, changed = a, jt, True; log(f'  {name} sweep {sweep}: accepted alpha[{k + 1}] x{f} -> {a[k]:.6f}, J = {Jc:.4f}')
            else:
                log(f'  {name} sweep {sweep}: alpha[{k + 1}] unchanged (J x{1 - step} = {trials[0][0]:.4f}, x{1 + step} = {trials[1][0]:.4f}, current {Jc:.4f})')
        if not changed: break
    return alpha, Jc, sweep
alpha, Jc, sw2 = descend(alpha, Jc, 0.05, 5, 'stage 2 (5%)')
log(f'stage 2 result: beta = {[round(b0 * x, 6) for x in alpha]}, J = {Jc:.4f}')
alpha, Jc, sweep = descend(alpha, Jc, 0.01, 20, 'stage 3 (1%)')
beta = [b0 * x for x in alpha]
json.dump(dict(beta0=b0, alpha=alpha, beta=beta, J=Jc, sweeps=sweep), open(os.path.join(HERE, f'best_widths_dft_{TAG}_max.json'), 'w'))
log(f'=== result: beta0 = {b0}, alpha = {[round(x, 6) for x in alpha]}, beta = {[round(x, 6) for x in beta]}, J = {Jc:.4f}')
