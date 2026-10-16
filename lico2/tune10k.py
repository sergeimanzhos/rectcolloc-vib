"""Parameter optimisation at fixed N on a FIXED point set (first --M fit points by index, first --nvalraw validation records),
judged only by the median validation residual over all nev states (val points with V <= vval). No reference levels.
Start: best of the iscan runs (p, Ec); then coordinate descent over the six Gaussian widths (x0.95 / x1.05, both tried
concurrently), --sweeps times. Log: tune10k.log; result: results/best10k.json. Finished trials are reused."""
import os, sys, json, subprocess, argparse, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument('--N', type=int, default=10000); ap.add_argument('--M', type=int, default=9000); ap.add_argument('--nev', type=int, default=40)
ap.add_argument('--workers', type=int, default=4); ap.add_argument('--vval', type=float, default=4500); ap.add_argument('--nvalraw', type=int, default=1920)
ap.add_argument('--sweeps', type=int, default=2); ap.add_argument('--tol', type=float, default=0.02); ap.add_argument('--cache', default='cache_adz')
a = ap.parse_args()
LOG = os.path.join(HERE, 'tune10k.log')
def log(s): print(s, flush=True); open(LOG, 'a').write(s + '\n')
def tag(P): return f"t10k_N{a.N}_M{a.M}_ec{P['ec']:g}_v{P['vmax']:g}_p{P['p']:g}_b" + '-'.join(f'{b:.4f}' for b in P['beta'])
def launch(P):
    fn = os.path.join(HERE, 'results', tag(P) + '.json')
    if os.path.exists(fn): return None
    cmd = [sys.executable, os.path.join(HERE, 'run_rc.py'), '--pes', 'dft', '--cache', a.cache, '--N', str(a.N), '--ec', str(P['ec']),
           '--beta', ','.join(f'{b:.4f}' for b in P['beta']), '--vmax', str(P['vmax']), '--p', str(P['p']), '--Ms', str(a.M), '--nev', str(a.nev),
           '--workers', str(a.workers), '--threads', '1', '--vval', str(a.vval), '--nvalraw', str(a.nvalraw), '--tag', tag(P)]
    return subprocess.Popen(cmd, cwd=HERE, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
def result(P, proc):
    if proc is not None:
        err = proc.communicate()[1]
        if proc.returncode: log(f'FAILED {tag(P)}: {err[-300:]}'); return np.inf
    r = json.load(open(os.path.join(HERE, 'results', tag(P) + '.json')))['recs'][-1]
    J = r['res_med'] if r['maxIm'] < 5 else np.inf
    log(f"{tag(P)} | M={r['M']} ZPE={r['zpe']:.2f} res med20 {r['res_med20']:.2f} all {r['res_med']:.2f} max {r['res_max']:.1f} | maxIm {r['maxIm']:.1e} cond {r['condS']:.1e} | {r['t']:.0f}s -> J={J:.3f}")
    return J
def trials(Ps):
    procs = [launch(P) for P in Ps]; return [result(P, pr) for P, pr in zip(Ps, procs)]
# start: p and Ec from the iscan runs (same N, M, validation set)
P = dict(ec=900.0, vmax=5500.0, p=1.0, beta=[0.99] * 6)
cands = [P, dict(P, p=0.0), dict(P, ec=1300.0)]
Js = trials(cands); best = Js[0]
log(f'=== tune10k N={a.N} M={a.M}: start candidates J = {np.round(Js, 2).tolist()}')
if Js[1] < best * (1 - a.tol): P, best = cands[1], Js[1]
if Js[2] < best * (1 - a.tol):
    Q = dict(P, ec=1300.0); J = trials([Q])[0]
    if J < best * (1 - a.tol): P, best = Q, J
log(f'start {P} J={best:.3f}')
for sweep in range(a.sweeps):
    changed = False
    for k in range(6):
        Qs = []
        for f in (0.95, 1.05):
            b = list(P['beta']); b[k] *= f; Qs.append(dict(P, beta=b))
        Js = trials(Qs); i = int(np.argmin(Js))
        if Js[i] < best * (1 - a.tol): P, best, changed = Qs[i], Js[i], True; log(f'  accepted beta[{k}] -> {P["beta"][k]:.4f}, J={best:.3f}')
    if not changed: break
P['J'] = best; json.dump(P, open(os.path.join(HERE, 'results', 'best10k.json'), 'w')); log(f'=== best10k: {P}')
