"""Basis tuning for the 100-state run at fixed N on a FIXED point set (first --M fit points by index, first --nvalraw
validation records), judged only by the median validation residual over all nev states. No reference levels.
Stage 1: centre energy cutoff Ec (list).  Stage 2: common scale of the widths (walk in steps of 5%).
Stage 3: one coordinate-descent sweep over the six widths (x0.95 / x1.05).  Two trials run concurrently.
Log: tune100.log; result: results/best100.json. Finished trials are reused."""
import os, sys, json, subprocess, argparse, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument('--N', type=int, default=10000); ap.add_argument('--M', type=int, default=12000); ap.add_argument('--nev', type=int, default=100)
ap.add_argument('--workers', type=int, default=5); ap.add_argument('--vval', type=float, default=5000); ap.add_argument('--nvalraw', type=int, default=3920)
ap.add_argument('--vmax', type=float, default=5500); ap.add_argument('--p', type=float, default=0.0); ap.add_argument('--tol', type=float, default=0.02)
ap.add_argument('--ecs', default='900,1200,1500,1800'); ap.add_argument('--beta0', default='0.8935,1.0395,0.99,0.99,1.0395,0.99'); ap.add_argument('--cache', default='cache_adz')
a = ap.parse_args()
LOG = os.path.join(HERE, 'tune100.log')
def log(s): print(s, flush=True); open(LOG, 'a').write(s + '\n')
def tag(P): return f"t100_N{a.N}_M{a.M}_ec{P['ec']:g}_v{P['vmax']:g}_p{P['p']:g}_b" + '-'.join(f'{b:.4f}' for b in P['beta'])
def launch(P):
    if os.path.exists(os.path.join(HERE, 'results', tag(P) + '.json')): return None
    cmd = [sys.executable, os.path.join(HERE, 'run_rc.py'), '--pes', 'dft', '--cache', a.cache, '--N', str(a.N), '--ec', str(P['ec']),
           '--beta', ','.join(f'{b:.4f}' for b in P['beta']), '--vmax', str(P['vmax']), '--p', str(P['p']), '--Ms', str(a.M), '--nev', str(a.nev),
           '--workers', str(a.workers), '--threads', '1', '--vval', str(a.vval), '--nvalraw', str(a.nvalraw), '--tag', tag(P)]
    return subprocess.Popen(cmd, cwd=HERE, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
def result(P, proc):
    if proc is not None:
        err = proc.communicate()[1]
        if proc.returncode: log(f'FAILED {tag(P)}: {err[-300:]}'); return np.inf
    r = json.load(open(os.path.join(HERE, 'results', tag(P) + '.json')))['recs'][-1]; res = np.array(r['res'])
    J = r['res_med'] if r['maxIm'] < 5 else np.inf
    log(f"{tag(P)} | M={r['M']} ZPE={r['zpe']:.2f} res med 1-40 {np.median(res[:40]):.2f} 41-100 {np.median(res[40:]):.2f} all {r['res_med']:.2f} max {r['res_max']:.1f} | maxIm {r['maxIm']:.1e} cond {r['condS']:.1e} | {r['t']:.0f}s -> J={J:.3f}")
    return J
def trials(Ps):
    out = []
    for i in range(0, len(Ps), 2):
        procs = [launch(P) for P in Ps[i:i + 2]]; out += [result(P, pr) for P, pr in zip(Ps[i:i + 2], procs)]
    return out
b0 = [float(x) for x in a.beta0.split(',')]
log(f'=== tune100 N={a.N} M={a.M} nev={a.nev} vmax={a.vmax} p={a.p} vval={a.vval}')
# stage 1: Ec
cands = [dict(ec=float(e), vmax=a.vmax, p=a.p, beta=b0) for e in a.ecs.split(',')]
Js = trials(cands); i = int(np.argmin(Js)); P, best = cands[i], Js[i]; log(f'stage 1 (Ec): J = {np.round(Js, 2).tolist()} -> Ec = {P["ec"]:g}')
# stage 2: common width scale
Js = trials([dict(P, beta=[b * f for b in P['beta']]) for f in (0.95, 1.05)]); d = (0.95, 1.05)[int(np.argmin(Js))]
if min(Js) < best * (1 - a.tol):
    P, best = dict(P, beta=[b * d for b in P['beta']]), min(Js)
    for _ in range(3):
        Q = dict(P, beta=[b * d for b in P['beta']]); J = trials([Q])[0]
        if J < best * (1 - a.tol): P, best = Q, J
        else: break
log(f'stage 2 (width scale): beta = {np.round(P["beta"], 4).tolist()} J = {best:.3f}')
# stage 3: one sweep over the six widths
for k in range(6):
    Qs = []
    for f in (0.95, 1.05):
        b = list(P['beta']); b[k] *= f; Qs.append(dict(P, beta=b))
    Js = trials(Qs); i = int(np.argmin(Js))
    if Js[i] < best * (1 - a.tol): P, best = Qs[i], Js[i]; log(f'  accepted beta[{k}] -> {P["beta"][k]:.4f}, J={best:.3f}')
P['J'] = best; json.dump(P, open(os.path.join(HERE, 'results', 'best100.json'), 'w')); log(f'=== best100: {P}')
