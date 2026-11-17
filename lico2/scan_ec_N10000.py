"""Ec scan for LiCO2 at N = 10000 (14 Nov 2026), after the N = 5000 scan was judged unreliable.
Each trial: run_rc.py --pes dft --cache data/cache_adz --N 10000 --Ms 12000 (first 12000 fit points by index with V <= 5500,
energy + gradient rows, weights (1 - V/Vmax)^p with p = 0) --vmax 5500 --p 0 --beta 1.0 x 6 --nev 100 --nvalraw 3920 --vval 5000,
Ec (centre cutoff, surrogate V <= Ec) in {700, 900, 1200, 1500}, then (user request) also {800, 1000, 1100, 1300}.  Two trials at a time, each 1 process, 1 thread (Ec 700, 1200); the Ec 900 trial was killed (out of memory, 3.2 GB each),
so the remaining trials (900, 1500, and 800, 1000, 1100, 1300) were run one at a time with 1 process, 2 threads.
Reported per Ec: max, 2nd largest and median of the 100 validation residuals, med20, maxIm, ZPE, highest level, time.
Results cached in results/<tag>.json; restartable."""
import numpy as np, json, os, subprocess
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, 'scan_ec_N10000.log')
N, M, NEV, NVALRAW, VVAL, CACHE, VMAX, P, BETA = 10000, 12000, 100, 3920, 5000, 'data/cache_adz', 5500, 0, [1.0] * 6
def log(s): print(s, flush=True); open(LOG, 'a').write(s + '\n')
def tag(ec): return f"tmax_N{N}_M{M}_ec{ec:g}_v{VMAX:g}_p{P:g}_b" + '-'.join(f"{b:.6f}" for b in BETA)
def trial(ec):
    t = tag(ec); fn = os.path.join(HERE, 'results', t + '.json')
    if not os.path.exists(fn):
        env = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2')
        cmd = ['python3', 'run_rc.py', '--pes', 'dft', '--cache', CACHE, '--N', str(N), '--ec', f"{ec:g}",
               '--beta', ','.join(repr(float(b)) for b in BETA), '--vmax', f"{VMAX:g}", '--p', f"{P:g}", '--Ms', str(M),
               '--nev', str(NEV), '--workers', '1', '--threads', '2', '--vval', str(VVAL), '--nvalraw', str(NVALRAW), '--batch', '1500', '--tag', t]
        subprocess.run(cmd, cwd=HERE, env=env, check=True, stdout=subprocess.DEVNULL)
        v = os.path.join(HERE, 'results', f'vec_{t}_M{M}.npz')
        if os.path.exists(v): os.remove(v)
    r = json.load(open(fn))['recs'][-1]; res = np.array(r['res']); tr = np.array(r['trans'])
    log(f"Ec={ec:g}: max {res.max():.2f} (level #{res.argmax()}, {tr[res.argmax()]:.1f}) | 2nd {np.sort(res)[-2]:.2f} | med {np.median(res):.2f} | "
        f"med20 {np.median(res[:20]):.2f} | n>150 {(res > 150).sum()} | maxIm {r['maxIm']:.2f} | ZPE {r['zpe']:.2f} | top level {tr[-1]:.1f} | t {r['t']:.0f}s")
log(f'=== Ec scan LiCO2 aDZ: N={N}, M={M}, Vmax={VMAX}, p={P}, beta={BETA}, nev={NEV}, validation first {NVALRAW} records with V <= {VVAL}')
with ThreadPoolExecutor(1) as ex: list(ex.map(trial, [900.0, 700.0, 1200.0, 1500.0, 800.0, 1000.0, 1100.0, 1300.0]))
log('=== done')
