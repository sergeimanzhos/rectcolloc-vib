"""Pre-generate and cache basis centres exactly as run_rc.py --pes dft does: python gen_centres.py <cache> <ec> <N>"""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import coords as c
from r2 import r2
from surr import polyfit
cache, ec, N = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]); H = c.H2CM
cfn = os.path.join(HERE, 'results', f"centres_dft_{os.path.basename(cache)}_ec{ec:g}_N{N}.npy")
if os.path.exists(cfn): sys.exit(0)
d = np.load(os.path.join(HERE, cache, 'fit.npz')); ok = d['conv'] & (d['E'] < 2.5 * 5500 / H)
Vs, rms = polyfit(d['q'][ok], d['E'][ok])
b = np.load(os.path.join(HERE, cache, 'box.npz')); lo, hi = b['lo'], b['hi']
Vc = lambda q: np.where(c.V_ff(q) <= 2.5 * ec / H, Vs(q), np.inf)
out, i, n = [], 20_000_000, 0
while n < N:
    q = lo + (hi - lo) * r2(np.arange(i, i + 200000), c.NQ); i += 200000; out.append(q[Vc(q) <= ec / H]); n += len(out[-1])
np.save(cfn + f'.{os.getpid()}.npy', np.concatenate(out)[:N]); os.replace(cfn + f'.{os.getpid()}.npy', cfn); print('saved', cfn, flush=True)
