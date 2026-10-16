"""Merge all dft_{tag}_cache*.jsonl files into dft_{tag}_c.npz (collocation) and dft_{tag}_v.npz (validation).
usage: python3 rebuild_npz.py <tag> <basis>"""
import numpy as np, json, sys, glob
tag, bas = sys.argv[1], sys.argv[2]
QC = np.load('pts_1_5000_15000.npz')['q']; QV = np.load('pts_10000000_5000_2000.npz')['q']
done = {}
for fn in sorted(glob.glob(f'dft_{tag}_cache*.jsonl')):
    for line in open(fn):
        try: r = json.loads(line); done[(r['s'], r['i'])] = r
        except Exception: pass
for s, Q in (('c', QC), ('v', QV)):
    n = len(Q); E = np.full(n, np.nan); g = np.full((n, 6), np.nan); mu = np.full((n, 3), np.nan)
    for (ss, i), r in done.items():
        if ss == s: E[i] = r['E']; g[i] = r['g']; mu[i] = r['mu']
    ok = ~np.isnan(E); pre = int(np.argmin(ok)) if not ok.all() else n
    np.savez(f'dft_{tag}_{s}.npz', q=Q, E=E, g=g, mu=mu, done=ok, basis=bas)
    print(f'{s}: {int(ok.sum())} points done, contiguous prefix {pre}')
