"""B3LYP energies, gradients (in collocation coordinates q) and dipoles at collocation and validation points.

usage:  NT=2 python3 colloc_dft.py <basis> <tag>                       # single process, all points in order
        NT=2 SHARD=k/n python3 colloc_dft.py <basis> <tag>             # process k of n (k = 0..n-1), round-robin split
env:    NT     threads per process (default 2)
        MAXC   number of collocation points to compute (default 15000; prefix of pts_1_5000_15000.npz)
        MAXV   number of validation points (default 1000; prefix of pts_10000000_5000_2000.npz)
        MAXMEM PySCF max_memory per process in MB (default 4000)
        XC     functional (default b3lyp); GRID PySCF grid level (default 3)

Order (single process): collocation 0-2500, validation 0-500, collocation 2500-5000, validation 500-1000, collocation 5000+.
Every finished point is appended to dft_{tag}_cache[_k].jsonl, so runs are resumable and shards never clash.
Merge all caches into dft_{tag}_c.npz / dft_{tag}_v.npz with:  python3 rebuild_npz.py <tag> <basis>
Points are the nested R2 sets used for the NN-PES runs (selected with V_NN <= 5000 cm-1).
Workstation example (48 cores, aTZ, M = 5000 + 1000 validation):
    for k in $(seq 0 23); do NT=2 SHARD=$k/24 MAXC=5000 MAXV=1000 nohup python3 colloc_dft.py aug-cc-pvtz atz > log_atz_$k.txt 2>&1 & done
"""
import numpy as np, json, os, sys, time, glob
from pyscf import gto, dft, lib
import coords6 as c
BAS, TAG = sys.argv[1], sys.argv[2]
lib.num_threads(int(os.environ.get('NT', 2)))
XC = os.environ.get('XC', 'b3lyp'); GRID = int(os.environ.get('GRID', 3))
MAXC = int(os.environ.get('MAXC', 15000)); MAXV = int(os.environ.get('MAXV', 1000))
sh = os.environ.get('SHARD'); K, NS = (map(int, sh.split('/')) if sh else (0, 1))
SYM = ['C', 'O', 'Cl', 'Cl']
QC = np.load('pts_1_5000_15000.npz')['q']; QV = np.load('pts_10000000_5000_2000.npz')['q']
order = [('c', i) for i in range(min(2500, MAXC))] + [('v', i) for i in range(min(500, MAXV))] \
      + [('c', i) for i in range(2500, min(5000, MAXC))] + [('v', i) for i in range(500, MAXV)] \
      + [('c', i) for i in range(5000, MAXC)]
order = order[K::NS]
CACHE = f'dft_{TAG}_cache.jsonl' if NS == 1 else f'dft_{TAG}_cache_{K}.jsonl'
done = {}
for fn in glob.glob(f'dft_{TAG}_cache*.jsonl'):                 # skip points done by any earlier run/shard
    for line in open(fn):
        try: r = json.loads(line); done[(r['s'], r['i'])] = r
        except Exception: pass
if os.path.exists(CACHE) and open(CACHE, 'rb').read()[-1:] not in (b'', b'\n'):
    open(CACHE, 'a').write('\n')                                # guard against a truncated last line after a crash
fh = open(CACHE, 'a'); dm = [None]; t0 = time.time(); m = 0
todo = [p for p in order if p not in done]
print(f'shard {K}/{NS}: {len(order)} points assigned, {len(todo)} to do', flush=True)
for s, i in todo:
    q = (QC if s == 'c' else QV)[i]; X = c.X_of_q(q)[0]
    mol = gto.M(atom=[(a, tuple(x)) for a, x in zip(SYM, X)], basis=BAS, unit='B', verbose=0, max_memory=int(os.environ.get('MAXMEM', 4000)))
    mf = dft.RKS(mol); mf.xc = XC; mf.grids.level = GRID; mf.conv_tol = 1e-10; mf.conv_tol_grad = 1e-6
    mf.kernel(dm0=dm[0])
    if not mf.converged: mf.kernel()
    if not mf.converged: print('NOT CONVERGED', s, i, flush=True); continue
    dm[0] = mf.make_rdm1()
    gr = mf.nuc_grad_method(); gr.grid_response = True; G = gr.kernel()
    r = dict(s=s, i=int(i), E=float(mf.e_tot), g=(G.ravel() @ c.A).tolist(), mu=list(mf.dip_moment(unit='AU', verbose=0)))
    fh.write(json.dumps(r) + '\n'); fh.flush(); m += 1
    if m % 25 == 0: print('%d/%d done  %.1f s/pt  last %s%d' % (m, len(todo), (time.time() - t0) / m, s, i), flush=True)
print('finished shard', K, flush=True)
