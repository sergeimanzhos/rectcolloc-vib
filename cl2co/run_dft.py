"""Collocation directly from DFT points (energies + gradients from colloc_dft.py).
usage: python3 run_dft.py tag N vmax p "M1,M2,..." outtag
V = E_DFT - E_min(DFT); gradient-row weights (1 - V/vmax)^p; validation residual on the DFT validation points done so far.
No reference levels exist for the DFT surface: records ZPE, transition energies, residuals."""
import numpy as np, json, sys, time
from collo6 import Basis, accumulate, solve, residual, H2CM
H = H2CM
EMIN = {'adz': -1033.7840144379, 'atz': -1033.8413302173}     # B3LYP minima from vp_dft optimisations
tag, N, vmax, p, Ms, out = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), [int(x) for x in sys.argv[5].split(',')], sys.argv[6]
BASEW = np.array([float(x) for x in __import__('os').environ['WIDTHS'].split(',')]) if __import__('os').environ.get('WIDTHS') else np.array([0.857375, 0.857375, 0.95, 0.9025, 0.95, 0.95])
dc = np.load(f'dft_{tag}_c.npz'); dv = np.load(f'dft_{tag}_v.npz')
nc = int(np.argmin(dc['done'])) if not dc['done'].all() else len(dc['done'])     # contiguous prefix done
assert max(Ms) <= nc, f'only {nc} collocation points available'
qf = dc['q'][:max(Ms)]; Vf = dc['E'][:max(Ms)] - EMIN[tag]; gf = dc['g'][:max(Ms)]
mv = dv['done']; qv = dv['q'][mv]; Vv = dv['E'][mv] - EMIN[tag]
if mv.sum() == 0:   # no validation points yet: use collocation points beyond max(Ms)
    qv = dc['q'][max(Ms):nc]; Vv = dc['E'][max(Ms):nc] - EMIN[tag]; mv = np.ones(len(qv), bool)
print(f'{tag}: {nc} colloc points, {mv.sum()} validation points; V range {Vf.min()*H:.0f}..{Vf.max()*H:.0f} cm-1', flush=True)
lam = np.clip(1 - Vf * H / vmax, 0, 1) ** p
cq = np.load(f'pts_20000000_1150_{N}.npz')['q']; bas = Basis(cq, BASEW)
S = Hm = None; done = 0; t0 = time.time(); recs = []
for M in Ms:
    S, Hm = accumulate(bas, qf[done:M], Vf[done:M], gf[done:M], lam=lam[done:M], S=S, H=Hm); done = M
    E, Cv = solve(S, Hm, 100, 2250 / H); o = np.argsort(E.real); E, Cv = E[o], Cv[:, o]
    r = residual(bas, qv, Vv, E, Cv) if len(qv) else np.full(len(E), np.nan); Er = E.real * H
    rec = dict(M=M, zpe=Er[0], E=(Er - Er[0]).tolist(), res=r.tolist(), res_med=float(np.median(r)), maxim=float(np.abs(E.imag).max() * H), t=time.time() - t0)
    recs.append(rec); json.dump(dict(tag=tag, N=N, vmax=vmax, p=p, nval=int(mv.sum()), recs=recs), open(f'{out}.json', 'w'))
    print(f"{out} M={M} ZPE {Er[0]:.2f} | nu: {np.round(Er[1:8]-Er[0],2)} | res med {rec['res_med']:.1f} | maxIm {rec['maxim']:.2f} | {rec['t']:.0f}s", flush=True)
np.save(f'{out}_C.npy', Cv); np.save(f'{out}_E.npy', E)
