"""Polynomial dipole surface mu(q) from DFT collocation points; checked on DFT validation points.
usage: python3 fit_dip.py tag deg  ->  dipfit_{tag}_d{deg}.npz (mono, coef) for ir6.py (DIP=...)"""
import numpy as np, itertools, sys
tag, deg = sys.argv[1], int(sys.argv[2])
dc = np.load(f'dft_{tag}_c.npz'); dv = np.load(f'dft_{tag}_v.npz')
q, mu = dc['q'][dc['done']], dc['mu'][dc['done']]; qv, muv = dv['q'][dv['done']], dv['mu'][dv['done']]
if len(qv) == 0: n = int(0.8 * len(q)); qv, muv, q, mu = q[n:], mu[n:], q[:n], mu[:n]
MONO = np.array([a for a in itertools.product(range(deg + 1), repeat=6) if sum(a) <= deg])
F = lambda x: np.prod(x[:, None, :] ** MONO[None], axis=2)
coef, *_ = np.linalg.lstsq(F(q), mu, rcond=None)
rt = np.sqrt(np.mean((F(q) @ coef - mu) ** 2, 0)); rv = np.sqrt(np.mean((F(qv) @ coef - muv) ** 2, 0))
lin = [coef[[i for i, a in enumerate(MONO) if a.sum() == 1 and a[k] == 1][0]] for k in range(6)]
print(f'deg {deg}: {len(MONO)} terms, {len(q)} fit pts, {len(qv)} val pts; rms fit {rt} val {rv} (au); |mu| spread {mu.std(0)}')
print('linear coefficients (dmu/dq):', np.round(lin, 5))
np.savez(f'dipfit_{tag}_d{deg}.npz', mono=MONO, coef=coef, rms_fit=rt, rms_val=rv)
