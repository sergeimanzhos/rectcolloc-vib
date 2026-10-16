"""Free CO2: cylindrically symmetric force field in rectilinear dimensionless normal coordinates (q1 sym. stretch, q3 asym. stretch,
b2 = qa^2 + qb^2 bend) by least squares on the B3LYP-D3(BJ)/aug-cc-pVTZ energies of the collocation points (pts_fit.npz).
Symmetry-allowed monomials q1^i (q3^2)^j (b2)^k up to total degree DEG; the quadratic, cubic and quartic coefficients give the
quartic force field for GVPT2 / VCI (vff_fit.json, same format as vff_cyl.json).   python fit_ff.py"""
import json, numpy as np
import co2c as c, co2ho as ho, co2_ffbench as fb
H = c.H2CM; old = json.load(open(f'{c.HERE}/vff_cyl.json'))
def terms(deg): return [(i, j, k) for i in range(deg + 1) for j in range(deg // 2 + 1) for k in range(deg // 2 + 1) if 1 <= i + 2 * j + 2 * k <= deg]
def fit(vmax, deg):
    x, E = ho.load('fit', vmax); q1, q3, b2 = fb.q_rect(c.geom(x)); T = terms(deg); A = np.array([q1 ** i * q3 ** (2 * j) * b2 ** k for i, j, k in T]).T
    co, *_ = np.linalg.lstsq(A, E, rcond=None); r = A @ co - E; return dict(zip(T, co)), np.sqrt(np.mean(r ** 2)) * H, np.abs(r).max() * H, len(E)
def ff(co):
    w0 = np.array(old['w']); cq = np.array([co[(0, 0, 1)], co[(2, 0, 0)], co[(0, 1, 0)]]); w = np.sqrt(2 * cq * w0); s = np.sqrt(w0 / w); sb, s1, s3 = s     # rescale q -> q sqrt(w / w0)
    g = lambda i, j, k: co.get((i, j, k), 0.0) * s1 ** i * s3 ** (2 * j) * sb ** (2 * k)
    return dict(w=w.tolist(), k111=6 * g(3, 0, 0), k133=2 * g(1, 1, 0), k1bb=2 * g(1, 0, 1), p1111=24 * g(4, 0, 0), p3333=24 * g(0, 2, 0), p1133=4 * g(2, 1, 0), p11bb=4 * g(2, 0, 1), p33bb=4 * g(0, 1, 1), pbbbb=24 * g(0, 0, 2)), co.get((1, 0, 0), 0.0)
keys = ['k111', 'k133', 'k1bb', 'p1111', 'p3333', 'p1133', 'p11bb', 'p33bb', 'pbbbb']
print('%-22s %s' % ('', '  wb      w1      w3    ' + ' '.join('%8s' % k for k in keys) + '   rms    max   npts  lin(q1)'))
print('%-22s %s' % ('finite-difference FF', ' '.join('%7.2f' % (v * H) for v in old['w']) + ' ' + ' '.join('%8.2f' % (old[k] * H) for k in keys)))
best = None
for vmax in (1500, 2500, 4000, 6000, 8000):
    for deg in (4, 6, 8):
        co, rms, mx, n = fit(vmax, deg)
        if n < 3 * len(co): continue
        F, lin = ff(co); print('%-22s %s' % ('Vmax %d deg %d' % (vmax, deg), ' '.join('%7.2f' % (v * H) for v in F['w']) + ' ' + ' '.join('%8.2f' % (F[k] * H) for k in keys) + ' %6.3f %6.2f %5d %7.3f' % (rms, mx, n, lin * H)))
        if (vmax, deg) == (4000, 8): best = F
json.dump(best, open(f'{c.HERE}/vff_fit.json', 'w'), indent=1); print('saved vff_fit.json (Vmax 4000, degree 8)')
