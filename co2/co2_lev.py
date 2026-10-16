"""Free CO2: many levels by rectangular collocation, blocks of vibrational angular momentum l = 0 ... 4.
Block l is the J = l block spanned by  F_l(X) phi_j(x1, x3, q),  F_l = Re or Im of (p_Z + i n_Z)^l  (l even: Re, l odd: Im), with
n = a^ x b^ and p = (a^ + b^) sin(theta/2) (a^, b^ unit C->O vectors; n, p are orthogonal, of equal length sin(theta), and both
perpendicular to the molecular axis), i.e. the ZZ..Z component of the rank-l harmonic tensor with |K| = l.  For l = 1 the
unnormalised (a x b)_Z of co2c.py is kept.  Band origin = E(J = l) - E0 - Be l   (rotational energy Be [J (J + 1) - l^2]).
Radial functions: Gaussians of co2_ffbench.py (GL0) or the pruned HO products of co2ho.py.  Energy rows, optionally + gradient rows.
  python co2_lev.py ff|dft <nproc>"""
import sys, json, time, numpy as np
from multiprocessing import Pool
import co2c as c, co2ho as ho, co2_ffbench as fb, co2_grad as gr
H = c.H2CM; PGR = 3; VMAX = 8000.0            # gradient rows weighted by lambda = (1 - V / Vmax)^p, p = 3
NEV = {0: 15, 1: 10, 2: 9, 3: 6, 4: 4}
def ang(X, l):
    a = X[:, 1] - X[:, 0]; b = X[:, 2] - X[:, 0]
    if l == 1: return a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    a = a / np.linalg.norm(a, axis=1)[:, None]; b = b / np.linalg.norm(b, axis=1)[:, None]; s = np.sqrt(np.clip((1 - (a * b).sum(1)) / 2, 0, None))
    z = ((a[:, 2] + b[:, 2]) * s + 1j * (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0])) ** l; return z.real if l % 2 == 0 else z.imag
class GL(fb.GL0):
    def __init__(s, C, beta, l): fb.GL0.__init__(s, C, beta, 0); s.l = l
    def __call__(s, X): B = fb.GL0.__call__(s, X); return B if s.l == 0 else ang(X, s.l)[:, None] * B
class HOl(ho.HO):
    def __init__(s, ecut, l): ho.HO.__init__(s, ecut, l); s.l = l; s.J = 0
    def __call__(s, X): B = ho.HO.__call__(s, X); return B if s.l == 0 else ang(X, s.l)[:, None] * B
def data(src, M, grad):
    if src == 'ff':
        x, V = fb.points(1, M, 8000.0); xv, Ev = fb.points(5_000_000, 3000, 7000.0); return x, V, (gr.gradV_ff(x) if grad else None), xv, Ev
    xv, Ev = ho.load('val', 7000)
    if grad: d = np.load(f'{c.HERE}/pts_grad.npz'); return d['x'][:M], d['E'][:M], d['gu'][:M], xv, Ev
    x, V = ho.load('fit', 8000); return x[:M], V[:M], None, xv, Ev
def job(a):
    ci, src, kind, par, M, grad, l = a; t0 = time.time(); x, V, GV, xv, Ev = data(src, M, grad)
    bas = GL(fb.centres(par[0], par[1]), par[2], l) if kind == 'g' else HOl(par[0], l)
    S, Hm = gr.accumulate(bas, x, V, GV, np.clip(1 - V * H / VMAX, 0, None) ** PGR); E, im, res, kept = gr.finish(bas, S, Hm, xv, Ev, NEV[l], 1e-12 if kind == 'g' else 1e-13)
    return dict(ci=ci, src=src, kind=kind, par=par, M=len(x), grad=grad, l=l, N=bas.N, kept=kept, E=E.tolist(), res=res.tolist(), im=im, t=time.time() - t0)
def run(src, cfgs, nproc, out):
    jobs = [(ci, src, k, p, M, g, l) for ci, (k, p, M, g) in enumerate(cfgs) for l in range(5)]; jobs.sort(key=lambda j: -(j[3][0] if j[2] == 'g' else j[3][0] / 20) * (4 if j[5] else 1))
    with Pool(nproc) as pool: R = pool.map(job, jobs, chunksize=1)
    runs = []
    for ci, (k, p, M, g) in enumerate(cfgs):
        rr = {r['l']: r for r in R if r['ci'] == ci}; E0 = rr[0]['E'][0]
        lev = {l: [e - E0 - (c.BE * H * l) for e in rr[l]['E']][(1 if l == 0 else 0):] for l in rr}; res = {l: rr[l]['res'][(1 if l == 0 else 0):] for l in rr}
        tag = ('Gauss N=%d Ec=%g beta=%s' % (p[0], p[1], list(p[2])) if k == 'g' else 'HO Ecut=%g' % p[0]) + ' M=%d %s' % (rr[0]['M'], 'E+grad p=%d' % PGR if g else 'E')
        runs.append(dict(tag=tag, kind=k, par=p, M=rr[0]['M'], grad=g, N=[rr[l]['N'] for l in range(5)], kept=[rr[l]['kept'] for l in range(5)], zpe=E0, lev=lev, res=res, im=max(rr[l]['im'] for l in rr), t=sum(rr[l]['t'] for l in rr)))
        print('%s | N %s kept %s | ZPE %.3f | maxIm %.1e | %.0f s' % (tag, runs[-1]['N'], runs[-1]['kept'], E0, runs[-1]['im'], runs[-1]['t']))
        for l in range(5): print('   l=%d: ' % l + ' '.join('%.2f(%.1f)' % (e, r) for e, r in zip(lev[l], res[l])), flush=True)
    json.dump(runs, open(f'{c.HERE}/{out}', 'w'), indent=1)
if __name__ == '__main__':
    mode, nproc = sys.argv[1], int(sys.argv[2]); b = (1.0, 1.0, 0.8)
    if mode == 'test': run('ff', [('g', (600, 4000.0, b), 7200, False), ('h', (17000.0,), 20000, False)], nproc, 'lev_test.json')
    elif mode == 'ff':
        run('ff', [('g', (600, 4000.0, b), 7200, False), ('g', (1200, 6000.0, b), 14400, False), ('g', (2400, 8000.0, b), 19200, False),
                   ('h', (25000.0,), 40000, False), ('h', (33000.0,), 40000, False), ('h', (41000.0,), 40000, False),
                   ('g', (600, 4000.0, b), 600, True), ('g', (600, 4000.0, b), 2400, True), ('g', (1200, 6000.0, b), 1200, True), ('g', (1200, 6000.0, b), 3600, True)], nproc, 'lev_ff.json')
    elif mode == 'dft':                      # Gaussian basis tuned on the fitted quartic field (co2_tune.py): widths 0.8 / 0.8 / 0.65, centres up to 8000 cm-1
        g = (0.8, 0.8, 0.65); run('dft', [('g', (n, 8000.0, g), 100000, False) for n in (600, 1200, 2400, 4000)] + [('g', (2400, 8000.0, g), m, False) for m in (2400, 4800)], nproc, 'lev_dft2.json')
    elif mode == 'dftgrad4k':
        g = (0.8, 0.8, 0.65); run('dft', [('g', (4000, 8000.0, g), 2400, True)], nproc, 'lev_dftgrad4k.json')
    elif mode == 'dftgrad':
        g = (0.8, 0.8, 0.65); run('dft', [('g', (n, 8000.0, g), m, True) for n, m in ((1200, 1200), (1200, 2400), (2400, 1200), (2400, 2400))], nproc, 'lev_dftgrad2.json')
