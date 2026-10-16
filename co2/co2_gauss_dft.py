"""Free CO2 on the B3LYP-D3(BJ)/aug-cc-pVTZ points with the Gaussian basis tuned on the quartic force field (co2_ffbench.py):
energy rows only, J = 0 and J = 1 f blocks.  python co2_gauss_dft.py"""
import json, numpy as np
from multiprocessing import Pool
import co2c as c, co2ho as ho, co2_ffbench as fb
H = c.H2CM
def run(a):
    N, ec, beta, vmax, M, thr = a; xf, Ef = ho.load('fit', vmax); xf, Ef = xf[:M], Ef[:M]; xv, Ev = ho.load('val', 6000); C = fb.centres(N, ec)
    E0, im0, r0, k0 = fb.solve(fb.GL0(C, beta, 0), xf, Ef, xv, Ev, 8, thr); E1, im1, r1, k1 = fb.solve(fb.GL0(C, beta, 1), xf, Ef, xv, Ev, 3, thr)
    return dict(N=N, ec=ec, beta=list(beta), vmax=vmax, M=len(xf), thr=thr, kept=(k0, k1), zpe=float(E0[0]), nu2=float(E1[0] - E0[0] - c.BE * H), J0=(E0[1:] - E0[0]).tolist(), res0=r0.tolist(), resF=r1.tolist(), im=max(im0, im1))
if __name__ == '__main__':
    b = (1.0, 1.0, 0.8); jobs = [(n, 4000.0, b, 8000.0, 100000, 1e-12) for n in (200, 300, 450, 600, 900, 1200)] + [(600, 4000.0, b, 8000.0, m, 1e-12) for m in (1800, 3600, 5400)] + \
        [(600, 4000.0, b, v, 100000, 1e-12) for v in (6500.0, 7250.0)] + [(600, ec, b, 8000.0, 100000, 1e-12) for ec in (3000.0, 5000.0)] + [(600, 4000.0, bb, 8000.0, 100000, 1e-12) for bb in ((0.9, 0.9, 0.8), (1.1, 1.1, 0.8), (1.0, 1.0, 0.7), (1.0, 1.0, 0.9))]
    with Pool(12) as pool: R = pool.map(run, jobs, chunksize=1)
    for r in R:
        print('N=%-5d Ec=%-5g beta=%s Vmax=%-5g M=%-5d kept=%s | ZPE %.3f nu2 %.3f dyad %.3f %.3f nu3 %.3f | res %s ; f %s | maxIm %.1e' % (r['N'], r['ec'], np.round(r['beta'], 2).tolist(), r['vmax'], r['M'], r['kept'], r['zpe'], r['nu2'],
              r['J0'][0], r['J0'][1], r['J0'][2], ' '.join('%.2f' % x for x in r['res0'][:4]), ' '.join('%.2f' % x for x in r['resF'][:1]), r['im']), flush=True)
    json.dump(R, open(f'{c.HERE}/co2_gauss_dft.json', 'w'), indent=1)
