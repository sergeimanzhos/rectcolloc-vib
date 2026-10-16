"""Free CO2, Gaussian basis: tuning of widths / centre cutoff / N / M on the fitted cylindrical quartic force field (vff_fit.json),
where VCI (vptlev_co2fit_K14.json) is the exact reference for all levels with l = 0 ... 4 up to 4350 cm-1.
  python co2_tune.py <stage> <nproc>      stage = widths | nm | grad"""
import sys, json, time, numpy as np
from multiprocessing import Pool
import co2c as c, co2_ffbench as fb, co2_grad as gr, co2_lev as lv
fb.FF = json.load(open(f'{c.HERE}/vff_fit.json')); H = c.H2CM; EMAX = 4350.0
def vci_ref():
    vp = json.load(open(f'{c.HERE}/vptlev_co2fit_K14.json'))['vci']; out = {l: [] for l in range(5)}
    for r in vp:
        l = int(round(np.sqrt(max(r['l2'], 0))))
        if l <= 4 and not any(abs(e - r['E']) < 0.02 for e in out[l]): out[l].append(r['E'])
    return {l: sorted(v) for l, v in out.items()}
REF = vci_ref()
def job(a):
    ci, (N, ec, beta, M, grad, thr), l = a; t0 = time.time(); x, V = fb.points(1, M, 8000.0); xv, Ev = fb.points(5_000_000, 3000, 7000.0); bas = lv.GL(fb.centres(N, ec), beta, l)
    S, Hm = gr.accumulate(bas, x, V, gr.gradV_ff(x) if grad else None, np.clip(1 - V * H / 8000.0, 0, None) ** 3); E, im, res, kept = gr.finish(bas, S, Hm, xv, Ev, lv.NEV[l], thr)
    return dict(ci=ci, l=l, E=E.tolist(), res=res.tolist(), kept=kept, im=im, t=time.time() - t0)
def run(cfgs, nproc, out):
    jobs = [(ci, cf, l) for ci, cf in enumerate(cfgs) for l in range(5)]; jobs.sort(key=lambda j: -j[1][0] * j[1][3] * (13 if j[1][4] else 1))
    with Pool(nproc) as pool: R = pool.map(job, jobs, chunksize=1)
    rows = []
    for ci, (N, ec, beta, M, grad, thr) in enumerate(cfgs):
        rr = {r['l']: r for r in R if r['ci'] == ci}; E0 = rr[0]['E'][0]; lev = {l: [e - E0 - c.BE * H * l for e in rr[l]['E']][(1 if l == 0 else 0):] for l in rr}
        err = [(lev[l][i] - e, e) for l in range(5) for i, e in enumerate(REF[l]) if e <= EMAX and i < len(lev[l])]
        e2, e3, e4 = [max(abs(d) for d, e in err if e <= lim) for lim in (2000, 3000, EMAX)]; rms = float(np.sqrt(np.mean([d * d for d, e in err])))
        rows.append(dict(N=N, ec=ec, beta=list(beta), M=M, grad=grad, thr=thr, kept=[rr[l]['kept'] for l in range(5)], zpe=E0, lev=lev, err2000=e2, err3000=e3, errall=e4, rms=rms, n=len(err), im=max(rr[l]['im'] for l in rr),
                         res={l: rr[l]['res'] for l in rr}, t=sum(rr[l]['t'] for l in rr)))
        print('N=%-5d Ec=%-5g beta=%s M=%-6d %-6s thr=%.0e kept=%s | ZPE %.3f | max|err| vs VCI: <2000 %.3f  <3000 %.3f  <4350 %.3f  rms %.3f (%d levels) | maxIm %.0e | %.0f s' % (N, ec, list(beta), M, 'E+grad' if grad else 'E', thr,
              rows[-1]['kept'], E0, e2, e3, e4, rms, len(err), rows[-1]['im'], rows[-1]['t']), flush=True)
    json.dump(rows, open(f'{c.HERE}/{out}', 'w'), indent=1)
if __name__ == '__main__':
    stage, nproc = sys.argv[1], int(sys.argv[2])
    if stage == 'widths':
        N = 2400; run([(N, ec, (b, b, bb), 4 * N, False, 1e-12) for ec in (6000.0, 8000.0) for b in (0.6, 0.8, 1.0) for bb in (0.5, 0.65, 0.8)] + [(N, 8000.0, (0.8, 0.8, 0.65), 4 * N, False, t) for t in (1e-10, 1e-14)], nproc, 'tune_widths.json')
    elif stage == 'nm':
        ec, b, bb, thr = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]); B = (b, b, bb)
        run([(n, ec, B, 4 * n, False, thr) for n in (600, 1200, 2400, 4000, 6000)] + [(2400, ec, B, m, False, thr) for m in (2400, 4800, 19200)], nproc, 'tune_nm.json')
    elif stage == 'grad':
        ec, b, bb, thr = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]); B = (b, b, bb)
        run([(n, ec, B, m, True, thr) for n, m in ((1200, 1200), (1200, 2400), (2400, 1200), (2400, 2400), (2400, 4800), (4000, 2400))], nproc, 'tune_grad.json')
