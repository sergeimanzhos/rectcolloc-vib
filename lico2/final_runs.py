"""Final LiCO2 collocation runs with the current best (max-residual, N = 10000) parameters, 17 Nov 2026:
Ec = 1200, Vmax = 5000, p = 0, beta = (0.932065, 0.99, 1.0, 1.018814, 1.0, 1.0) (state of tune_li_max10k.py when stopped, sweep 7).
Same set of runs as for the manuscript version (f100_*): per basis N = 14000 and N = 10000 on the first point set (cache_adz / cache_atz)
and on the second, independent point set (cache_adz7 / cache_atz7); same --Ms and --nvalraw as the f100_* runs; Vmax = 5000 for all.
Intensities (intens.py, region V <= 5000) for the first-set N = 14000 runs, as for the manuscript version.
The aDZ first-set N = 14000 run (fmax_adz_N14000) was started separately with the same arguments.
One step at a time (1 process, 2 threads; two N = 10000 runs at once exceeded the memory). Restartable: finished steps are skipped.
Exits (code 3) instead of starting a step that would not finish inside the 7000 s window (estimates below); relaunch to continue."""
import os, sys, json, time, subprocess
HERE = os.path.dirname(os.path.abspath(__file__)); T0 = time.time()
BETA = '0.932065,0.99,1.0,1.018814,1.0,1.0'
ENV = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2')
STEPS = [  # ('run', tag, cache, N, Ms, nvalraw, est s)  /  ('int', tag, cache, est s)
    ('run', 'fmax_adz_N14000', 'cache_adz', 14000, '1000000', 3920, 5400),
    ('run', 'fmax_atz_N14000', 'cache_atz', 14000, '1000000', 0, 6000),
    ('int', 'fmax_adz_N14000', 'cache_adz', 3600),
    ('int', 'fmax_atz_N14000', 'cache_atz', 3600),
    ('run', 'fmax_adz_N10000', 'cache_adz', 10000, '1000000', 3920, 2700),
    ('run', 'fmax_atz_N10000', 'cache_atz', 10000, '1000000', 0, 4500),
    ('run', 'fmax_adz7_N10000', 'cache_adz7', 10000, '1000000', 6900, 4500),
    ('run', 'fmax_adz7_N14000', 'cache_adz7', 14000, '1000000', 6900, 6000),
    ('run', 'fmax_atz7_N10000', 'cache_atz7', 10000, '16000', 1735, 4500),
    ('run', 'fmax_atz7_N14000', 'cache_atz7', 14000, '16000', 1735, 6000)]
for st in STEPS:
    kind, tag, cache, est = st[0], st[1], st[2], st[-1]
    out = os.path.join(HERE, 'results', (tag if kind == 'run' else f'intens_{tag}') + '.json')
    if os.path.exists(out): continue
    if time.time() - T0 + est > 7000: print('window nearly used up; exiting before', kind, tag, flush=True); sys.exit(3)
    print(time.strftime('%H:%M'), 'start', kind, tag, flush=True)
    if kind == 'run':
        N, Ms, nval = st[3], st[4], st[5]
        cmd = ['python3', 'run_rc.py', '--pes', 'dft', '--cache', f'data/{cache}', '--N', str(N), '--ec', '1200', '--beta', BETA, '--vmax', '5000',
               '--p', '0', '--Ms', Ms, '--nev', '100', '--workers', '1', '--threads', '2', '--vval', '5000', '--nvalraw', str(nval), '--batch', '1500', '--tag', tag]
        log = tag + '.out'
    else:
        M = json.load(open(os.path.join(HERE, 'results', tag + '.json')))['recs'][-1]['M']
        cmd = ['python3', 'intens.py', f'data/{cache}', f'results/vec_{tag}_M{M}.npz', tag, '5000']; log = f'intens_{tag}.out'
    with open(os.path.join(HERE, log), 'w') as f: subprocess.run(cmd, cwd=HERE, env=ENV, check=True, stdout=f, stderr=subprocess.STDOUT)
    if kind == 'run':
        r = json.load(open(out))['recs'][-1]
        print(time.strftime('%H:%M'), 'done', tag, 'M', r['M'], 'ZPE %.2f' % r['zpe'], 'res med %.1f max %.1f' % (r['res_med'], r['res_max']), 't %.0f' % r['t'], flush=True)
    else: print(time.strftime('%H:%M'), 'done intensities', tag, flush=True)
print('all done', flush=True)
