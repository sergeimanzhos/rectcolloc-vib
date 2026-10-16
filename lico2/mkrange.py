"""Add the convergence range of each level to a collocation intensity file: python mkrange.py ir_colloc_<tag>.json out.json run1 run2 ...
(run1 must be the run the intensity file was made from; Elo/Ehi = min/max of each level over the listed runs, matched by energy order)."""
import sys, json, numpy as np
src, out, runs = sys.argv[1], sys.argv[2], sys.argv[3:]
T = np.array([json.load(open(f'results/{t}.json'))['recs'][-1]['trans'] for t in runs]); j = json.load(open(src))
for k, r in enumerate(j['rows'], 1):
    assert abs(r['E'] - T[0, k]) < 0.01; r['Elo'], r['Ehi'] = float(T[:, k].min()), float(T[:, k].max())
json.dump(j, open(out, 'w'), indent=1)
w = np.array([r['Ehi'] - r['Elo'] for r in j['rows']]); print('range over %d runs, width (cm-1): median %.2f, 90th percentile %.2f, max %.2f' % (len(runs), np.median(w), np.percentile(w, 90), w.max()))
