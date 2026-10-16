"""Image of LiCO2 (C2v, bidentate) at an equilibrium geometry with its bond lengths and angles.  usage: python mol_lico2.py atz|adz out.png"""
import sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc
import coords as c
B = 0.529177210903; tag = sys.argv[1]
def geom(X):
    C, O1, O2, Li = X * B; d = lambda a, b: np.linalg.norm(a - b)
    ang = lambda a, b, cc: np.degrees(np.arccos(np.dot(a - b, cc - b) / d(a, b) / d(cc, b)))
    return dict(rCO=(d(C, O1) + d(C, O2)) / 2, rLiO=(d(Li, O1) + d(Li, O2)) / 2, rLiC=d(Li, C), aOCO=ang(O1, C, O2), aOLiO=ang(O1, Li, O2), asym=abs(d(C, O1) - d(C, O2)) + abs(d(Li, O1) - d(Li, O2)))
G = {t: c.X_of_q(np.load(f'cache_{t}/qmin.npy')[None])[0] for t in ('adz', 'atz')}
for t, X in G.items():
    g = geom(X); print('%s: r(C-O) %.4f A  r(Li-O) %.4f A  r(Li-C) %.4f A  O-C-O %.2f deg  O-Li-O %.2f deg  (asymmetry %.1e)' % (t, g['rCO'], g['rLiO'], g['rLiC'], g['aOCO'], g['aOLiO'], g['asym']))
g = geom(G[tag])
# planar drawing: C at the top, Li at the bottom on the C2 axis, O atoms left and right
h = np.radians(g['aOCO'] / 2)
P = {'C': np.array([0.0, 0.0]), 'O1': g['rCO'] * np.array([-np.sin(h), -np.cos(h)]), 'O2': g['rCO'] * np.array([np.sin(h), -np.cos(h)]), 'Li': np.array([0.0, -g['rLiC']])}
col = {'C': '#555555', 'O': '#d62728', 'Li': '#9467bd'}; rad = {'C': 0.30, 'O': 0.28, 'Li': 0.34}
plt.rcParams.update({'font.size': 20}); fig, ax = plt.subplots(figsize=(8, 7.6))
for a in ('O1', 'O2'):
    ax.plot(*zip(P['C'], P[a]), color='k', lw=5, zorder=1, solid_capstyle='round')
    ax.plot(*zip(P['Li'], P[a]), color='k', lw=3, ls=(0, (3, 2)), zorder=1)
for a, p in P.items():
    el = a.rstrip('12'); ax.add_patch(Circle(p, rad[el], facecolor=col[el], edgecolor='k', lw=1.5, zorder=2)); ax.text(p[0], p[1], el, ha='center', va='center', color='w', fontsize=22, fontweight='bold', zorder=3)
m = 0.5 * (P['C'] + P['O2']); ax.text(m[0] + 0.17, m[1] + 0.13, '%.3f Å' % g['rCO'], ha='left', va='center')
m = 0.5 * (P['Li'] + P['O2']); ax.text(m[0] + 0.17, m[1] - 0.10, '%.3f Å' % g['rLiO'], ha='left', va='center')
ra = 1.15; ax.add_patch(Arc(P['C'], ra, ra, theta1=270 - g['aOCO'] / 2, theta2=270 + g['aOCO'] / 2, lw=2, color='0.3')); ax.text(0, -ra / 2 - 0.06, '%.1f°' % g['aOCO'], ha='center', va='top')
rb = 1.25; ax.add_patch(Arc(P['Li'], rb, rb, theta1=90 - g['aOLiO'] / 2, theta2=90 + g['aOLiO'] / 2, lw=2, color='0.3')); ax.text(0, P['Li'][1] + rb / 2 + 0.06, '%.1f°' % g['aOLiO'], ha='center', va='bottom')
ax.set_aspect('equal'); ax.set_xlim(-2.1, 2.5); ax.set_ylim(-g['rLiC'] - 0.6, 0.55); ax.axis('off')
fig.tight_layout(); fig.savefig(sys.argv[2], dpi=200, bbox_inches='tight', pad_inches=0.2); print('saved', sys.argv[2])
