"""Image of Cl2CO at an equilibrium geometry with its bond lengths and angles.  usage: python mol_cl2co.py out.png"""
import sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc
B = 0.529177210903
def geom(X):
    C, O, Cl1, Cl2 = X * B; d = lambda a, b: np.linalg.norm(a - b)
    ang = lambda a, b, c: np.degrees(np.arccos(np.dot(a - b, c - b) / d(a, b) / d(c, b)))
    return dict(rCO=d(C, O), rCCl=(d(C, Cl1) + d(C, Cl2)) / 2, dCCl=abs(d(C, Cl1) - d(C, Cl2)), aClCCl=ang(Cl1, C, Cl2), aOCCl=(ang(O, C, Cl1) + ang(O, C, Cl2)) / 2,
                planar=abs(np.linalg.det(np.array([O - C, Cl1 - C, Cl2 - C]))))
G = {'B3LYP/aug-cc-pVTZ': np.load('vp_atz.npz')['Xeq'], 'B3LYP/aug-cc-pVDZ': np.load('vp_adz.npz')['Xeq'], 'NN PES': np.load('eq.npz')['X0']}
for k, X in G.items():
    g = geom(X); print('%-20s r(C=O) %.4f A  r(C-Cl) %.4f A (asym %.1e)  Cl-C-Cl %.2f deg  O-C-Cl %.2f deg  planarity %.1e' % (k, g['rCO'], g['rCCl'], g['dCCl'], g['aClCCl'], g['aOCCl'], g['planar']))
LEVEL = 'NN PES'; g = geom(G[LEVEL])
# planar drawing: C at origin, C=O up (C2 axis), Cl atoms below
h = np.radians(g['aClCCl'] / 2)
P = {'C': np.array([0.0, 0.0]), 'O': np.array([0.0, g['rCO']]), 'Cl1': g['rCCl'] * np.array([-np.sin(h), -np.cos(h)]), 'Cl2': g['rCCl'] * np.array([np.sin(h), -np.cos(h)])}
col = {'C': '#555555', 'O': '#d62728', 'Cl': '#2ca02c'}; rad = {'C': 0.30, 'O': 0.28, 'Cl': 0.40}
plt.rcParams.update({'font.size': 20}); fig, ax = plt.subplots(figsize=(8, 6.9))
for a in ('Cl1', 'Cl2'): ax.plot(*zip(P['C'], P[a]), color='k', lw=5, zorder=1, solid_capstyle='round')
for dx in (-0.06, 0.06): ax.plot([dx, dx], [0, g['rCO']], color='k', lw=4, zorder=1)                      # C=O double bond
for a, p in P.items():
    el = a.rstrip('12'); ax.add_patch(Circle(p, rad[el], facecolor=col[el], edgecolor='k', lw=1.5, zorder=2)); ax.text(p[0], p[1], el, ha='center', va='center', color='w', fontsize=22, fontweight='bold', zorder=3)
ax.text(0.16, 0.55 * g['rCO'], '%.3f Å' % g['rCO'], ha='left', va='center')
m = 0.55 * P['Cl2']; ax.text(m[0] + 0.14, m[1] + 0.10, '%.3f Å' % g['rCCl'], ha='left', va='center')
ra = 1.15; ax.add_patch(Arc((0, 0), ra, ra, theta1=270 - g['aClCCl'] / 2, theta2=270 + g['aClCCl'] / 2, lw=2, color='0.3'))
ax.text(0, -ra / 2 - 0.22, '%.1f°' % g['aClCCl'], ha='center', va='top')
rb = 0.85; ax.add_patch(Arc((0, 0), rb, rb, theta1=90, theta2=90 + g['aOCCl'], lw=2, color='0.3'))
t = np.radians(90 + g['aOCCl'] / 2); ax.text(1.02 * rb / 2 * np.cos(t) - 0.08, 1.02 * rb / 2 * np.sin(t) + 0.05, '%.1f°' % g['aOCCl'], ha='right', va='center')
ax.set_aspect('equal'); ax.set_xlim(-2.25, 2.25); ax.set_ylim(-1.75, 1.65); ax.axis('off')
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=200); print('saved', sys.argv[1])
