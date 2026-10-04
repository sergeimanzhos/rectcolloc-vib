"""Mirror stick spectra: first two sets up, last two down (negative y).
usage: EXP=exp.csv python3 plot_ir_mirror.py outstem f1:lab1:col f2:... f3:... f4:..."""
import json, sys, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
plt.rcParams.update({'font.size': 20})
EXP = os.environ.get('EXP'); out = sys.argv[1]; sets = [(a.split(':') + ['-'])[:4] for a in sys.argv[2:]]
D = [{r['n']: r for r in json.load(open(f))['rows']} for f, _, _, _ in sets]
up, dn = [0, 1], [2, 3]; off = [-2.5, 2.5]; LW = 3.0
lab = {7: r'$\nu_5$', 67: r'$\nu_1$', 68: r'$\nu_2\!+\!\nu_5\!+\!\nu_6$', 3: r'$\nu_2$', 4: r'$\nu_4$',
       48: r'$2\nu_5$', 12: r'$\nu_2\!+\!\nu_6$', 1: r'$\nu_3$', 2: r'$\nu_6$', 28: r'$\nu_2\!+\!\nu_5$', 19: r'$2\nu_4$'}
absfmt = FuncFormatter(lambda v, p: f'{abs(v):g}')
XL = (200, 2150)

def sticks(ax, idx, sign, thr, base):
    for k, o in zip(idx, off):
        f, l, c, ls = sets[k]; rows = [r for r in D[k].values() if r['A'] > thr]
        ax.vlines([r['E'] + o for r in rows], sign * base, [sign * r['A'] for r in rows], color=c, lw=LW, linestyles=ls)
        ax.plot([], [], ls, color=c, lw=3, label=l)

# ---------- linear ----------
fig, ax = plt.subplots(figsize=(12, 9))
sticks(ax, up, +1, 0.05, 0); hu, lu = ax.get_legend_handles_labels()
sticks(ax, dn, -1, 0.05, 0); h, l = ax.get_legend_handles_labels(); hd, ld = h[len(hu):], l[len(lu):]
ax.axhline(0, color='k', lw=0.8)
if EXP:
    ex = np.loadtxt(EXP, delimiter=','); ax2 = ax.twinx()
    for s in (1, -1): eh, = ax2.plot(ex[:, 0], s * ex[:, 1], color='k', lw=1.5, alpha=0.75)
    ax2.set_ylim(-0.36, 0.36); ax2.set_ylabel('absorbance'); ax2.yaxis.set_major_formatter(absfmt)
    ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
    hu.append(eh); lu.append('experiment (absorbance, right axis)')
ax.set_ylim(-520, 520); ax.yaxis.set_major_formatter(absfmt); ax.set_xlim(*XL)
ax.set_xlabel(r'wavenumber (cm$^{-1}$)'); ax.set_ylabel('intensity (km/mol)'); ax.grid(alpha=0.3, axis='y')
leg1 = ax.legend(hu, lu, fontsize=15, frameon=False, loc='upper center', bbox_to_anchor=(0.6, 1.0)); ax.add_artist(leg1)
ax.legend(hd, ld, fontsize=15, frameon=False, loc='lower center', bbox_to_anchor=(0.6, 0.0))
pos = {7: (-15, -40, 'right'), 67: (-50 if EXP else -15, -40, 'right'), 68: (12, 5, 'left'), 3: (-25, 12, 'right'),
       12: (0, 12, 'center'), 48: (0, None, 'center')}
for n, (dx, dy, ha) in pos.items():
    E = D[0][n]['E']; A = max(D[k][n]['A'] for k in up)
    y = (66 if EXP else A + 12) if dy is None else A + dy
    ax.text(E + dx, y, lab[n], ha=ha, fontsize=17)
fig.tight_layout(); fig.savefig(out + '.png', dpi=200); plt.close(fig)

# ---------- log (two mirrored panels) ----------
fig, (a1, a2) = plt.subplots(2, 1, sharex=True, figsize=(12, 10), gridspec_kw=dict(hspace=0))
for a, idx in ((a1, up), (a2, dn)):
    sticks(a, idx, +1, 1e-4, 1e-4); a.set_yscale('log'); a.set_ylim(1e-4, 1e5); a.grid(alpha=0.3, axis='y')
    a.legend(fontsize=15, frameon=False, loc='upper center' if a is a1 else 'lower center', ncol=2)
a2.invert_yaxis(); a2.set_xlim(*XL); a2.set_xlabel(r'wavenumber (cm$^{-1}$)')
a1.set_yticks([1e-4, 1e-2, 1, 1e2, 1e4]); a2.set_yticks([1e-2, 1, 1e2, 1e4])
fig.supylabel('intensity (km/mol)', x=0.02)
for n, t in lab.items():
    E = D[0][n]['E']; A = max(D[k][n]['A'] for k in up)
    dx, ha = {3: (-8, 'right'), 4: (8, 'left'), 67: (-12, 'right'), 68: (12, 'left')}.get(n, (0, 'center'))
    a1.text(E + dx, A * 2.2, t, ha=ha, fontsize=17)
fig.tight_layout(); fig.subplots_adjust(hspace=0); fig.savefig(out + '_log.png', dpi=200); plt.close(fig)
