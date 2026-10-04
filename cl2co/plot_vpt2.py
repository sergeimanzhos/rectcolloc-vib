"""Mirror plot: VPT2 sticks up, collocation sticks down, experiment on both sides (linear); log version in two panels.
usage: python3 plot_vpt2.py out vpt2.json:label:color colloc.json:label:color"""
import json, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
plt.rcParams.update({'font.size': 20})
out = sys.argv[1]; (fu, lu, cu) = sys.argv[2].split(':'); DS = [(a.split(':') + ['-'])[:4] for a in sys.argv[3:]]
U = json.load(open(fu))['rows']; DR = [json.load(open(f))['rows'] for f, _, _, _ in DS]
OFF = [0.0] if len(DS) == 1 else list((np.arange(len(DS)) - (len(DS) - 1) / 2) * 5.0)
lab = {'000010': r'$\nu_5$', '000001': r'$\nu_1$', '011010': r'$\nu_2\!+\!\nu_5\!+\!\nu_6$', '001000': r'$\nu_2$',
       '000020': r'$2\nu_5$', '011000': r'$\nu_2\!+\!\nu_6$'}
labl = dict(lab, **{'000100': r'$\nu_4$', '100000': r'$\nu_3$', '010000': r'$\nu_6$', '001010': r'$\nu_2\!+\!\nu_5$', '000200': r'$2\nu_4$'})
absfmt = FuncFormatter(lambda v, p: f'{abs(v):g}'); XL = (200, 2150)
ex = np.loadtxt('exp_digitized.csv', delimiter=',')
def get(R, l): return [r for r in R if r['label'] == l]
# linear
fig, ax = plt.subplots(figsize=(12, 9))
for R, s, c, l, o, ls in [(U, 1, cu, lu, 0.0, '-')] + [(R_, -1, c_, l_, o_, ls_) for R_, (_, l_, c_, ls_), o_ in zip(DR, DS, OFF)]:
    rr = [r for r in R if r['A'] > 0.05]; ax.vlines([r['E'] + o for r in rr], 0, [s * r['A'] for r in rr], color=c, lw=3, linestyles=ls)
    ax.plot([], [], ls, color=c, lw=3, label=l)
h, l = ax.get_legend_handles_labels()
ax.axhline(0, color='k', lw=0.8); ax2 = ax.twinx()
for s in (1, -1): eh, = ax2.plot(ex[:, 0], s * ex[:, 1], color='k', lw=1.5, alpha=0.75)
YL = 520 if len(DS) == 1 else 610; ax2.set_ylim(-0.36 * YL / 520, 0.36 * YL / 520); ax2.set_ylabel('absorbance'); ax2.yaxis.set_major_formatter(absfmt)
ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
ax.set_ylim(-YL, YL); ax.yaxis.set_major_formatter(absfmt); ax.set_xlim(*XL)
ax.set_xlabel(r'wavenumber (cm$^{-1}$)'); ax.set_ylabel('intensity (km/mol)'); ax.grid(alpha=0.3, axis='y')
leg1 = ax.legend([h[0], eh], [l[0], 'experiment (absorbance, right axis)'], fontsize=15, frameon=False, loc='upper center', bbox_to_anchor=(0.6, 1.0)); ax.add_artist(leg1)
ax.legend(h[1:], l[1:], fontsize=15, frameon=False, loc='lower center', bbox_to_anchor=(0.6, 0.0))
pos = {'000010': (-15, -40, 'right'), '000001': (-50, -40, 'right'), '011010': (12, 5, 'left'), '001000': (-25, 12, 'right'),
       '011000': (0, 12, 'center'), '000020': (0, None, 'center')}
for k, (dx, dy, ha) in pos.items():
    r = get(U, k)
    if not r: continue
    E, A = r[0]['E'], r[0]['A']; y = 66 if dy is None else A + dy
    if k == '011010' and A < 50: dx, y, ha = 10, 105, 'right'
    if k == '000020': y = A + 12
    ax.text(E + dx, y, lab[k], ha=ha, fontsize=17)
fig.tight_layout(); fig.savefig(out + '.png', dpi=200); plt.close(fig)
# log
fig, (a1, a2) = plt.subplots(2, 1, sharex=True, figsize=(12, 10), gridspec_kw=dict(hspace=0))
for a, R, c, l, o, ls in [(a1, U, cu, lu, 0.0, '-')] + [(a2, R_, c_, l_, o_, ls_) for R_, (_, l_, c_, ls_), o_ in zip(DR, DS, OFF)]:
    rr = [r for r in R if r['A'] > 1e-4]; a.vlines([r['E'] + o for r in rr], 1e-4, [r['A'] for r in rr], color=c, lw=3, linestyles=ls)
    a.plot([], [], ls, color=c, lw=3, label=l)
for a in (a1, a2):
    a.set_yscale('log'); a.set_ylim(1e-4, 1e5); a.grid(alpha=0.3, axis='y')
    a.legend(fontsize=15, frameon=False, loc='upper center' if a is a1 else 'lower center')
a2.invert_yaxis(); a2.set_xlim(*XL); a2.set_xlabel(r'wavenumber (cm$^{-1}$)')
a1.set_yticks([1e-4, 1e-2, 1, 1e2, 1e4]); a2.set_yticks([1e-2, 1, 1e2, 1e4]); fig.supylabel('intensity (km/mol)', x=0.02)
for k, t in labl.items():
    r = get(U, k)
    if not r: continue
    dx, ha = {'001000': (-8, 'right'), '000100': (8, 'left'), '000001': (-12, 'right'), '011010': (12, 'left')}.get(k, (0, 'center'))
    a1.text(r[0]['E'] + dx, r[0]['A'] * 2.2, t, ha=ha, fontsize=17)
fig.tight_layout(); fig.subplots_adjust(hspace=0); fig.savefig(out + '_log.png', dpi=200); plt.close(fig)
