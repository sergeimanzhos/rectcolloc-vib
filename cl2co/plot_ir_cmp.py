"""Stick spectra comparing dipole sources: python3 plot_ir_cmp.py outstem file1:label1:color ..."""
import json, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
import os; EXP = os.environ.get('EXP'); out = sys.argv[1]; sets = [a.split(':') for a in sys.argv[2:]]
D = [{r['n']: r for r in json.load(open(f))['rows']} for f, _, _ in sets]
K = len(sets); off = (np.arange(K) - (K - 1) / 2) * 5.0; w = 3.0 if K <= 2 else 2.4
MK = ['o', 's', '^', 'D', 'v', 'P']
lab = {7: r'$\nu_5$', 67: r'$\nu_1$', 68: r'$\nu_2\!+\!\nu_5\!+\!\nu_6$', 3: r'$\nu_2$', 4: r'$\nu_4$',
       48: r'$2\nu_5$', 12: r'$\nu_2\!+\!\nu_6$', 1: r'$\nu_3$', 2: r'$\nu_6$', 28: r'$\nu_2\!+\!\nu_5$',
       19: r'$2\nu_4$'}
for log in (False, True):
    fig, ax = plt.subplots(figsize=(12, 6.3))
    thr = 1e-4 if log else 0.05
    for i, ((f, l, c), d, o) in enumerate(zip(sets, D, off)):
        rows = [r for r in d.values() if r['A'] > thr]; x = [r['E'] + o for r in rows]; y = [r['A'] for r in rows]
        ax.vlines(x, 1e-4 if log else 0, y, color=c, lw=w)
        ax.plot([], [], '-', color=c, lw=3, label=l)
    for n, t in lab.items():
        E = D[0][n]['E']; A = max(d[n]['A'] for d in D)
        if log:
            dx, ha = {3: (-8, 'right'), 4: (8, 'left'), 67: (-12, 'right'), 68: (12, 'left')}.get(n, (0, 'center'))
            ax.text(E + dx, A * 2.2, t, ha=ha, fontsize=17)
        elif n == 67 and EXP:
            ax.text(E - 50, A - 40, t, ha='right', fontsize=17)
        elif n == 48 and EXP:
            ax.text(E, 66, t, ha='center', fontsize=17)
        elif n in (7, 67):
            ax.text(E - 15, A - 40, t, ha='right', fontsize=17)
        elif n in (68,):
            ax.text(E + 12, A + 5, t, ha='left', fontsize=17)
        elif n in (3, 12, 48):
            ax.text(E - (25 if n == 3 else 0), A + 12, t, ha='right' if n == 3 else 'center', fontsize=17)
    ax.set_xlabel(r'wavenumber (cm$^{-1}$)'); ax.set_ylabel('intensity (km/mol)'); ax.set_xlim(200, 2150)
    if log:
        ax.set_yscale('log'); ax.set_ylim(1e-4, 1e5 if K <= 3 else 1e6)
        ax.legend(fontsize=15 if K <= 3 else 14, frameon=False, loc='upper center', ncol=K if K <= 2 else (3 if K == 3 else 2))
    else:
        ax.set_ylim(0, 520); h, l = ax.get_legend_handles_labels()
        if EXP:
            ex = np.loadtxt(EXP, delimiter=','); ax2 = ax.twinx()
            h += ax2.plot(ex[:, 0], ex[:, 1], color='k', lw=1.6, alpha=0.8, label='experiment (absorbance, right axis)')
            l.append(h[-1].get_label()); ax2.set_ylim(0, 0.36); ax2.set_ylabel('absorbance'); ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
        ax.legend(h, l, fontsize=15, frameon=False, loc='upper center', bbox_to_anchor=(0.6, 1.0))
    ax.grid(alpha=0.3, axis='y'); fig.tight_layout()
    fig.savefig(out + ('_log' if log else '') + '.png', dpi=200); plt.close(fig)
# csv
with open(out + '.csv', 'w') as fh:
    fh.write('n,label,E_cm-1,' + ','.join(l.replace(',', ';') for _, l, _ in sets) + '\n')
    for n, r in D[0].items():
        fh.write(f"{n},{r['label']},{r['E']:.3f}," + ','.join(f"{d[n]['A']:.6g}" for d in D) + '\n')
