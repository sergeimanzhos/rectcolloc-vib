"""Summary mirror plot for LiCO2: collocation sticks up, VPT2 sticks down; basis sets distinguished by colour and line style.
Collocation sticks are narrow with a pale band of the same colour over the convergence range of each level (rows with Elo, Ehi).
Linear and log versions.  usage: XMIN=.. XMAX=.. FIGW=.. FIGH=.. FONT=.. python plot_summary.py out up1 [up2] "|" down1 [down2] [band=text]
each set: file:label:color:linestyle"""
import json, sys, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch
FS = (float(os.environ.get('FIGW', 9.69)), float(os.environ.get('FIGH', 6.27))); FONT = float(os.environ.get('FONT', 16)); plt.rcParams.update({'font.size': FONT})
out = sys.argv[1]; A_ = sys.argv[2:]
BANDL = ([a[5:] for a in A_ if a[:5] == 'band='] or [''])[0]; A_ = [a for a in A_ if a[:5] != 'band=']; k_ = A_.index('|')
US, DS = [(a.split(':') + ['-'])[:4] for a in A_[:k_]], [(a.split(':') + ['-'])[:4] for a in A_[k_ + 1:]]
XL = (float(os.environ.get('XMIN', 200)), float(os.environ.get('XMAX', 2300)))
inx = lambda R: [r for r in R if XL[0] <= r['E'] <= XL[1]]
UR = [inx(json.load(open(f))['rows']) for f, _, _, _ in US]; DR = [inx(json.load(open(f))['rows']) for f, _, _, _ in DS]
lab = {'100000': r'$\nu_1$', '010000': r'$\nu_2$', '001000': r'$\nu_3$', '000100': r'$\nu_4$', '000010': r'$\nu_5$', '000001': r'$\nu_6$', '020000': r'$2\nu_2$'}
absfmt = FuncFormatter(lambda v, p: f'{abs(v):g}')
CMP = os.environ.get('COMPACT') == '1'     # one-row legends, band entry with the lower legend, less headroom (short figures)
YL = (1.3 if CMP else 1.75) * max([r['A'] for R in UR + DR for r in R] + [1e-3])
LEG = dict(fontsize=0.72 * FONT, frameon=True, facecolor='w', edgecolor='none', framealpha=0.85)
if CMP: LEG.update(fontsize=0.62 * FONT, columnspacing=1.0, handlelength=1.6, borderpad=0.2)
def labels():           # one label per state: the strongest stick carrying that label among the upper sets
    best = {}
    for R in UR:
        for r in R:
            if r['label'] in lab and r['A'] > best.get(r['label'], (0, 0))[1]: best[r['label']] = (r['E'], r['A'])
    return best
for kind in ('lin', 'log'):
    if kind == 'lin':
        fig, ax = plt.subplots(figsize=FS); aU = aD = ax; s = -1; y0 = 0.0; thr = 0.002 * YL
    else:
        fig, (aU, aD) = plt.subplots(2, 1, sharex=True, figsize=FS, gridspec_kw=dict(hspace=0)); s = 1; y0 = thr = 1e-4
    hU, lU, hD, lD = [], [], [], []
    for R, (_, l, c, ls) in zip(UR, US):
        rr = [r for r in R if r['A'] > thr]
        for r in rr:
            if 'Elo' in r: aU.fill_between([r['Elo'], r['Ehi']], [y0, y0], [r['A'], r['A']], color=c, alpha=0.25, lw=0)
        aU.vlines([r['E'] for r in rr], y0, [r['A'] for r in rr], color=c, lw=1.8, linestyles=ls); hU.append(aU.plot([], [], ls, color=c, lw=1.8)[0]); lU.append(l)
    if BANDL: hU.append(Patch(color='0.5', alpha=0.3, lw=0)); lU.append(BANDL)
    for R, (_, l, c, ls) in zip(DR, DS):
        rr = [r for r in R if r['A'] > thr]
        aD.vlines([r['E'] for r in rr], y0, [s * r['A'] for r in rr], color=c, lw=1.8, linestyles=ls); hD.append(aD.plot([], [], ls, color=c, lw=1.8)[0]); lD.append(l)
    if kind == 'lin':
        ax.axhline(0, color='k', lw=0.8); ax.set_ylim(-YL, YL); ax.yaxis.set_major_formatter(absfmt); ax.set_ylabel('intensity (km/mol)'); ax.grid(alpha=0.3, axis='y')
        if CMP: hD += hU[2:]; lD += lU[2:]; hU, lU = hU[:2], lU[:2]
        l1 = ax.legend(hU, lU, loc='upper center', ncol=2 if CMP else 1, **LEG); ax.add_artist(l1); ax.legend(hD, lD, loc='lower center', ncol=3 if CMP else 1, **LEG)
        for k, (E, A) in labels().items():
            if A > thr:
                big = A > 0.55 * YL
                ax.text(E - (0.012 * (XL[1] - XL[0]) if big else 0), 0.93 * A if big else A + 0.02 * YL, lab[k], ha='right' if big else 'center', fontsize=0.85 * FONT)
        aD.set_xlabel(r'wavenumber (cm$^{-1}$)'); fig.tight_layout()
    else:
        for a in (aU, aD): a.set_yscale('log'); a.set_ylim(1e-4, 1e6 if CMP else 1e8); a.grid(alpha=0.3, axis='y')
        if CMP: hD += hU[2:]; lD += lU[2:]; hU, lU = hU[:2], lU[:2]
        aU.legend(hU, lU, loc='upper center', ncol=2 if CMP else 1, **LEG); aD.legend(hD, lD, loc='lower center', ncol=3 if CMP else 1, **LEG)
        aD.invert_yaxis(); aD.set_xlabel(r'wavenumber (cm$^{-1}$)'); aU.set_yticks([1e-4, 1e-2, 1, 1e2]); aD.set_yticks([1e-2, 1, 1e2])
        fig.supylabel('intensity (km/mol)', x=0.02)
        for k, (E, A) in labels().items():
            if A > thr: aU.text(E, A * 2.2, lab[k], ha='center', fontsize=0.85 * FONT)
        fig.tight_layout(); fig.subplots_adjust(hspace=0)
    aD.set_xlim(*XL); fig.savefig(out + ('.png' if kind == 'lin' else '_log.png'), dpi=200, bbox_inches='tight', pad_inches=0.08); plt.close(fig)
print('saved', out + '.png', out + '_log.png')
