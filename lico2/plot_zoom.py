"""Zoomed mirror plots for LiCO2 (format of cl2co/plot_vpt2.py): VPT2 sticks up, collocation sticks down, one digitized
experimental spectrum per half (absorbance, right axis), linear and log versions.  Collocation sticks are narrow, with a
pale band of the same colour spanning the convergence range of each level (rows with Elo, Ehi).
usage: XMIN=.. XMAX=.. python plot_zoom.py out up:label:color "|" down:label:color expU=csv:label expD=csv:label band="legend text" """
import json, sys, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch
FS = (float(os.environ.get('FIGW', 12)), float(os.environ.get('FIGH', 9))); FONT = float(os.environ.get('FONT', 20)); plt.rcParams.update({'font.size': FONT})
out = sys.argv[1]; A_ = sys.argv[2:]
EX = {a[:4]: a[5:].split(':', 1) for a in A_ if a[:5] in ('expU=', 'expD=')}
BANDL = ([a[5:] for a in A_ if a[:5] == 'band='] or [''])[0]
A_ = [a for a in A_ if a[:5] not in ('expU=', 'expD=', 'band=')]; k_ = A_.index('|')
(fu, lu, cu), (fd, ld, cd) = A_[0].split(':')[:3], A_[k_ + 1].split(':')[:3]
XL = (float(os.environ.get('XMIN', 200)), float(os.environ.get('XMAX', 1700)))
inx = lambda R: [r for r in R if XL[0] <= r['E'] <= XL[1]]
U, D = inx(json.load(open(fu))['rows']), inx(json.load(open(fd))['rows'])
lab = {'100000': r'$\nu_1$', '010000': r'$\nu_2$', '001000': r'$\nu_3$', '000100': r'$\nu_4$', '000010': r'$\nu_5$', '000001': r'$\nu_6$', '020000': r'$2\nu_2$'}
absfmt = FuncFormatter(lambda v, p: f'{abs(v):g}')
YL = 1.25 * max([r['A'] for r in U + D] + [1e-3])
def absorbance(fn):
    d = np.loadtxt(fn, delimiter=',', skiprows=1); o = np.argsort(d[:, 1]); x, T = d[o, 1], d[o, 2]
    xg = np.arange(np.ceil(x.min()), np.floor(x.max()) + 1, 1.0); Tg = np.interp(xg, x, T)
    base = np.array([Tg[max(0, i - 75):i + 76].max() for i in range(len(xg))]); k = np.ones(151) / 151
    base = np.convolve(np.pad(base, 75, mode='edge'), k, mode='valid'); base = np.maximum(base, Tg)
    m = (xg >= XL[0]) & (xg <= XL[1]); return xg[m], -np.log10(Tg / base)[m]
EXD = {k: absorbance(v[0]) + (v[1],) for k, v in EX.items()}
AL = 1.25 * max([a.max() for x, a, _ in EXD.values() if len(a)] + [1e-3])
def curve(axx, key, sgn):
    x, a, _ = EXD[key]; return axx.plot(x, sgn * a, color='k', lw=1.5, alpha=0.75)[0]
LEG = dict(fontsize=0.75 * FONT, frameon=True, facecolor='w', edgecolor='none', framealpha=0.85)
for kind in ('lin', 'log'):
    if kind == 'lin':
        fig, ax = plt.subplots(figsize=FS); aU = aD = ax; s = -1; y0 = 0.0; thr = 0.002 * YL
    else:
        fig, (aU, aD) = plt.subplots(2, 1, sharex=True, figsize=(FS[0], FS[1] + (1 if 'FIGH' not in os.environ else 0)), gridspec_kw=dict(hspace=0)); s = 1; y0 = thr = 1e-4
    rr = [r for r in U if r['A'] > thr]; aU.vlines([r['E'] for r in rr], y0, [r['A'] for r in rr], color=cu, lw=3)
    hU, lU = [aU.plot([], [], '-', color=cu, lw=3)[0]], [lu]
    rr = [r for r in D if r['A'] > thr]
    for r in rr:
        if 'Elo' in r: aD.fill_between([r['Elo'], r['Ehi']], [y0, y0], [s * r['A'], s * r['A']], color=cd, alpha=0.3, lw=0)
    aD.vlines([r['E'] for r in rr], y0, [s * r['A'] for r in rr], color=cd, lw=1.5)
    hD, lD = [aD.plot([], [], '-', color=cd, lw=1.5)[0]], [ld]
    if BANDL: hD.append(Patch(color=cd, alpha=0.3, lw=0)); lD.append(BANDL)
    if kind == 'lin':
        ax.axhline(0, color='k', lw=0.8); ax.set_ylim(-YL, YL); ax.yaxis.set_major_formatter(absfmt); ax.set_ylabel('intensity (km/mol)'); ax.grid(alpha=0.3, axis='y')
        if EXD:
            ax2 = ax.twinx(); ax2.set_ylim(-AL, AL); ax2.set_ylabel('absorbance'); ax2.yaxis.set_major_formatter(absfmt)
            if 'expU' in EXD: hU.append(curve(ax2, 'expU', 1)); lU.append(EXD['expU'][2])
            if 'expD' in EXD: hD.append(curve(ax2, 'expD', -1)); lD.append(EXD['expD'][2])
            ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
        l1 = ax.legend(hU, lU, loc='upper center', **LEG); ax.add_artist(l1); ax.legend(hD, lD, loc='lower center', **LEG)
        for r in U:
            if r['label'] in lab and r['A'] > thr:
                big = r['A'] > 0.6 * YL      # tall stick: label beside it, not above (the legend is there)
                ax.text(r['E'] - (0.012 * (XL[1] - XL[0]) if big else 0), 0.93 * r['A'] if big else r['A'] + 0.02 * YL, lab[r['label']], ha='right' if big else 'center', fontsize=0.85 * FONT)
        aD.set_xlabel(r'wavenumber (cm$^{-1}$)'); fig.tight_layout()
    else:
        for a, key, hh, ll in ((aU, 'expU', hU, lU), (aD, 'expD', hD, lD)):
            a.set_yscale('log'); a.set_ylim(1e-4, 1e6); a.grid(alpha=0.3, axis='y')
            if key in EXD:
                b = a.twinx(); b.set_ylim(0, AL); hh.append(curve(b, key, 1)); ll.append(EXD[key][2])
                if a is aD: b.invert_yaxis()
                a.set_zorder(b.get_zorder() + 1); a.patch.set_visible(False); b.set_yticks([t for t in b.get_yticks() if 0 < t < AL])
            a.legend(hh, ll, loc='upper center' if a is aU else 'lower center', **LEG)
        aD.invert_yaxis(); aD.set_xlabel(r'wavenumber (cm$^{-1}$)'); aU.set_yticks([1e-4, 1e-2, 1, 1e2]); aD.set_yticks([1e-2, 1, 1e2])
        fig.supylabel('intensity (km/mol)', x=0.02)
        if EXD: fig.text(0.985, 0.5, 'absorbance', rotation=90, va='center', ha='right')
        for r in U:
            if r['label'] in lab and r['A'] > thr: aU.text(r['E'], r['A'] * 2.2, lab[r['label']], ha='center', fontsize=0.85 * FONT)
        fig.tight_layout(); fig.subplots_adjust(hspace=0, right=0.9 if EXD else None)
    aD.set_xlim(*XL); fig.savefig(out + ('.png' if kind == 'lin' else '_log.png'), dpi=200); plt.close(fig)
print('saved', out + '.png', out + '_log.png')
