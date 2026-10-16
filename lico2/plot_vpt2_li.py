"""Mirror plot for LiCO2 in the format of cl2co/plot_vpt2.py: VPT2 sticks up, collocation sticks down (linear), and a log
version in two panels.  Optional experiment: one digitized spectrum per half (absorbance, right axis).
usage: python plot_vpt2_li.py out up1 [up2] "|" down1 [down2] [expU=csv:label] [expD=csv:label]   (sticks: file:label:color:linestyle)
Experimental csv: columns wavenumber_axis, wavenumber_band_corrected, transmittance_percent, clipped_flag (handoff format).
Absorbance = -log10(T / baseline), baseline = running maximum of T over +-75 cm-1, smoothed."""
import json, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
plt.rcParams.update({'font.size': 20})
out = sys.argv[1]; A_ = sys.argv[2:]
EX = {a[:4]: a[5:].split(':', 1) for a in A_ if a[:5] in ('expU=', 'expD=')}; A_ = [a for a in A_ if a[:5] not in ('expU=', 'expD=')]; k_ = A_.index('|')
US, DS = [(a.split(':') + ['-'])[:4] for a in A_[:k_]], [(a.split(':') + ['-'])[:4] for a in A_[k_ + 1:]]
UR = [json.load(open(f))['rows'] for f, _, _, _ in US]; DR = [json.load(open(f))['rows'] for f, _, _, _ in DS]
offs = lambda n: [0.0] if n == 1 else list((np.arange(n) - (n - 1) / 2) * 5.0)
OFFU, OFF = offs(len(US)), offs(len(DS))
lab = {'100000': r'$\nu_1$', '010000': r'$\nu_2$', '001000': r'$\nu_3$', '000001': r'$\nu_6$'}
labl = dict(lab, **{'000100': r'$\nu_4$', '000010': r'$\nu_5$', '020000': r'$2\nu_2$'})
import os
absfmt = FuncFormatter(lambda v, p: f'{abs(v):g}'); XL = (200, float(os.environ.get('XMAX', 1700)))
def get(R, l): return [r for r in R if r['label'] == l]
inx = lambda R: [r for r in R if XL[0] <= r['E'] <= XL[1]]
UR = [inx(R) for R in UR]; DR = [inx(R) for R in DR]; U = UR[0]
YL = 1.18 * max(r['A'] for R in UR + DR for r in R)
def absorbance(fn):
    d = np.loadtxt(fn, delimiter=',', skiprows=1); o = np.argsort(d[:, 1]); x, T = d[o, 1], d[o, 2]
    xg = np.arange(np.ceil(x.min()), np.floor(x.max()) + 1, 1.0); Tg = np.interp(xg, x, T)
    base = np.array([Tg[max(0, i - 75):i + 76].max() for i in range(len(xg))]); k = np.ones(151) / 151
    base = np.convolve(np.pad(base, 75, mode='edge'), k, mode='valid'); base = np.maximum(base, Tg)
    return xg, -np.log10(Tg / base)
EXD = {k: absorbance(v[0]) + (v[1],) for k, v in EX.items()}
# right-axis scale: the 1569 cm-1 band of the lower (or only) experimental trace reaches the height of the strongest stick
if EXD:
    xr, ar, _ = EXD.get('expD', list(EXD.values())[0]); AL = 1.18 * ar[(xr > 1540) & (xr < 1600)].max()
def expcurve(axx, key, sgn):
    x, a, _ = EXD[key]; m = (x >= XL[0]) & (x <= XL[1]); return axx.plot(x[m], sgn * a[m], color='k', lw=1.5, alpha=0.75)[0]
# ---------------- linear
fig, ax = plt.subplots(figsize=(12, 9))
for R, s, c, l, o, ls in [(R_, 1, c_, l_, o_, ls_) for R_, (_, l_, c_, ls_), o_ in zip(UR, US, OFFU)] + [(R_, -1, c_, l_, o_, ls_) for R_, (_, l_, c_, ls_), o_ in zip(DR, DS, OFF)]:
    rr = [r for r in R if r['A'] > 0.05]; ax.vlines([r['E'] + o for r in rr], 0, [s * r['A'] for r in rr], color=c, lw=3, linestyles=ls)
    ax.plot([], [], ls, color=c, lw=3, label=l)
h, l = ax.get_legend_handles_labels(); nU = len(US); hU, lU, hD, lD = h[:nU], l[:nU], h[nU:], l[nU:]
ax.axhline(0, color='k', lw=0.8); ax.set_ylim(-YL, YL); ax.yaxis.set_major_formatter(absfmt); ax.set_xlim(*XL)
if EXD:
    ax2 = ax.twinx(); ax2.set_ylim(-AL, AL); ax2.set_ylabel('absorbance'); ax2.yaxis.set_major_formatter(absfmt)
    if 'expU' in EXD: hU.append(expcurve(ax2, 'expU', 1)); lU.append(EXD['expU'][2])
    if 'expD' in EXD: hD.append(expcurve(ax2, 'expD', -1)); lD.append(EXD['expD'][2])
    ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
ax.set_xlabel(r'wavenumber (cm$^{-1}$)'); ax.set_ylabel('intensity (km/mol)'); ax.grid(alpha=0.3, axis='y')
leg1 = ax.legend(hU, lU, fontsize=15, frameon=True, facecolor='w', edgecolor='none', framealpha=1, loc='upper center', bbox_to_anchor=(0.42, 1.0)); ax.add_artist(leg1)
ax.legend(hD, lD, fontsize=15, frameon=True, facecolor='w', edgecolor='none', framealpha=1, loc='lower center', bbox_to_anchor=(0.42, 0.0))
for k, t in lab.items():
    r = get(U, k)
    if r: ax.text(r[0]['E'] - (22 if k in ('001000', '000001') else 0), max(x['A'] for R_ in UR for x in get(R_, k)) + 0.02 * YL, t, ha='right' if k in ('001000', '000001') else 'center', fontsize=17)
fig.tight_layout(); fig.savefig(out + '.png', dpi=200); plt.close(fig)
# ---------------- log
fig, (a1, a2) = plt.subplots(2, 1, sharex=True, figsize=(12, 10), gridspec_kw=dict(hspace=0))
for a, R, c, l, o, ls in [(a1, R_, c_, l_, o_, ls_) for R_, (_, l_, c_, ls_), o_ in zip(UR, US, OFFU)] + [(a2, R_, c_, l_, o_, ls_) for R_, (_, l_, c_, ls_), o_ in zip(DR, DS, OFF)]:
    rr = [r for r in R if r['A'] > 1e-4]; a.vlines([r['E'] + o for r in rr], 1e-4, [r['A'] for r in rr], color=c, lw=3, linestyles=ls)
    a.plot([], [], ls, color=c, lw=3, label=l)
for a, key in ((a1, 'expU'), (a2, 'expD')):
    a.set_yscale('log'); a.set_ylim(1e-4, 1e5); a.grid(alpha=0.3, axis='y')
    hh, ll = a.get_legend_handles_labels()
    if key in EXD:
        b = a.twinx(); b.set_ylim(0, AL); hh.append(expcurve(b, key, 1)); ll.append(EXD[key][2])
        if a is a2: b.invert_yaxis()
        a.set_zorder(b.get_zorder() + 1); a.patch.set_visible(False)
        b.set_yticks([t for t in b.get_yticks() if 0 < t < AL])
    a.legend(hh, ll, fontsize=15, frameon=True, facecolor='w', edgecolor='none', framealpha=1, loc='upper center' if a is a1 else 'lower center', bbox_to_anchor=(0.42, 1.0) if a is a1 else (0.42, 0.0))
a2.invert_yaxis(); a2.set_xlim(*XL); a2.set_xlabel(r'wavenumber (cm$^{-1}$)')
a1.set_yticks([1e-4, 1e-2, 1, 1e2, 1e4]); a2.set_yticks([1e-2, 1, 1e2, 1e4]); fig.supylabel('intensity (km/mol)', x=0.02)
if EXD: fig.text(0.985, 0.5, 'absorbance', rotation=90, va='center', ha='right')
for k, t in labl.items():
    r = get(U, k)
    if r: a1.text(r[0]['E'], max(x['A'] for R_ in UR for x in get(R_, k)) * 2.2, t, ha='center', fontsize=17)
fig.tight_layout(); fig.subplots_adjust(hspace=0, right=0.9 if EXD else None); fig.savefig(out + '_log.png', dpi=200); plt.close(fig)
print('saved', out + '.png', out + '_log.png')
