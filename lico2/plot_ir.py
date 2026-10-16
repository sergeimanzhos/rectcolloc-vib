"""IR spectrum plot from results/intens_<tag>.json: python plot_ir.py <tag> <title>"""
import sys, json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
tag, title = sys.argv[1], sys.argv[2]
j = json.load(open(f'results/intens_{tag}.json')); R = j['rows']; Hm = j['harm']
col = {'A1': '#1f77b4', 'B1': '#2ca02c', 'B2': '#d62728', 'A2': '#7f7f7f'}
plt.rcParams.update({'font.size': 16})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(12, 9.5), sharex=True)
x = np.linspace(200, 1700, 6001); g = 5.0                      # Lorentzian HWHM 5 cm-1
env = sum(r['I_b'] * (g / np.pi) / ((x - r['nu']) ** 2 + g ** 2) for r in R)
envh = sum(h['I'] * (g / np.pi) / ((x - h['nu']) ** 2 + g ** 2) for h in Hm)
a1.plot(x, envh, color='0.6', lw=1.5, ls='--', label='double harmonic')
a1.plot(x, env, color='k', lw=1.8, label='collocation')
a1.set_ylabel('Absorption (km mol$^{-1}$ / cm$^{-1}$)'); a1.legend(frameon=False); a1.set_title(title, fontsize=16)
for s in ('A1', 'B1', 'B2'):
    rr = [r for r in R if r['sym'] == s and r['I_b'] > 1e-3]
    a2.vlines([r['nu'] for r in rr], 1e-3, [r['I_b'] for r in rr], color=col[s], lw=2.5, label=f'collocation, {s}')
a2.plot([h['nu'] for h in Hm], [h['I'] for h in Hm], 'kv', ms=10, mfc='none', mew=1.8, label='double harmonic')
a2.set_yscale('log'); a2.set_ylim(1e-3, 2e3); a2.set_xlim(200, 1700)
a2.set_xlabel('Wavenumber (cm$^{-1}$)'); a2.set_ylabel('Band intensity (km mol$^{-1}$)'); a2.legend(frameon=False, ncol=2, fontsize=13, loc='upper left')
fig.tight_layout(); fig.savefig(f'results/ir_{tag}.png', dpi=130); print('saved', f'results/ir_{tag}.png')
