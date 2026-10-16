"""Final free-CO2 level table (Gaussian basis only): harmonic | GVPT2 | VCI (fitted quartic FF) | collocation E (N = 4000, M = 8157)
with N and M convergence | collocation E + grad (N = 2400, M = 2400 points = 9600 rows).  python lev_final.py"""
import json, numpy as np
EMAX = 4350.0
vp = json.load(open('vptlev_co2fit_K14.json')); w = vp['w']; wb, w1, w3 = 0.5 * (w[0] + w[1]), w[2], w[3]
def bylist(L):
    out = {l: [] for l in range(5)}
    for r in L:
        l = int(round(np.sqrt(max(r['l2'], 0))))
        if l <= 4 and not any(abs(e - r['E']) < 0.02 for e in out[l]): out[l].append(r['E'])
    return {l: sorted(v) for l, v in out.items()}
vci, gv = bylist(vp['vci']), bylist(vp['gvpt2'])
zero = {l: sorted([(v1 * w1 + v2 * wb + v3 * w3, (v1, v2, v3)) for v1 in range(6) for v3 in range(4) for v2 in range(l, 12, 2) if v1 + v2 + v3 > 0]) for l in range(5)}
E = json.load(open('lev_dft2.json')); G = json.load(open('lev_dftgrad2.json'))
def one(runs, N, M): r = [r for r in runs if r['N'][0] == N and r['M'] == M][0]; return {int(l): v for l, v in r['lev'].items()}, r
eN4, rN4 = one(E, 4000, 8157); eN2, rN2 = one(E, 2400, 8157); eM48, rM48 = one(E, 2400, 4800); gg, rg = one(G, 2400, 2400); g12, _ = one(G, 1200, 2400); g21, _ = one(G, 2400, 1200)
rows = sorted((e, l, i) for l in range(5) for i, e in enumerate(vci[l]) if e <= EMAX)
SUP = str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹')
lab = lambda l, i: '(%d %d %d)' % zero[l][i][1]
g = lambda d, l, i: d[l][i] if l in d and i < len(d[l]) else None
f = lambda v, fmt='%.2f': fmt % v if v is not None else '-'
hdr = ['#', '(v1 v2 v3)', 'l', 'harmonic', 'GVPT2', 'VCI', 'colloc. E', 'ΔN', 'ΔM', 'colloc. E+grad', 'E+grad − E']
out = ['| ' + ' | '.join(hdr) + ' |', '|' + '---|' * len(hdr)]; csv = [','.join(hdr)]; dN, dM, dG, dGN, dGM = [], [], [], [], []
a, b, m, x = rN4['zpe'], rN2['zpe'], rM48['zpe'], rg['zpe']
c = ['0', 'ZPE (absolute)', '', '%.2f' % (wb + 0.5 * (w1 + w3)), '%.2f' % vp['zpe_gvpt2'], '%.2f' % vp['zpe_vci'], '%.2f' % a, '%+.2f' % (a - b), '%+.2f' % (b - m), '%.2f' % x, '%+.2f' % (x - a)]
out.append('| ' + ' | '.join(c) + ' |'); csv.append(','.join(c))
for n, (e, l, i) in enumerate(rows, 1):
    a, b, m, x = g(eN4, l, i), g(eN2, l, i), g(eM48, l, i), g(gg, l, i); dN.append(a - b); dM.append(b - m); dG.append(x - a); dGN.append(x - g(g12, l, i)); dGM.append(x - g(g21, l, i))
    c = ['%d' % n, lab(l, i), '%d' % l, '%.2f' % zero[l][i][0], f(g(gv, l, i)), f(g(vci, l, i)), f(a), f(a - b, '%+.2f'), f(b - m, '%+.2f'), f(x), f(x - a, '%+.2f')]
    out.append('| ' + ' | '.join(c) + ' |'); csv.append(','.join(c))
print('\n'.join(out)); open('co2_levels_final.md', 'w').write('\n'.join(out) + '\n'); open('co2_levels_final.csv', 'w', encoding='utf-8').write('\n'.join(csv) + '\n')
lo = [k for k, (e, l, i) in enumerate(rows) if e <= 3000]
for name, d in (('E: N 4000 vs 2400 (M 8157)', dN), ('E: M 8157 vs 4800 (N 2400)', dM), ('E+grad: N 2400 vs 1200 (M 2400)', dGN), ('E+grad: M 2400 vs 1200 (N 2400)', dGM), ('E+grad vs E', dG)):
    d = np.abs(d); print('%-34s max all %.2f, max below 3000 %.2f, rms %.2f' % (name, d.max(), d[lo].max(), np.sqrt(np.mean(d ** 2))))
print('ZPE: E N4000 %.3f, E+grad %.3f, GVPT2 %.2f, VCI %.2f' % (rN4['zpe'], rg['zpe'], vp['zpe_gvpt2'], vp['zpe_vci']))
