"""Free-CO2 level table: harmonic | GVPT2 | VCI (fitted quartic FF) | collocation with energy + gradient rows (Gaussian basis).
Uses N = 4000 if lev_dftgrad4k.json exists, else N = 2400 (M = 2400 points).  Prints the N / M convergence below the table."""
import os, json, numpy as np
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
G = json.load(open('lev_dftgrad2.json')) + (json.load(open('lev_dftgrad4k.json')) if os.path.exists('lev_dftgrad4k.json') else [])
def one(N, M):
    r = [r for r in G if r['N'][0] == N and r['M'] == M]; return ({int(l): v for l, v in r[0]['lev'].items()}, r[0]) if r else (None, None)
NB = 4000 if one(4000, 2400)[0] else 2400
best, rb = one(NB, 2400); prevN, rpN = one(2400 if NB == 4000 else 1200, 2400); prevM, rpM = one(NB, 1200)
rows = sorted((e, l, i) for l in range(5) for i, e in enumerate(vci[l]) if e <= EMAX)
g = lambda d, l, i: d[l][i] if d and l in d and i < len(d[l]) else None
f = lambda v: '%.2f' % v if v is not None else '-'
hdr = ['#', '(v1 v2 v3)', 'l', 'harmonic', 'GVPT2', 'VCI', 'collocation']
out = ['| ' + ' | '.join(hdr) + ' |', '|' + '---|' * len(hdr)]; csv = [','.join(hdr)]
c = ['0', 'ZPE (absolute)', '', '%.2f' % (wb + 0.5 * (w1 + w3)), '%.2f' % vp['zpe_gvpt2'], '%.2f' % vp['zpe_vci'], '%.2f' % rb['zpe']]; out.append('| ' + ' | '.join(c) + ' |'); csv.append(','.join(c))
dN, dM = [abs(rb['zpe'] - rpN['zpe'])], [abs(rb['zpe'] - rpM['zpe'])] if rpM else [np.nan]
for n, (e, l, i) in enumerate(rows, 1):
    x = g(best, l, i); c = ['%d' % n, '(%d %d %d)' % zero[l][i][1], '%d' % l, '%.2f' % zero[l][i][0], f(g(gv, l, i)), f(g(vci, l, i)), f(x)]
    out.append('| ' + ' | '.join(c) + ' |'); csv.append(','.join(c)); dN.append(abs(x - g(prevN, l, i)))
    if prevM: dM.append(abs(x - g(prevM, l, i)))
print('\n'.join(out)); open('co2_levels_final.md', 'w').write('\n'.join(out) + '\n'); open('co2_levels_final.csv', 'w', encoding='utf-8').write('\n'.join(csv) + '\n')
lo = [0] + [k + 1 for k, (e, l, i) in enumerate(rows) if e <= 3000]
for name, d in (('N %d -> %d (M 2400)' % (rpN['N'][0], NB), dN), ('M 1200 -> 2400 (N %d)' % NB, dM)):
    d = np.array(d); print('%-26s max all %.2f, max below 3000 %.2f, rms %.2f' % (name, np.nanmax(d), np.nanmax(d[lo]), np.sqrt(np.nanmean(d ** 2))))
