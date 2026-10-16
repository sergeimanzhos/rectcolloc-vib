"""Assemble the free-CO2 level table: harmonic | GVPT2 | VCI (quartic FF) | collocation on the quartic FF | collocation on the DFT points.
Levels are paired within each l block by energy rank (zero-order states sorted by harmonic energy).  python lev_table.py [EMAX]"""
import sys, os, json, numpy as np
EMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 4350.0; VT = sys.argv[2] if len(sys.argv) > 2 else 'co2fit'; SETS = sys.argv[3].split(',') if len(sys.argv) > 3 else ['DFT']
vp = json.load(open('vptlev_%s_K14.json' % VT)); vp12 = json.load(open('vptlev_%s_K12.json' % VT)); w = vp['w']; wb, w1, w3 = 0.5 * (w[0] + w[1]), w[2], w[3]
def bylist(L):
    out = {l: [] for l in range(5)}
    for r in L:
        l = int(round(np.sqrt(max(r['l2'], 0))))
        if l <= 4 and not any(abs(e - r['E']) < 0.02 for e in out[l]): out[l].append(r['E'])
    return {l: sorted(v) for l, v in out.items()}
vci, vci12, gv = bylist(vp['vci']), bylist(vp12['vci']), bylist(vp['gvpt2'])
zero = {l: sorted([(v1 * w1 + v2 * wb + v3 * w3, '%d %d^%d %d' % (v1, v2, l, v3)) for v1 in range(6) for v3 in range(4) for v2 in range(l, 12, 2) if v1 + v2 + v3 > 0]) for l in range(5)}
def load(f): return json.load(open(f)) if os.path.exists(f) else []
ff, dft, dg = load('lev_ff.json'), load('lev_dft.json'), load('lev_dftgrad.json')
def pick(runs, kind, grad):
    c = [r for r in runs if r['kind'] == kind and bool(r['grad']) == grad]; return sorted(c, key=lambda r: (r['N'][0], r['M']))
cols = [('GVPT2', gv, None), ('VCI', vci, vci12)]
for name, runs in [x for x in (('QFF', ff), ('DFT', dft + dg)) if x[0] in SETS]:
    for kind, grad, lab in (('h', False, 'HO E'), ('g', False, 'Gauss E'), ('g', True, 'Gauss E+grad')):
        c = pick(runs, kind, grad)
        if c:
            lev = lambda r: {int(l): v for l, v in r['lev'].items()}
            cols.append(('%s %s' % (name, lab), lev(c[-1]), lev(c[-2]) if len(c) > 1 else None)); print('column %-22s: %s   [previous: %s]' % (cols[-1][0], c[-1]['tag'], c[-2]['tag'] if len(c) > 1 else '-'))
            print('     ZPE %.3f' % c[-1]['zpe'] + ('  (previous %.3f)' % c[-2]['zpe'] if len(c) > 1 else ''))
rows = []
for l in range(5):
    for i, e in enumerate(vci[l]):
        if e > EMAX: break
        rows.append((e, l, i))
rows.sort()
get = lambda d, l, i: d[l][i] if d and l in d and i < len(d[l]) else None
f = lambda v: '%.2f' % v if v is not None else '-'
hdr = ['#', 'l', 'zero-order state', 'harmonic'] + [c[0] for c in cols]; out = ['| ' + ' | '.join(hdr) + ' |', '|' + '---|' * len(hdr)]; csv = [','.join(hdr)]
for n, (e, l, i) in enumerate(rows, 1):
    cells = ['%d' % n, '%d' % l, zero[l][i][1].replace('^', '^').replace(' ', ' '), '%.2f' % zero[l][i][0]] + [f(get(c[1], l, i)) for c in cols]; out.append('| ' + ' | '.join(cells) + ' |'); csv.append(','.join(cells))
print('\n'.join(out))
print('\nconvergence: max |last - previous| over the %d tabulated levels' % len(rows))
for name, cur, prev in cols:
    if prev is None: continue
    d = [abs(get(cur, l, i) - get(prev, l, i)) for e, l, i in rows if get(cur, l, i) is not None and get(prev, l, i) is not None]
    d2 = [abs(get(cur, l, i) - get(prev, l, i)) for e, l, i in rows if e <= 3000 and get(cur, l, i) is not None and get(prev, l, i) is not None]
    print('  %-22s all: %.2f   levels below 3000: %.2f' % (name, max(d), max(d2)))
print('\ndeviation from VCI (quartic-FF columns) / from DFT HO column (DFT columns): max |diff| all, below 3000')
for name, cur, prev in cols:
    ref = vci if (name.startswith('QFF') or name == 'GVPT2') else ([c for c in cols if c[0] == 'DFT HO E'] or [None])[0]
    if ref is None or name in ('VCI', 'DFT HO E'): continue
    ref = ref if isinstance(ref, dict) else ref[1]
    d = [(abs(get(cur, l, i) - get(ref, l, i)), e) for e, l, i in rows if get(cur, l, i) is not None and get(ref, l, i) is not None]
    print('  %-22s all: %.2f   levels below 3000: %.2f' % (name, max(x[0] for x in d), max(x[0] for x in d if x[1] <= 3000)))
open('co2_levels_%s.md' % '_'.join(SETS), 'w').write('\n'.join(out) + '\n'); open('co2_levels_%s.csv' % '_'.join(SETS), 'w').write('\n'.join(csv) + '\n')
