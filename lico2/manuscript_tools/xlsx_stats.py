"""Statistics quoted in the manuscript (18 Nov 2026), OLD (manuscript) and NEW (max-residual) collocation:
colloc - GVPT2 (levels matched by label, colloc HO weight >= 0.4, E <= 2200), aTZ - aDZ (both weights >= 0.4), convergence range widths.
usage: python xlsx_stats.py workbook.xlsx   (the GVPT2 sheets and the matched-label set are read from it; 17 Nov version)"""
import openpyxl, numpy as np, sys
import os; HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); from xlsx_gen import coll, gvpt2, matched
R = os.path.join(HERE, '..', 'results')
wb = openpyxl.load_workbook(sys.argv[1])
G = {T: gvpt2(wb[f'GVPT2 {T[0]}{T[1:].upper()}']) for T in ('adz', 'atz')}
ws = wb['Matched by label']; keep = {(ws.cell(r, 1).value, ws.cell(r, 2).value) for r in range(5, ws.max_row + 1)}
def stats(tag, fa, fra, ft, frt):
    C = {'adz': coll(f'{R}/{fa}', f'{R}/{fra}'), 'atz': coll(f'{R}/{ft}', f'{R}/{frt}')}
    rows = matched(C['adz']['rows'], G['adz'], C['atz']['rows'], G['atz'], keep)
    print('=====', tag, 'ZPE', C['adz']['zpe'], C['atz']['zpe'])
    for b in ('aDZ', 'aTZ'):
        d = np.array([r['d'][f'{b} collocation']['E'] - r['d'][f'{b} GVPT2']['E'] for r in rows if f'{b} collocation' in r['d'] and f'{b} GVPT2' in r['d']
                      and r['d'][f'{b} collocation']['w'] >= 0.4 and r['d'][f'{b} collocation']['E'] <= 2200])
        print(f'{b} colloc-GVPT2 (w>=0.4, E<=2200): n {len(d)} MAE {abs(d).mean():.2f} RMS {np.sqrt((d**2).mean()):.2f} max {abs(d).max():.2f}')
    d = np.array([r['d']['aTZ collocation']['E'] - r['d']['aDZ collocation']['E'] for r in rows if 'aDZ collocation' in r['d'] and 'aTZ collocation' in r['d']
                  and r['d']['aDZ collocation']['w'] >= 0.4 and r['d']['aTZ collocation']['w'] >= 0.4 and r['d']['aDZ collocation']['E'] <= 2200])
    print(f'aTZ-aDZ colloc (both w>=0.4, aDZ E<=2200): n {len(d)} MAE {abs(d).mean():.2f} RMS {np.sqrt((d**2).mean()):.2f} max {abs(d).max():.2f}')
    for b in ('adz', 'atz'):
        x = np.array([r['Ehi'] - r['Elo'] for r in C[b]['rows']]); E = np.array([r['E'] for r in C[b]['rows']])
        print(f'{b} range width over 4 runs: n<1 {(x<1).sum()} n<2 {(x<2).sum()} of {len(x)}; first k levels <1: {np.argmax(x>=1)}; below 1000 max {x[E<1000].max():.2f} (n {(E<1000).sum()}); '
              f'RMS {np.sqrt((x**2).mean()):.2f} max {x.max():.2f} (level at {E[x.argmax()]:.1f}); median {np.median(x):.2f}')
    return C, rows
stats('OLD (manuscript)', 'intens_adz100_old_N14000.json', 'ir_colloc_adz100_old_N14000_range.json', 'intens_atz100_N14000.json', 'ir_colloc_atz100_N14000_range.json')
stats('NEW (max residual)', 'intens_fmax_adz_N14000.json', 'ir_colloc_fmax_adz_N14000_range.json', 'intens_fmax_atz_N14000.json', 'ir_colloc_fmax_atz_N14000_range.json')
