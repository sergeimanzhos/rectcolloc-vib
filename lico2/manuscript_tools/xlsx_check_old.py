"""Check that xlsx_gen.py regenerates results/LiCO2_levels_aDZ_aTZ.xlsx (manuscript version) exactly from the f100 result files."""
import os
import openpyxl, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from xlsx_gen import coll, gvpt2, matched
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
wb = openpyxl.load_workbook(os.path.join(R, 'LiCO2_levels_aDZ_aTZ.xlsx'))
C = {'adz': coll(f'{R}/intens_adz100_old_N14000.json', f'{R}/ir_colloc_adz100_old_N14000_range.json'),
     'atz': coll(f'{R}/intens_atz100_N14000.json', f'{R}/ir_colloc_atz100_N14000_range.json')}
bad = 0
for T, sh in (('adz', 'Collocation aDZ'), ('atz', 'Collocation aTZ')):
    ws = wb[sh]; c = C[T]
    if abs(ws['D2'].value - c['zpe']) > 1e-9: print('zpe', T)
    for i, r in enumerate(c['rows']):
        row = 5 + i
        exp = {1: r['n'], 2: r['E'], 3: r['Elo'], 4: r['Ehi'], 6: r['sym'], 7: r['label'], 8: r['state'], 9: r['w'], 10: r['A'], 11: r['Aa'], 12: r['res'], 13: r['flag']}
        for col, e in exp.items():
            v = ws.cell(row, col).value
            ok = (v is not None and abs(v - e) < 1e-9) if isinstance(e, float) else v == e
            if not ok:
                bad += 1
                if bad < 10: print(sh, row, col, repr(v), repr(e))
print('collocation sheets mismatches', bad)
G = {T: gvpt2(wb[f'GVPT2 {T[0]}{T[1:].upper()}']) for T in ('adz', 'atz')}
ws = wb['Matched by label']
keep = {(ws.cell(r, 1).value, ws.cell(r, 2).value) for r in range(5, ws.max_row + 1)}
rows = matched(C['adz']['rows'], G['adz'], C['atz']['rows'], G['atz'], keep)
old = [(ws.cell(r, 1).value, ws.cell(r, 2).value) for r in range(5, ws.max_row + 1)]
new = [(r['label'], r['state']) for r in rows]
print(len(new), len(old)); print('extra', [x for x in new if x not in old]); print('missing', [x for x in old if x not in new]); print('order same', new == old)
bad = 0
if new == old:
    for i, r in enumerate(rows):
        d = r['d']; g = lambda n, f: d[n][f] if n in d else None
        exp = {1: r['label'], 2: r['state'], 4: g('aDZ collocation', 'E'), 5: g('aDZ collocation', 'w'), 6: g('aDZ GVPT2', 'E'), 11: g('aTZ collocation', 'E'),
               12: g('aTZ collocation', 'w'), 13: g('aTZ GVPT2', 'E'), 20: g('aDZ collocation', 'A'), 21: g('aDZ GVPT2', 'A'), 22: g('aTZ collocation', 'A'),
               23: g('aTZ GVPT2', 'A'), 24: r['note']}
        for col, e in exp.items():
            v = ws.cell(5 + i, col).value
            ok = (v is not None and e is not None and abs(v - e) < 1e-9) if isinstance(e, float) else v == e
            if not ok:
                bad += 1
                if bad < 12: print('diff', 5 + i, col, repr(v)[:200], '|', repr(e)[:200])
print('matched mismatches', bad)
# levels sheet
ws = wb['Levels (energy order)']; bad = 0
for T, c0 in (('adz', 2), ('atz', 8)):
    for i, r in enumerate(C[T]['rows']):
        row = 8 + i
        for col, e in ((c0, r['E']), (c0 + 1, r['label']), (c0 + 2, r['A'])):
            v = ws.cell(row, col).value
            ok = abs(v - e) < 1e-9 if isinstance(e, float) else v == e
            if not ok: bad += 1; print('lev', row, col, v, e) if bad < 5 else None
print('levels mismatches', bad)
