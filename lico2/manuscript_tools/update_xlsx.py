"""Update LiCO2_levels_aDZ_aTZ.xlsx (version sent by Sergei, 17 Nov 2026) with the max-residual collocation results (fmax_* runs, 17-18 Nov 2026).
GVPT2 sheets unchanged. Generator rules (xlsx_gen.py) checked by regenerating the existing sheets exactly from the old results (check_old.py).
Usage: python update_xlsx.py in.xlsx out.xlsx   (in.xlsx = the 17 Nov version, git show 54d46d0:lico2/results/LiCO2_levels_aDZ_aTZ.xlsx)"""
import openpyxl, json, sys, copy, re
import os; HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); from xlsx_gen import coll, gvpt2, matched
R = os.path.join(HERE, '..', 'results')
src, dst = sys.argv[1], sys.argv[2]
wb = openpyxl.load_workbook(src)
C = {T: coll(f'{R}/intens_fmax_{T}_N14000.json', f'{R}/ir_colloc_fmax_{T}_N14000_range.json') for T in ('adz', 'atz')}
M = {T: json.load(open(f'{R}/fmax_{T}_N14000.json'))['recs'][-1]['M'] for T in ('adz', 'atz')}
# --- Collocation sheets
for T, sh in (('adz', 'Collocation aDZ'), ('atz', 'Collocation aTZ')):
    ws = wb[sh]; c = C[T]
    old = ws['A1'].value; assert 'Vmax = 5500 cm-1)' in old
    ws['A1'] = old.replace('Vmax = 5500 cm-1)', 'Vmax = 5000 cm-1; parameters tuned by the maximum validation residual)')
    ws['D2'] = float(c['zpe'])
    for i, r in enumerate(c['rows']):
        row = 5 + i; assert ws.cell(row, 1).value == r['n']
        for col, v in zip((2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13),
                          (r['E'], r['Elo'], r['Ehi'], r['sym'], r['label'], r['state'], r['w'], r['A'], r['Aa'], r['res'], r['flag'])):
            ws.cell(row, col).value = float(v) if isinstance(v, float) else v
# --- Levels sheet (collocation columns B-D aDZ, H-J aTZ)
ws = wb['Levels (energy order)']
for T, c0 in (('adz', 2), ('atz', 8)):
    ws.cell(6, c0).value = float(C[T]['zpe'])
    for i, r in enumerate(C[T]['rows']):
        row = 8 + i; assert ws.cell(row, 1).value == r['n']
        ws.cell(row, c0).value = float(r['E']); ws.cell(row, c0 + 1).value = r['label']; ws.cell(row, c0 + 2).value = float(r['A'])
# --- Matched by label (rows regenerated; GVPT2-only rows: same selection as the existing sheet; formulas copied from row 5 pattern)
ws = wb['Matched by label']
keep = {(ws.cell(r, 1).value, ws.cell(r, 2).value) for r in range(5, ws.max_row + 1)}
G = {T: gvpt2(wb[f'GVPT2 {T[0]}{T[1:].upper()}']) for T in ('adz', 'atz')}
rows = matched(C['adz']['rows'], G['adz'], C['atz']['rows'], G['atz'], keep)
NC = 24
style = {col: copy.copy(ws.cell(5, col)._style) for col in range(1, NC + 1)}
FCOL = [3, 7, 8, 9, 10, 14, 15, 16, 17, 18, 19]
form = {col: ws.cell(5, col).value for col in FCOL}
for col in FCOL: assert isinstance(form[col], str) and form[col].startswith('=')
old_n = ws.max_row - 4
if old_n > len(rows): ws.delete_rows(5 + len(rows), old_n - len(rows))
for i, r in enumerate(rows):
    row = 5 + i; d = r['d']; g = lambda n, f: float(d[n][f]) if n in d else None
    vals = {1: r['label'], 2: r['state'], 4: g('aDZ collocation', 'E'), 5: g('aDZ collocation', 'w'), 6: g('aDZ GVPT2', 'E'), 11: g('aTZ collocation', 'E'),
            12: g('aTZ collocation', 'w'), 13: g('aTZ GVPT2', 'E'), 20: g('aDZ collocation', 'A'), 21: g('aDZ GVPT2', 'A'), 22: g('aTZ collocation', 'A'),
            23: g('aTZ GVPT2', 'A'), 24: r['note']}
    for col in range(1, NC + 1):
        cell = ws.cell(row, col); cell._style = copy.copy(style[col])
        cell.value = re.sub(r'(?<![0-9$])([A-Z]{1,2})5(?![0-9])', lambda m: f'{m.group(1)}{row}', form[col].replace('$A5', f'$A{row}')) if col in FCOL else vals[col]
print('matched rows', old_n, '->', len(rows))
# --- Notes
ws = wb['Notes']
a8 = ws['A8'].value
old8 = 'Gaussian basis, N = 14000, centres with V <= 900 cm-1, widths 0.8935/1.0395/0.99/0.99/1.0395/0.99; energy and gradient rows with constant weight; Vmax = 5500 cm-1; 16150 (aDZ) and 13831 (aTZ) DFT points.'
assert old8 in a8
ws['A8'] = a8.replace(old8, 'Gaussian basis, N = 14000, centres with V <= 1200 cm-1, widths 0.9321/0.99/1.00/1.0188/1.00/1.00; energy and gradient rows with '
                      f'constant weight (p = 0); Vmax = 5000 cm-1; {M["adz"]} (aDZ) and {M["atz"]} (aTZ) DFT points. Centre cutoff, Vmax, p and widths were '
                      'tuned on the aDZ data at N = 10000, M = 12000 by minimizing the maximum validation residual of the 100 levels '
                      '(1275 validation points with V <= 5000 cm-1); the width tuning was stopped in the seventh 1% sweep (17 November 2026).')
a9 = ws['A9'].value; assert 'over four runs: N = 10000 and 14000, each on two independent DFT point sets.' in a9
ws['A9'] = a9.replace('over four runs: N = 10000 and 14000, each on two independent DFT point sets.',
                      'over four runs: N = 10000 and 14000, each on two independent DFT point sets (second set: 16462 aDZ and 16000 aTZ points), all with '
                      'the parameters of note 8.')
ws['A13'] = ('Convergence of the collocation levels (range width over the four runs): the 11 levels below 1000 cm-1 to within 0.3 cm-1 and the '
             'lowest 19 levels to within 1 cm-1; 57 (aDZ) and 54 (aTZ) of the 99 levels within 1 cm-1, 70 within 2 cm-1; median width 0.8 cm-1; '
             'largest widths 12.4 (aDZ) and 7.6 (aTZ) cm-1 near 2100 cm-1. See the range columns.')
a14 = ws['A14'].value
A = {T: [r for r in C[T]['rows'] if r['label'] == '000001'] for T in ('adz', 'atz')}
s = sum(r['A'] for r in A['adz']); P = [r for r in C['atz']['rows'] if r['label'] == '012000' and r['sym'] == 'B2' and 1540 < r['E'] < 1600][0]
ws['A14'] = (f"aDZ: the asymmetric CO stretch v6 shares its intensity with a second B2 state {A['adz'][1]['E'] - A['adz'][0]['E']:.0f} cm-1 above it "
             f"(both carry label 000001; {A['adz'][0]['A']:.0f} and {A['adz'][1]['A']:.0f} km/mol); the split ratio is not converged, the sum "
             f"(about {s:.0f} km/mol) is. At aTZ the partner (label 012000, {P['E']:.1f} cm-1) keeps only {P['A']:.0f} km/mol.")
ws['A16'] = ('Source of the numbers: collocation (updated 18 November 2026, parameters tuned by the maximum validation residual): rectcolloc-vib/lico2/'
             'results intens_fmax_adz_N14000.json, intens_fmax_atz_N14000.json, ir_colloc_fmax_adz_N14000_range.json, ir_colloc_fmax_atz_N14000_range.json; '
             'ranges from fmax_<basis>_N14000/_N10000 and fmax_<basis>7_N14000/_N10000 (calculations of 17-18 November 2026 on the workstation). '
             'GVPT2 (unchanged): vpt2_adz_hi.json, vpt2_atz_hi.json (aDZ recomputed 14 October 2026, see the previous version of this note). '
             'Workbook generated 13 October 2026; collocation columns updated 18 November 2026.')
wb.save(dst); print('saved', dst)
print('old note 14:', a14)
