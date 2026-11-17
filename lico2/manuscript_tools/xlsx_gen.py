"""Generate the collocation-dependent content of LiCO2_levels_aDZ_aTZ.xlsx (rules inferred from the existing workbook and checked by
regenerating the existing sheets from the old result files, check_old.py)."""
import json, numpy as np
def state(label):
    t = [f"{'' if int(q) == 1 else q}ν{m}" for q, m in zip(label, '123456') if int(q) > 0]
    return ' + '.join(t) if t else 'ground'
def coll(fint, frange):
    """fint: results/intens_<tag>.json (intens.py); frange: results/ir_colloc_<tag>_range.json (mkrange.py)"""
    d = json.load(open(fint)); rg = json.load(open(frange))['rows']
    out = []
    for r, g in zip(d['rows'], rg):
        assert abs(r['nu'] - g['E']) < 1e-9 and r['assign'] == g['label']
        out.append(dict(n=r['level'], E=r['nu'], Elo=g['Elo'], Ehi=g['Ehi'], sym=r['sym'], label=r['assign'], state=state(r['assign']), w=r['weight'],
                        A=r['I_b'], Aa=r['I_a'], res=r['res'], flag='mixed' if r['weight'] < 0.4 else None))
    return dict(zpe=float(d['zpe']), rows=out)
def gvpt2(ws):
    out = []
    for r in range(5, ws.max_row + 1):
        if ws.cell(r, 2).value is None: continue
        out.append(dict(label=ws.cell(r, 2).value, state=ws.cell(r, 3).value, E=ws.cell(r, 4).value, A=ws.cell(r, 5).value, polyad=ws.cell(r, 7).value))
    return out
def matched(cadz, gadz, catz, gatz, keep_gvpt2_only=None):
    """rows: one per (label, k): k-th level carrying the label in each of the four lists"""
    src = [('aDZ collocation', cadz), ('aDZ GVPT2', gadz), ('aTZ collocation', catz), ('aTZ GVPT2', gatz)]
    occ, ncount = {}, {name: {} for name, _ in src}
    for name, lst in src:
        cnt = {}
        for r in lst:
            k = cnt.get(r['label'], 0) + 1; cnt[r['label']] = k; ncount[name][r['label']] = k
            occ.setdefault((r['label'], k), {})[name] = r
    rows = []
    for (lab, k), d in occ.items():
        order = ['aDZ GVPT2', 'aDZ collocation', 'aTZ GVPT2', 'aTZ collocation']
        key = next(d[n]['E'] for n in order if n in d)
        st = next(d[n]['state'] for n in order if n in d)
        notes = []
        for pre in ('aDZ', 'aTZ'):
            c = d.get(f'{pre} collocation'); g = d.get(f'{pre} GVPT2')
            if c is not None:
                if c['w'] < 0.4: notes.append(f'{pre} collocation state mixed (weight {c["w"]:.2f})')
                if ncount[f'{pre} collocation'][lab] > 1: notes.append(f'{ncount[f"{pre} collocation"][lab]} {pre} collocation levels carry this label')
            if g is not None and g.get('polyad'): notes.append(f'{pre} GVPT2 polyad: {g["polyad"]}')
        rows.append(dict(key=key, label=lab, state=st + (f' (#{k})' if k > 1 else ''), d=d, note='; '.join(notes) or None))
    if keep_gvpt2_only is not None:
        rows = [r for r in rows if any('collocation' in n for n in r['d']) or (r['label'], r['state']) in keep_gvpt2_only]
    rows.sort(key=lambda r: r['key'])
    return rows
