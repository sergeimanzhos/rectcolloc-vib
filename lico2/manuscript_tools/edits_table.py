"""Table 3 (LiCO2): replace the collocation DZ/TZ cells (columns 4 and 5) with the new values, UNTRACKED.
usage: python edits_table.py word/document.xml   (unpacked manuscript docx; 18 Nov 2026 update)"""
import re, json, sys
X = sys.argv[1]; x = open(X, encoding='utf8').read()
import os; new = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'table3_new.json')))
i = x.find('Li out'); a = x.rfind('<w:tbl>', 0, i); b = x.find('</w:tbl>', i); tbl = x[a:b]
txt = lambda c: ''.join(re.findall(r'<w:t(?: [^>]*)?>([^<]*)</w:t>', c))
def setcell(c, s):
    ts = list(re.finditer(r'(<w:t(?: [^>]*)?>)([^<]*)(</w:t>)', c)); assert ts
    out = c[:ts[0].start()] + f'<w:t xml:space="preserve">{s}</w:t>'; pos = ts[0].end()
    for m in ts[1:]: out += c[pos:m.start()] + '<w:t xml:space="preserve"></w:t>'; pos = m.end()
    return out + c[pos:]
rows = re.findall(r'<w:tr[ >].*?</w:tr>', tbl, flags=re.S); n = 0; nt = tbl
for row in rows:
    cells = re.findall(r'<w:tc>.*?</w:tc>', row, flags=re.S)
    band = txt(cells[0]).strip()
    ks = [k for k in new if band == k or band.startswith(k + ' ')]; key = max(ks, key=len) if ks else None
    if key is None or len(cells) < 8: continue
    r2 = row
    for ci, s in ((3, new[key][0]), (4, new[key][1])):
        print(f'{band[:22]:24s} col{ci}: {txt(cells[ci])!r:32s} -> {s!r}  runs={len(re.findall("<w:t[ >]", cells[ci]))}')
        r2 = r2.replace(cells[ci], setcell(cells[ci], s), 1); n += 1
    nt = nt.replace(row, r2, 1)
assert n == 22, n
open(X, 'w', encoding='utf8').write(x[:a] + nt + x[b:]); print('cells replaced', n)
