"""Tracked text edits (author Claude) of the LiCO2 section for the max-residual results (18 Nov 2026).
usage: python edits_text.py word/document.xml   (after edits_table.py; text anchors are those of the 17 Nov manuscript version)"""
import re, sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import track as T
X = sys.argv[1]; x = open(X, encoding='utf8').read()
P1 = 'Collocation calculations of the lowest 100 vibrational levels were done'

# 1. p = 1 -> p = 0 inside the OMath (tracked math runs)
a, b = T.para_bounds(x, P1); p = x[a:b]
om = '<m:oMath><m:r><w:rPr><w:rFonts w:ascii="Cambria Math" w:hAnsi="Cambria Math"/></w:rPr><m:t>p=1</m:t></m:r></m:oMath>'
assert p.count(om) == 1
rp = '<w:rPr><w:rFonts w:ascii="Cambria Math" w:hAnsi="Cambria Math"/></w:rPr>'
om2 = (f'<m:oMath><m:r>{rp}<m:t>p=</m:t></m:r>'
       f'<w:del w:id="{T.nid()}" w:author="{T.AUTH}" w:date="{T.DATE}"><m:r>{rp}<m:t>1</m:t></m:r></w:del>'
       f'<w:ins w:id="{T.nid()}" w:author="{T.AUTH}" w:date="{T.DATE}"><m:r>{rp}<m:t>0</m:t></m:r></w:ins></m:oMath>')
x = x[:a] + p.replace(om, om2) + x[b:]

# 2. tuning procedure sentence
x = T.replace_span(x, ' was found to be best for this floppier molecule. Basis parameters', 'Basis parameters, energy cutoffs, ',
                   ' were first pre-tuned using a harmonic PES', 'tuned further using DFT collocation points.',
                   [('The centre cutoff ', ''), ('V', 'i'), ('max', 'sub'), ('′', 'sup'), (', ', ''), ('V', 'i'), ('max', 'sub'), (', ', ''), ('p', 'i'),
                    (' and the Gaussian widths were tuned on the UB3LYP-D3(BJ)/aug-cc-pVDZ collocation points with the procedure of Section 2.3 (', ''),
                    ('J', 'i'), (' = ', ''), ('r', 'i'), ('max', 'sub'), (') at ', ''), ('N', 'i'), (' = 10000 and ', ''), ('M', 'i'),
                    (' = 12000, with 1275 validation points: ', ''), ('V', 'i'), ('max', 'sub'), ('′', 'sup'),
                    (' was scanned from 700 to 1500 cm⁻¹, then ', ''), ('V', 'i'), ('max', 'sub'), (' = 5000 and 5500 cm⁻¹ and ', ''), ('p', 'i'),
                    (' = 0, 1 and 3 were compared, the common width was scanned in steps of 0.05, and the widths were then optimized by coordinate '
                     'descent with 5% and then 1% steps (stopped during the seventh 1% sweep). The same parameters were used for aug-cc-pVTZ.', '')], ctx=P1)

# 3. Vmax, V'max values; widths
x = T.replace_span(x, ' were 5500 and 900 cm', '5500 and 900', ' were 5500 and 900 cm', '5500 and 900', [('5000 and 1200', '')], ctx=P1)
x = T.replace_span(x, 'were 0.8935, 1.0395, 0.99, 0.99, 1.0395, and 0.99, respectively', '0.8935, 1.0395, 0.99, 0.99, 1.0395, and 0.99',
                   'were 0.8935, 1.0395, 0.99, 0.99, 1.0395, and 0.99, respectively', '0.8935, 1.0395, 0.99, 0.99, 1.0395, and 0.99',
                   [('0.9321, 0.99, 1.00, 1.0188, 1.00, and 1.00', '')], ctx=P1)

# 4. N, M and convergence statements
R = ' of about 14000 (the results shown in the figure) were sufficient to converge the first 25 levels to better than 1 cm'
x = T.replace_span(x, R, 'of about 14000', R, 'of about 14000', [('= 14000 and about 12000 (12248 points at aug-cc-pVDZ, 11690 at aug-cc-pVTZ), respectively', '')], ctx=P1)
R = ' (the results shown in the figure) were sufficient to converge the first 25 levels to better than 1 cm'
x = T.replace_span(x, R, '25', R, '25', [('19', '')], ctx=P1)
R = '; Of all levels, 61 levels were converged within 1 cm⁻¹, and 78 within 2 cm⁻¹.'
x = T.replace_span(x, R, '61 levels', R, '61 levels', [('57 (aDZ) and 54 (aTZ) levels', '')], ctx=P1)
R = ' were converged within 1 cm⁻¹, and 78 within 2 cm⁻¹.'
x = T.replace_span(x, R, '78', R, '78', [('70', '')], ctx=P1)
R = ', are converged to 0.2 cm⁻¹ or better. In '
x = T.replace_span(x, R, '0.2', R, '0.2', [('0.3', '')], ctx=P1)
# convergence definition (asked for in comment 25) + deletion of the Figure 6 convergence-bar sentence
R = ' cm⁻¹ or better. In '
x = T.replace_span(x, R, 'In ', R, 'In ',
                   [('The convergence was estimated as the range of each level over four calculations: ', ''), ('N', 'i'),
                    (' = 10000 and 14000, each with two independent sets of DFT points. ', '')], ctx=P1)
# delete the runs of "Figure 6, for the collocation calculations, ... for most lines" (incl. the REF field) and the ")." that follows
a, b = T.para_bounds(x, P1); p = x[a:b]
s = p.find('<w:commentRangeStart w:id="25"/>'); e = p.find('<w:commentRangeEnd w:id="25"/>'); assert 0 < s < e
def delrun(m):
    r = m.group(0)
    if '<w:del ' in r or '<w:ins ' in r: return r
    head, rpr, body, ts = T._parts(r)
    body2 = re.sub(r'<w:t( [^>]*)?>', '<w:delText xml:space="preserve">', body).replace('</w:t>', '</w:delText>')
    body2 = re.sub(r'<w:instrText( [^>]*)?>', r'<w:delInstrText\1>', body2).replace('</w:instrText>', '</w:delInstrText>')
    return f'<w:del w:id="{T.nid()}" w:author="{T.AUTH}" w:date="{T.DATE}"><w:r>{rpr}{body2}</w:r></w:del>'
mid = T.RUN.sub(delrun, p[s:e])
p = p[:s] + mid + p[e:]
x = x[:a] + p + x[b:]
R = '). All levels are given in Supporting Information. '
x = T.replace_span(x, R, ').', R, ').', [], ctx=P1)

# 5. agreement paragraph: convergence numbers
R = 'the collocation convergence error (0.2 cm⁻¹ below 1000 cm⁻¹, 1.9 cm⁻¹ RMS over all 99 levels)'
x = T.replace_span(x, R, '0.2 cm⁻¹ below 1000 cm⁻¹, 1.9 cm⁻¹ RMS', R, '0.2 cm⁻¹ below 1000 cm⁻¹, 1.9 cm⁻¹ RMS',
                   [('0.3 cm⁻¹ below 1000 cm⁻¹, 2.4–2.9 cm⁻¹ RMS', '')])

# 6. Figure 6 caption: remove the convergence-range sentence
R = ' For collocation calculations, the convergence range with respect to N and M is shown in pale color (invisibly small for most lines). '
x = T.replace_span(x, R, 'For collocation', R, 'most lines). ', [])

# 7. basis-set paragraph
P7 = 'Increasing the DFT basis size from aug-cc-pVDZ'
R = 'maximum about 25 cm⁻¹). The shift is systematic'
x = T.replace_span(x, R, 'about 25', R, 'about 25', [('about 27', '')], ctx=P7)
R = 'ν4, ν5 and ν6 increase by 9–10 cm⁻¹'
x = T.replace_span(x, R, '9–10', R, '9–10', [('9–11', '')], ctx=P7)
R = 'the first 25 levels to 1 cm⁻¹ or better, and all 99 levels to an RMS of 1.9 cm⁻¹ (maximum 9 cm⁻¹ at -pVDZ, 7 cm⁻¹ at -pVTZ)'
x = T.replace_span(x, R, '25 levels', R, '25 levels', [('19 levels', '')], ctx=P7)
R = ' to 1 cm⁻¹ or better, and all 99 levels to an RMS of 1.9 cm⁻¹ (maximum 9 cm⁻¹ at -pVDZ, 7 cm⁻¹ at -pVTZ)'
x = T.replace_span(x, R, '1.9 cm⁻¹ (maximum 9 cm⁻¹ at -pVDZ, 7 cm⁻¹', R, '1.9 cm⁻¹ (maximum 9 cm⁻¹ at -pVDZ, 7 cm⁻¹',
                   [('2.9 cm⁻¹ at -pVDZ and 2.4 cm⁻¹ at -pVTZ (maximum 12 cm⁻¹ at -pVDZ, 8 cm⁻¹', '')], ctx=P7)
R = 'its intensity splits 435:147 km/mol between levels about 5 cm⁻¹ apart'
x = T.replace_span(x, R, '435:147', R, '435:147', [('354:224', '')], ctx=P7)
R = ' km/mol between levels about 5 cm⁻¹ apart'
x = T.replace_span(x, R, '5 cm⁻¹', R, '5 cm⁻¹', [('4 cm⁻¹', '')], ctx=P7)
R = 'the partner keeps only 13 km/mol, giving a single ν6 band at 1570.3 cm⁻¹'
x = T.replace_span(x, R, '13', R, '13', [('11', '')], ctx=P7)
R = ' km/mol, giving a single ν6 band at 1570.3 cm⁻¹'
x = T.replace_span(x, R, '1570.3', R, '1570.3', [('1570.1', '')], ctx=P7)
open(X, 'w', encoding='utf8').write(x); print('text edits done')
