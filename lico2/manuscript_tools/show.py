"""Show the runs of the paragraph containing a phrase: python show.py word/document.xml phrase ..."""
import re,sys
x=open(sys.argv[1],encoding='utf8').read()
def para_at(i):
    a=x.rfind('<w:p ',0,i); a2=x.rfind('<w:p>',0,i); a=max(a,a2); b=x.find('</w:p>',i)+6; return a,b
def show(phrase):
    i=x.find(phrase); assert i>=0, phrase
    a,b=para_at(i); p=x[a:b]
    toks=re.findall(r'(<m:oMath>.*?</m:oMath>|<w:r>.*?</w:r>|<w:r .*?</w:r>|<w:fldSimple.*?</w:fldSimple>|<w:commentRangeStart[^>]*>|<w:commentRangeEnd[^>]*>)',p,flags=re.S)
    print(f'--- para [{a}:{b}] tokens={len(toks)}')
    for t in toks:
        if t.startswith('<m:oMath'): print('  [MATH]', ''.join(re.findall(r'<m:t[^>]*>([^<]*)</m:t>',t)))
        elif t.startswith('<w:comment'): print('  ',t[:40])
        else:
            rpr=re.findall(r'<w:rPr>(.*?)</w:rPr>',t,flags=re.S); txt=''.join(re.findall(r'<w:t(?: [^>]*)?>([^<]*)</w:t>',t)); it=re.findall(r'<w:instrText[^>]*>([^<]*)',t); fc=re.findall(r'w:fldCharType="(\w+)"',t)
            print('  R', repr(txt), 'instr:'+it[0][:40] if it else '', fc if fc else '', ('rPr:'+rpr[0][:80]) if rpr else '')
for ph in sys.argv[2:]: show(ph)
