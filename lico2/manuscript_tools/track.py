"""Tracked-change helper for word/document.xml (raw XML string editing; nothing outside the edited runs is touched)."""
import re
from xml.sax.saxutils import escape
AUTH = 'Claude'; DATE = '2026-11-14T00:00:00Z'
_id = [91000]
def nid():
    _id[0] += 1; return _id[0]
RUN = re.compile(r'<w:r(?: [^>]*)?>(?:(?!</w:r>).)*</w:r>', re.S)

def _run_at(x, pos):
    """the <w:r> element that contains position pos"""
    st = max(x.rfind('<w:r>', 0, pos), x.rfind('<w:r ', 0, pos))
    m = RUN.match(x, st); assert m and m.start() <= pos < m.end(), ('no run at', pos)
    return m.start(), m.end()

def _parts(run):
    m = re.match(r'(<w:r(?: [^>]*)?>)((?:<w:rPr>.*?</w:rPr>)?)(.*)</w:r>$', run, re.S)
    head, rpr, body = m.groups()
    ts = re.findall(r'<w:t(?: [^>]*)?>([^<]*)</w:t>', body)
    return head, rpr, body, ts

def _split(run):
    """rPr, body before the single w:t, its text, body after it (fields, breaks etc. are preserved)"""
    head, rpr, body, ts = _parts(run)
    assert len(ts) == 1, ('run must have exactly one w:t', ts)
    m = re.search(r'<w:t(?: [^>]*)?>([^<]*)</w:t>', body)
    return rpr, body[:m.start()], m.group(1), body[m.end():]

def _keep(rpr, pre, text, post):
    if not (pre or text or post): return ''
    t = f'<w:t xml:space="preserve">{text}</w:t>' if text else ''
    return f'<w:r>{rpr}{pre}{t}{post}</w:r>'

def _mk(rpr, text, kind, style=''):
    if style:
        add = {'i': '<w:i/><w:iCs/>', 'sub': '<w:vertAlign w:val="subscript"/>', 'sup': '<w:vertAlign w:val="superscript"/>',
               'plain': ''}[style]
        inner = re.sub(r'</?w:rPr>', '', rpr)
        inner = re.sub(r'<w:i/>|<w:iCs/>|<w:iCs w:val="0"/>|<w:vertAlign [^>]*/>', '', inner) if style in ('plain', 'i', 'sub', 'sup') else inner
        rpr = f'<w:rPr>{inner}{add}</w:rPr>' if (inner + add) else ''
    if kind == 'keep':
        return f'<w:r>{rpr}<w:t xml:space="preserve">{text}</w:t></w:r>' if text else ''
    if kind == 'del':
        return f'<w:del w:id="{nid()}" w:author="{AUTH}" w:date="{DATE}"><w:r>{rpr}<w:delText xml:space="preserve">{text}</w:delText></w:r></w:del>' if text else ''
    if kind == 'ins':
        return f'<w:ins w:id="{nid()}" w:author="{AUTH}" w:date="{DATE}"><w:r>{rpr}<w:t xml:space="preserve">{escape(text)}</w:t></w:r></w:ins>' if text else ''

def _del_whole(run):
    head, rpr, body, ts = _parts(run)
    assert not re.search(r'fldChar|instrText|<w:drawing|<w:object', body), 'refusing to delete field/drawing run'
    body2 = re.sub(r'<w:t( [^>]*)?>', '<w:delText xml:space="preserve">', body).replace('</w:t>', '</w:delText>')
    return f'<w:del w:id="{nid()}" w:author="{AUTH}" w:date="{DATE}"><w:r>{rpr}{body2}</w:r></w:del>'

def para_bounds(x, ctx):
    assert x.count(ctx) == 1, ('ctx not unique', ctx, x.count(ctx))
    i = x.find(ctx); a = max(x.rfind('<w:p ', 0, i), x.rfind('<w:p>', 0, i)); b = x.find('</w:p>', i) + 6
    return a, b

def replace_span(x, start_anchor, start_sub, end_anchor, end_sub, new_parts, ins_style_base=None, ctx=None):
    if ctx is not None:
        a, b = para_bounds(x, ctx)
        return x[:a] + replace_span(x[a:b], start_anchor, start_sub, end_anchor, end_sub, new_parts, ins_style_base) + x[b:]
    """Delete text from the first occurrence of start_sub inside the w:t containing start_anchor (unique in x) through the end of
    end_sub inside the w:t containing end_anchor (unique, at or after the start), and insert new_parts = [(text, style)] there.
    style: '' = formatting of the start run; 'i', 'sub', 'sup', 'plain' modify it. Text in new_parts is plain (escaped here)."""
    assert x.count(start_anchor) == 1, ('start anchor not unique', start_anchor, x.count(start_anchor))
    assert x.count(end_anchor) == 1, ('end anchor not unique', end_anchor, x.count(end_anchor))
    pa, pb = x.find(start_anchor), x.find(end_anchor)
    a0, a1 = _run_at(x, pa); b0, b1 = _run_at(x, pb)
    assert b0 >= a0
    rA, preA, tA, postA = _split(x[a0:a1])
    ia = tA.find(start_sub); assert ia >= 0 and tA.count(start_sub) == 1, ('start_sub', start_sub, tA)
    base = rA if ins_style_base is None else ins_style_base
    ins = ''.join(_mk(base, t, 'ins', s) for t, s in new_parts)
    if a0 == b0:
        ie = tA.find(end_sub, ia); assert ie >= 0, ('end_sub', end_sub)
        ie += len(end_sub)
        new = _keep(rA, preA, tA[:ia], '') + _mk(rA, tA[ia:ie], 'del') + ins + _keep(rA, '', tA[ie:], postA)
        return x[:a0] + new + x[a1:]
    assert not postA, 'start run has content after its text'
    rB, preB, tB, postB = _split(x[b0:b1]); assert not preB, 'end run has content before its text'
    ie = tB.find(end_sub); assert ie >= 0 and tB.count(end_sub) == 1, ('end_sub', end_sub, tB); ie += len(end_sub)
    mid = x[a1:b0]
    assert '<m:oMath' not in mid and 'fldChar' not in mid and 'instrText' not in mid and '<w:p ' not in mid and '<w:p>' not in mid, 'span crosses math/field/paragraph'
    mid2 = RUN.sub(lambda m: _del_whole(m.group(0)), mid)
    new = _keep(rA, preA, tA[:ia], '') + _mk(rA, tA[ia:], 'del') + mid2 + _mk(rB, tB[:ie], 'del') + ins + _keep(rB, '', tB[ie:], postB)
    return x[:a0] + new + x[b1:]

def insert_after(x, anchor, after_sub, new_parts):
    """tracked insertion right after after_sub inside the w:t containing anchor"""
    return replace_span(x, anchor, after_sub, anchor, after_sub, [(after_sub, '__keep__')] + new_parts) if False else _ins_after(x, anchor, after_sub, new_parts)

def ins_end_of_para(x, ctx, new_parts, base_anchor):
    """tracked insertion at the end of the paragraph containing ctx, formatted like the run containing base_anchor"""
    a, b = para_bounds(x, ctx); p = x[a:b]
    r0, r1 = _run_at(p, p.find(base_anchor)); rA = _parts(p[r0:r1])[1]
    ins = ''.join(_mk(rA, t, 'ins', s) for t, s in new_parts)
    return x[:b - 6] + ins + x[b - 6:]

def _ins_after(x, anchor, after_sub, new_parts):
    assert x.count(anchor) == 1, ('anchor not unique', anchor)
    a0, a1 = _run_at(x, x.find(anchor)); rA, preA, tA, postA = _split(x[a0:a1])
    i = tA.find(after_sub); assert i >= 0 and tA.count(after_sub) == 1, ('after_sub', after_sub); i += len(after_sub)
    ins = ''.join(_mk(rA, t, 'ins', s) for t, s in new_parts)
    return x[:a0] + _keep(rA, preA, tA[:i], '') + ins + _keep(rA, '', tA[i:], postA) + x[a1:]

def ins_before(x, anchor, before_sub, new_parts):
    """tracked insertion right before before_sub inside the w:t containing anchor, formatted like that run"""
    assert x.count(anchor) == 1, ('anchor not unique', anchor)
    a0, a1 = _run_at(x, x.find(anchor)); rA, preA, tA, postA = _split(x[a0:a1])
    i = tA.find(before_sub); assert i >= 0 and tA.count(before_sub) == 1, ('before_sub', before_sub)
    ins = ''.join(_mk(rA, t, 'ins', s) for t, s in new_parts)
    return x[:a0] + _keep(rA, preA, tA[:i], '') + ins + _keep(rA, '', tA[i:], postA) + x[a1:]
