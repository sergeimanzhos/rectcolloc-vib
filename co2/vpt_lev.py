"""GVPT2 and VCI level lists (all states up to EMAX above the ground state) on the quartic force field of vpt2_cyl.py, with the
vibrational angular momentum l of every level from <l_z^2>, l_z = q_a p_b - q_b p_a (bend components a, b = modes 0, 1).
  KMAX=12 python vpt_lev.py co2cyl [EMAX=4400]      -> vptlev_<tag>_K<K>.json"""
import sys, os, json, numpy as np
src = open('vpt2_cyl.py').read().split('# ---------- GVPT2 ----------')[0].replace("json.dump(C, open('vff_cyl.json', 'w'), indent=1)", 'pass')
# FFJSON=<file>: use the cylindrical quartic force field of that file (format of vff_cyl.json) instead of the finite-difference one
src = src.replace("    cub = {(0, 0, 2): kbb1,", "    if 'FFJSON' in os.environ: C = json.load(open(os.environ['FFJSON'])); kbb1, pb, p11b, p33b = C['k1bb'], C['pbbbb'], C['p11bb'], C['p33bb']" + chr(10) + "    cub = {(0, 0, 2): kbb1,")
assert 'FFJSON' in src
exec(src)
EMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 4400.0
Rz = op({0: qpow[1], 1: dq1}) - op({1: qpow[1], 0: dq1}); L2 = -(Rz @ Rz); g0 = int(np.argmin(E0)); lab = lambda i: ''.join('%d' % v for v in S[i])
# ---------- GVPT2: targets = all states with harmonic energy <= EMAX, polyad partners up to EMAX + 300 ----------
targets = [i for i in range(NS) if (E0[i] - E0[g0]) * HC <= EMAX]; cand = [i for i in range(NS) if (E0[i] - E0[g0]) * HC <= EMAX + 300]; H1sq = H1 @ H1
parent = {i: i for i in cand}
def find(i):
    while parent[i] != i: parent[i] = parent[parent[i]]; i = parent[i]
    return i
for a in targets:
    for b in cand:
        if b != a and abs(E0[a] - E0[b]) < THR and (abs(H1[a, b]) > 1e-9 or abs(H2[a, b]) > 1e-9 or abs(H1sq[a, b]) > 1e-12): parent[find(a)] = find(b)
groups = {}
for i in cand: groups.setdefault(find(i), []).append(i)
gv = []
for g in groups.values():
    if not any(t in g for t in targets): continue
    P = list(g); Ps = set(P); Q = np.array([i for i in range(NS) if i not in Ps]); He = np.zeros((len(P), len(P)))
    for x, a in enumerate(P):
        for y, b in enumerate(P):
            He[x, y] = H1[a, b] + H2[a, b] + (E0[a] if a == b else 0) + 0.5 * np.sum(H1[a, Q] * H1[Q, b] * (1 / (E0[a] - E0[Q]) + 1 / (E0[b] - E0[Q])))
    ev, V = np.linalg.eigh(0.5 * (He + He.T)); Lp = L2[np.ix_(P, P)]
    for k in range(len(P)):
        gv.append(dict(E=float(ev[k]), l2=float(V[:, k] @ Lp @ V[:, k]), dom=lab(P[int(np.argmax(V[:, k] ** 2))]), wdom=float((V[:, k] ** 2).max()), harm=float((E0[P[int(np.argmax(V[:, k] ** 2))]] - E0[g0]) * HC),
                       v3=float(V[:, k] ** 2 @ S[P, 3]), quanta=int(S[P].sum(1).max()), npoly=len(P)))
Ezp = [r['E'] for r in gv if r['dom'] == lab(g0)][0]
for r in gv: r['E'] = (r['E'] - Ezp) * HC
gv = sorted([r for r in gv if 1 < r['E'] <= EMAX + 200], key=lambda r: r['E'])
# ---------- VCI ----------
evv, Vv = np.linalg.eigh(np.diag(E0) + H1 + H2); vci = []
for k in range(1, NS):
    e = (evv[k] - evv[0]) * HC
    if e > EMAX + 200: break
    v = Vv[:, k]; d = int(np.argmax(v ** 2)); vci.append(dict(E=float(e), l2=float(v @ L2 @ v), dom=lab(d), wdom=float(v[d] ** 2), v3=float(v ** 2 @ S[:, 3])))
print('K = %d, %d basis states; GVPT2 ZPE %.2f, VCI ZPE %.2f cm-1; harmonic %s' % (K, NS, Ezp * HC, evv[0] * HC, np.round(w * HC, 2)))
print('GVPT2 levels: %d, VCI levels: %d (up to %.0f cm-1)' % (len(gv), len(vci), EMAX + 200))
for name, L in (('GVPT2', gv), ('VCI', vci)):
    print(name); print(' '.join('%.1f(l2=%.1f,%s)' % (r['E'], r['l2'], r['dom']) for r in L[:60]))
json.dump(dict(K=K, w=(w * HC).tolist(), zpe_gvpt2=float(Ezp * HC), zpe_vci=float(evv[0] * HC), gvpt2=gv, vci=vci), open('vptlev_%s_K%d.json' % (os.environ.get('OUT', TAG), K), 'w'), indent=1)
