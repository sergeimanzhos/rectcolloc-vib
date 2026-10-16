"""Eckart-frame dimensionless normal coordinates for LiCO2 (atom order C, O, O, Li), atomic units.
Reference: UB3LYP-D3(BJ)/aug-cc-pVDZ (grid 8) C2v minimum and harmonic modes from ff_lico2.npz.
Also provides the quartic force field (same as GVPT2/VCI) as a cheap surrogate V_ff(q)."""
import numpy as np, math, os
H2CM = 219474.6313705; AMU = 1822.888486209
_d = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'vpt', 'ff_lico2.npz'))
X0 = _d['X']; m = _d['m'] * AMU; NA = len(m); mm = np.repeat(m, 3)
omega = _d['omega']; L = _d['L'].copy(); NQ = len(omega); EL = [str(e) for e in _d['el']]
# project translations/rotations out of L exactly (stored modes satisfy Eckart only to ~1e-5)
_TR = []
for _a in range(3):
    _t = np.zeros((NA, 3)); _t[:, _a] = 1; _TR.append((_t * np.sqrt(mm.reshape(NA, 3))).ravel())
    _e = np.zeros(3); _e[_a] = 1; _TR.append((np.cross(_e, X0) * np.sqrt(mm.reshape(NA, 3))).ravel())
_U = np.linalg.qr(np.array(_TR).T)[0]; L = L - _U @ (_U.T @ L)
_w, _v = np.linalg.eigh(L.T @ L); L = L @ _v @ np.diag(_w ** -0.5) @ _v.T
for k in range(NQ):                                    # fix sign convention
    if L[np.argmax(np.abs(L[:, k])), k] < 0: L[:, k] *= -1
Q = np.sqrt(omega)[:, None] * (L.T * np.sqrt(mm)[None, :])   # q = Q (x - X0)
A = (L / np.sqrt(mm)[:, None]) / np.sqrt(omega)[None, :]     # x = X0 + A q
def X_of_q(q):
    q = np.atleast_2d(q); return X0[None] + (q @ A.T).reshape(-1, NA, 3)
def eckart_rotate(X):
    com = (m[None, :, None] * X).sum(1) / m.sum(); x = X - com[:, None, :]
    Hc = np.einsum('i,nia,ib->nab', m, x, X0); U, S, Vt = np.linalg.svd(Hc)
    d = np.sign(np.linalg.det(np.einsum('nba,ncb->nac', Vt, U)))
    D = np.ones((len(X), 3)); D[:, 2] = d
    R = np.einsum('nba,nb,ncb->nac', Vt, D, U)
    return np.einsum('nab,nib->nia', R, x), S
def q_of_X(X, eckart=True):
    X = np.asarray(X).reshape(-1, NA, 3); x = eckart_rotate(X)[0] if eckart else X
    return (x - X0[None]).reshape(len(X), 3 * NA) @ Q.T
# ---- quartic force field surrogate (dimensionless q), identical to the one used for GVPT2/VCI
def _ff():
    Hq0, Hqp, Hqm, dQ = _d['Hq0'], _d['Hqp'], _d['Hqm'], float(_d['dQ']); acc = {}
    # NB: _d L may differ in sign from L above; rebuild constants in the sign-fixed basis
    sg = np.sign(np.einsum('ik,ik->k', _d['L'], L))
    for k in range(NQ):
        D1 = (Hqp[k] - Hqm[k]) / (2 * dQ) * sg[k]; D2 = (Hqp[k] + Hqm[k] - 2 * Hq0) / dQ ** 2
        for i in range(NQ):
            for j in range(NQ):
                s = sg[i] * sg[j]
                acc.setdefault(tuple(sorted((i, j, k))), []).append(D1[i, j] * s)
                acc.setdefault(tuple(sorted((i, j, k, k))), []).append(D2[i, j] * s)
    terms = []
    for key, v in acc.items():
        phi = np.mean(v) / np.sqrt(np.prod(omega[list(key)]))
        c = 1.0
        for a in set(key): c *= math.factorial(key.count(a))
        if abs(phi) > 1e-12: terms.append((np.array([key.count(a) for a in range(NQ)]), phi / c))
    return terms
TERMS = _ff()
def V_ff(q):
    q = np.atleast_2d(q); v = 0.5 * (q ** 2) @ omega
    for p, c in TERMS: v = v + c * np.prod(q ** p[None, :], axis=1)
    return v
def V_harm(q):
    q = np.atleast_2d(q); return 0.5 * (q ** 2) @ omega
