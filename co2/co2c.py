"""Free CO2: coordinates, basis and exact Cartesian finite-difference KEO for rectangular collocation (skill sections 2-4).
Basis coordinates (rotation- and translation-invariant, smooth at linearity):
  x1 = (r1 + r2 - 2 re) / (sqrt2 s1)      symmetric stretch, dimensionless
  x3 = (r1 - r2) / (sqrt2 s3)             antisymmetric stretch, dimensionless
  xb = (1 + cos theta) * kappa            bend; xb = q_bend^2 / 2 near linearity (q_bend = dimensionless 2D bend radius)
Harmonic limit: V = w1 x1^2/2 + w3 x3^2/2 + w2 xb, ground state exp(-x1^2/2 - x3^2/2 - xb).
J = 0 (l = 0) states: basis phi_j(x).  J = 1 f-parity states (lowest = bend fundamental, l = 1): basis (a x b)_Z phi_j(x),
where a, b are the two C->O vectors; the Cartesian Laplacian is applied to these functions directly."""
import os, numpy as np, numexpr as ne
HERE = os.path.dirname(os.path.abspath(__file__))
H2CM = 219474.6313705; AMU = 1822.888486209; BOHR = 0.529177210903
mC, mO = 12.0 * AMU, 15.994914620 * AMU; M3 = np.array([mC, mO, mO]); MM = np.repeat(M3, 3)
_c = np.load(os.path.join(HERE, 'ckpt_co2.npz')); _X = _c['X']; _w = np.sort(_c['omega'])
RE = float(np.mean([np.linalg.norm(_X[1] - _X[0]), np.linalg.norm(_X[2] - _X[0])]))
W2, W1, W3 = float(_w[:2].mean()), float(_w[2]), float(_w[3])          # bend (x2), symmetric stretch, antisymmetric stretch (Eh)
S1 = 1 / np.sqrt(W1 * mO); G33 = 1 / mO + 2 / mC; S3 = np.sqrt(G33 / W3); GTH = 2 * (1 / mO + 2 / mC) / RE ** 2; KAP = W2 / GTH
BE = 1 / (2 * 2 * mO * RE ** 2)                                         # equilibrium rotational constant (Eh)
FDW = np.array([-1, 16, -30, 16, -1]) / 12.0; FDS = np.array([-2, -1, 0, 1, 2])
def geom(x):
    """x (n,3) -> Cartesian geometries (n,3,3) in the XY plane: C at the origin, O1 on +x, O2 at angle theta."""
    x = np.atleast_2d(x); r1 = RE + (S1 * x[:, 0] + S3 * x[:, 1]) / np.sqrt(2); r2 = RE + (S1 * x[:, 0] - S3 * x[:, 1]) / np.sqrt(2)
    ct = x[:, 2] / KAP - 1; st = np.sqrt(np.clip(1 - ct ** 2, 0, None)); X = np.zeros((len(x), 3, 3)); X[:, 1, 0] = r1; X[:, 2, 0] = r2 * ct; X[:, 2, 1] = r2 * st; return X
def internal(X):
    a = X[:, 1] - X[:, 0]; b = X[:, 2] - X[:, 0]; r1 = np.linalg.norm(a, axis=1); r2 = np.linalg.norm(b, axis=1); ct = (a * b).sum(1) / (r1 * r2)
    x = np.stack([(r1 + r2 - 2 * RE) / (np.sqrt(2) * S1), (r1 - r2) / (np.sqrt(2) * S3), (1 + ct) * KAP], 1); return x, a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
def V_harm(x): x = np.atleast_2d(x); return 0.5 * W1 * x[:, 0] ** 2 + 0.5 * W3 * x[:, 1] ** 2 + W2 * x[:, 2]
class Basis:
    def __init__(s, C, beta, J=0): s.beta = np.asarray(beta, float); s.Cs = np.asarray(C) / s.beta; s.c2 = (s.Cs ** 2).sum(1); s.J = J
    def __call__(s, X):
        x, axb = internal(X); xs = x / s.beta; G = xs @ s.Cs.T; x2 = (xs ** 2).sum(1)[:, None]; c2 = s.c2[None, :]; B = ne.evaluate('exp(-0.5*(x2+c2-2*G))')
        return B if s.J == 0 else axb[:, None] * B
def T_basis(bas, X, h=1e-3):
    n = len(X); Xf = X.reshape(n, 9); out = FDW[2] * (1 / MM).sum() * bas(X)
    for i in range(9):
        for sft, w in zip(FDS, FDW):
            if sft == 0: continue
            Xd = Xf.copy(); Xd[:, i] += sft * h; out += (w / MM[i]) * bas(Xd.reshape(n, 3, 3))
    return -0.5 * out / h ** 2
def r2seq(idx, d=3):
    phi = 2.0
    for _ in range(60): phi = (1 + phi) ** (1.0 / (d + 1))
    return (0.5 + np.outer(np.asarray(idx, float), 1 / phi ** np.arange(1, d + 1))) % 1.0
def box(vbox):
    """harmonic turning points at energy vbox, widened on the soft (stretched) side of the symmetric stretch"""
    a1 = np.sqrt(2 * vbox / W1); a3 = np.sqrt(2 * vbox / W3); return np.array([-a1, -1.25 * a3, 0.0]), np.array([1.7 * a1, 1.25 * a3, 1.35 * vbox / W2])
