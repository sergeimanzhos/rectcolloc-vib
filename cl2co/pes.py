"""Python port of Cl2COpes.m (NN ensemble PES for Cl2CO, symmetrised over Cl1<->Cl2).
Input Rij in bohr, columns [CO, CCl1, CCl2, OCl1, OCl2, Cl1Cl2]; output cm-1."""
import numpy as np
BOHR_PER_A = 1.8897261254535
REQ = np.array([1.1816, 1.7518, 1.7518, 2.6047, 2.6047, 2.9021]) * BOHR_PER_A
class NNPES:
    def __init__(self, fname):
        rows = [list(map(float, l.split())) for l in open(fname) if l.strip()]
        h = rows[0]; self.nd, self.Nterms, self.nNNs = int(h[0]), int(h[1]), h[2]; self.info = h[3]
        nd = self.nd
        self.minp = np.array(rows[1][:nd]); self.maxp = np.array(rows[2][:nd])
        self.mint, self.maxt = rows[1][nd], rows[1][nd + 1]
        W = np.array(rows[3:3 + self.Nterms]); assert W.shape == (self.Nterms, nd + 3)
        self.w, self.c, self.b, self.d = W[:, :nd], W[:, nd], W[:, nd + 1], W[:, nd + 2]
    def _import(self, x):
        xn = 2 * (x - self.minp) / (self.maxp - self.minp) - 1
        E = np.tanh(xn @ self.w.T + self.b) @ self.c + self.d.sum()
        E /= self.nNNs
        return self.mint + 0.5 * (E + 1) * (self.maxt - self.mint)
    def __call__(self, Rij):
        R = np.atleast_2d(Rij); x = np.exp(-5.0 * (R - REQ) / REQ); out = np.empty(len(x))
        for i in range(0, len(x), 20000):
            s = slice(i, i + 20000)
            out[s] = 0.5 * (self._import(x[s]) + self._import(x[s][:, [0, 2, 1, 4, 3, 5]]))
        return out
def rij(X):
    """X (..., 4, 3) bohr, atom order C, O, Cl1, Cl2 -> Rij bohr"""
    X = np.asarray(X); d = lambda a, b: np.linalg.norm(X[..., a, :] - X[..., b, :], axis=-1)
    return np.stack([d(0, 1), d(0, 2), d(0, 3), d(1, 2), d(1, 3), d(2, 3)], -1)
