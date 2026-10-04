"""Polynomial force field V(q) and dipole mu(q) in dimensionless normal coords from the vp_dft stencil.
V: energies + gradients (gradient rows scaled by delta); monomials limited to those the stencil identifies.
mu: values only."""
import numpy as np, itertools
H2CM = 219474.6313705
QUAD = (5, 2, 4, 1)          # nu1, nu2, nu5, nu6 (label order nu3,nu6,nu2,nu4,nu5,nu1)

def candidates(kind):
    out = []
    for a in itertools.product(range(7), repeat=6):
        deg = sum(a); sup = [i for i in range(6) if a[i]]; ns = len(sup)
        if deg == 0: continue
        if kind == 'V':
            ok = (ns == 1 and deg <= 6) or (ns == 2 and deg <= 4) or (ns == 3 and deg <= 4 and max(a) <= 2) \
                 or (ns == 4 and deg == 4 and set(sup) == set(QUAD))
        else:
            ok = (ns == 1 and deg <= 4) or (ns == 2 and (deg <= 3 or a.count(2) == 2)) \
                 or (ns == 3 and deg == 3 and set(sup) <= set(QUAD)) or (ns == 4 and deg == 4 and set(sup) == set(QUAD))
        if ok: out.append(a)
    return out

def mono(q, a): return np.prod(q ** np.array(a), 1)
def dmono(q, a, k):
    if a[k] == 0: return np.zeros(len(q))
    b = list(a); b[k] -= 1; return a[k] * np.prod(q ** np.array(b), 1)

def fit(npz):
    d = np.load(npz); q = d['q']; dl = float(d['delta']); V = d['E'] - d['E'][0]; g = d['g']; mu = d['mu'] - d['mu'][0]
    MV = candidates('V')
    A = np.vstack([np.array([mono(q, a) for a in MV]).T] + [dl * np.array([dmono(q, a, k) for a in MV]).T for k in range(6)])
    y = np.concatenate([V] + [dl * g[:, k] for k in range(6)])
    s = np.linalg.svd(A / np.linalg.norm(A, axis=0), compute_uv=False)
    cV, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = A @ cV - y
    MM = candidates('mu'); B = np.array([mono(q, a) for a in MM]).T
    sm = np.linalg.svd(B / np.linalg.norm(B, axis=0), compute_uv=False)
    cM, *_ = np.linalg.lstsq(B, mu, rcond=None)
    info = dict(nV=len(MV), rankV_cond=s[0] / s[-1], rmsV_cm=np.sqrt(np.mean(resid[:len(q)] ** 2)) * H2CM,
                nM=len(MM), rankM_cond=sm[0] / sm[-1], rmsM=np.sqrt(np.mean((B @ cM - mu) ** 2)), npts=len(q))
    cVd = {a: v for a, v in zip(MV, cV)}; cMd = {a: v for a, v in zip(MM, cM)}
    return cVd, cMd, info, d

if __name__ == '__main__':
    import sys
    cV, cM, info, d = fit(sys.argv[1])
    print(info)
    om = d['omega'] * H2CM
    print('omega Hessian:', np.round(om, 2))
    print('2*c_ii      :', np.round([2 * cV[tuple(2 * (np.arange(6) == k))] * H2CM for k in range(6)], 2))
    print('linear terms (cm-1):', np.round([cV[tuple(1 * (np.arange(6) == k))] * H2CM for k in range(6)], 4))
