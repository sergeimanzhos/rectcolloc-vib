"""DFT stencil data for numerical GVPT2 of LiCO2 (adapted from cl2co/vp_dft.py).  usage: python vp_dft_li.py basis tag [nproc]
Reference geometry, normal modes and dimensionless normal coordinates are those of coords.py (the collocation coordinates),
and the DFT settings are those of the collocation points (UB3LYP-D3(BJ), grid 4, grid-response gradients), so VPT2 and
collocation refer to the same surface.  Stencil (delta = 0.25): +-d, +-2d along each mode; +-d+-d for every pair.
E, gradient (in q) and dipole at every point -> vp_<tag>.npz (format of cl2co/vp_dft.py)."""
import sys, os, itertools, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import coords as c
BAS, TAG = sys.argv[1], sys.argv[2]; NP = int(sys.argv[3]) if len(sys.argv) > 3 else 12
d = float(os.environ.get('DELTA', 0.25)); CACHE = os.path.join(HERE, f'cache_{TAG}')
def init():
    global DM0
    from pyscf import lib; lib.num_threads(1); DM0 = np.load(os.path.join(CACHE, 'dm0.npy'))
def calc(q):
    from pyscf import gto, dft as D
    X = c.X_of_q(np.asarray(q)[None])[0]
    mol = gto.M(atom=[[e, x] for e, x in zip(c.EL, X)], basis=BAS, unit='Bohr', spin=1, symmetry=False, verbose=0)
    mf = D.UKS(mol).set(xc='b3lyp', conv_tol=1e-11, max_cycle=200); mf.grids.level = 4; mf.disp = 'd3bj'
    mf.kernel(dm0=DM0); assert mf.converged
    g = mf.nuc_grad_method(); g.grid_response = True
    return mf.e_tot, g.kernel().ravel() @ c.A, mf.dip_moment(unit='AU', verbose=0)
if __name__ == '__main__':
    pts = [np.zeros(6)]
    for i in range(6):
        for s in (-2, -1, 1, 2): v = np.zeros(6); v[i] = s * d; pts.append(v)
    for i, j in itertools.combinations(range(6), 2):
        for si, sj in itertools.product((-1, 1), repeat=2): v = np.zeros(6); v[i] = si * d; v[j] = sj * d; pts.append(v)
    pts = np.array(pts)
    from multiprocessing import Pool
    # other levels of theory: centre the stencil at that level's minimum (qmin from dft2.py ref), still in the aDZ normal coordinates;
    # the small off-diagonal quadratic terms this leaves are part of the fitted force field (H2 in vp2)
    qf = os.path.join(CACHE, 'qmin.npy'); Q0 = np.load(qf) if (TAG != 'adz' and os.path.exists(qf)) else np.zeros(6)
    with Pool(NP, init) as pool: R = pool.map(calc, list(pts + Q0[None]), chunksize=1)
    X0 = c.X_of_q(Q0[None])[0]; X0 = X0 - (c.m[:, None] * X0).sum(0) / c.m.sum(); print('stencil centre q0 =', np.round(Q0, 4))
    np.savez(os.path.join(HERE, f'vp_{TAG}.npz'), q=pts, E=np.array([r[0] for r in R]), g=np.array([r[1] for r in R]), mu=np.array([r[2] for r in R]),
             omega=c.omega, L=c.L, Dq=c.A, Xeq=X0, m=c.m, delta=d, basis=BAS)
    print('stencil done:', len(pts), 'points; |g| at centre (a.u.)', np.abs(R[0][1]).max())
