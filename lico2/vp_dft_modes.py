"""GVPT2 stencil for LiCO2 in the normal coordinates of the level of theory itself.  usage: python vp_dft_modes.py basis tag [nproc]
Step 1 (if cache_<tag>/modes.npz is missing): take the quadratic force field of the existing stencil vp_<tag>.npz (done in the aDZ
normal coordinates, centred at this level's minimum), diagonalise it and store this level's normal modes.
Step 2: stencil (delta = 0.25: +-d, +-2d per mode, +-d+-d per pair) in those modes -> vp_<tag>.npz (old file kept as vp_<tag>_adzmodes.npz)."""
import sys, os, itertools, shutil, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import coords as c
BAS, TAG = sys.argv[1], sys.argv[2]; NP = int(sys.argv[3]) if len(sys.argv) > 3 else 8
d = float(os.environ.get('DELTA', 0.25)); CACHE = os.path.join(HERE, f'cache_{TAG}'); MF = os.path.join(CACHE, 'modes.npz'); H = c.H2CM
mm = np.repeat(c.m, 3)
if not os.path.exists(MF):
    from vp_fit import fit
    cV, cM, info, dd = fit(os.path.join(HERE, f'vp_{TAG}.npz')); w0 = c.omega; F = np.zeros((6, 6))
    for i in range(6):
        for j in range(i, 6):
            a = [0] * 6; a[i] += 1; a[j] += 1; v = cV.get(tuple(a), 0.0)
            if i == j: F[i, i] = 2 * v * w0[i]
            else: F[i, j] = F[j, i] = v * np.sqrt(w0[i] * w0[j])
    w2, U = np.linalg.eigh(F); om = np.sqrt(w2); L = c.L @ U
    for k in range(6):
        if L[np.argmax(np.abs(L[:, k])), k] < 0: L[:, k] *= -1; U[:, k] *= -1
    A = (L / np.sqrt(mm)[:, None]) / np.sqrt(om)[None, :]
    Q0 = np.load(os.path.join(CACHE, 'qmin.npy')); Xmin = c.X_of_q(Q0[None])[0]
    np.savez(MF, omega=om, L=L, A=A, Xmin=Xmin, U=U)
    print('harmonic frequencies in the aDZ modes (diagonal only):', np.round(np.sqrt(np.diag(F)) * H, 2))
    print('harmonic frequencies of this level                  :', np.round(om * H, 2))
    print('mode mixing |U|:'); print(np.round(np.abs(U), 3))
M = np.load(MF); om, L, A, Xmin = M['omega'], M['L'], M['A'], M['Xmin']
def init():
    global DM0
    from pyscf import lib; lib.num_threads(1); DM0 = np.load(os.path.join(CACHE, 'dm0.npy'))
def calc(q):
    from pyscf import gto, dft as D
    X = Xmin + (A @ np.asarray(q)).reshape(-1, 3)
    mol = gto.M(atom=[[e, x] for e, x in zip(c.EL, X)], basis=BAS, unit='Bohr', spin=1, symmetry=False, verbose=0)
    mf = D.UKS(mol).set(xc='b3lyp', conv_tol=1e-11, max_cycle=200); mf.grids.level = 4; mf.disp = 'd3bj'
    mf.kernel(dm0=DM0); assert mf.converged
    g = mf.nuc_grad_method(); g.grid_response = True
    return mf.e_tot, g.kernel().ravel() @ A, mf.dip_moment(unit='AU', verbose=0)
if __name__ == '__main__':
    pts = [np.zeros(6)]
    for i in range(6):
        for s in (-2, -1, 1, 2): v = np.zeros(6); v[i] = s * d; pts.append(v)
    for i, j in itertools.combinations(range(6), 2):
        for si, sj in itertools.product((-1, 1), repeat=2): v = np.zeros(6); v[i] = si * d; v[j] = sj * d; pts.append(v)
    pts = np.array(pts)
    from multiprocessing import Pool
    with Pool(NP, init) as pool: R = pool.map(calc, list(pts), chunksize=1)
    out = os.path.join(HERE, f'vp_{TAG}.npz')
    if os.path.exists(out) and not os.path.exists(out.replace('.npz', '_adzmodes.npz')): shutil.copy(out, out.replace('.npz', '_adzmodes.npz'))
    X0 = Xmin - (c.m[:, None] * Xmin).sum(0) / c.m.sum()
    np.savez(out, q=pts, E=np.array([r[0] for r in R]), g=np.array([r[1] for r in R]), mu=np.array([r[2] for r in R]),
             omega=om, L=L, Dq=A, Xeq=X0, m=c.m, delta=d, basis=BAS)
    print('stencil done:', len(pts), 'points; |g| at centre (a.u.)', np.abs(R[0][1]).max())
