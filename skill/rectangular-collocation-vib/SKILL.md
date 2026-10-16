---
name: "rectangular-collocation-vib"
description: "Use to compute molecular vibrational energy levels (and optionally IR intensities) by rectangular collocation with an exact Cartesian finite-difference KEO, accumulated square matrices, optional potential-gradient rows, and only the requested number of eigenpairs."
---

# Rectangular collocation for vibrational levels (exact Cartesian FD KEO, accumulated matrices)

Computes the lowest `nev` J=0 vibrational levels of a molecule. Method: rectangular collocation (M points > N basis functions), exact kinetic energy operator applied by finite differences in Cartesian coordinates, and accumulation of the square N x N matrices S = sum B^T B and H = sum B^T (T+V) B over batches of points. The rectangular B matrix is never stored. Options: potential-gradient (derivative) collocation rows (section 5b) and IR intensities from a dipole surface (section 8b).

## 0. Before starting

**Check the machine.** Read the core count and RAM (`nproc`/`free` on Linux; `$env:NUMBER_OF_PROCESSORS` and `Get-CimInstance Win32_ComputerSystem` on Windows).
- Many cores (for example the user's 48-core / 223 GB Windows PC): run full-size jobs. Parallelize as in section 7, run long jobs in the background with a log file, and checkpoint.
- Small cloud workspace (2-4 cores): say so and use reduced test sizes unless the user insists on full size.

**Reuse existing code.** If the project folder already has tested scripts (e.g. `accum.py`, `collo.py`, a `run_*.py`, a ported PES), read them first and reuse or extend them. Only write new code where it is missing. Read the project's CLAUDE.md for molecule-specific facts. Worked examples with tested code are in the repository https://github.com/sergeimanzhos/rectcolloc-vib (`cl2co/`: 6D, NN PES and DFT points; `lico2/`: 6D doublet radical, DFT points, 100 states, IR intensities, GVPT2).

**Inputs.** Ask only for missing essentials; use the defaults for everything else and state them once at the start.

| Input | Default | Notes |
|---|---|---|
| Molecule: atoms, isotopes | none, required | Masses from isotope table, in electron masses (AMU2ME = 1822.888486209) |
| Energy source | none, required | (a) provided PES subroutine, or (b) ab initio method/basis/program |
| `nev`, number of eigenpairs | none, required | Only these are computed |
| Initial basis size N | none, required | Points are added in batches of N |
| Initial Vmax (cm^-1) | none, required | Point region V <= Vmax; optimized only if the user asks |
| Coordinates for the basis | Eckart-frame dimensionless normal coordinates | Alternatives: bond/valence internals, curvilinear normal coordinates, user-defined |
| Basis type | Gaussian | Alternatives: inverse multiquadric (IMQ), harmonic-oscillator products, user-defined |
| FD stencil size | 5 | 3, 7, 9, ... on request |
| FD step h | 1e-3 bohr for 5-point | Adjust with stencil (section 3) |
| Max points M / stopping rule | until levels 1..nev change < 0.1 cm^-1 between batches, or user-given M_max | |
| Potential-gradient rows | off | On request, or offer it when energies are ab initio and gradients are cheap (section 5b) |
| IR intensities | off | On request (section 8b). With ab initio energies, always save dipoles at the points anyway |
| Reference levels | none | If given (often ZPE followed by transitions), report errors against them |

Units: atomic units internally, cm^-1 in output (HARTREE2CM = 219474.6313705).

## 1. Energy source

**Provided PES subroutine** (.m, .f/.f90, .py, ...):
- Port it to vectorized NumPy, or wrap Fortran with f2py/ctypes. Keep the original's coordinate convention and atom order; write the convention down in the code.
- Verify the port before using it. Compare against values the user supplies, or against the original run in Octave/MATLAB if available, to about 1e-10 relative. Also check that the minimum and harmonic frequencies are sensible.
- Check for holes: sample the point region densely and flag any V far below the minimum or any NaN. Polynomial force fields often turn over at large displacements.
- For gradient rows, get dV/dx analytically or by autodiff, and check against FD of energies.

**Ab initio:**
- Energies are needed only at the collocation points, not at the FD-displaced geometries: V is multiplicative, and only the basis is differentiated. That makes ab initio energies the dominant cost.
- Cache every result keyed by its R2-sequence index (`cache/pts_<tag>.npz`): energy, Cartesian gradient (if used), dipole vector, and the SCF convergence flag. Never recompute a cached point; changing N, M, Ec, widths or the basis reuses the same data. Save after every ~100 candidate points so that progress is visible and the run can restart.
- Always store the dipole from each converged density (it is essentially free), so that intensities can be computed later without rerunning anything.
- Run many single-point calculations in parallel (processes x threads about equal to the core count; set OMP_NUM_THREADS / lib.num_threads per process to avoid oversubscription). Start the SCF from the equilibrium density, and log any unconverged points.
- DFT settings (PySCF): grid level 3 with conv_tol 1e-9 costs about half as much as level 4 / 1e-11 and changes energies by < 0.1 cm^-1; confirm this once on a few points. Compute the reference minimum energy E0 with exactly the same settings. Use gradients with grid response (`grad.grid_response = True`) so they are consistent with the energies.
- Geometry frame: pass geometries in the Eckart frame and stop the program from reorienting them (PySCF: symmetry=False; Gaussian: `nosymm`), so gradients and dipoles are body-fixed.
- Pre-screening: compute only the energy first, and the gradient only if V <= Vmax. After a few hundred candidates, fit a low-order polynomial surrogate of V(q) to all computed energies and skip candidates with surrogate V > ~1.35 Vmax. Check on the unscreened early points that the screen would not have rejected any in-region point.
- PySCF does not run on native Windows: use WSL2 there, or the user's quantum chemistry program (e.g. Gaussian) through generated inputs and parsed outputs.

## 2. Reference geometry, coordinates, point region

1. **Minimum and harmonic analysis.** Optimize the geometry and compute the Hessian (analytic, or FD of energies or gradients). Project translations and rotations out of the mass-weighted Hessian before diagonalizing, so the vibrational modes satisfy the Eckart conditions exactly. Get the mass-weighted normal modes, harmonic frequencies, the harmonic ZPE, the map q = diag(sqrt(omega)) L^T sqrt(m) d to dimensionless normal coordinates, and its inverse d = A q.
2. **Basis coordinates must be rotation- and translation-invariant functions of the Cartesian geometry.** This is required for the Cartesian Laplacian to give the exact J=0 KEO. For normal coordinates: remove the centre of mass, rotate to the Eckart frame relative to the equilibrium geometry (a mass-weighted Kabsch/SVD rotation with the determinant sign fix satisfies the Eckart conditions exactly), then q = Q (R x - Xeq). Internal coordinates are invariant already and need no Eckart step. For curvilinear normal coordinates: do a GF analysis in internals, and treat torsions as periodic.
3. **Region.** Find the 1D turning points along each basis coordinate and take the bounding box. Use a box energy a little above Vmax, so that Vmax can later be lowered without regenerating points. Generate points with the R2 (Kronecker / plastic-number) quasi-random sequence in the box and keep those with V <= Vmax. Prefer R2 over Sobol: it has a better minimum distance and gives a lower cond(B^T B).
4. **Coordinate singularities - check before production.**
   - If Vmax exceeds a barrier to linearity (e.g. the H2O bend) or reaches another configuration where the Eckart frame is ill-defined, two things go wrong. The Eckart rotation becomes nearly singular: the second singular value of sum m x Xeq^T drops towards 0. And rectilinear normal coordinates become double-valued past linearity.
   - The FD KEO is then garbage at those geometries, and a handful of such points poisons the whole fit with spurious states and huge residuals.
   - Detect it by testing the h-sensitivity of T*phi (h vs 2h) on the validation points, and inspect the worst geometries.
   - Remedies: exclude points and Gaussian centres beyond a geometric cut (e.g. bond angle > 165 deg), switch to valence coordinates, or lower Vmax.
   - Flag states that lie above such a barrier as coordinate-limited.
5. **Linear molecules (e.g. CO2).** J = 0 contains only l = 0 states. Get the l > 0 vibrational states from J = l blocks:
   - Use rotation-invariant internal coordinates (symmetric/antisymmetric stretch and a bend variable smooth at linearity, e.g. xb = (1 + cos theta) kappa, which is q_bend^2/2 near linearity).
   - Block l: basis F_l(X) phi_j(internal), with F_l = Re (l even) or Im (l odd) of (p_Z + i n_Z)^l, where n = a^ x b^ and p = (a^ + b^) sin(theta/2) (a^, b^ unit bond vectors). This is the ZZ..Z component of the rank-l harmonic tensor with |K| = l. For l = 1, (a x b)_Z alone is enough. Apply the Cartesian FD Laplacian to the full product.
   - Band origin = E(J = l) - E0 - B_e l (rotational energy B [J(J+1) - l^2]). Coriolis coupling to other |K| is of order B^2/dE (about 1e-3 cm^-1) and is neglected.
   - The bend radial function must be the l = 0 projection of a 2D Gaussian, exp(-(q - c)^2 / (2 beta^2)) I0(q c / beta^2) with q = sqrt(2 xb) (use scipy.special.ive for stability). Plain Gaussians in xb do not converge (residuals of 50-300 cm^-1).
   - Comparison tables then list l for every level; (v1 v2 v3) labels repeat for the different l components.
6. **Polynomial force fields are not valid references far from the minimum.** Check a quartic field for holes before using it as a benchmark (LiCO2: a hole 1316 cm^-1 above the minimum along the soft Li mode). For a degenerate mode, check cylindrical symmetry (phi_aaaa = 3 phi_aabb); a numerically differentiated CO2 field violated it by 30%. A least-squares fit of a symmetry-adapted polynomial (degree 8) to the computed energies below ~4000 cm^-1 gave stable quartic constants.
7. **Validation set.** Reserve a fixed, separate stretch of the R2 sequence (never used in fitting) for residual checks, typically 1000 points.

## 3. Exact KEO by Cartesian finite differences

T = -1/2 sum_i (1/m_i) d^2/dx_i^2 over all 3*N_atoms Cartesian coordinates. For each collocation point and each Cartesian coordinate, displace that coordinate by the stencil offsets k*h, convert every displaced geometry to the basis coordinates (including the Eckart step if used), evaluate all basis functions, and combine with the central second-derivative weights.

- Generate the weights with Fornberg's algorithm for any stencil size s (odd). For reference: 3-point [1,-2,1]/h^2; 5-point [-1,16,-30,16,-1]/(12h^2); 7-point [2,-27,270,-490,270,-27,2]/(180h^2).
- Step size: about 1e-4 bohr for s=3, 1e-3 for s=5, about 5e-3 for s>=7. If the user changes s, rescan h quickly and pick the plateau.
- Basis evaluations per point: 3*N_atoms*(s-1) + 1 (the centre is shared).
- **Mandatory checks before production.** Stop and fix the code if any of them fails.
  - In rectilinear normal coordinates with a Gaussian basis (no Eckart rotation), compare the FD T*B against the analytic -1/2 sum_k omega_k d^2/dq_k^2 at about 40 points. Expect a relative error of about 1e-8 or better for s=5.
  - With the Eckart step, check |q_eckart(geom(q)) - q| at about 1e-12, and that q is invariant under rotation and translation.
  - Exact-KEO check: apply the FD KEO (with the Eckart step) to a function of one bond length f(r_AB) and compare with -1/(2 mu_AB) (f'' + 2 f'/r). Expect about 1e-8.

## 4. Basis functions

- **Gaussian (default).** phi_j(x) = exp(-sum_k (x_k - c_jk)^2 / (2 beta_k^2)) with anisotropic per-coordinate widths beta_k.
  - Centres: N points drawn from the R2 sequence with V <= Ec, always with Ec < Vmax. With ab initio energies, use the polynomial surrogate of V to select the centres only, and apply the same geometric cut as for the points.
  - Widths: start with beta_k about 1 in dimensionless normal coordinates. For other coordinates, scale each by its harmonic ground-state width so that beta = 1 means the same thing. Keep symmetry-equivalent coordinates equal.
- **Choosing Ec and widths - this depends on the dimension, on N, and on how high the levels go.**
  - In higher dimension with many levels close to the ZPE, a low Ec worked best: about 0.5-0.6 x ZPE, or about 0.3 x E_top. Examples: H2CO and Cl2CO, 6D, N = 5000-10000, top level about 2 x ZPE.
  - In low dimension, or when the levels reach far above the ZPE, a low Ec packs the centres far too densely (spacing << beta). S then becomes hopelessly over-complete, with spurious states and erratic levels. Example: H2O, 3D, N = 1000, 30 levels reaching about 3.4 x ZPE above the minimum; there Ec of about 0.6-0.8 x E_top (12000 cm^-1 with Vmax = 18000) and beta about 0.85-0.9 worked best.
  - Scale reasoning: keep the ratio of width to centre spacing similar to what worked before. Watch how many vectors canonical orthogonalization keeps: many dropped means the centres are too dense or too wide.
  - **The kept count is the real basis size.** CO2 (3D, 34 levels to 4350 cm^-1): with beta = 1.0/1.0/0.8 only about 470 of N = 2400 Gaussians survived, and levels near 4000 cm^-1 stayed 5-7 cm^-1 off however large N was. Narrower widths (0.8/0.8/0.65) with centres up to Ec = Vmax = 8000 cm^-1 kept about 800 and converged all 34 levels to 0.2 cm^-1. When the top levels converge slowly in N, check the kept count before adding functions.
  - **Tune against an exact reference when one exists.** Fit a quartic field to the computed energies, get exact levels from VCI on it, and tune N, M, Ec and widths by collocation on that same field (all levels, not only the fundamentals). Then apply the settings to the ab initio points. Selecting centres or points with an analytic surrogate is fine, but say which surface was used.
  - **Harmonic sanity test (cheap; use it to pick N, Ec and beta).** Replace V by the harmonic potential at the same points and use the rectilinear KEO (eckart=False), so the exact levels are known. A good basis gives median errors of about 0.01 cm^-1 and validation residuals of about 0 for the lowest levels. If the harmonic test is poor, the real run will be worse.
  - Then tune on the real surface. Scan Ec and an isotropic beta first, then refine beta_k per coordinate in 0.95/1.05 steps (coordinate descent), judging by validation residual and level stability vs M (section 6).
- **Inverse multiquadric.** phi_j(x) = (1 + sum_k ((x_k - c_jk)/beta_k)^2)^(-1/2). Same centre and width logic as the Gaussians; the slower decay usually needs larger M/N.
- **Harmonic-oscillator products.** Products of 1D HO functions in the chosen coordinates, pruned by harmonic excitation energy <= E_cut, with E_cut chosen to give about N functions (report the actual N). These are global and oscillating: expect M of about 5-6N before the levels stabilize, against about 2N for Gaussians. Spurious eigenvalues (even negative ZPEs) appear below about 3N.
- **User-defined basis.** Implement it as a function `basis(x) -> (npts, N)`, plus `grad(x) -> (npts, N, ncoord)` if gradient rows are used, so that sections 3 and 5 are unchanged.

## 5. Accumulation and eigensolver

For each batch of points (and its FD-displaced geometries):
```
B  = basis(points)                       # (N_b, N)
TB = FD-KEO applied to basis              # (N_b, N)
HB = TB + V[:,None] * B
S += B.T @ B ;  H += B.T @ HB             # H is generally non-symmetric
```
- Checkpoint S, H, M, the point indices and all inputs after every batch (`cache/SH_<tag>_N<N>_M<M>.npz`). Adding points resumes from the last checkpoint, and nothing rectangular is stored.
- **Canonical orthogonalization.** Eigendecompose S, keep eigenvalues > thr*lambda_max (default thr = 1e-12; check 1e-10 and 1e-14 once - in CO2, 1e-14 admitted spurious states with errors > 100 cm^-1, while 1e-10 and 1e-12 agreed to 0.05 cm^-1), and form X = U_k Lambda_k^(-1/2). Report the number kept and cond(S). Raising thr does not cure an over-complete basis and can even create spurious states; fix Ec, N or beta instead.
- **Only nev eigenpairs.** A = X^T H X is non-symmetric. Use shift-invert Arnoldi (`scipy.sparse.linalg.eigs(A, k=nev, sigma=s0)`) with s0 just below the expected ZPE (for example 0.9 x harmonic ZPE), so the eigenvalues nearest s0 are the lowest. Sort by real part and report max|Im|. Never run a full dense eig of A unless N is small or the user asks.
- For large N, orthogonalizing S dominates the cost. Report the timing per stage so the user can see where time goes.

## 5b. Option: potential-gradient (derivative) collocation rows

Use this when the energy source gives gradients cheaply: an ab initio gradient costs about 1.5-2 energies, and a gradient from an analytic or autodiff PES is nearly free. Each point then gives 1 + n_q equations instead of 1.

**Equations.** For a direction D = sum_i a_i d/dx_i whose vector a is constant in space, [T, D] = 0 exactly, so D R = (T + V - E) D psi + (D V) psi = 0.
- Choose a_k = dX/dq_k, the columns of A. For Eckart normal coordinates generated as X(q) = Xeq + A q, q(X(q) + eps a_k) = q + eps e_k exactly, so D_k is simply d/dq_k at the point.
- Rows per point and coordinate k:
  - B_k = d phi_j / d q_k (analytic).
  - A_k = T(d phi_j/dq_k) + V d phi_j/dq_k + V_k phi_j, with V_k = g_cart . a_k.
- Computing T(d phi_j/dq_k): use d/dq_k (T phi_j), from a 5-point first-derivative FD in q_k (eps = 1e-2) of the exact FD-KEO rows evaluated at q +/- eps e_k and q +/- 2 eps e_k. This keeps the KEO exact without differentiating the Eckart map. Check it against the analytic rectilinear result (expect about 1e-8).
- For internal or other coordinates, take a_k as a right inverse of dq/dX (the mass-weighted pseudo-inverse) at each point. Keep a_k constant (a straight line X + eps a_k), so D_k = a_k . grad still commutes exactly with T; take D_k phi and D_k (T phi) by 5-point FD along that line (eps = 1e-2). Check a_k . grad u_j = delta_jk (about 1e-9) and the commutation T D phi = D T phi (about 1e-9). With a bend radius q = sqrt(2 xb), use grad q = grad xb / q.
- Weight the gradient rows by lambda(V) = (1 - V/Vmax)^p. Best p differs: p = 3 for Cl2CO (NN PES), p = 0 (constant) for LiCO2; CO2 used p = 3. Scan p in {0, 1, 3} once.

**Accumulation.** S += B^T B + lambda sum_k B_k^T B_k ; H += B^T A + lambda sum_k B_k^T A_k.
- Use lambda about 1 in dimensionless coordinates (it scales like beta^2). Results are insensitive between 0.25 and 4, so scan it once.
- Nothing else changes: orthogonalization, eigs and diagnostics are as in sections 5 and 6.

**Cost.** Basis-side work grows by about 4 n_q x; for a triatomic that is 13x more KEO evaluations per point, which is negligible against ab initio energies. Gradient rows also improve the conditioning of S.

**What to expect and report.**
- In H2O at B3LYP/aug-cc-pVDZ with N = 1000, 500 points with gradients matched 1500-2000 energy-only points. In CO2 at B3LYP-D3(BJ)/aug-cc-pVTZ (34 levels, N = 4000), 2400 points with gradients matched 8157 energy-only points to 0.25 cm^-1. At equal DFT cost that is about a 2x saving at low M/N, and both variants converge to the same limit (median 0.04 cm^-1).
- Report the convergence of both variants against their own all-points limits, the difference between the two limits, and the DFT cost ratio (gradient time / energy time) measured on the actual machine.

## 6. Convergence diagnostics (after every batch)

- Levels: ZPE and transitions E_i - E_0 for i = 1..nev-1; change since the previous batch (median and max, for levels 1-20 and for all levels).
- Validation residual: for each state, |(H - E_i) psi_i| at the validation points divided by |psi_i|, in cm^-1, median over the lowest 20 and over all nev. In a healthy run this is about 0.1-2 cm^-1 for the low states. Residuals of hundreds of cm^-1 or more that do not drop with M signal a basis or coordinate problem (sections 2.4 and 4), not too few points.
- max|Im(E)|, the number of kept orthogonalized vectors, and cond(S).
- If reference levels are given: errors vs the reference (median and max, for levels 1-20 and for all levels). Use them only to report, never to tune widths, Vmax or Ec, unless the user explicitly says to.
- Flag spurious states: eigenvalues below the ZPE region, or with large Im or residual.
- **Final convergence statement (required wording).** Give N and M used, then per level the largest deviation from the other runs (larger and smaller N, and a second independent point set), and summarize as: "with these N and M, the first n levels are converged to 1 cm^-1 or better; over all levels the RMS (and MAE) deviation is x and the maximum y cm^-1". Say which statistic is meant. Matching levels across runs by energy order inflates the maxima where levels are dense; say so.

## 7. Running and parallelizing

- The accumulation is embarrassingly parallel over point chunks. Use joblib/multiprocessing workers, each returning partial (S, H) that are summed, with BLAS threads per worker set so that workers x threads is about the core count (e.g. 12 x 4 on 48 cores). The final dense steps (eigh of S, products, eigs) use all threads in one process.
- Memory: S and H are 8*N^2 bytes each (N = 20000 needs 3.2 GB per matrix). Size the chunks so B and TB for one chunk fit comfortably; with gradient rows, chunk over the shifted q sets as well.
- Run long jobs in the background with a log file and keep printed output terse (one line per batch). Report results when done instead of streaming output.
- To stop a background job, kill it by PID. Never run `pkill -f <script>` from a shell whose own command line contains the script name: it kills that shell too.
- File layout: `config.json` (all inputs and defaults used), `cache/` (energies, gradients, dipoles, checkpoints), `results/levels_<tag>_N<N>_M<M>.csv`, `run.log`.

## 8. Options the user may invoke later

- **Optimize Vmax (only on request).** Scan Vmax (and Ec for Gaussian/IMQ centres) around the initial value at fixed N and M/N. Judge by validation residual and level stability. Extending Vmax upward reuses cached energies for points already computed. Report the chosen value and why.
- **Change N.** Rebuild the basis with the new N, regenerate the points in batches of the new N (reusing cached energies by R2 index), restart the accumulation, and compare the spectra across N at matched M/N.
- **Change coordinates, basis type or stencil.** Re-run the section 3 checks, re-tune widths, and compare with the previous run in one table.
- **Gradient rows.** See section 5b.
- **IR intensities.** See section 8b.

## 8b. Option: IR intensities (from ab initio dipoles at the points)

1. **Wavefunctions.** psi_f = B C_f from the eigenvectors. Make them real with the phase that makes sum C^2 real, and check that the leftover imaginary part is about 0.
2. **Volume element - essential.**
   - The Cartesian Laplacian restricted to rotation-invariant functions of q is Hermitian with the measure induced by d^3N X, not with dq.
   - That measure is rho(q) = |det[dX/dt_a, dX/domega_a (= e_a x X), dX/dq_k]| at X(q) in the Eckart frame, i.e. translations, infinitesimal rotations and vibrational coordinates.
   - Use rho in every integral: norms, overlaps and matrix elements.
   - Check: with rho, the collocation states are orthogonal (H2O: median overlap 2e-4, max 0.05, against max 0.3 with plain dq). Without it, the permanent dipole leaks into the A1 bands; the H2O nu1 band came out at 576 instead of 3 km/mol.
3. **Transition dipoles.** Use mu_0f = <0| mu - <0|mu|0> |f> with rho-weighted normalization. Subtracting <0|mu|0> removes the error from any residual non-orthogonality. Report <0|f> to show that it is small.
4. **Two integration routes, as a cross-check.**
   - (a) Quasi-Monte-Carlo over the DFT points themselves (fit plus validation points, uniform in the region), with weights rho.
   - (b) Fit each dipole component with a polynomial in q (degree about 6), checking its RMS error on the validation points against the range of the dipole. Then integrate on a dense R2 grid (2e5+ points) in the same region, with weights rho.
   - The two should agree for all bands above about 0.1 km/mol. Prefer (b) for the weak bands.
5. **Units.** The integrated band intensity is A(km/mol) = 8 pi^3 N_A nu |mu_0f|^2 / (3 h c 4 pi eps0) / 1000, with nu in m^-1 and mu in C m (converted from atomic units via e a0).
   - Check the constant against the double-harmonic formula (42.256 km/mol per (D/A)^2/amu).
   - Compare the fundamentals with double-harmonic intensities from FD dipole derivatives at the minimum.
6. **Assignment and symmetry.** Get the symmetry from the sign of psi under the symmetry operation expressed in q (e.g. q_antisym -> -q_antisym). Get approximate (v1 v2 ...) labels from the largest overlap with harmonic-oscillator product functions, and report that overlap weight. Flag weights < 0.4 as mixed.
7. **Output.** A CSV of levels with assignment, symmetry, transition energy, intensity from both routes, residual and flags. A stick-spectrum plot, with a log intensity axis when the intensities span many decades, and the harmonic fundamentals as markers if useful.

## 8c. Comparison with GVPT2 / VCI (when requested)

- Use ONE harmonic spectrum per surface: the harmonic frequencies of the surface the anharmonic methods actually use (same DFT settings as the collocation points). If the GVPT2 stencil coordinates are scaled by frequencies w0 from another calculation, the fitted diagonal quadratic coefficient c_ii gives the true frequency w = sqrt(2 c_ii w0), not 2 c_ii; rescale q' = q sqrt(w/w0) (and all polynomial coefficients) before setting up H0. Otherwise every GVPT2 level is off by (w - w0)/2 per quantum (LiCO2 aDZ: up to 0.7 cm^-1 per quantum).
- Do the GVPT2 stencil in the level of theory's own normal coordinates; a stencil in another basis set's coordinates gave wrong aTZ GVPT2 levels.
- Extend GVPT2 targets to all states up to the top collocation level (not just <= 2 quanta), with polyad partners up to ~300 cm^-1 higher, and converge VCI in its basis size (K).
- Match levels by harmonic-oscillator label only for clearly assigned states (weight >= 0.4 in the collocation state); mixed states give spurious label-matched differences of hundreds of cm^-1. Report the number of pairs used.
- For complexes, compare shifts relative to the free fragment (e.g. LiCO2 vs free CO2) level by level and per method; Fermi-resonant fragment levels make harmonic shifts misleading.

## 9. Report

Start with: method, coordinates, basis, N, M, Vmax, Ec, widths, stencil and h, energy source, and any geometric cut. Then give:
- the ZPE and the first transitions;
- comparison tables with the harmonic values AND every anharmonic method side by side for the same 20-50 levels (not just fundamentals), the ZPE as level 0, a state label and (for linear molecules) l; convergence information as text under the table, not as extra columns unless asked;
- a convergence table vs M (and vs N if run): ZPE, level changes, validation residual, max|Im|, errors vs reference if given; with gradient rows, the energies-only vs gradient comparison at matched DFT cost;
- the full list of nev levels saved as CSV (with intensities if computed);
- run time per stage and cores used.

State plainly what is not converged or coordinate-limited and what would most likely improve it (larger N, different coordinates, Ec/Vmax, widths). Any figures use large fonts and only the panels asked for, no figure or panel titles, legends that state level of theory, N and M, legends that do not cover data, and no wasted white space. Mirror stick spectra (method 1 up, method 2 down) come in a linear and a log-intensity version; convergence ranges can be drawn as a pale band around each stick. For a manuscript page, size figures to the printable area (A4 landscape with 1-inch margins: 9.69 x 6.27 in; two stacked plots about 9.5 x 2.5 in each).