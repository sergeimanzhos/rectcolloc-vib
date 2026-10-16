# Rectangular collocation for molecular vibrations with potential-gradient rows

Rectangular (M > N) collocation for vibrational levels with an exact Cartesian finite-difference
kinetic-energy operator, optional potential-derivative (gradient) rows, and IR intensities.
Applications: Cl2CO (NN PES and DFT points), LiCO2 (DFT points, 100 states, IR spectra, comparison with GVPT2 and Ar-matrix data)
and free CO2 (linear molecule, vibrational angular momentum blocks, comparison with GVPT2 and VCI).

## Contents

| Folder | Content |
|---|---|
| `skill/rectangular-collocation-vib/SKILL.md` | The procedure as a Claude skill (updated October 2026 with the LiCO2 / CO2 lessons) |
| `cl2co/` | Cl2CO: NN-PES collocation, DFT-point collocation (aug-cc-pVDZ / aug-cc-pVTZ), GVPT2, IR intensities |
| `lico2/` | LiCO2 (C2v, doublet): UB3LYP-D3(BJ) collocation for 100 states, IR intensities, GVPT2, figures, Excel summary |
| `co2/` | Free CO2: B3LYP-D3(BJ)/aug-cc-pVTZ collocation (Gaussian basis, energy and gradient rows), fitted quartic field, GVPT2, VCI |

## Cl2CO (`cl2co/`)

| File | Purpose |
|---|---|
| `pes.py`, `fitNN100NNs4.dat`, `Cl2COpes.m` | NN PES (Python port of the MATLAB evaluator); other NN fits `fitN*.dat` |
| `eq.py`, `coords6.py`, `r2.py` | Equilibrium and Hessian; Eckart-frame dimensionless normal coordinates, volume element; R2 sequence |
| `collo6.py`, `run.py` | Gaussian basis, exact FD KEO, gradient rows, accumulated S/H, solver; driver `python3 run.py N vmax cfg p fac ec "M1,..." tag` |
| `colloc_dft.py`, `run_dft.py`, `workflow_dft_collocation.sh` | Collocation on B3LYP aug-cc-pVDZ / aug-cc-pVTZ points (`dft_*_c.npz` fit points, `dft_*_v.npz` validation points) |
| `ir6.py`, `fit_dip.py`, `dipder*.py` | IR intensities; dipole-surface fits (`dipfit_*.npz`); linear dipole derivatives (B3LYP, CCSD(T)) |
| `vp_dft.py`, `vp_fit.py`, `vp2.py`, `vp_nn.py` | Numerical 6-mode GVPT2 (stencils `vp_*.npz`, results `vpt2_*.json`) |
| `results/` | Level and intensity JSON/CSV files, Excel summaries (`Cl2CO_PES_collocation_vs_reference.xlsx`, `Cl2CO_aDZ_aTZ_collocation_vs_VPT2.xlsx`) |
| `figures/` | Spectra and geometry figures |

Best settings (NN PES): widths (0.857375, 0.857375, 0.95, 0.9025, 0.95, 0.95); Ec = 1150 cm-1; Vmax = 5000 (N <= 5000) or 5500 cm-1 (N = 10000);
gradient rows weighted by (1 - V/Vmax)^3. Reminder: the nu1 / nu2+nu5+nu6 Fermi resonance affects the GVPT2 nu1.

## LiCO2 (`lico2/`)

UB3LYP-D3(BJ), aug-cc-pVDZ (aDZ) and aug-cc-pVTZ (aTZ), PySCF (grid 4, gradients with grid response).

| File | Purpose |
|---|---|
| `coords.py`, `collo.py`, `run_rc.py`, `r2.py`, `gen_centres.py` | Eckart normal coordinates, Gaussian basis + FD KEO + gradient rows, driver |
| `dft_points.py`, `dft2.py` | DFT energies, gradients and dipoles at R2 points (cached in `data/cache_*`) |
| `tune10k.py`, `tune100.py`, `f100.params` | Width / centre tuning; final parameters (Ec = 900 cm-1, widths 0.8935, 1.0395, 0.99, 0.99, 1.0395, 0.99) |
| `intens.py`, `mkrange.py` | IR intensities (volume element, degree-6 dipole fit, complex-pair handling); convergence ranges |
| `vp_dft_li.py`, `vp_dft_modes.py`, `vp_fit.py`, `vp2_li_hi.py` | GVPT2: stencil, fit, solver (zeroth-order frequencies rescaled to the fitted surface, Oct 2026) |
| `plot_vpt2_li.py`, `plot_zoom.py`, `plot_summary.py` | Mirror plots (linear and log), convergence bands, A4 layouts |
| `mol_lico2.py`, `orbitals_lico2.py`, `orb_inplane.py` | Geometry and orbital figures |
| `vpt_ff/` | Finite-difference force fields (LiCO2, CO2) and early VPT2 analysis |
| `data/cache_{adz,adz7,atz,atz7}/` | DFT points: `fit.npz`, `val.npz` (energies, gradients, dipoles), E0, reference density, box |
| `results/` | Level / intensity JSON and CSV, `LiCO2_levels_aDZ_aTZ.xlsx` (100 levels x 4 methods), `LiCO2_summary.md` |
| `figures/`, `exp/` | Figures; digitized Ar-matrix spectra and band table (Kafafi et al., JACS 105, 3886 (1983)) |

Final runs: N = 14000, Vmax = 5500 cm-1, M = 16150 (aDZ) / 13831 (aTZ), energy + gradient rows (p = 0).
Convergence (largest deviation over N = 10000/14000 and two independent point sets): first 25 levels within 1 cm-1;
RMS 1.9 cm-1 and max 8.8 (aDZ) / 7.2 (aTZ) cm-1 over all 99 levels. Wavefunction files `vec_*.npz` (23 MB each) are not in the repository.

## Free CO2 (`co2/`)

B3LYP-D3(BJ)/aug-cc-pVTZ, re = 1.1604 A.

| File | Purpose |
|---|---|
| `co2c.py` | Internal coordinates (x1, x3, xb), geometry, FD KEO |
| `co2_points.py` | DFT energies at R2 points (`pts_fit.npz` 8157 points <= 8000 cm-1, `pts_val.npz`) |
| `co2_ffbench.py` | Gaussian basis with I0-projected bend function; quartic-field benchmark |
| `co2_lev.py` | l = 0..4 blocks (J = l harmonic-tensor factor), energy and gradient rows |
| `co2_grad.py` | Gradient rows in internal coordinates (constant Cartesian directions); DFT gradients (`pts_grad.npz`, 2400 points) |
| `co2_tune.py` | Width / N / M tuning against VCI on the fitted quartic field |
| `fit_ff.py`, `vff_fit.json` | Symmetry-adapted polynomial fit to the DFT energies -> quartic force field |
| `vpt2_cyl.py`, `vpt_lev.py` | GVPT2 and VCI on the quartic field, with vibrational angular momentum l |
| `lev_final2.py`, `results/co2_levels_final.*` | Final table: harmonic, GVPT2, VCI, collocation for 34 levels (l = 0-4) up to 4350 cm-1 |
| `co2ho.py` | Harmonic-oscillator product basis (used only as an internal check; imported by `co2_lev.py`) |

Final collocation: Gaussians, widths 0.8/0.8/0.65, centres and points with V <= 8000 cm-1, N = 4000, M = 2400 points with
energy + gradient rows (weight (1 - V/Vmax)^3). Convergence: N 2400 -> 4000 changes levels by <= 0.24 cm-1 (0.05 below 3000 cm-1);
M 1200 -> 2400 by <= 0.31 cm-1 (0.06 below 3000 cm-1).

Requirements: Python 3, numpy, scipy, numexpr, matplotlib, openpyxl (Excel); PySCF (with dftd3) for the DFT/CCSD(T) scripts.
On Windows, PySCF runs in WSL.
