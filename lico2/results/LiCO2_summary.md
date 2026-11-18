# LiCO₂ (C₂ᵥ): rectangular collocation vs GVPT2, aug-cc-pVDZ vs aug-cc-pVTZ

Status: 13 October 2026. All calculations finished; nothing is running.

**Superseded for the collocation numbers (18 Nov 2026):** the collocation results were redone with parameters tuned by the maximum validation
residual (Ec = 1200, Vmax = 5000, p = 0, widths 0.9321 / 0.99 / 1.00 / 1.0188 / 1.00 / 1.00); see `../MAXRES_HANDOFF.md` section 4b and
`LiCO2_levels_aDZ_aTZ.xlsx`. The tables below are the October (Ec = 900) version.

## Method

- **Surface:** UB3LYP-D3(BJ) with aug-cc-pVDZ (aDZ) and aug-cc-pVTZ (aTZ), PySCF, grid 4, gradients with grid response.
- **Collocation:** exact Cartesian finite-difference KEO, Eckart-frame dimensionless normal coordinates (aDZ reference),
  Gaussian basis, N = 14000, centres with V ≤ 900 cm⁻¹, widths 0.8935 / 1.0395 / 0.99 / 0.99 / 1.0395 / 0.99,
  energy plus gradient rows with constant weight (p = 0), Vmax = 5500 cm⁻¹, 100 states (to about 2165 cm⁻¹ above the ZPE).
- **GVPT2:** numerical GVPT2 code of the Cl₂CO repository, 85-point stencil (δ = 0.25) with the same DFT settings,
  each basis set in its own normal coordinates, intensities from the full polynomial dipole.
- **Intensities (collocation):** degree-6 fitted DFT dipole surface, ρ-weighted dense grid, Löwdin orthonormalisation.

## Fundamentals

Position in cm⁻¹; intensity in km/mol in brackets.

| Mode | aDZ collocation | aDZ GVPT2 | aTZ collocation | aTZ GVPT2 | Ar matrix (Kafafi 1983) |
|---|---|---|---|---|---|
| ν1 (B1, Li out-of-plane) | 311.6 (31.3) | 312.2 (31.4) | 303.6 (30.5) | 305.2 (30.6) | not measured |
| ν2 (B2, Li in-plane) | 366.2 (24.9) | 364.0 (24.6) | 364.5 (24.9) | 364.8 (24.6) | not observed |
| ν3 (A1, Li–O₂ stretch) | 555.0 (118.7) | 550.1 (117.3) | 555.0 (115.6) | 553.0 (115.1) | 533.4, 526.0, 520.0 |
| ν4 (A1, OCO bend) | 780.9 (0.77) | 781.5 (1.05) | 790.2 (0.53) | 791.2 (0.87) | 798.7 |
| ν5 (A1, sym. CO stretch) | 1294.2 (0.015) | 1295.4 (0.048) | 1303.5 (0.035) | 1304.8 (0.070) | 1329.9, 1325.0 |
| ν6 (B2, asym. CO stretch) | 1559.9 (435) + 1564.7 (147) | 1562.0 (596) | 1570.3 (585) | 1571.0 (611) | 1568.6, 1569.9 |
| Zero-point energy | 2478.3 | 2478.4 | 2487.9 | 2488.3 | |

Harmonic frequencies of the same DFT surfaces (cm⁻¹), ν1 … ν6: aDZ 311.1, 388.2, 581.4, 786.6, 1322.3, 1600.9;
aTZ 300.5, 387.2, 581.7, 796.8, 1332.3, 1609.7.

*Correction, 14 October 2026:* the aDZ GVPT2 values were recomputed. The earlier run took its zeroth-order frequencies as
2c_ii in coordinates scaled by the reference (integration grid 8) frequencies, which put the aDZ GVPT2 levels too high by
0.05–0.7 cm⁻¹ per quantum. aTZ was not affected.

## Overtones and combinations with visible intensity

| Band | aDZ collocation | aDZ GVPT2 | aTZ collocation | aTZ GVPT2 |
|---|---|---|---|---|
| 2ν1 | 625.0 (0.16) | 627.4 (0.09) | 611.1 (0.25) | 616.7 (0.39) |
| 2ν2 | 728.9 (9.2) | 724.9 (10.5) | 726.3 (9.9) | 728.5 (10.4) |
| 2ν3 | 1098.0 (0.79) | 1079.8 (0.42) | 1097.5 (0.77) | 1087.8 (0.38) |
| ν6 partner (mixed B2 state) | see ν6 | — | 1557.4 (13.4) | — |
| 2ν4 | 1561.4 (1.17) | 1563.1 (0.61) | 1580.2 (0.74) | 1582.4 (0.63) |
| ν2 + ν5 region, lower | 1654.0 (10.3) | 1647.6 (3.9) | 1660.8 (11.1) | 1661.1 (6.1) |
| ν2 + ν5 region, upper | 1669.7 (4.2) | 1663.5 (7.4) | 1677.3 (5.5) | 1676.1 (6.4) |
| ν4 + ν5 | 2066.3 (2.7) | 2068.6 (2.0) | 2084.9 (3.0) | 2087.9 (2.3) |
| ν3 + ν6 | 2116.7 (0.95) | 2109.8 (0.94) | 2125.2 (0.93) | 2122.3 (1.0) |

## Main findings

1. **Collocation and GVPT2 agree on the fundamentals** to 0.3–2 cm⁻¹ at aTZ and 0.5–5.0 cm⁻¹ at aDZ, and to a few percent in the
   intensities of the strong bands.
2. **The basis-set effect is larger than the difference between the two vibrational methods:** about 8–10 cm⁻¹ on ν1, ν4, ν5 and ν6.
3. **The ν6 doublet seen at aDZ is a basis-set effect.** At aDZ ν6 shares its intensity 435 : 147 with a neighbouring B2 state; at aTZ
   the partner keeps only 13 km/mol and collocation agrees with GVPT2 (585 vs 611 km/mol, 0.7 cm⁻¹ apart).
4. **Where the methods differ:** 2ν3 (collocation 18 cm⁻¹ higher at aDZ, 10 at aTZ, and twice as intense), 2ν1 (3–6 cm⁻¹ lower),
   and the intensity pattern of the 1650–1680 cm⁻¹ pair (collocation about 2 : 1, GVPT2 about 1 : 1 at aTZ and reversed at aDZ).
5. **Against the Ar-matrix values, aTZ is clearly better than aDZ** for ν4 (790 vs 799), ν5 (1304 vs 1325–1330) and ν6 (1570 vs 1569–1570).
   ν3 is 22–35 cm⁻¹ too high with both basis sets and both methods, and ν5 is computed far weaker than its observed weak band suggests.

## Convergence of the 100 collocation levels

Maximum (median) change of the levels, cm⁻¹.

| Check | Basis | Levels 1–40 | Levels 41–70 | Levels 71–99 |
|---|---|---|---|---|
| N = 10000 vs 14000 | aDZ | 1.0 (0.1) | 4.5 (0.8) | 8.8 (2.1) |
| N = 10000 vs 14000 | aTZ | 1.1 (0.06) | 3.9 (0.6) | 7.2 (1.6) |
| Two independent point sets, N = 14000 | aDZ | 0.5 (0.1) | 1.9 (0.2) | 2.1 (0.2) |
| Two independent point sets, N = 14000 | aTZ | 0.9 (0.09) | 1.4 (0.3) | 1.4 (0.2) |
| Vmax 5000 vs 5500, N = 10000 | aDZ | 0.4 (0.05) | 0.5 (0.09) | 0.8 (0.14) |
| Vmax 5000 vs 5500, N = 10000 | aTZ | 0.2 (0.01) | 0.5 (0.06) | 0.8 (0.11) |
| M = 12000 vs 16150, N = 10000 | aDZ | 0.3 (0.03) | 1.4 (0.08) | 1.2 (0.2) |

Width of the uncertainty band (range of each level over N = 10000, 14000 and the two point sets), cm⁻¹:

| Region above ZPE | aDZ levels | aDZ median / max | aTZ levels | aTZ median / max |
|---|---|---|---|---|
| below 1600 cm⁻¹ | 39 | 0.20 / 1.4 | 39 | 0.26 / 2.0 |
| 1600–1900 cm⁻¹ | 26 | 1.1 / 4.7 | 28 | 1.1 / 4.7 |
| 1900–2300 cm⁻¹ | 34 | 2.6 / 9.5 | 32 | 2.5 / 7.6 |

Runs behind the tables:

| Run | Points | Equations per basis function | Median validation residual (cm⁻¹) |
|---|---|---|---|
| aDZ, N = 10000, first point set | 16150 | 11.3 | 23.8 |
| aDZ, N = 14000, first point set (plotted) | 16150 | 8.1 | 16.7 |
| aDZ, N = 10000, second point set | 11396 | 8.0 | 26.0 |
| aDZ, N = 14000, second point set | 12606 | 6.3 | 18.9 |
| aTZ, N = 10000, first point set | 13831 | 9.7 | 25.1 |
| aTZ, N = 14000, first point set (plotted) | 13831 | 6.9 | 21.3 |
| aTZ, N = 10000, second point set | 14547 | 10.2 | 30.1 |
| aTZ, N = 14000, second point set | 14730 | 7.4 | 20.9 |

Other checks: harmonic sanity test for 100 states at N = 10000 gives a median error of 0.003 cm⁻¹ (max 0.7); the coordinate-singularity
test (h vs 2h) is clean to 10⁻⁴ cm⁻¹ and the sampled region stays bidentate; the centre cutoff and the six widths were tuned on the
validation residual only.

## What is not converged

- **Levels above about 1600 cm⁻¹ in N.** They still move by 1–2 cm⁻¹ typically and up to 7–9 cm⁻¹ between N = 10000 and 14000; no larger
  basis was run.
- **The aDZ ν6 intensity split** (435 : 147 here; between 384 : 200 and 487 : 95 in other runs). The sum, about 580 km/mol, is stable.
- **Labels and weak intensities of a few near-degenerate groups above 1900 cm⁻¹**, where states overlap strongly.
- **Validation residuals** of the upper states (tens of cm⁻¹) are above what the skill calls healthy, although the levels are stable in M and Vmax.

## DFT data available

| Set | Points computed | Below 5500 cm⁻¹ | Validation points |
|---|---|---|---|
| aDZ, first (original box) | 36240 | 16150 | 3920 computed |
| aDZ, second (wider box, to 6500 cm⁻¹) | about 110000 | 40184 below 6500 | 6900 computed |
| aTZ, first | 14463 | 13831 | 758 |
| aTZ, second (wider box) | 27376 | 22119 | 1735 |

The second aTZ set was only partly used (14730 of its points); the rest is available for a larger basis.

## Figures

| File | Content |
|---|---|
| `lico2_IR_summary_adz_atz.png`, `_log.png` | Summary: collocation up, GVPT2 down; aDZ blue solid, aTZ red dashed; uncertainty bands; A4 landscape |
| `lico2_IR_summary_adz_atz_zoom_1480_1720.png`, `_log.png` | Same, zoomed on the ν6 region |
| `lico2_IR_adz100_A4.png`, `_log.png` | aDZ only: GVPT2 up, collocation down, with bands |
| `lico2_IR_atz100_N14000_A4.png`, `_log.png` | aTZ only: GVPT2 up, collocation down, with bands |
| `lico2_IR_vpt2_adz_colloc_adz100_old_N14000.png`, `_log.png` | aDZ with the two Ar-matrix spectra |
| `lico2_IR_vpt2_atz_colloc_atz100_N14000.png`, `_log.png` | aTZ with the two Ar-matrix spectra |
| `lico2_IR_adz100_zoom_*.png` | aDZ zooms with bands and experiment |

Level and intensity tables: `intens_adz100_old_N14000.csv`, `intens_atz100_N14000.csv`. Every plotted spectrum is recorded in
`plotted_spectra_manifest.json` in the WSL project folder (`~/lico2_work/lico2_project/rc`).
