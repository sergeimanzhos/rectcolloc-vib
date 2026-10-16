# LiCO2 experimental IR spectra — handoff

Goal: experimental vibrational reference data for LiCO2 to compare against computed spectra.

## Bottom line
- **No gas-phase vibrational spectrum of LiCO2 (neutral or anion) exists** in the literature (searched Google Scholar, reference lists, citing articles). Jestilä et al. (2020) explicitly could not record cryogenic ion-trap spectra of MCO2⁻ (signal too low).
- Best available data: **inert-gas matrix IR** — Kafafi et al. 1983 (Ar, 15 K; also Kr, Xe, N2, CO2) and Manceron et al. 1985 (Kr, 9 K). Two isomers: **C2v** (bidentate, Li bound to both O) and **Cs** (Li bound to one O).
- **Use tabulated band positions** (`data/LiCO2_matrix_IR_band_positions.csv`) for frequency comparisons. Use digitized curves only for band shapes / relative intensities.

## Key frequencies (cm⁻¹)

| Isomer | Mode | Ar 15 K (Kafafi) | Kr 9 K (Manceron, nat. Li) |
|---|---|---|---|
| C2v | ν4(B2) asym CO2 str | 1568.6, 1569.9 (sh 1561.9, 1564.4) | 1568.0 / 1566.0 |
| C2v | ν1(A1) sym CO2 str | 1329.9, 1325.0 | 1326 / 1304 |
| C2v | ν2(A1) δ(OCO) | 798.7 | 795 |
| C2v | ν3(A1) sym Li–O str | 533.4, 526.0, 520.0 | 522 / 509 |
| C2v | ν5(B2) asym Li–O str | not obs. (force-field fit 400) | — |
| Cs | ν1(A′) C=O str | 1750.9, 1755.7 | 1750 / 1745 |
| Cs | ν2(A′) C–O str | 1221.4, 1208.7 | 1222 |
| Cs | ν3(A′) δ(OCO) | 739.5 | 740 |
| Cs | ν4(A′) Li–O str | 611.6, 599.3 | 588 |
| Cs | ν5(A′) δ(LiOC) | not obs. (force-field fit 200) | — |

Paired values = matrix-site splittings. Kafafi "calcd" values (Tables II, VI) are force-field fits, not ab initio. Kafafi also tabulate ⁶Li, ¹³C, ¹⁸O isotopomer frequencies (Tables II, VI) — useful for checking computed isotope shifts.

Assignment disagreements: ~1304 cm⁻¹ = C2v doublet partner (Manceron, Kr) vs LiCO2·CO2 adduct (Kafafi, Ar). 509 cm⁻¹ in Kr is listed for both C2v LiCO2 and Li2(CO2)2.

Co-existing products in the spectra (don't misassign): LiCO2·CO2 (1574.9, 1304.9, 807.0), LiC2O4 (1968.7, 1159.8), Li2C2O4 (1662.2, 1315.5, 807.0, 497.7), Li2CO2 (1447.9, 984.2, 565.9), CO2 bend ~662 (saturated).

## Files
- `data/LiCO2_matrix_IR_band_positions.csv` — band table from Kafafi 1983 Table I/II/VI and Manceron 1985 Table 1 (all species).
- `data/Kafafi1983_Fig1A_LiCO2_Ar_matrix.csv`, `...Fig1B...` — digitized Fig. 1 of Kafafi 1983 (Ar, 15 K). A = as deposited; B = after 5 h IR (Nernst glower) exposure (Cs → C2v conversion). Columns:
  - `wavenumber_axis_cm-1` — from printed axis ticks (up to ~20 cm⁻¹ low in places).
  - `wavenumber_band_corrected_cm-1` — piecewise-linear correction anchored on tabulated bands (A: 1968.7, 1750.9, 1662.2, 1159.8, 662.5; B: 1968.7, 1662.2, 1568.6, 1325.0, 1159.8, 662.5). Good to ~±5 cm⁻¹ near anchors; residual errors up to ~10–20 cm⁻¹ remain in 750–850 and 1400–1450 cm⁻¹.
  - `transmittance_percent` — on each trace's own printed scale (A: 100–50 %, B: 80–10 %; B baseline slightly >100 % as printed).
  - `clipped_overlap_flag` — 1 where trace A's deepest bands (≈1575, ≈662) run into trace B and were clipped (depth unreliable).
- `figures/LiCO2_Ar_matrix_summary.png` — both traces + stick rows for all assigned species (Ar) and LiCO2 isomers (Kr).
- `figures/Kafafi1983_Fig1_overlay_check.png` — digitized traces over the source figure (QC).
- `figures/Kafafi1983_Fig1_source_capture.png` — screen capture of the original figure (800×800 px, ~2.75 cm⁻¹/px).
- `scripts/` — run from inside `scripts/`: `python digitize.py && python finalize.py && python peaks.py && python summary_plot.py` (numpy, scipy, pillow, matplotlib). Reproduces `data/` exactly.

## Caveats for comparison with computation
- Matrix shifts relative to gas phase: typically a few to tens of cm⁻¹ (Ar usually smallest; compare Ar vs Kr columns).
- LiCO2 is strongly ionic (Li⁺CO2⁻); environment effects may be larger than for neutral closed-shell molecules.
- Experimental values are anharmonic fundamentals; Kafafi estimate anharmonicity x ≈ −60 cm⁻¹ for the ν4-type CO2⁻ stretch in Ar.

## Draft text for a manuscript

No gas-phase vibrational spectrum of LiCO2 has been reported, to the best of our knowledge; even cryogenic ion-trap spectroscopy of the related anions MCO2⁻ did not yield sufficient signal [3]. We therefore compare with infrared spectra of LiCO2 isolated in inert-gas matrices. Kafafi et al. [1] co-deposited Li atoms with CO2 in excess argon at 15 K and identified two geometrical isomers: a symmetric (C2v) form, with bands at 1568.6, 1325.0, 798.7 and 533.4 cm⁻¹, and an asymmetric (Cs) form, with bands at about 1751/1756, 1209/1221, 739.5 and 599/612 cm⁻¹; the assignments were confirmed by ⁶Li, ¹³C and ¹⁸O substitution. Manceron et al. [2] independently observed both isomers in a krypton matrix at 9 K at very similar frequencies. Matrix environments typically shift vibrational frequencies by a few to several tens of cm⁻¹ relative to the free molecule, and the spectra also contain bands of the LiCO2·CO2 adduct, LiC2O4, Li2C2O4 and Li2CO2, which should be kept in mind when comparing with computed spectra.

[1] Z. H. Kafafi, R. H. Hauge, W. E. Billups, J. L. Margrave, "Carbon dioxide activation by lithium metal. 1. Infrared spectra of Li⁺CO2⁻, Li⁺C2O4⁻, and Li2²⁺CO2²⁻ in inert-gas matrices," J. Am. Chem. Soc. 105, 3886–3893 (1983). https://doi.org/10.1021/ja00350a025

[2] L. Manceron, A. Loutellier, J. P. Perchard, "Reduction of carbon dioxide to oxalate by lithium atoms: a matrix isolation study of the intermediate steps," J. Mol. Struct. 129, 115–124 (1985). https://doi.org/10.1016/0022-2860(85)80197-6

[3] J. S. Jestilä, J. K. Denton, E. H. Perez, T. Khuu, E. Aprà, S. S. Xantheas, M. A. Johnson, E. Uggerud, "Characterization of the alkali metal oxalates (MC2O4⁻) and their formation by CO2 reduction via the alkali metal carbonites (MCO2⁻)," Phys. Chem. Chem. Phys. 22, 7460–7473 (2020). https://doi.org/10.1039/D0CP00547A

(Kafafi et al., Inorg. Chem. 23, 177 (1984) covers Na/K/Cs only — not needed for Li.)

## Possible next steps
- Digitize Kafafi Fig. 3 (¹⁸O isotopomers, Ar) and Fig. 4 (Xe matrix) if isotope-resolved or Xe spectra are needed.
- Compare computed harmonic/anharmonic frequencies and IR intensities for both isomers against the table; check ⁶Li/¹³C/¹⁸O shifts against Kafafi Tables II & VI.
