#!/bin/bash
# Figure 6 of the manuscript (LiCO2 spectra, 18 Nov 2026 version): collocation with the max-residual parameters (N = 14000, first point set)
# up, GVPT2 down; no convergence band (the ranges are in results/LiCO2_levels_aDZ_aTZ.xlsx). DejaVu Sans font.
# Output: figures/lico2_IR_summary_fmax.png and _log.png (reproduces the committed files byte for byte in the cloud session's matplotlib).
set -e
cd "$(dirname "$0")"
export MPLCONFIGDIR=$(mktemp -d); echo "font.family: DejaVu Sans" > $MPLCONFIGDIR/matplotlibrc
COMPACT=1 XMIN=250 XMAX=2200 FIGW=9.75 FIGH=2.75 FONT=12 python3 plot_summary.py figures/lico2_IR_summary_fmax \
  "results/ir_colloc_fmax_adz_N14000.json:collocation, B3LYP-D3(BJ)/aug-cc-pVDZ, N = 14000, M = 12248:C0:-" \
  "results/ir_colloc_fmax_atz_N14000.json:collocation, B3LYP-D3(BJ)/aug-cc-pVTZ, N = 14000, M = 11690:C3:-." "|" \
  "vpt2_adz_hi.json:GVPT2, B3LYP-D3(BJ)/aug-cc-pVDZ:C0:-" "vpt2_atz_hi.json:GVPT2, B3LYP-D3(BJ)/aug-cc-pVTZ:C3:-."
