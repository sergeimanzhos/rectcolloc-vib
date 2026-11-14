#!/bin/bash
# (run in cl2co/; the outputs were then moved to results/)
# 0) check: old widths (BASEW default) must reproduce results/dft_adz_M5000.json of the repository
# 1) option 3: DFT collocation with the PES max-residual widths, otherwise exactly as workflow_dft_collocation.sh steps 3 and 5
set -e
OMP_NUM_THREADS=2 python3 run_dft.py adz 5000 5000 3 "2500,5000" check_oldw_dft_adz_M5000
export WIDTHS=$(python3 -c "import json;print(','.join(repr(b) for b in json.load(open('results/best_widths_pes_max.json'))['beta']))")
echo "WIDTHS=$WIDTHS"
for T in adz atz; do
  OMP_NUM_THREADS=2 python3 run_dft.py $T 5000 5000 3 "2500,5000" pesw_dft_${T}_M5000
  DIP=dipfit_${T}_d4.npz OMP_NUM_THREADS=2 python3 ir6.py pesw_dft_${T}_M5000 5000 5000
done
echo finished
