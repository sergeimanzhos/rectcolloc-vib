#!/bin/bash
# (run in cl2co/; the outputs were then moved to results/)
# option 4: DFT collocation with the widths tuned on the aDZ DFT data (best_widths_dft_adz_max.json), for both aDZ and aTZ,
# otherwise exactly as run_pesw.sh (workflow_dft_collocation.sh steps 3 and 5)
set -e
export WIDTHS=$(python3 -c "import json;print(','.join(repr(b) for b in json.load(open('best_widths_dft_adz_max.json'))['beta']))")
echo "WIDTHS=$WIDTHS"
for T in adz atz; do
  OMP_NUM_THREADS=2 python3 run_dft.py $T 5000 5000 3 "2500,5000" dftw_dft_${T}_M5000
  DIP=dipfit_${T}_d4.npz OMP_NUM_THREADS=2 python3 ir6.py dftw_dft_${T}_M5000 5000 5000
done
echo finished
