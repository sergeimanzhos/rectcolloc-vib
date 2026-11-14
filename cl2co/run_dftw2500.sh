#!/bin/bash
# (run in cl2co/; the outputs were then moved to results/)
# option 4, M = 2500 intensities (for the Excel M = 2500 columns): same as run_dftw.sh but M list "2500"
set -e
export WIDTHS=$(python3 -c "import json;print(','.join(repr(b) for b in json.load(open('best_widths_dft_adz_max.json'))['beta']))")
echo "WIDTHS=$WIDTHS"
for T in adz atz; do
  OMP_NUM_THREADS=2 python3 run_dft.py $T 5000 5000 3 "2500" dftw_dft_${T}_M2500
  DIP=dipfit_${T}_d4.npz OMP_NUM_THREADS=2 python3 ir6.py dftw_dft_${T}_M2500 5000 5000
done
echo finished
