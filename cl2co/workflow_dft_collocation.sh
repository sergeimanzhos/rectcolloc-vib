#!/bin/bash
# Collocation directly from DFT points, from scratch to spectrum.  usage: bash workflow_dft_collocation.sh <basis> <tag> <nproc>
#   e.g. bash workflow_dft_collocation.sh aug-cc-pvtz atz 24      (24 processes x 2 threads = 48 cores)
# Steps 1-2 are resumable: re-running skips points already in dft_<tag>_cache*.jsonl.
set -e
BAS=${1:-aug-cc-pvtz}; TAG=${2:-atz}; NP=${3:-24}
# 1. DFT energies + gradients + dipoles at M=5000 collocation points and 1000 validation points
for k in $(seq 0 $((NP-1))); do
  NT=2 SHARD=$k/$NP MAXC=5000 MAXV=1000 nohup python3 colloc_dft.py $BAS $TAG > log_${TAG}_$k.txt 2>&1 &
done
wait
# 2. merge caches
python3 rebuild_npz.py $TAG $BAS
# 3. collocation, N=5000, Vmax=5000 cm-1, gradient-row fade p=3, M = 2500 then 5000 (needs EMIN[tag] in run_dft.py)
OMP_NUM_THREADS=8 python3 run_dft.py $TAG 5000 5000 3 "2500" dft_${TAG}_M2500
OMP_NUM_THREADS=8 python3 run_dft.py $TAG 5000 5000 3 "2500,5000" dft_${TAG}_M5000
# 4. degree-4 polynomial dipole surface from the DFT dipoles
python3 fit_dip.py $TAG 4
# 5. intensities: full DFT dipole surface, and the linear B3LYP dipole used for the fitted-PES spectra
for M in 2500 5000; do
  DIP=dipfit_${TAG}_d4.npz OMP_NUM_THREADS=8 python3 ir6.py dft_${TAG}_M$M 5000 5000
done
echo "done: ir_dft_${TAG}_M*_dipfit_${TAG}_d4.json"
