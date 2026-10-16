#!/bin/bash
# after the DFT points: collocation parameter scans and convergence series for CO2 (J = 0 and J = 1 f blocks). Log: colloc.log
cd ~/lico2_work/lico2_project/co2_atz; . ../rc/env.sh
while ! grep -q -E "pts done" pts.log; do sleep 30; done
export OPENBLAS_NUM_THREADS=24 NUMEXPR_NUM_THREADS=24 NUMEXPR_MAX_THREADS=24 OMP_NUM_THREADS=24
python co2_colloc.py scan
echo "colloc done $? $(date)"
