#!/bin/bash
# free CO2, B3LYP-D3(BJ)/aug-cc-pVTZ, grid 8: optimisation, harmonic analysis, semi-quartic force field, GVPT2 + VCI. Log: vpt.log
cd ~/lico2_work/lico2_project/co2_atz; . ../rc/env.sh
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 MKL_NUM_THREADS=32
python ff.py co2 && python vpt2.py co2
echo "vpt done $? $(date)"
