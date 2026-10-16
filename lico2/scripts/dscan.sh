#!/bin/bash
# aDZ DFT surface: gradient-row fade p and centre energy Ec at N = 5000 (beta 0.99 from the tuner so far), M = half and all points,
# then N = 10000 with all points. Judged by validation residual (val points with V <= 4500) and stability. Log: dscan.log
cd ~/lico2_work/lico2_project/rc; . ./env.sh
run() { python run_rc.py --pes dft --cache cache_adz --N $1 --ec $2 --beta $4 --vmax 5000 --p $3 --Ms 2666,1000000 --nev 40 --workers 6 --threads 1 --vval 4500 --tag dscan_N$1_ec$2_p$3_b$4; }
run 5000 900 3 0.99; run 5000 900 1 0.99; run 5000 900 0 0.99; run 5000 1400 3 0.99; run 5000 1400 1 0.99
run 10000 900 1 0.9; run 10000 900 1 0.99
echo "dscan done $(date)"
