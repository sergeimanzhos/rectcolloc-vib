#!/bin/bash
# N = 10000 on the aDZ DFT surface with a fixed set of the first 9000 fit points (V <= 5500): width, Ec and p. Log: iscan.log
cd ~/lico2_work/lico2_project/rc; . ./env.sh
run() { python run_rc.py --pes dft --cache cache_adz --N 10000 --ec $1 --beta $3 --vmax 5500 --p $2 --Ms 9000 --nev 40 --workers 6 --threads 1 --vval 4500 --tag iscan_N10000_ec$1_p$2_b$3; }
run 900 1 0.99; run 900 1 1.05; run 900 0 0.99
while [ ! -f results/centres_dft_cache_adz_ec1300_N10000.npy ]; do sleep 30; done
run 1300 1 0.99
echo "iscan done $(date)"
