#!/bin/bash
# quartic-FF scan at N = 5000: gradient-row fade exponent p, centre energy Ec, and M up to 2N; log: pscan.log
cd ~/lico2_work/lico2_project/rc; . ./env.sh
run() { python run_rc.py --pes ff --N 5000 --ec $1 --beta 0.9 --vmax 5000 --p $2 --Ms 5000,10000 --nev 40 --workers 6 --threads 1 --ref ../vpt/vpt2_lico2_big.json --tag ffscan_ec$1_p$2; }
run 900 3; run 900 1; run 900 0; run 1400 3; run 2000 3; run 1400 1
echo "pscan done $(date)"
