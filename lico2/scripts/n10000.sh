#!/bin/bash
# N = 10000: harmonic sanity test, then quartic-FF validation against VCI; log: n10000.log
cd ~/lico2_work/lico2_project/rc; . ./env.sh
python run_rc.py --pes harm --N 10000 --ec 900 --beta 0.9 --vmax 5000 --p 3 --Ms 5000,10000 --nev 30 --workers 4 --threads 1 --tag harm_N10000
python run_rc.py --pes ff --N 10000 --ec 900 --beta 0.9 --vmax 5000 --p 3 --Ms 5000,10000 --nev 40 --workers 4 --threads 1 --ref ../vpt/vpt2_lico2_big.json --tag ffval_N10000
echo "n10000 done $(date)"
