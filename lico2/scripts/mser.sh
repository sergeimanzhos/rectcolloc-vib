#!/bin/bash
# M convergence at N = 10000 with the widths from the first optimisation pass (p = 0, Ec 900); fixed validation set. Log: mser.log
cd ~/lico2_work/lico2_project/rc; . ./env.sh
B=0.9405,1.0395,0.99,0.99,0.99,0.99
python run_rc.py --pes dft --cache cache_adz --N 10000 --ec 900 --beta $B --vmax 5500 --p 0 --Ms 6000,9000,12000,1000000 --nev 40 --workers 12 --threads 1 --vval 4500 --nvalraw 1920 --tag mser_v5500
python run_rc.py --pes dft --cache cache_adz --N 10000 --ec 900 --beta $B --vmax 5000 --p 0 --Ms 6000,9000,1000000 --nev 40 --workers 12 --threads 1 --vval 4500 --nvalraw 1920 --tag mser_v5000
echo "mser done $(date)"
