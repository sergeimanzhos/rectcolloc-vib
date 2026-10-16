#!/bin/bash
cd ~/lico2_work/lico2_project/rc; . ./env.sh
for EC in 900 1500; do
python run_rc.py --pes harm --N 10000 --ec $EC --beta 0.95 --vmax 5500 --p 0 --Ms 10000 --nev 100 --workers 4 --threads 1 --tag harm100_N10000_ec$EC
done
echo "harm100 done $(date)"
