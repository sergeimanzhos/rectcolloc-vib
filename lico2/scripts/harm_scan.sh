#!/bin/bash
# harmonic sanity test at N = 5000 (extends the HANDOFF table beyond N = 3000); log: harm_N5000.log
cd ~/lico2_work/lico2_project/rc; . ./env.sh
python run_rc.py --pes harm --N 5000 --ec 900 --beta 0.9 --vmax 5000 --p 3 --Ms 2500,5000 --nev 30 --workers 4 --threads 1 --tag harm_N5000
echo "harm done $(date)"
