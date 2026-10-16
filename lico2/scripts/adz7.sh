#!/bin/bash
# extra aDZ points for the 100-state run: wide box (cache_adz7/box.npz), region V <= 6500 cm-1. Resumable. Log: night_adz7.log
cd ~/lico2_work/lico2_project/rc; . ./env.sh
{ echo "adz7 start $(date)"
python dft2.py fit --tag adz7 --basis aug-cc-pvdz --n 24000 --vmax 6500 --nproc 36
python dft2.py val --tag adz7 --basis aug-cc-pvdz --n 2500 --vmax 6500 --nproc 36
echo "adz7 exit $? $(date)"; } >> night_adz7.log 2>&1
