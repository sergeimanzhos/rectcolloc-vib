#!/bin/bash
# more aDZ points for N and M convergence (resumes cache_adz): fit to 16000 in region, val to 1600. Log: night.log
cd ~/lico2_work/lico2_project/rc; . ./env.sh
while pgrep -f "dft2.py val --tag atz" > /dev/null; do sleep 20; done
{
echo "adz_more start $(date)"
python dft2.py fit --tag adz --basis aug-cc-pvdz --n 16000 --vmax 5500 --nproc 40
python dft2.py val --tag adz --basis aug-cc-pvdz --n 1600 --vmax 5500 --nproc 40
echo "adz_more exit $? $(date)"
} >> night.log 2>&1
touch adz_more.done
