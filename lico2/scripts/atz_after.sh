#!/bin/bash
# after the aTZ points: aTZ collocation with the aDZ-optimised parameters (N = 10000), 10 workers. Log: final_atz.log
cd ~/lico2_work/lico2_project/rc; . ./env.sh
while pgrep -f "dft2.py (fit|val) --tag atz" > /dev/null; do sleep 60; done
read EC VMAX P BETA < final_adz.params
{ echo "final aTZ: ec=$EC p=$P beta=$BETA  $(date)"
python run_rc.py --pes dft --cache cache_atz --N 10000 --ec $EC --beta $BETA --vmax 5000 --p $P --Ms 6000,9000,1000000 --nev 40 --workers 10 --threads 1 --vval 4500 --tag final_atz_N10000
echo "final aTZ exit $? $(date)"; } >> final_atz.log 2>&1
