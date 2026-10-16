#!/bin/bash
cd ~/lico2_work/lico2_project/rc; . ./env.sh
while ! grep -q "iscan done" iscan.log; do sleep 30; done
# the iscan runs are the start candidates: give them the tuner's tags (same N, M=9000, validation set)
python - <<'PY'
import shutil,os
m={'iscan_N10000_ec900_p1_b0.99':'t10k_N10000_M9000_ec900_v5500_p1_b'+'-'.join(['0.9900']*6),
   'iscan_N10000_ec900_p0_b0.99':'t10k_N10000_M9000_ec900_v5500_p0_b'+'-'.join(['0.9900']*6),
   'iscan_N10000_ec1300_p1_b0.99':'t10k_N10000_M9000_ec1300_v5500_p1_b'+'-'.join(['0.9900']*6)}
for a,b in m.items():
    if os.path.exists(f'results/{a}.json') and not os.path.exists(f'results/{b}.json'): shutil.copy(f'results/{a}.json',f'results/{b}.json')
PY
python tune10k.py --workers 4
echo "tune10k done $(date)"
