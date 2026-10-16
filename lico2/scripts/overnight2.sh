#!/bin/bash
# Overnight chain for the 100-state aDZ run (2026-10-08), self-contained. Master log: rc/overnight2.log
#  A. when the aTZ point job ends: restart the extra aDZ point job on 64 processes (target 40000 points, V <= 6500) and
#     compute its validation points in parallel (10 processes)                                         -> night_adz7.log, night_adz7val.log
#  B. when the tuner ends: 100-state runs on the existing points (M series, N = 7000/10000/14000, Vmax 5000/5500) -> f100_old.log
#  C. when >= 22000 new points are in: 100-state runs on the new points (Vmax 5500/6000/6500, N = 10000/14000)    -> f100_new.log
#  D. intensities, GVPT2 over the same range, mirror plots with the experimental spectra                -> post100.log
cd ~/lico2_work/lico2_project/rc; . ./env.sh
exec >> overnight2.log 2>&1
echo "overnight2 start $(date)"
H2CM=219474.6313705
(   # ---- A
while pgrep -f "dft2.py (fit|val) --tag atz" > /dev/null; do sleep 60; done
echo "A: aTZ points finished, restarting adz7 on 64 processes $(date)"
for p in $(pgrep -f "adz7.sh"); do kill $p 2>/dev/null; done; sleep 1
for p in $(pgrep -f "dft2.py fit --tag adz7"); do kill $p 2>/dev/null; done; sleep 5; rm -f cache_adz7/*.tmp.npz
( echo "adz7 val start $(date)"; python dft2.py val --tag adz7 --basis aug-cc-pvdz --n 2500 --vmax 6500 --nproc 10; echo "adz7 val exit $? $(date)" ) >> night_adz7val.log 2>&1 &
{ echo "adz7 restart (64 procs) $(date)"; python dft2.py fit --tag adz7 --basis aug-cc-pvdz --n 40000 --vmax 6500 --nproc 64; echo "adz7 exit $? $(date)"; } >> night_adz7.log 2>&1
wait
) &
# ---- B
while pgrep -f tune100.py > /dev/null; do sleep 30; done
read EC BETA <<< $(python -c "
import json,os
P=json.load(open('results/best100.json')) if os.path.exists('results/best100.json') else dict(ec=900.0,beta=[0.8935,1.0395,0.99,0.99,1.0395,0.99])
print(P['ec'],','.join('%.4f'%b for b in P['beta']))")
echo "B: parameters ec=$EC beta=$BETA $(date)"; echo "$EC $BETA" > f100.params
run() { python run_rc.py --pes dft --ec $EC --beta $BETA --p 0 --nev 100 --workers 10 --threads 1 --vval 5000 "$@"; }
{
run --cache cache_adz --N 10000 --vmax 5500 --Ms 9000,12000,1000000 --nvalraw 3920 --tag f100_old_N10000
run --cache cache_adz --N 7000  --vmax 5500 --Ms 1000000 --nvalraw 3920 --tag f100_old_N7000
run --cache cache_adz --N 14000 --vmax 5500 --Ms 1000000 --nvalraw 3920 --tag f100_old_N14000
run --cache cache_adz --N 10000 --vmax 5000 --Ms 1000000 --nvalraw 3920 --tag f100_old_N10000_v5000
echo "B done $(date)"
} >> f100_old.log 2>&1
# ---- C
count() { python -c "
import numpy as np,sys
try:
    d=np.load('cache_adz7/$1.npz'); print(int((d['conv']&(d['E']<=6500/$H2CM)).sum()), len(d['E']))
except Exception: print(0, 0)"; }
while true; do read NF _ <<< $(count fit); read NV NVRAW <<< $(count val); [ "$NF" -ge 15000 ] && [ "$NV" -ge 800 ] && break
  HR=$(date +%-H); [ "$HR" -ge 4 ] && [ "$HR" -lt 12 ] && [ "$NF" -ge 10000 ] && [ "$NV" -ge 500 ] && break     # do not wait past 04:00
  pgrep -f "dft2.py fit --tag adz7" > /dev/null || { sleep 180; pgrep -f "dft2.py fit --tag adz7" > /dev/null || { [ "$NF" -ge 12000 ] && break; }; }; sleep 120; done
echo "C: new points in region: fit $NF, val $NV (raw $NVRAW) $(date)"
{
run --cache cache_adz7 --N 10000 --vmax 6000 --Ms 8000,12000,22000 --nvalraw $NVRAW --tag f100_new_N10000_v6000
run --cache cache_adz7 --N 10000 --vmax 5500 --Ms 16000 --nvalraw $NVRAW --tag f100_new_N10000_v5500
run --cache cache_adz7 --N 10000 --vmax 6500 --Ms 11000,22000 --nvalraw $NVRAW --tag f100_new_N10000_v6500
run --cache cache_adz7 --N 14000 --vmax 6000 --Ms 22000 --nvalraw $NVRAW --tag f100_new_N14000_v6000
run --cache cache_adz7 --N 7000  --vmax 6000 --Ms 22000 --nvalraw $NVRAW --tag f100_new_N7000_v6000
echo "C done $(date)"
} >> f100_new.log 2>&1
# ---- D
{
echo "D start $(date)"
export OPENBLAS_NUM_THREADS=16 NUMEXPR_NUM_THREADS=16 NUMEXPR_MAX_THREADS=16
python - <<'PY'
s = open('vp2_li.py').read()
for old, new in (("T = [i for i in range(NS) if S[i].sum() <= 5 and (E0[i] - ZPE0) * H2CM < 1900]", "T = [i for i in range(NS) if S[i].sum() <= 6 and (E0[i] - ZPE0) * H2CM < 2600]"),
                 ("    if nu > 1750: continue\n", "    if nu > 2400: continue\n"), ("open(f'vpt2_{TAG}.json', 'w')", "open(f'vpt2_{TAG}_hi.json', 'w')")):
    assert s.count(old) == 1, old; s = s.replace(old, new)
open('vp2_li_hi.py', 'w').write(s)
PY
python vp2_li_hi.py adz 10 > vp2_adz_hi.log 2>&1; tail -3 vp2_adz_hi.log
EXPU="expU=exp/Kafafi1983_Fig1A_LiCO2_Ar_matrix.csv:experiment, Ar matrix, as deposited (absorbance, right axis)"
EXPD="expD=exp/Kafafi1983_Fig1B_LiCO2_Ar_matrix.csv:experiment, Ar matrix, after 5 h IR exposure (absorbance, right axis)"
post() {   # cache vec tag vmax label
  [ -f "$2" ] || { echo "missing $2"; return; }
  python intens.py $1 $2 $3 $4 > intens_$3.log 2>&1; grep -E "overlaps|largest|dipole fit" intens_$3.log
  XMAX=2300 python plot_vpt2_li.py results/lico2_IR_vpt2_adz_colloc_$3 "vpt2_adz_hi.json:GVPT2, B3LYP-D3(BJ)/aug-cc-pVDZ:#1f77b4:-" "|" "results/ir_colloc_$3.json:$5:#2ca02c:-" "$EXPU" "$EXPD"
}
post cache_adz results/vec_f100_old_N10000_M16150.npz adz100_old 5500 "collocation, B3LYP-D3(BJ)/aug-cc-pVDZ points, N = 10000, M = 16150"
V=$(ls results/vec_f100_new_N10000_v6000_M*.npz 2>/dev/null | sort -t M -k 2 -n | tail -1); M=$(echo $V | sed 's/.*_M\([0-9]*\).npz/\1/')
post cache_adz7 "$V" adz100_new 6000 "collocation, B3LYP-D3(BJ)/aug-cc-pVDZ points, N = 10000, M = $M"
echo "D done $(date)"
} >> post100.log 2>&1
echo "overnight2 chain done $(date)"; wait; echo "overnight2 all done $(date)"; touch overnight2.done
