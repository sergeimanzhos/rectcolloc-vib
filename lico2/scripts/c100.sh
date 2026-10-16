#!/bin/bash
# 100-state aDZ runs on the new (wider box, V <= 6500) points, then intensities / GVPT2 / plots.  Logs: f100_new.log, post100.log
# The basis centres are the ones used with the old point set (same surface, same Ec), so the two point sets share one basis.
cd ~/lico2_work/lico2_project/rc; . ./env.sh
read EC BETA < f100.params
NVRAW=$(python -c "import numpy as np; print(len(np.load('cache_adz7/val.npz')['E']))")
run() { python run_rc.py --pes dft --ec $EC --beta $BETA --p 0 --nev 100 --workers 12 --threads 1 --vval 5000 "$@"; }
{
echo "C start $(date), val raw $NVRAW"
run --cache cache_adz7 --N 10000 --vmax 6000 --Ms 8000,1000000 --nvalraw $NVRAW --tag f100_new_N10000_v6000
run --cache cache_adz7 --N 10000 --vmax 6500 --Ms 1000000 --nvalraw $NVRAW --tag f100_new_N10000_v6500
run --cache cache_adz7 --N 10000 --vmax 5500 --Ms 1000000 --nvalraw $NVRAW --tag f100_new_N10000_v5500
run --cache cache_adz7 --N 14000 --vmax 6000 --Ms 1000000 --nvalraw $NVRAW --tag f100_new_N14000_v6000
echo "C done $(date)"
} >> f100_new.log 2>&1
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
[ -f vpt2_adz_hi.json ] || python vp2_li_hi.py adz 10 > vp2_adz_hi.log 2>&1; tail -2 vp2_adz_hi.log
EXPU="expU=exp/Kafafi1983_Fig1A_LiCO2_Ar_matrix.csv:experiment, Ar matrix, as deposited (absorbance, right axis)"
EXPD="expD=exp/Kafafi1983_Fig1B_LiCO2_Ar_matrix.csv:experiment, Ar matrix, after 5 h IR exposure (absorbance, right axis)"
post() {   # cache vec tag vmax label
  [ -f "$2" ] || { echo "missing $2"; return; }
  python intens.py $1 $2 $3 $4 > intens_$3.log 2>&1; grep -E "overlaps|largest|dipole fit|Traceback|Error" intens_$3.log
  XMAX=2300 python plot_vpt2_li.py results/lico2_IR_vpt2_adz_colloc_$3 "vpt2_adz_hi.json:GVPT2, B3LYP-D3(BJ)/aug-cc-pVDZ:#1f77b4:-" "|" "results/ir_colloc_$3.json:$5:#2ca02c:-" "$EXPU" "$EXPD"
}
post cache_adz results/vec_f100_old_N14000_M16150.npz adz100_old_N14000 5500 "collocation, B3LYP-D3(BJ)/aug-cc-pVDZ points, N = 14000, M = 16150"
V=$(ls results/vec_f100_new_N14000_v6000_M*.npz 2>/dev/null | tail -1); M=$(echo $V | sed 's/.*_M\([0-9]*\).npz/\1/')
post cache_adz7 "$V" adz100_new_N14000 6000 "collocation, B3LYP-D3(BJ)/aug-cc-pVDZ points, N = 14000, M = $M"
echo "D done $(date)"
} >> post100.log 2>&1
touch c100.done
