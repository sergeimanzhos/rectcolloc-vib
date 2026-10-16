#!/bin/bash
# aTZ, 100 states (2026-10-10), self-contained.  Master log: rc/atz_chain.log
#  A. second, independent aTZ point set at the wider-box aDZ points (cache_atz7): first 18000 points with V(aDZ) <= 5800, then its
#     validation points, then the remaining points so the CPUs keep working                        -> night_atz7.log
#  B. on the existing aTZ points (cache_atz): GVPT2 stencil + GVPT2, collocation N = 10000 / 14000 (100 states), intensities, plots -> atz100.log
#  D. when the first 18000 new points and the validation set are in: collocation on cache_atz7, convergence ranges, plots          -> atz100.log
cd ~/lico2_work/lico2_project/rc; . ./env.sh
exec >> atz_chain.log 2>&1
echo "atz_chain start $(date)"
BAS=aug-cc-pvtz
mkdir -p cache_atz7; cp -n cache_atz/E0.npy cache_atz/dm0.npy cache_atz/dip0.npy cache_atz/qmin.npy cache_atz7/; cp -n cache_adz7/box.npz cache_atz7/
for N in 10000 14000; do cp -n results/centres_dft_cache_adz_ec900_N$N.npy results/centres_dft_cache_atz7_ec900_N$N.npy; cp -n results/centres_dft_cache_adz_ec900_N$N.npy results/centres_dft_cache_atz_ec900_N$N.npy; done
rm -f atzB.done atzA1.done atzD.done
D2() { python dft2.py $1 --tag atz7 --basis $BAS --src adz7 --vmax 6000 --srcmax 5800 "${@:2}"; }
(   # ---- A
echo "A start $(date)"
L=3500; D2 fit --limit $L --nproc 44
while [ ! -f atzB.done ] && [ $L -lt 17000 ]; do L=$((L+600)); D2 fit --limit $L --nproc 44; done     # 48 processes while B needs memory
D2 fit --limit 18000 --nproc 66
D2 val --nproc 66
touch atzA1.done; echo "A1 (18000 fit + val) done $(date)"
D2 fit --nproc 44
echo "A exit $? $(date)"
) >> night_atz7.log 2>&1 &
read EC BETA < f100.params
run() { python run_rc.py --pes dft --ec $EC --beta $BETA --p 0 --nev 100 --workers 8 --threads 1 --vval 5000 "$@"; }
LVL="B3LYP-D3(BJ)/aug-cc-pVTZ"
EXPU="expU=exp/Kafafi1983_Fig1A_LiCO2_Ar_matrix.csv:experiment, Ar matrix, as deposited (absorbance, right axis)"
EXPD="expD=exp/Kafafi1983_Fig1B_LiCO2_Ar_matrix.csv:experiment, Ar matrix, after 5 h IR exposure (absorbance, right axis)"
plots() {   # tag  M  band-legend  range-runs...
  T=$1; M=$2; BL=$3; shift 3
  python mkrange.py results/ir_colloc_$T.json results/ir_colloc_${T}_range.json "$@"
  XMIN=200 XMAX=2300 FIGW=9.69 FIGH=6.27 FONT=16 python plot_zoom.py results/lico2_IR_${T}_A4 "vpt2_atz_hi.json:GVPT2, $LVL:#1f77b4" "|" "results/ir_colloc_${T}_range.json:collocation, $LVL points, N = 14000, M = $M:#2ca02c" "band=$BL"
  XMAX=2300 python plot_vpt2_li.py results/lico2_IR_vpt2_atz_colloc_$T "vpt2_atz_hi.json:GVPT2, $LVL:#1f77b4:-" "|" "results/ir_colloc_$T.json:collocation, $LVL points, N = 14000, M = $M:#2ca02c:-" "$EXPU" "$EXPD"
}
{   # ---- B
echo "B start $(date)"
[ -f vp_atz.npz ] || python vp_dft_li.py $BAS atz 20
python vp_fit.py vp_atz.npz; python vp2_li_hi.py atz 10 > vp2_atz_hi.log 2>&1; tail -2 vp2_atz_hi.log
run --cache cache_atz --N 10000 --vmax 5500 --Ms 9000,1000000 --tag f100_atz_N10000
run --cache cache_atz --N 14000 --vmax 5500 --Ms 1000000 --tag f100_atz_N14000
run --cache cache_atz --N 10000 --vmax 5000 --Ms 1000000 --tag f100_atz_N10000_v5000
export OPENBLAS_NUM_THREADS=16 NUMEXPR_NUM_THREADS=16 NUMEXPR_MAX_THREADS=16
V=$(ls results/vec_f100_atz_N14000_M*.npz | tail -1); M=$(echo $V | sed 's/.*_M\([0-9]*\).npz/\1/')
python intens.py cache_atz $V atz100_N14000 5500 > intens_atz100_N14000.log 2>&1; grep -E "complex pairs|overlaps|largest|dipole fit|Traceback|Error" intens_atz100_N14000.log
plots atz100_N14000 $M "range over N = 10000 and 14000" f100_atz_N14000 f100_atz_N10000
unset OPENBLAS_NUM_THREADS NUMEXPR_NUM_THREADS NUMEXPR_MAX_THREADS; . ./env.sh
echo "B done $(date)"
} >> atz100.log 2>&1
touch atzB.done
while [ ! -f atzA1.done ]; do sleep 120; done
{   # ---- D
echo "D start $(date)"
NVRAW=$(python -c "import numpy as np; print(len(np.load('cache_atz7/val.npz')['E']))")
run --cache cache_atz7 --N 10000 --vmax 5500 --Ms 16000 --nvalraw $NVRAW --tag f100_atz7_N10000
run --cache cache_atz7 --N 14000 --vmax 5500 --Ms 16000 --nvalraw $NVRAW --tag f100_atz7_N14000
export OPENBLAS_NUM_THREADS=16 NUMEXPR_NUM_THREADS=16 NUMEXPR_MAX_THREADS=16
V=$(ls results/vec_f100_atz_N14000_M*.npz | tail -1); M=$(echo $V | sed 's/.*_M\([0-9]*\).npz/\1/')
plots atz100_N14000 $M "range over N = 10000, 14000 and two point sets" f100_atz_N14000 f100_atz_N10000 f100_atz7_N14000 f100_atz7_N10000
echo "D done $(date)"
} >> atz100.log 2>&1
touch atzD.done; echo "chain B, D done $(date)"; wait; echo "atz_chain all done $(date)"; touch atz_chain.done
