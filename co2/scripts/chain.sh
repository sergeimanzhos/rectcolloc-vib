cd ~/lico2_work/lico2_project/co2_atz; . ../rc/env.sh
export OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 OMP_NUM_THREADS=1
while pgrep -f "co2_grad.py dftpts" > /dev/null; do sleep 15; done
[ -f pts_grad.npz ] && python co2_lev.py dftgrad 20 > lev_dftgrad.log 2>&1
while pgrep -f "co2_lev.py ff" > /dev/null; do sleep 15; done
python lev_table.py 4350 co2fit DFT > lev_table.log 2>&1; python lev_table.py 4350 co2cyl QFF > lev_table_qff.log 2>&1
touch chain.done
