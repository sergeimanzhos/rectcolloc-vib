# LiCO2 max-residual tuning at N = 10000 -- state when stopped (17 Nov 2026, ~09:10 JST, by user request)
Script: tune_li_max10k.py (restartable: every finished trial is cached in results/tmax_N10000_M12000_*.json; the script
replays cached trials and continues with the first missing one). Log: tune_li_max10k.log.
Restart command (from this directory):  timeout 7150 python3 tune_li_max10k.py > tune_li_max10k.out 2>&1
Position when stopped: stage 3 (1% descent), sweep 7, beta[6] x0.99 trial was running and is lost (will be recomputed).
Current best (sweep 7, after beta[1..5]): Ec=1200, Vmax=5000, p=0, beta=(0.932065, 0.99, 1.0, 1.018814, 1.0, 1.0), J=113.9164 cm-1.
Sweep 7 has accepted beta[1] -> 0.932065, so a sweep 8 follows after beta[6].
