# LiCO2 max-residual redo — handoff for a new Claude session (state 17 Nov 2026, 18:55 JST)

Read this whole file before doing anything. It replaces the conversation history.

## 0. Standing rules from Sergei (apply to everything)
- No improvisation: do exactly what is asked and report exactly what and how was done ("always state exactly how the calculation was done").
- "Can you ...?" is a capability question: answer it, do not start calculations until told.
- Obey "stop" immediately.
- Do not update Word/Excel or commit unless asked. Commits go to the **main** branch of github.com/sergeimanzhos/rectcolloc-vib.
- Word edits: numbers in tables **untracked**; text changes **tracked** (author "Claude").
- When results change a figure in the manuscript, regenerate and replace the figure without being asked.
- Report spectrum/table updates as **changes from the manuscript values**, not only absolute values.
- Large fonts and lean figures.

## 1. What was asked (chronology)
1. "Redo the calculations for the section LiCO2 with max residual. Tell me if they deviate a lot from before and if you think that term
   helping stability with M is needed" (the term = Δ(M, M2), median level change between two values of M, from the old objective J = r̄ + Δ).
   Ab initio data are NOT redone (DFT point sets data/cache_adz, cache_atz, cache_adz7, cache_atz7 are reused).
2. Tuning N: "do option 2 first" (N = 5000), scope "Ec, Vmax, p, then widths", widths with 5% then 1% steps as for Cl2CO.
3. N = 5000 tuning chose Ec = 700; Sergei did not trust it -> aborted. Ec scan at N = 10000 (M = 12000) over 700..1500 -> "1200 be it".
4. Remaining tuning at N = 10000 for all (Vmax, p, beta0, 5% and 1% width descent): "keep going" until a sweep accepts nothing.
5. 17 Nov 09:08: "1) Stop now keeping all information needed to restart 2) update the manuscript with the current best results (Word and
   Excel attached; tables without tracking, text with tracking) 3) continue tuning".
   Convergence ranges: "do same as before but do not include them on the graph. include them in the Excel instead".
6. 09:38: complete handoff package + update the GitHub repository (this file and the commit are the result).

## 2. Definitions used
- Collocation: run_rc.py (exact Cartesian FD KEO, Gaussian basis, energy + gradient rows, weights (1 - V/Vmax)^p, shift-invert, nev = 100).
- J = maximum over the 100 levels of the unweighted validation residual r_i = ||(H - E_i)psi_i|| / ||psi_i|| on a fixed validation set
  (first 3920 validation records of data/cache_adz with V <= 5000 cm-1; collo.residual).
- Ec = centre cutoff (V'max): centres are R2 points with degree-4 polynomial surrogate V <= Ec (and quartic FF <= 2.5 Ec), run_rc.py line 67.
- surr.py was MISSING from the repository and was reconstructed (degree-4 polynomial fit, as in dft2.py). Consequence: centre sets are not
  bit-identical to the original runs (check chk_final_adz_N5000: ZPE -0.007, level |change| median 0.08, max 1.0 cm-1).
- Manuscript says p = 1 but the original runs (f100_*) used p = 0.

## 3. Results so far
N = 5000, M = 6000 (tune_li_max.py, log tune_li_max.log): beta0 = 1.00, Ec 700 (from 700/900/1200/1500), Vmax 5500, p 0 -> abandoned.
Ec scan N = 10000, M = 12000, Vmax 5500, p 0, beta 1 (scan_ec_N10000.py, log scan_ec_N10000.log), J (max residual):
  700 161.4 | 800 151.9 | 900 141.8 | 1000 137.8 | 1100 129.8 | **1200 126.8** | 1300 127.4 | 1500 137.4.
N = 10000 tuning (tune_li_max10k.py, log tune_li_max10k.log), Ec = 1200 fixed:
  Vmax 5000 (J 118.05) vs 5500 (126.79) -> 5000 accepted; p: 0 118.05, 1 165.72, 3 244.26 -> 0; beta0: 0.95 141.7, 1.00 118.0, 1.05 135.3.
  5%: beta4 -> 1.05 (116.35); second 5% sweep accepted nothing.
  1% sweeps: beta1 x0.99 in every sweep (1 -> 0.932065, gains 0.19, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03); beta2 -> 0.99 (sweep 1);
  beta4 1.05 -> 1.0395 -> 1.0291 -> 1.0188 (sweeps 1-3). Stopped in sweep 7 after beta5 (beta6 x0.99 trial lost).
  **Current best: Ec 1200, Vmax 5000, p 0, beta = (0.932065, 0.99, 1.0, 1.018814, 1.0, 1.0), J = 113.9164.**
Level changes per 1% sweep are hundredths of cm-1 (levels 1-39) to a few tenths (upper levels); see the conversation summary in section 6.
Fundamentals at N = 10000 with the 5%-stage widths vs the manuscript (N = 14000, old parameters): within 0.1 cm-1 except the nu6 pair
(-0.7, -1.1 cm-1); mean |change| levels 1-39 0.37, 40-69 1.07, 70-99 1.09 cm-1 (mixes N and parameter changes).

## 4a. UPDATE 17 Nov 2026 18:55 JST — where things stand now (read this before section 4)
- In the cloud container, fmax_adz_N14000 was killed by the OOM killer (6.1 GB of 7 GB). The cloud container was also reclaimed twice
  while long jobs were running (overnight, and ~14:18), killing the running job (finished, cached results survive on disk).
  => No long calculations in the cloud container.
- The session is linked to Sergei's Windows workstation (desktop-dj0n4p7); the Linux VM reachable through device_bash has only 4 GB / 2 cores
  (not usable for the runs; it cannot reach WSL). Folder granted: D:\Claude_use (Sergei's designated work folder).
- ALL final runs (steps 1-10 of final_runs.py: 4 x N = 10000, 4 x N = 14000, 2 intensities) are being run by **Claude Code on the
  workstation in WSL**, from D:\Claude_use\lico2_final_cc_20261117.tar.gz unpacked to D:\Claude_use\lico2_final, following
  lico2_final/CLAUDE_CODE_TASK.md. Command: N14000=1 WINDOW=1e9 THREADS=<cores> WORKERS=<k> python3 final_runs.py > final_runs_ws.out.
  Sergei (18:50): Claude Code expects to finish on the morning of 18 Nov. He will say when it is done.
- final_runs.py now takes THREADS, WORKERS, WINDOW, N14000 from the environment (defaults: 2, 1, 7000 s, N = 14000 steps refused).
- NEXT, when Sergei says Claude Code is done:
  1. Read D:\Claude_use\lico2_final\redo_lico2\final_runs_ws.out; stage results/fmax_*.json, intens_fmax_*.json/.csv,
     ir_colloc_fmax_*.json, levels_fmax_*.csv into the cloud container (device_stage_files) and copy them into redo_lico2/results.
  2. Steps B-F of section 4 (ranges with mkrange.py over the 4 runs per basis; Excel; Figure 6 without bars; Word Table 3 untracked,
     text tracked; report with a table of CHANGES from the manuscript; answer on the Delta(M) term).
  3. Step 3 of Sergei's 09:08 request, "continue tuning": proposed to hand the tuning to Claude Code on the workstation as well
     (package + one-paragraph instruction; resumes from the 115 cached trials) — confirm with Sergei before doing it.
- Promised (not yet done): prepare and test the Excel/Word update scripts on the manuscript data before the results arrive.

## 4. In progress / pending (do in this order)
A. Final runs: final_runs.py (restartable, skips finished steps; in the cloud session launch as
   `timeout 7150 python3 final_runs.py > final_runs.out 2>&1` and relaunch after each exit code 3). Steps: fmax_adz_N14000 (was started by
   hand 09:11 17 Nov), fmax_atz_N14000, intensities (intens.py, VMAX 5000) for both, then N = 10000 set 1 and N = 10000/14000 on set 2.
   One step at a time (memory: N = 10000 needs 3.2 GB; 7 GB container). N = 10000 run ~33 min with 2 threads; N = 14000 unknown.
B. Ranges: mkrange.py ir_colloc_fmax_<b>_N14000.json -> _range.json over fmax_<b>_N14000, fmax_<b>_N10000, fmax_<b>7_N14000, fmax_<b>7_N10000.
C. Excel (input: the workbook Sergei sent, LiCO2_levels_aDZ_aTZ.xlsx): update with gen.py (upd_li/xl/gen.py; check_old.py proves it
   regenerates the existing workbook exactly). Follow the Cl2CO update_xlsx.py pattern (rectcolloc-vib commit eb7e134 /
   upd_opt4/xl/update_xlsx.py): Collocation aDZ/aTZ sheets (E, range low/high, sym, label, state, weight, A fitted, A DFT points, residual,
   flag), title rows (Vmax 5000, widths, Ec), ZPE, Levels sheet columns B-D / H-J, Matched-by-label regenerated (GVPT2-only rows: keep the
   existing selection), Notes (collocation line: N, Ec 1200, widths, Vmax 5000, M values; range definition; convergence line; source line).
   Formulas in the matched sheet: copy the existing pattern (harmonic, differences).
D. Figure 6 (word/media/image11.png linear, image12.png log; 1904x505 and 1880x505 px): plot_summary.py with MPLCONFIGDIR containing
   matplotlibrc "font.family: DejaVu Sans", COMPACT=1 XMIN=250 XMAX=2200 FIGW=9.75 FIGH=2.75 FONT=13; upper: new ir_colloc aDZ (C0 '-'),
   aTZ (C3 '-.'), labels "collocation, B3LYP-D3(BJ)/aug-cc-pVDZ, N = 14000, M = <M>"; lower: vpt2_adz_hi.json, vpt2_atz_hi.json
   "GVPT2, B3LYP-D3(BJ)/aug-cc-pV?Z"; NO band (use the json without Elo/Ehi, no band= argument).
E. Word (input: the docx Sergei sent): Table 3 collocation columns (DZ, TZ: "E (A)") untracked; text tracked:
   parameters paragraph ("p = 1 was found to be best", Vmax 5500 and V'max 900, widths 0.8935...), convergence sentences (first 25 levels,
   61/78 counts — not reproducible from the old range files, use range width over the 4 runs and tell Sergei), method comparison
   statistics (clearly assigned = colloc HO weight >= 0.4, matched k-th by label, E <= 2200: old data give aDZ 64 states MAE 5.2 RMS 7.1,
   aTZ 62 MAE 7.5 RMS 10.8; aTZ-aDZ: both weights >= 0.4 -> MAE 10.4 RMS 12.6 max 25.5), nu6 split numbers (435:147), 2nu3/2nu1 deviations,
   nu4/nu5/nu6 vs Ar matrix, Figure 6 caption (remove the convergence-range sentence) and the text sentence about pale convergence bars.
   Use upd_li/tools/track.py (replace_span, ins_before, ins_end_of_para).
F. Report to Sergei: exactly how everything was done; table of CHANGES from the manuscript values; whether the results deviate a lot;
   answer on the Delta(M) term (evidence: level changes with M, e.g. earlier M 12000 vs 16150 at N = 10000: median <= 0.2, max 1.4 cm-1);
   caveats (surr.py reconstructed, p = 1 in text vs p = 0 used, old 61/78 counts).
G. Then restart the tuning: `timeout 7150 python3 tune_li_max10k.py > tune_li_max10k.out 2>&1` (relaunch on each exit code 3).
   Sergei asked "is this ever gonna end??" — offered options (stop now / 0.05 cm-1 threshold / direct beta1 trials 0.90, 0.85, 0.80);
   he chose to stop, update the manuscript, then continue tuning as specified.

## 5. Environment notes (cloud session)
- Background commands are killed after ~2 h: drivers exit with code 3 before a step that would not fit; relaunch by hand.
- One N = 10000 job at a time (two at once were OOM-killed). Centre generation for a new (cache, Ec, N) takes 20-45 min and is cached
  in results/centres_dft_<cache>_ec<Ec>_N<N>.npy.
- coords.py needs ../vpt/ff_lico2.npz: `ln -s <dir>/vpt_ff vpt` next to the lico2 working directory.
- Restart package: lico2_restart_20261117.tar.gz (sent to Sergei; contains redo_lico2/ and upd_li/ incl. his Word/Excel inputs).

## 6. Files
Scripts: tune_li_max.py, scan_ec_N10000.py, tune_li_max10k.py, final_runs.py, surr.py (reconstructed), mkrange.py, intens.py, plot_summary.py.
Logs: tune_li_max.log, scan_ec_N10000.log, tune_li_max10k.log. State: TUNING_STATE.md. Trials: results/tmax_*.json (+ levels_tmax_*.csv).
Excel/Word tools: in the repository lico2/manuscript_tools/ (xlsx_gen.py, xlsx_check_old.py, track.py, show.py); in the restart package
as upd_li/xl/gen.py, check_old.py and upd_li/tools/track.py, show.py, with Sergei's input files upd_li/in/ms.docx and li.xlsx (17 Nov 2026).
Repository: results/LiCO2_levels_aDZ_aTZ.xlsx was replaced on 17 Nov 2026 by the version Sergei sent (the manuscript version, without the
free-CO2 column); the Word manuscript is not in the repository (package only).
To work from the repository instead of the package: copy lico2/ to a working directory, `ln -s <workdir>/vpt_ff <workdir>/../vpt`,
and run the scripts there (they write to ./results).
