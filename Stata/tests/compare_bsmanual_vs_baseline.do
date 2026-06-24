*===============================================================================
* compare_bsmanual_vs_baseline.do
*
* Cross-check: SE from `bsmanual()` on the NEW ado vs. a hand-rolled
* manual bootstrap on the OLD baseline ado using the SAME per-rep
* cluster samples (produced by did_multiplegt_dyn_bootsample).
*
* Logic:
*  - Sampler is deterministic given the seed -- both runs see identical
*    cluster-count vectors per rep.
*  - On test_sample.dta, ac_uq_id is nested in unique_small_grid_id, so the
*    NEW ado's reghdfe absorbing ac_uq_id is a no-op vs. unit FE -- slopes
*    match what the OLD ado computes via matrix accum.
*  - Therefore per-rep point estimates match, the SDs across reps match,
*    and the bootstrap SEs should match to numerical precision.
*===============================================================================

clear all
set more off
set linesize 140
cd "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata"

local data_path "test_sample.dta"
local old_ado   "did_multiplegt_dyn_baseline.ado"
local new_ado   "did_multiplegt_dyn.ado"
local log_path  "tests/log_bsmanual_vs_baseline.log"
local csv_path  "tests/bs_samples.csv"
local csv_dta   "tests/bs_samples.dta"

local Y    "countk"
local G    "unique_small_grid_id"
local T    "monthyear"
local D    "downup_ac"
local CTRL "av_wind_speed wind_direction"
local CLUS "ac_uq_id"
local HDFE "ac_uq_id"
local NEFF 3
local NPL  3
local BREP 5
local BSED 1234

capture log close
log using "`log_path'", text replace

display as text _newline "{hline 78}"
display as text "compare_bsmanual_vs_baseline.do"
display as text "Comparing bsmanual() on NEW ado against a manual loop on OLD ado,"
display as text "both using identical cluster samples (CSV from the sampler)."
display as text "reps=`BREP'  seed=`BSED'"
display as text "{hline 78}"

*===============================================================================
* 1. Generate the CSV of per-rep cluster counts (same seed as bsmanual)
*===============================================================================
display as text _newline(2) "Step 1: generate cluster-resample CSV via did_multiplegt_dyn_bootsample"
use "`data_path'", clear
do "`new_ado'"
adopath ++ "."

did_multiplegt_dyn_bootsample, cluster(`CLUS') reps(`BREP') seed(`BSED') saving("`csv_path'")

* Read the CSV into a tempfile dataset and rename cluster_id to match the
* user's cluster variable so we can merge it onto the analysis data below.
qui import delimited using "`csv_path'", clear
qui rename cluster_id `CLUS'
qui save "`csv_dta'", replace
display as text "CSV written to `csv_path' and reshaped at `csv_dta'."

*===============================================================================
* 2. Run NEW ado with bsmanual() -- captures bsmanual bootstrap SEs.
*    bsmanual internally generates the same CSV (deterministic from seed),
*    so the per-rep cluster counts are identical to those in `csv_path'.
*===============================================================================
display as text _newline(2) "{hline 78}"
display as text "Step 2: NEW ado with bsmanual(`BREP', `BSED')"
display as text "{hline 78}"

use "`data_path'", clear
do "`new_ado'"

did_multiplegt_dyn `Y' `G' `T' `D',                ///
    effects(`NEFF') placebo(`NPL')                 ///
    controls(`CTRL')                               ///
    hdfe_controls(`HDFE')                          ///
    cluster(`CLUS')                                ///
    bsmanual(`BREP', `BSED')                       ///
    graph_off _no_updates

forvalues i = 1/`NEFF' {
    scalar SE_NEW_E_`i' = e(se_effect_`i')
}
forvalues i = 1/`NPL' {
    scalar SE_NEW_P_`i' = e(se_placebo_`i')
}
scalar SE_NEW_A = e(se_avg_total_effect)

*===============================================================================
* 3. Manual bootstrap loop on the OLD baseline ado, using the SAME CSV.
*===============================================================================
display as text _newline(2) "{hline 78}"
display as text "Step 3: manual loop on OLD baseline ado (same per-rep cluster counts)"
display as text "{hline 78}"

* Load OLD ado.
do "`old_ado'"

* Storage for per-rep OLD estimates: BREP rows x (NEFF + NPL + 1) cols.
local NPARAM = `NEFF' + `NPL' + 1
matrix OLD_reps = J(`BREP', `NPARAM', .)
local n_failed = 0

* Snapshot original data so we can reload it cleanly per rep.
use "`data_path'", clear
tempfile orig_data
qui save "`orig_data'", replace

forvalues b = 1/`BREP' {

    * (i) reload original data
    qui use "`orig_data'", clear

    * (ii) merge in this rep's per-cluster counts
    qui merge m:1 `CLUS' using "`csv_dta'", keepusing(rep_`b') keep(match master) nogen
    qui replace rep_`b' = 0 if missing(rep_`b')

    * (iii) build the single bootstrap weight (count is the only weight,
    *       since the test does not pass a user weight() to either run)
    capture drop _boot_w_XX
    qui gen double _boot_w_XX = rep_`b'

    * (iv) recursive call to OLD ado with weight = bootstrap count.
    *      No hdfe_controls (the OLD ado doesn't accept it).
    capture qui did_multiplegt_dyn `Y' `G' `T' `D',          ///
        effects(`NEFF') placebo(`NPL')                       ///
        controls(`CTRL')                                     ///
        cluster(`CLUS')                                      ///
        weight(_boot_w_XX)                                   ///
        graph_off _no_updates

    if _rc != 0 {
        local n_failed = `n_failed' + 1
        display as text "  rep " %3.0f `b' "  FAILED (rc=" _rc ")"
        continue
    }

    * (v) capture this rep's estimates
    forvalues i = 1/`NEFF' {
        capture local _v = e(Effect_`i')
        if _rc == 0 {
            matrix OLD_reps[`b', `i'] = `_v'
        }
    }
    forvalues i = 1/`NPL' {
        capture local _v = e(Placebo_`i')
        if _rc == 0 {
            matrix OLD_reps[`b', `NEFF' + `i'] = `_v'
        }
    }
    capture local _v = e(Av_tot_effect)
    if _rc == 0 {
        matrix OLD_reps[`b', `NPARAM'] = `_v'
    }
}

display as text "Manual loop done.  failed reps: " `n_failed' " / `BREP'"

*===============================================================================
* 4. Compute SE = sample SD across surviving reps via Mata.
*===============================================================================
mata: M = st_matrix("OLD_reps")
mata: keep = (rowmissing(M) :== 0)
mata: M = select(M, keep)
mata: st_local("B_eff", strofreal(rows(M)))
mata: m = mean(M)
mata: d = M :- m
mata: covM = (d' * d) / (rows(M) - 1)
mata: seM = sqrt(diagonal(covM))
mata: st_matrix("OLD_se", seM)

display as text "Effective surviving reps for OLD: " `B_eff'

forvalues i = 1/`NEFF' {
    scalar SE_OLD_E_`i' = OLD_se[`i', 1]
}
forvalues i = 1/`NPL' {
    scalar SE_OLD_P_`i' = OLD_se[`NEFF' + `i', 1]
}
scalar SE_OLD_A = OLD_se[`NPARAM', 1]

*===============================================================================
* 5. Side-by-side comparison of bsmanual SE (NEW) vs manual loop SE (OLD).
*===============================================================================
display as text _newline(2) "{hline 90}"
display as text "Side-by-side: bsmanual SE (NEW ado) vs manual-loop SE (OLD ado)"
display as text "{hline 90}"

display as text _column(4) %-14s "quantity"           ///
       _column(22) %16s "SE_NEW (bsmanual)" ///
       _column(46) %18s "SE_OLD (manual)"   ///
       _column(70) %12s "diff"

forvalues i = 1/`NEFF' {
    display as text _column(4) %-14s "Effect_`i'"                       ///
        _column(22) %16.10f SE_NEW_E_`i'                                ///
        _column(46) %18.10f SE_OLD_E_`i'                                ///
        _column(70) %12.4e SE_NEW_E_`i' - SE_OLD_E_`i'
}
forvalues i = 1/`NPL' {
    display as text _column(4) %-14s "Placebo_`i'"                      ///
        _column(22) %16.10f SE_NEW_P_`i'                                ///
        _column(46) %18.10f SE_OLD_P_`i'                                ///
        _column(70) %12.4e SE_NEW_P_`i' - SE_OLD_P_`i'
}
display as text _column(4) %-14s "Av_tot_effect"                        ///
    _column(22) %16.10f SE_NEW_A                                        ///
    _column(46) %18.10f SE_OLD_A                                        ///
    _column(70) %12.4e SE_NEW_A - SE_OLD_A

*===============================================================================
* 6. Verdict (pass if max abs SE diff < 1e-6)
*===============================================================================
scalar max_diff = 0
forvalues i = 1/`NEFF' {
    scalar _d = abs(SE_NEW_E_`i' - SE_OLD_E_`i')
    if (_d > max_diff) {
        scalar max_diff = _d
    }
}
forvalues i = 1/`NPL' {
    scalar _d = abs(SE_NEW_P_`i' - SE_OLD_P_`i')
    if (_d > max_diff) {
        scalar max_diff = _d
    }
}
scalar _d = abs(SE_NEW_A - SE_OLD_A)
if (_d > max_diff) {
    scalar max_diff = _d
}

display as text _newline "max abs SE diff = " %12.4e max_diff
if (max_diff < 1e-6) {
    display as result "bsmanual matches manual baseline-ado bootstrap (max diff < 1e-6)."
}
else {
    display as error "bsmanual SE differs from manual baseline bootstrap (max diff = " max_diff ")."
}

log close
