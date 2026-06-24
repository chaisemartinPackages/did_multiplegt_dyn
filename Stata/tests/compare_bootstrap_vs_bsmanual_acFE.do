*===============================================================================
* compare_bootstrap_vs_bsmanual_acFE.do
*
* Compares the OLD baseline ado's built-in bootstrap() against the NEW ado's
* bsmanual() under matched-seed cluster resampling. Both should produce
* IDENTICAL effects, placebos, average total effect, coefficients, and SEs.
*
* FE structure tested: per-ac fixed effects.
*   OLD path: controls(av_wind_speed wind_direction + ac dummies)
*             with ac == last ac level dropped as reference.
*             The ac dummies are time-invariant -> no-op in the first-diff
*             residualization, but included for explicit symmetry with NEW.
*   NEW path: controls(av_wind_speed wind_direction) hdfe_controls(i.ac_uq_id).
*             ac_uq_id is nested in unique_small_grid_id, so absorbing it on
*             top of unit FE is a no-op in the level fit -> equivalent
*             residualization to the OLD path.
*
* RNG mechanics:
*   OLD: bootstrap(5, 1234) -> Stata's built-in cluster bootstrap (uses
*        `bsample, cluster()` internally with set seed 1234).
*   NEW: bsmanual(5, 1234)  -> our sampler also uses `bsample, cluster()`
*        with set seed 1234, drawing the same per-rep cluster sequence.
*===============================================================================

clear all
set more off
set linesize 140
set maxvar 50000
cd "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata"

local data_path "test_sample.dta"
local old_ado   "did_multiplegt_dyn_baseline.ado"
local new_ado   "did_multiplegt_dyn.ado"
local log_path  "tests/log_compare_bootstrap_vs_bsmanual_acFE.log"

local Y    "countk"
local G    "unique_small_grid_id"
local T    "monthyear"
local D    "downup_ac"
local CTRL_BASE "av_wind_speed wind_direction"
local CLUS "ac_uq_id"
local HDFE "i.ac_uq_id"
local NEFF 3
local NPL  3
local BREP 2
local BSED 1234
local TOL  1e-10

capture log close
log using "`log_path'", text replace

display as text _newline "{hline 78}"
display as text "compare_bootstrap_vs_bsmanual_acFE.do"
display as text "OLD bootstrap(`BREP', `BSED') vs NEW bsmanual(`BREP', `BSED')"
display as text "FE: per-ac fixed effects (OLD via dummies, NEW via hdfe_controls)"
display as text "Tolerance: `TOL'"
display as text "{hline 78}"

use "`data_path'", clear
quietly count
display as text "Observations: " r(N)

*===============================================================================
* 1. Build ac dummies for OLD path (drop last ac as reference).
*===============================================================================
quietly levelsof `CLUS', local(acs_kept)
local n_acs : word count `acs_kept'
local ref_ac : word `n_acs' of `acs_kept'
display as text _newline "ac levels: `acs_kept'"
display as text "Reference ac (excluded from dummies): " `ref_ac'

local ac_dummies ""
foreach a of local acs_kept {
    if `a' == `ref_ac' continue
    capture drop _acFE_`a'_XX
    qui gen byte _acFE_`a'_XX = (`CLUS' == `a')
    local ac_dummies "`ac_dummies' _acFE_`a'_XX"
}
display as text "OLD path ac dummies: `ac_dummies'"

* Save the data with the ac dummies prepared so both runs see the same data.
tempfile data_with_dummies
qui save "`data_with_dummies'", replace

*===============================================================================
* 2. OLD path: did_multiplegt_dyn_baseline.ado with bootstrap(`BREP', `BSED')
*===============================================================================
display as text _newline(2) "{hline 78}"
display as text "OLD PATH: baseline ado + bootstrap(`BREP', `BSED')"
display as text "{hline 78}"

use "`data_with_dummies'", clear
do "`old_ado'"
adopath ++ "."

did_multiplegt_dyn `Y' `G' `T' `D',                       ///
    effects(`NEFF') placebo(`NPL')                        ///
    controls(`CTRL_BASE' `ac_dummies')                    ///
    cluster(`CLUS')                                       ///
    bootstrap(`BREP', `BSED')                             ///
    graph_off _no_updates save_sample
est store reg1
matrix b_OLD = e(b)
matrix V_OLD = e(V)

forvalues i = 1/`NEFF' {
    scalar E_OLD_`i'    = e(Effect_`i')
    scalar SE_E_OLD_`i' = e(se_effect_`i')
}
forvalues i = 1/`NPL' {
    scalar P_OLD_`i'    = e(Placebo_`i')
    scalar SE_P_OLD_`i' = e(se_placebo_`i')
}
scalar A_OLD    = e(Av_tot_effect)
scalar SE_A_OLD = e(se_avg_total_effect)
scalar PJE_OLD  = e(p_jointeffects)
scalar PJP_OLD  = e(p_jointplacebo)

*===============================================================================
* 3. NEW path: did_multiplegt_dyn.ado with bsmanual(`BREP', `BSED')
*===============================================================================
display as text _newline(2) "{hline 78}"
display as text "NEW PATH: new ado + bsmanual(`BREP', `BSED')"
display as text "{hline 78}"

use "`data_with_dummies'", clear
do "`new_ado'"

did_multiplegt_dyn `Y' `G' `T' `D',                       ///
    effects(`NEFF') placebo(`NPL')                        ///
    controls(`CTRL_BASE')                                 ///
    hdfe_controls(`HDFE')                                 ///
    cluster(`CLUS')                                       ///
    bsmanual(`BREP', `BSED')                              ///
    graph_off _no_updates
est store reg2
matrix b_NEW = e(b)
matrix V_NEW = e(V)

forvalues i = 1/`NEFF' {
    scalar E_NEW_`i'    = e(Effect_`i')
    scalar SE_E_NEW_`i' = e(se_effect_`i')
}
forvalues i = 1/`NPL' {
    scalar P_NEW_`i'    = e(Placebo_`i')
    scalar SE_P_NEW_`i' = e(se_placebo_`i')
}
scalar A_NEW    = e(Av_tot_effect)
scalar SE_A_NEW = e(se_avg_total_effect)
scalar PJE_NEW  = e(p_jointeffects)
scalar PJP_NEW  = e(p_jointplacebo)

*===============================================================================
* 4. Side-by-side comparison
*===============================================================================
display as text _newline(2) "{hline 120}"
display as text "Side-by-side: OLD bootstrap() vs NEW bsmanual() under matched seed"
display as text "{hline 120}"

display as text _column(4) %-16s "quantity" ///
       _column(22) %16s "OLD estimate" ///
       _column(40) %16s "NEW estimate" ///
       _column(58) %12s "est diff"     ///
       _column(74) %16s "OLD SE"       ///
       _column(92) %16s "NEW SE"       ///
       _column(110) %12s "SE diff"

forvalues i = 1/`NEFF' {
    display as text _column(4) %-16s "Effect_`i'"                          ///
        _column(22) %16.10f E_OLD_`i'                                      ///
        _column(40) %16.10f E_NEW_`i'                                      ///
        _column(58) %12.2e E_OLD_`i' - E_NEW_`i'                           ///
        _column(74) %16.10f SE_E_OLD_`i'                                   ///
        _column(92) %16.10f SE_E_NEW_`i'                                   ///
        _column(110) %12.2e SE_E_OLD_`i' - SE_E_NEW_`i'
}
forvalues i = 1/`NPL' {
    display as text _column(4) %-16s "Placebo_`i'"                         ///
        _column(22) %16.10f P_OLD_`i'                                      ///
        _column(40) %16.10f P_NEW_`i'                                      ///
        _column(58) %12.2e P_OLD_`i' - P_NEW_`i'                           ///
        _column(74) %16.10f SE_P_OLD_`i'                                   ///
        _column(92) %16.10f SE_P_NEW_`i'                                   ///
        _column(110) %12.2e SE_P_OLD_`i' - SE_P_NEW_`i'
}
display as text _column(4) %-16s "Av_tot_effect"                           ///
    _column(22) %16.10f A_OLD                                              ///
    _column(40) %16.10f A_NEW                                              ///
    _column(58) %12.2e A_OLD - A_NEW                                       ///
    _column(74) %16.10f SE_A_OLD                                           ///
    _column(92) %16.10f SE_A_NEW                                           ///
    _column(110) %12.2e SE_A_OLD - SE_A_NEW

display as text _newline "Joint significance p-values (bootstrap covariance):"
display as text _column(4) %-20s "p_jointeffects:" ///
    _column(28) %16.10f PJE_OLD ///
    _column(46) %16.10f PJE_NEW ///
    _column(64) %12.2e PJE_OLD - PJE_NEW
display as text _column(4) %-20s "p_jointplacebo:" ///
    _column(28) %16.10f PJP_OLD ///
    _column(46) %16.10f PJP_NEW ///
    _column(64) %12.2e PJP_OLD - PJP_NEW

display as text _newline "mreldif(b_OLD, b_NEW) = " %12.4e mreldif(b_OLD, b_NEW)
display as text "mreldif(V_OLD, V_NEW) = " %12.4e mreldif(V_OLD, V_NEW)

*===============================================================================
* 5. Verdict (pass if max abs diff < TOL across estimates, SEs, p-values)
*===============================================================================
scalar max_est_diff = 0
scalar max_se_diff  = 0
scalar max_p_diff   = 0

forvalues i = 1/`NEFF' {
    scalar _d = abs(E_OLD_`i' - E_NEW_`i')
    if (_d > max_est_diff) {
        scalar max_est_diff = _d
    }
    scalar _d = abs(SE_E_OLD_`i' - SE_E_NEW_`i')
    if (_d > max_se_diff) {
        scalar max_se_diff = _d
    }
}
forvalues i = 1/`NPL' {
    scalar _d = abs(P_OLD_`i' - P_NEW_`i')
    if (_d > max_est_diff) {
        scalar max_est_diff = _d
    }
    scalar _d = abs(SE_P_OLD_`i' - SE_P_NEW_`i')
    if (_d > max_se_diff) {
        scalar max_se_diff = _d
    }
}
scalar _d = abs(A_OLD - A_NEW)
if (_d > max_est_diff) {
    scalar max_est_diff = _d
}
scalar _d = abs(SE_A_OLD - SE_A_NEW)
if (_d > max_se_diff) {
    scalar max_se_diff = _d
}
scalar _d = abs(PJE_OLD - PJE_NEW)
if (_d > max_p_diff) {
    scalar max_p_diff = _d
}
scalar _d = abs(PJP_OLD - PJP_NEW)
if (_d > max_p_diff) {
    scalar max_p_diff = _d
}

display as text _newline "max abs diff (point estimates) = " %12.4e max_est_diff
display as text "max abs diff (bootstrap SEs)   = " %12.4e max_se_diff
display as text "max abs diff (joint p-values)  = " %12.4e max_p_diff

if (max_est_diff < `TOL') & (max_se_diff < `TOL') & (max_p_diff < `TOL') {
    display as result _newline "OLD bootstrap() == NEW bsmanual()  --  all diffs under tolerance `TOL'."
}
else {
    display as error _newline "DIVERGENCE between OLD bootstrap() and NEW bsmanual()."
    display as error "est=" max_est_diff "  SE=" max_se_diff "  p=" max_p_diff "  tol=`TOL'"
}

log close
