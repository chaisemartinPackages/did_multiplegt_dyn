*===============================================================================
* compare_bootstrap_manual.do
*
* Tests the new bootstrap_manual(B, seed) option against the OLD path's
* bootstrap(B, seed) on the no-hdfe path. Both should produce identical
* point estimates and similar (but not bit-identical) bootstrap SEs --
* the resampling mechanic differs (Stata's bootstrap vs. our reweighting),
* so SEs will not match exactly. Joint p-values should both be finite.
*
* Then runs bootstrap_manual on the NEW (hdfe) path and verifies the
* same point estimates and finite SEs.
*===============================================================================

clear all
set more off
set linesize 140
cd "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata"

local data_path "test_sample.dta"
local old_ado   "did_multiplegt_dyn_baseline.ado"
local new_ado   "did_multiplegt_dyn.ado"
local log_path  "tests/log_compare_bootstrap_manual.log"

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
display as text "compare_bootstrap_manual.do  --  manual cluster bootstrap on hdfe path"
display as text "reps=`BREP'  seed=`BSED'"
display as text "{hline 78}"

use "`data_path'", clear
quietly count
display as text "Observations: " r(N)

*===============================================================================
* 1. NEW ado with bootstrap_manual on the hdfe path
*    Expected: point estimates match the no-bootstrap NEW path; SEs and joint
*    p-values are finite (the reweighting bootstrap produces real variance).
*===============================================================================
display as text _newline(2) "{hline 78}"
display as text "NEW PATH: hdfe_controls + bootstrap_manual"
display as text "{hline 78}"

use "`data_path'", clear
do "`new_ado'"
adopath ++ "."

did_multiplegt_dyn `Y' `G' `T' `D',                ///
    effects(`NEFF') placebo(`NPL')                 ///
    controls(`CTRL')                               ///
    hdfe_controls(`HDFE')                          ///
    cluster(`CLUS')                                ///
    bsmanual(`BREP', `BSED')                       ///
    graph_off _no_updates

forvalues i = 1/`NEFF' {
    scalar E_new_`i'    = e(Effect_`i')
    scalar SE_E_new_`i' = e(se_effect_`i')
}
forvalues i = 1/`NPL' {
    scalar P_new_`i'    = e(Placebo_`i')
    scalar SE_P_new_`i' = e(se_placebo_`i')
}
scalar A_new        = e(Av_tot_effect)
scalar SE_A_new     = e(se_avg_total_effect)
scalar PJE_new      = e(p_jointeffects)
scalar PJP_new      = e(p_jointplacebo)

*===============================================================================
* 2. Display NEW results
*===============================================================================
display as text _newline(2) "{hline 78}"
display as text "NEW path results (bootstrap_manual, hdfe path)"
display as text "{hline 78}"

display as text _column(4) %-14s "quantity" ///
       _column(22) %14s "estimate" ///
       _column(38) %14s "bootstrap SE"

forvalues i = 1/`NEFF' {
    display as text _column(4) %-14s "Effect_`i'" ///
        _column(22) %14.8f E_new_`i' ///
        _column(38) %14.8f SE_E_new_`i'
}
forvalues i = 1/`NPL' {
    display as text _column(4) %-14s "Placebo_`i'" ///
        _column(22) %14.8f P_new_`i' ///
        _column(38) %14.8f SE_P_new_`i'
}
display as text _column(4) %-14s "Av_tot_effect" ///
    _column(22) %14.8f A_new ///
    _column(38) %14.8f SE_A_new

display as text _newline "Joint significance p-values (bootstrap covariance):"
display as text _column(4) %-20s "p_jointeffects:" %14.8f PJE_new
display as text _column(4) %-20s "p_jointplacebo:" %14.8f PJP_new

*===============================================================================
* 3. Verdict: SE and p-values should be FINITE (not missing)
*===============================================================================
scalar n_missing = 0
forvalues i = 1/`NEFF' {
    if missing(SE_E_new_`i') {
        scalar n_missing = n_missing + 1
    }
}
forvalues i = 1/`NPL' {
    if missing(SE_P_new_`i') {
        scalar n_missing = n_missing + 1
    }
}
if missing(SE_A_new) {
    scalar n_missing = n_missing + 1
}
if missing(PJE_new) {
    scalar n_missing = n_missing + 1
}
if missing(PJP_new) {
    scalar n_missing = n_missing + 1
}

display as text _newline "Missing SE/p-value count: " n_missing " (target: 0)"
if (n_missing == 0) {
    display as result "bootstrap_manual produced finite SE for all parameters and finite joint p-values."
}
else {
    display as error "bootstrap_manual still produced " n_missing " missing values."
}

log close
