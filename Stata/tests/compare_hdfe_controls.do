*===============================================================================
* compare_hdfe_controls.do
*
* Parity test for the new hdfe_controls() option on did_multiplegt_dyn.ado.
*
* OLD path: original ado, controls(av_wind_speed wind_direction)
* NEW path: modified ado, controls(av_wind_speed wind_direction)
*                       + hdfe_controls(ac_uq_id)
*
* Because ac_uq_id is nested inside the auto-absorbed group FE
* (unique_small_grid_id), hdfe_controls(ac_uq_id) adds no identifying
* information, and the level fit on Y_level_XX is equivalent to the
* original first-difference fit. Effects, placebos, and average total
* effect should match to ~1e-10.
*
* No bootstrap on this run: we only check point estimates.
*===============================================================================

clear all
set more off
set linesize 140
cd "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata"
local data_path "test_sample.dta"
local old_ado   "did_multiplegt_dyn_baseline.ado"
local new_ado   "did_multiplegt_dyn.ado"
local log_path  "tests/log_compare_hdfe_controls.log"

local Y    "countk"
local G    "unique_small_grid_id"
local T    "monthyear"
local D    "downup_ac"
local CTRL "av_wind_speed wind_direction"
local CLUS "ac_uq_id"
local HDFE "ac_uq_id"
local NEFF 3
local NPL  3

capture log close
log using "`log_path'", text replace

display as text _newline "{hline 78}"
display as text "compare_hdfe_controls.do  --  OLD vs NEW point-estimate parity"
display as text "{hline 78}"

*-------------------------------------------------------------------------------
* 1. Quick sanity on data
*-------------------------------------------------------------------------------
capture confirm file "`data_path'"
if _rc {
    display as error "test_sample.dta not found at `data_path'."
    log close
    exit 198
}
use "`data_path'", clear
quietly count
display as text "Observations: " r(N)
quietly levelsof `G'
display as text "Distinct groups: " wordcount(`"`r(levels)'"')
quietly levelsof `CLUS'
display as text "Distinct ac_uq_id (cluster / hdfe): " wordcount(`"`r(levels)'"')

*===============================================================================
* 2. OLD path -- baseline ado, no hdfe_controls
*===============================================================================
display as text _newline(2) "{hline 78}"
display as text "OLD PATH"
display as text "{hline 78}"

use "`data_path'", clear
do "`old_ado'"

did_multiplegt_dyn `Y' `G' `T' `D',                ///
    effects(`NEFF') placebo(`NPL')                 ///
    controls(`CTRL')                               ///
    cluster(`CLUS') graph_off _no_updates

matrix b_old = e(b)
forvalues i = 1/`NEFF' {
    scalar E_old_`i' = e(Effect_`i')
}
forvalues i = 1/`NPL' {
    scalar P_old_`i' = e(Placebo_`i')
}
scalar A_old = e(Av_tot_effect)

*===============================================================================
* 3. NEW path -- modified ado, with hdfe_controls(ac_uq_id)
*===============================================================================
display as text _newline(2) "{hline 78}"
display as text "NEW PATH"
display as text "{hline 78}"

use "`data_path'", clear
do "`new_ado'"

did_multiplegt_dyn `Y' `G' `T' `D',                ///
    effects(`NEFF') placebo(`NPL')                 ///
    controls(`CTRL')                               ///
    hdfe_controls(`HDFE')                          ///
    cluster(`CLUS') graph_off _no_updates

matrix b_new = e(b)
forvalues i = 1/`NEFF' {
    scalar E_new_`i' = e(Effect_`i')
}
forvalues i = 1/`NPL' {
    scalar P_new_`i' = e(Placebo_`i')
}
scalar A_new = e(Av_tot_effect)

*===============================================================================
* 4. Side-by-side comparison
*===============================================================================
display as text _newline(2) "{hline 78}"
display as text "Side-by-side: effects, placebos, average total effect"
display as text "{hline 78}"

display as text _column(4) %-14s "quantity" ///
       _column(22) %14s "OLD" ///
       _column(40) %14s "NEW" ///
       _column(58) %14s "old-new"

forvalues i = 1/`NEFF' {
    display as text _column(4) %-14s "Effect_`i'" ///
        _column(22) %14.8f E_old_`i' ///
        _column(40) %14.8f E_new_`i' ///
        _column(58) %14.2e E_old_`i' - E_new_`i'
}
forvalues i = 1/`NPL' {
    display as text _column(4) %-14s "Placebo_`i'" ///
        _column(22) %14.8f P_old_`i' ///
        _column(40) %14.8f P_new_`i' ///
        _column(58) %14.2e P_old_`i' - P_new_`i'
}
display as text _column(4) %-14s "Av_tot_effect" ///
    _column(22) %14.8f A_old ///
    _column(40) %14.8f A_new ///
    _column(58) %14.2e A_old - A_new

display as text _newline "mreldif(b_old, b_new) = " %12.4e mreldif(b_old, b_new)

*===============================================================================
* 5. Verdict (pass if max abs diff < 1e-8)
*===============================================================================
scalar max_diff = 0
forvalues i = 1/`NEFF' {
    scalar _d = abs(E_old_`i' - E_new_`i')
    if (_d > max_diff) {
        scalar max_diff = _d
    }
}
forvalues i = 1/`NPL' {
    scalar _d = abs(P_old_`i' - P_new_`i')
    if (_d > max_diff) {
        scalar max_diff = _d
    }
}
scalar _d = abs(A_old - A_new)
if (_d > max_diff) {
    scalar max_diff = _d
}

display as text _newline "max abs diff across reported point estimates = " %12.4e max_diff
if max_diff < 1e-8 {
    display as result "PARITY CHECK PASSED (max diff < 1e-8)."
}
else {
    display as error "PARITY CHECK FAILED (max diff = " max_diff ", threshold 1e-8)."
}

log close
