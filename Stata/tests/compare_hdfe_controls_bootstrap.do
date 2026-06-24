*===============================================================================
* compare_hdfe_controls_bootstrap.do
*
* Bootstrap-on parity test for hdfe_controls() vs the baseline ado.
*
* OLD path: original ado, controls(av_wind_speed wind_direction)
* NEW path: modified ado, controls(av_wind_speed wind_direction)
*                       + hdfe_controls(ac_uq_id)
*
* Both runs use bootstrap(5, 1234). With the same seed, the bootstrap
* command draws identical cluster resamples; the point estimates on
* each resample are identical (proven by the no-bootstrap parity test);
* therefore the bootstrap SE should be identical.
*
* The seed is printed before each invocation so both runs can be
* visually confirmed to use the same value.
*===============================================================================

clear all
set more off
set linesize 140
cd "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata"

local data_path "test_sample.dta"
local old_ado   "did_multiplegt_dyn_baseline.ado"
local new_ado   "did_multiplegt_dyn.ado"
local log_path  "tests/log_compare_hdfe_controls_bootstrap.log"

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
display as text "compare_hdfe_controls_bootstrap.do  --  OLD vs NEW point estimates AND bootstrap SE"
display as text "Both runs use bootstrap(`BREP', `BSED'); seed is printed before each call."
display as text "{hline 78}"

*-------------------------------------------------------------------------------
* 1. Sanity on data
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

*===============================================================================
* 2. OLD path -- baseline ado, no hdfe_controls
*===============================================================================
display as text _newline(2) "{hline 78}"
display as text "OLD PATH"
display as text "{hline 78}"

use "`data_path'", clear
do "`old_ado'"

display as input _newline ">>> SEED FOR OLD RUN: `BSED'   (reps = `BREP')"
display as input ">>> Stata's c(seed) before OLD call: " c(seed)

did_multiplegt_dyn `Y' `G' `T' `D',                ///
    effects(`NEFF') placebo(`NPL')                 ///
    controls(`CTRL')                               ///
    cluster(`CLUS')                                ///
    bootstrap(`BREP', `BSED')                      ///
    graph_off _no_updates

matrix b_old = e(b)
matrix V_old = e(V)
forvalues i = 1/`NEFF' {
    scalar E_old_`i'    = e(Effect_`i')
    scalar SE_E_old_`i' = e(se_effect_`i')
}
forvalues i = 1/`NPL' {
    scalar P_old_`i'    = e(Placebo_`i')
    scalar SE_P_old_`i' = e(se_placebo_`i')
}
scalar A_old        = e(Av_tot_effect)
scalar SE_A_old     = e(se_avg_total_effect)
scalar PJE_old      = e(p_jointeffects)
scalar PJP_old      = e(p_jointplacebo)

*===============================================================================
* 3. NEW path -- modified ado, with hdfe_controls(ac_uq_id)
*===============================================================================
display as text _newline(2) "{hline 78}"
display as text "NEW PATH"
display as text "{hline 78}"

use "`data_path'", clear
do "`new_ado'"

display as input _newline ">>> SEED FOR NEW RUN: `BSED'   (reps = `BREP')"
display as input ">>> Stata's c(seed) before NEW call: " c(seed)

did_multiplegt_dyn `Y' `G' `T' `D',                ///
    effects(`NEFF') placebo(`NPL')                 ///
    controls(`CTRL')                               ///
    hdfe_controls(`HDFE')                          ///
    cluster(`CLUS')                                ///
    bootstrap(`BREP', `BSED')                      ///
    graph_off _no_updates

matrix b_new = e(b)
matrix V_new = e(V)
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
* 4. Side-by-side: estimates AND bootstrap SE
*===============================================================================
display as text _newline(2) "{hline 100}"
display as text "Side-by-side: point estimates AND bootstrap SE (reps=`BREP', seed=`BSED')"
display as text "{hline 100}"

display as text _column(4) %-14s "quantity" ///
       _column(22) %14s "OLD est" ///
       _column(38) %14s "NEW est" ///
       _column(54) %12s "est diff" ///
       _column(70) %14s "OLD SE" ///
       _column(86) %14s "NEW SE" ///
       _column(102) %12s "SE diff"

forvalues i = 1/`NEFF' {
    display as text _column(4) %-14s "Effect_`i'" ///
        _column(22) %14.8f E_old_`i' ///
        _column(38) %14.8f E_new_`i' ///
        _column(54) %12.2e E_old_`i' - E_new_`i' ///
        _column(70) %14.8f SE_E_old_`i' ///
        _column(86) %14.8f SE_E_new_`i' ///
        _column(102) %12.2e SE_E_old_`i' - SE_E_new_`i'
}
forvalues i = 1/`NPL' {
    display as text _column(4) %-14s "Placebo_`i'" ///
        _column(22) %14.8f P_old_`i' ///
        _column(38) %14.8f P_new_`i' ///
        _column(54) %12.2e P_old_`i' - P_new_`i' ///
        _column(70) %14.8f SE_P_old_`i' ///
        _column(86) %14.8f SE_P_new_`i' ///
        _column(102) %12.2e SE_P_old_`i' - SE_P_new_`i'
}
display as text _column(4) %-14s "Av_tot_effect" ///
    _column(22) %14.8f A_old ///
    _column(38) %14.8f A_new ///
    _column(54) %12.2e A_old - A_new ///
    _column(70) %14.8f SE_A_old ///
    _column(86) %14.8f SE_A_new ///
    _column(102) %12.2e SE_A_old - SE_A_new

display as text _newline "mreldif(b_old, b_new) = " %12.4e mreldif(b_old, b_new)
display as text "mreldif(V_old, V_new) = " %12.4e mreldif(V_old, V_new)

display as text _newline "Joint significance p-values (built from bootstrap vcov):"
display as text _column(4) %-20s "p_jointeffects:" %14.8f PJE_old _column(50) %14.8f PJE_new _column(70) %12.2e PJE_old - PJE_new
display as text _column(4) %-20s "p_jointplacebo:" %14.8f PJP_old _column(50) %14.8f PJP_new _column(70) %12.2e PJP_old - PJP_new

*===============================================================================
* 5. Verdict on BOTH point estimates and SE (pass if max abs diff < 1e-8)
*===============================================================================
scalar max_est_diff = 0
scalar max_se_diff  = 0

forvalues i = 1/`NEFF' {
    scalar _de = abs(E_old_`i'    - E_new_`i')
    scalar _ds = abs(SE_E_old_`i' - SE_E_new_`i')
    if (_de > max_est_diff) {
        scalar max_est_diff = _de
    }
    if (_ds > max_se_diff) {
        scalar max_se_diff = _ds
    }
}
forvalues i = 1/`NPL' {
    scalar _de = abs(P_old_`i'    - P_new_`i')
    scalar _ds = abs(SE_P_old_`i' - SE_P_new_`i')
    if (_de > max_est_diff) {
        scalar max_est_diff = _de
    }
    if (_ds > max_se_diff) {
        scalar max_se_diff = _ds
    }
}
scalar _de = abs(A_old    - A_new)
scalar _ds = abs(SE_A_old - SE_A_new)
if (_de > max_est_diff) {
    scalar max_est_diff = _de
}
if (_ds > max_se_diff) {
    scalar max_se_diff = _ds
}

display as text _newline "max abs diff (point estimates) = " %12.4e max_est_diff
display as text "max abs diff (bootstrap SE)    = " %12.4e max_se_diff
if (max_est_diff < 1e-8) & (max_se_diff < 1e-8) {
    display as result "BOOTSTRAP PARITY PASSED (estimates AND SE match within 1e-8)."
}
else if (max_est_diff < 1e-8) & (max_se_diff >= 1e-8) {
    display as error "POINT ESTIMATES MATCH but BOOTSTRAP SE DIFFER (max diff = " max_se_diff ")."
}
else {
    display as error "PARITY FAILED (est diff = " max_est_diff ", SE diff = " max_se_diff ")."
}

log close
