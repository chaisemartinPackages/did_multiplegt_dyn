*===============================================================================
* test_hdfe_controls.do
*
* Tests the NEW hdfe_controls(varlist) option against the OLD help-file
* workaround (manual ac_uq_id-by-monthyear linear-trend controls), using the
* pre-built sample at test_sample.dta.
*
* Run prepare_test_sample.do FIRST to build that sample.
*
* This test uses effects(2) placebo(2) and the 12-month slice to keep runtime
* short. By construction the two paths should produce identical e(b) and e(V).
*===============================================================================

clear all
set more off
set linesize 120

*-------------------------------------------------------------------------------
* 0. Paths
*-------------------------------------------------------------------------------

local ado_path    "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata"
local sample_path "`ado_path'/test_sample.dta"
local log_path    "`ado_path'/test_hdfe_controls.log"
local Y_var       "countk"

capture log close
log using "`log_path'", text replace

adopath ++ "`ado_path'"
which did_multiplegt_dyn

*-------------------------------------------------------------------------------
* 1. Load the prepared sample
*-------------------------------------------------------------------------------

capture confirm file "`sample_path'"
if _rc {
    display as error "test_sample.dta not found at `sample_path'."
    display as error "Run prepare_test_sample.do first."
    log close
    exit 198
}

use "`sample_path'", clear

quietly count
display as text "Observations loaded: " r(N)
quietly levelsof unique_small_grid_id
display as text "Distinct groups:     " wordcount(`"`r(levels)'"')

*-------------------------------------------------------------------------------
* 2. Build OLD workaround controls: ac_uq_id dummies * monthyear (linear trends)
*    Drop the LAST ac_uq_id as reference. With 5 ac values that gives 4 controls.
*-------------------------------------------------------------------------------

quietly levelsof ac_uq_id, local(acs_kept)
local n_acs : word count `acs_kept'
local ref_ac : word `n_acs' of `acs_kept'

local old_workaround_controls ""
foreach a of local acs_kept {
    if `a' == `ref_ac' continue
    capture drop _DXt_ac_`a'
    gen double _DXt_ac_`a' = (ac_uq_id == `a') * monthyear
    local old_workaround_controls "`old_workaround_controls' _DXt_ac_`a'"
}
display as text _newline "OLD workaround controls: `old_workaround_controls'"

*===============================================================================
* 3. OLD PATH
*===============================================================================

display as text _newline(2) "{hline 78}"
display as text "OLD PATH: controls(av_wind_speed wind_direction `old_workaround_controls')"
display as text "{hline 78}"

did_multiplegt_dyn `Y_var' unique_small_grid_id monthyear downup_ac, ///
    effects(2) placebo(2)                                            ///
    controls(av_wind_speed wind_direction `old_workaround_controls') ///
    cluster(ac_uq_id) graph_off _no_updates

matrix b_old = e(b)
matrix V_old = e(V)
local colnames_old : colfullnames b_old

*===============================================================================
* 4. NEW PATH
*===============================================================================

display as text _newline(2) "{hline 78}"
display as text "NEW PATH: controls(av_wind_speed wind_direction)"
display as text "          hdfe_controls(i.ac_uq_id)   <- reghdfe-style; auto multiplies by T"
display as text "{hline 78}"

did_multiplegt_dyn `Y_var' unique_small_grid_id monthyear downup_ac, ///
    effects(2) placebo(2)                                            ///
    controls(av_wind_speed wind_direction)                           ///
    hdfe_controls(i.ac_uq_id)                                        ///
    cluster(ac_uq_id) graph_off _no_updates

matrix b_new = e(b)
matrix V_new = e(V)
local colnames_new : colfullnames b_new

*===============================================================================
* 5. Comparison
*===============================================================================

display as text _newline(2) "{hline 78}"
display as text "Comparison: OLD vs NEW"
display as text "{hline 78}"

display as text "OLD coefficient names: `colnames_old'"
display as text "NEW coefficient names: `colnames_new'"

capture matrix diff_b = b_old - b_new
if _rc == 0 {
    matrix list diff_b, title("Coefficient differences (old - new)")
}
else {
    display as error "Coefficient vectors have different shapes; manual inspection needed."
    matrix list b_old, title("OLD coefficients")
    matrix list b_new, title("NEW coefficients")
}

capture matrix diff_V = V_old - V_new
if _rc == 0 {
    display as text _newline "Diagonal of V_old vs V_new:"
    forvalues i = 1/`= colsof(b_old)' {
        display as text "  " %-30s "`: word `i' of `colnames_old''" ///
            " var_old = " %12.6g V_old[`i', `i'] ///
            "   var_new = " %12.6g V_new[`i', `i']
    }
}
else {
    display as error "V matrices have different shapes; manual inspection needed."
}

display as text _newline "Summary metrics (small means a match):"
display as text "  mreldif(b_old, b_new) = " mreldif(b_old, b_new)
capture display as text "  mreldif(V_old, V_new) = " mreldif(V_old, V_new)

display as text _newline(2) "{hline 78}"
display as text "Expected: mreldif == 0 and diff_b all zeros, since hdfe_controls(ac_uq_id)"
display as text "auto-expands into the same D_a * monthyear controls the OLD path types"
display as text "by hand (with the same a == last-level reference dropped)."
display as text "{hline 78}"

log close
