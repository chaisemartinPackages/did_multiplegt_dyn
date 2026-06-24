*===============================================================================
* prepare_test_sample.do
*
* Builds a small, deterministic test sample from data_test.csv and saves it as
* test_sample.dta so that test_hdfe_controls.do and test_hdfe_cells.do can
* reuse it directly without re-importing the 14M-row CSV.
*
* The slice:
*   - 5 smallest ac_uq_id values
*   - first 12 distinct monthyear values
*   - rows non-missing in outcome, treatment, the two controls, ac_uq_id
*
* Validation checks before saving:
*   1. At least 5 distinct groups whose treatment changes over time (switchers).
*   2. Within each baseline treatment level, at least 2 distinct first-switch
*      dates F_g (otherwise the did_multiplegt_dyn baseline cell is dropped).
*
* If either check fails, the script aborts before writing the .dta so you do
* not silently feed a degenerate sample into the regressions.
*===============================================================================

clear all
set more off
set linesize 120

local data_path "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata/data_test.csv"
local save_path "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata/test_sample.dta"
local log_path  "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata/prepare_test_sample.log"

capture log close
log using "`log_path'", text replace

*-------------------------------------------------------------------------------
* 1. Load and slice
*-------------------------------------------------------------------------------

import delimited using "`data_path'", clear varnames(1)

* Build the outcome variable.
gen countk = count * 1000

* Sanity: required columns present.
confirm variable unique_small_grid_id
confirm variable monthyear
confirm variable downup_ac
confirm variable ac_uq_id
confirm variable av_wind_speed
confirm variable wind_direction
confirm variable countk

* Pick the 5 smallest ac_uq_id values.
quietly levelsof ac_uq_id, local(all_acs)
local first5 ""
local k = 0
foreach a of local all_acs {
    local ++k
    if `k' <= 5 local first5 "`first5' `a'"
}
display as text "Keeping ac_uq_id values:`first5'"

gen byte _keep_ac = 0
foreach a of local first5 {
    replace _keep_ac = 1 if ac_uq_id == `a'
}
keep if _keep_ac == 1
drop _keep_ac

* Keep first 12 distinct monthyear values.
sort monthyear
egen _month_rank = group(monthyear)
keep if _month_rank >= 1 & _month_rank <= 12
drop _month_rank

* Drop rows missing key variables.
keep if !missing(countk, downup_ac, av_wind_speed, wind_direction, ac_uq_id, ///
                 unique_small_grid_id, monthyear)

display as text _newline "Sample sizes after filtering:"
quietly count
display as text "    Observations:        " r(N)
quietly levelsof unique_small_grid_id
display as text "    Distinct groups:     " wordcount(`"`r(levels)'"')
quietly levelsof ac_uq_id
display as text "    Distinct ac_uq_id:   " wordcount(`"`r(levels)'"')
quietly levelsof monthyear
display as text "    Distinct months:     " wordcount(`"`r(levels)'"')

*-------------------------------------------------------------------------------
* 2. Validation: are there switching units, and is switching-date variation
*    present within each baseline-treatment cell?
*-------------------------------------------------------------------------------

* Mirror the ado's definitions of baseline treatment (D_g,1) and F_g
* (first time period where treatment differs from baseline).
sort unique_small_grid_id monthyear
by unique_small_grid_id: gen byte _is_first_obs = (_n == 1)
gen double _D_baseline_tmp = downup_ac if _is_first_obs == 1
bys unique_small_grid_id: egen double _D_baseline = mean(_D_baseline_tmp)
drop _D_baseline_tmp _is_first_obs

gen byte _diff_from_baseline = (downup_ac != _D_baseline)
bys unique_small_grid_id (monthyear): gen double _F_g_tmp = monthyear if _diff_from_baseline == 1 & _diff_from_baseline[_n-1] != 1
bys unique_small_grid_id: egen double _F_g = min(_F_g_tmp)
drop _diff_from_baseline _F_g_tmp

* Mark groups whose treatment ever changes (switchers).
gen byte _is_switcher = !missing(_F_g)

* Count switchers.
quietly levelsof unique_small_grid_id if _is_switcher == 1, local(_switchers_)
local n_switching_groups : word count `_switchers_'
display as text _newline "Validation:"
display as text "    Switching groups in slice: " `n_switching_groups'

if `n_switching_groups' < 5 {
    display as error "Too few switching groups in slice (< 5). Aborting prep."
    log close
    exit 198
}

* Within each baseline level, count distinct F_g values among switchers.
preserve
    keep unique_small_grid_id _D_baseline _F_g _is_switcher
    duplicates drop unique_small_grid_id, force
    keep if _is_switcher == 1
    bys _D_baseline: gegen long _n_F_g_in_d = nunique(_F_g)
    display as text _newline "Distinct switching dates F_g by baseline treatment:"
    table _D_baseline, statistic(mean _n_F_g_in_d) statistic(count _F_g)
    quietly sum _n_F_g_in_d
    local min_n_F_g = r(min)
    display as text "    Minimum distinct F_g across baseline cells: " `min_n_F_g'
    if `min_n_F_g' < 2 {
        display as error "Some baseline cells have <2 distinct F_g values. The"
        display as error "did_multiplegt_dyn estimator will drop those cells."
        display as error "The sample is still usable, but expect fewer effective groups."
    }
restore

* Drop helper variables before saving so the test do-files see a clean dataset.
drop _D_baseline _F_g _is_switcher

*-------------------------------------------------------------------------------
* 3. Save
*-------------------------------------------------------------------------------

save "`save_path'", replace
display as text _newline "Saved test sample to: `save_path'"
display as text "Variables in saved file:"
describe, simple

log close
