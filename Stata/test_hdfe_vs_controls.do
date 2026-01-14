********************************************************************************
* Comprehensive comparison: hdfe() vs controls() with FE dummies
* This validates that hdfe() produces equivalent results to controls()
********************************************************************************

clear all
set seed 12345
set more off

* Load local version
cap program drop did_multiplegt_dyn
cap program drop did_multiplegt_dyn_core_new
adopath ++ "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata"
which did_multiplegt_dyn

********************************************************************************
* Create test data with industry-specific trends
********************************************************************************

local N_groups = 300
local T_periods = 10
local N_industries = 15

set obs `=`N_groups' * `T_periods''
gen group = ceil(_n / `T_periods')
bysort group: gen time = _n

* Industry - time invariant within groups
bysort group: gen industry = ceil(runiform() * `N_industries') if _n == 1
bysort group: replace industry = industry[1]

* Treatment timing (staggered)
bysort group: gen treat_time = ceil(runiform() * (`T_periods' + 2)) if _n == 1
bysort group: replace treat_time = treat_time[1]
replace treat_time = . if treat_time > `T_periods'

gen D = (time >= treat_time) & treat_time != .

* Generate industry-specific LINEAR TRENDS
* Each industry i has trend coefficient = (i - 8) * 0.4
gen industry_trend = (industry - 8) * 0.4

* Group FE
bysort group: gen group_fe = rnormal(0, 2) if _n == 1
bysort group: replace group_fe = group_fe[1]

* True ATT
local true_att = 3.0

* Y = group_fe + industry_trend * time + ATT * D + noise
gen Y = group_fe + industry_trend * time + `true_att' * D + rnormal(0, 1)

di ""
di "=============================================="
di "DATA SUMMARY"
di "=============================================="
di "True ATT: `true_att'"
di "Groups: `N_groups', Periods: `T_periods', Industries: `N_industries'"
di ""
tab D

********************************************************************************
* Create industry×time dummies for controls() approach
********************************************************************************

* Method: Create industry dummies interacted with time
* This is equivalent to including industry-specific linear trends

tab industry, gen(ind_)

* Create industry × time interactions
forvalues i = 1/`N_industries' {
    gen ind_time_`i' = ind_`i' * time
}

* Build controls varlist
local controls_list ""
forvalues i = 1/`N_industries' {
    local controls_list "`controls_list' ind_time_`i'"
}

di "Created `N_industries' industry × time interaction variables for controls()"

********************************************************************************
* TEST 1: BASELINE (no adjustments)
********************************************************************************

di ""
di "=============================================="
di "TEST 1: BASELINE (no controls, no hdfe)"
di "=============================================="

timer clear 1
timer on 1

did_multiplegt_dyn Y group time D, effects(3) placebo(2) cluster(group) graph_off

timer off 1

matrix results_baseline = J(7,1,.)
matrix results_baseline[1,1] = e(effect_1)
matrix results_baseline[2,1] = e(se_effect_1)
matrix results_baseline[3,1] = e(effect_2)
matrix results_baseline[4,1] = e(effect_3)
matrix results_baseline[5,1] = e(placebo_1)
matrix results_baseline[6,1] = e(placebo_2)
matrix results_baseline[7,1] = e(N_effect_1)

********************************************************************************
* TEST 2: CONTROLS() with industry×time dummies
********************************************************************************

di ""
di "=============================================="
di "TEST 2: CONTROLS() with industry×time dummies"
di "=============================================="

timer clear 2
timer on 2

did_multiplegt_dyn Y group time D, effects(3) placebo(2) controls(`controls_list') cluster(group) graph_off

timer off 2

matrix results_controls = J(7,1,.)
matrix results_controls[1,1] = e(effect_1)
matrix results_controls[2,1] = e(se_effect_1)
matrix results_controls[3,1] = e(effect_2)
matrix results_controls[4,1] = e(effect_3)
matrix results_controls[5,1] = e(placebo_1)
matrix results_controls[6,1] = e(placebo_2)
matrix results_controls[7,1] = e(N_effect_1)

********************************************************************************
* TEST 3: HDFE(industry)
********************************************************************************

di ""
di "=============================================="
di "TEST 3: HDFE(industry)"
di "=============================================="

timer clear 3
timer on 3

did_multiplegt_dyn Y group time D, effects(3) placebo(2) hdfe(industry) cluster(group) graph_off

timer off 3

matrix results_hdfe = J(7,1,.)
matrix results_hdfe[1,1] = e(effect_1)
matrix results_hdfe[2,1] = e(se_effect_1)
matrix results_hdfe[3,1] = e(effect_2)
matrix results_hdfe[4,1] = e(effect_3)
matrix results_hdfe[5,1] = e(placebo_1)
matrix results_hdfe[6,1] = e(placebo_2)
matrix results_hdfe[7,1] = e(N_effect_1)

********************************************************************************
* TEST 4: trends_nonparam(industry) for comparison
********************************************************************************

di ""
di "=============================================="
di "TEST 4: trends_nonparam(industry)"
di "=============================================="

timer clear 4
timer on 4

cap noi did_multiplegt_dyn Y group time D, effects(3) placebo(2) trends_nonparam(industry) cluster(group) graph_off

timer off 4

matrix results_trends = J(7,1,.)
if _rc == 0 {
    matrix results_trends[1,1] = e(effect_1)
    matrix results_trends[2,1] = e(se_effect_1)
    matrix results_trends[3,1] = e(effect_2)
    matrix results_trends[4,1] = e(effect_3)
    matrix results_trends[5,1] = e(placebo_1)
    matrix results_trends[6,1] = e(placebo_2)
    matrix results_trends[7,1] = e(N_effect_1)
}

********************************************************************************
* RESULTS COMPARISON
********************************************************************************

di ""
di "=============================================="
di "RESULTS COMPARISON"
di "=============================================="
di ""
di "True ATT: `true_att'"
di ""
di "                        Effect_1     SE      Effect_2   Effect_3   Placebo_1  Placebo_2     N"
di "--------------------------------------------------------------------------------------------"

di "Baseline         " %12.4f results_baseline[1,1] %9.4f results_baseline[2,1] ///
   %11.4f results_baseline[3,1] %11.4f results_baseline[4,1] ///
   %11.4f results_baseline[5,1] %11.4f results_baseline[6,1] ///
   %8.0f results_baseline[7,1]

di "Controls(ind*t)  " %12.4f results_controls[1,1] %9.4f results_controls[2,1] ///
   %11.4f results_controls[3,1] %11.4f results_controls[4,1] ///
   %11.4f results_controls[5,1] %11.4f results_controls[6,1] ///
   %8.0f results_controls[7,1]

di "HDFE(industry)   " %12.4f results_hdfe[1,1] %9.4f results_hdfe[2,1] ///
   %11.4f results_hdfe[3,1] %11.4f results_hdfe[4,1] ///
   %11.4f results_hdfe[5,1] %11.4f results_hdfe[6,1] ///
   %8.0f results_hdfe[7,1]

if results_trends[1,1] != . {
    di "trends_nonparam  " %12.4f results_trends[1,1] %9.4f results_trends[2,1] ///
       %11.4f results_trends[3,1] %11.4f results_trends[4,1] ///
       %11.4f results_trends[5,1] %11.4f results_trends[6,1] ///
       %8.0f results_trends[7,1]
}

di "--------------------------------------------------------------------------------------------"
di ""

* Check if HDFE and Controls produce similar results
local diff_effect1 = abs(results_hdfe[1,1] - results_controls[1,1])
local diff_placebo1 = abs(results_hdfe[5,1] - results_controls[5,1])

di "Difference between HDFE and Controls:"
di "  Effect_1 difference: " %8.4f `diff_effect1'
di "  Placebo_1 difference: " %8.4f `diff_placebo1'
di ""

di "Timing comparison:"
timer list

di ""
di "=============================================="
di "TEST COMPLETED"
di "=============================================="
