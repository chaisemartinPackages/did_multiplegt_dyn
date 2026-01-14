********************************************************************************
* Performance comparison: hdfe() vs controls() with large dataset
* 100,000 observations, 50 FE categories
********************************************************************************

clear all
set seed 12345
set more off
set maxvar 10000

* Load local version
cap program drop did_multiplegt_dyn
cap program drop did_multiplegt_dyn_core_new
adopath ++ "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata"
which did_multiplegt_dyn

********************************************************************************
* Configuration
********************************************************************************

local N_obs = 100000
local N_industries = 50
local T_periods = 10
local N_groups = ceil(`N_obs' / `T_periods')

di ""
di "=============================================="
di "CONFIGURATION"
di "=============================================="
di "Target observations: `N_obs'"
di "Groups: `N_groups'"
di "Time periods: `T_periods'"
di "FE categories (industries): `N_industries'"
di "=============================================="

********************************************************************************
* Create test data with industry-specific trends
********************************************************************************

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
* Each industry i has trend coefficient = (i - 25) * 0.3
gen industry_trend = (industry - 25) * 0.3

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
di "Actual observations: " _N
di ""
tab D
di ""
tab industry if time == 1, summarize(industry_trend)

********************************************************************************
* Create industry×time dummies for controls() approach
********************************************************************************

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

di ""
di "Created `N_industries' industry × time interaction variables for controls()"
di ""

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
timer list 1
local time_baseline = r(t1)

matrix results_baseline = J(12,1,.)
matrix results_baseline[1,1] = e(effect_1)
matrix results_baseline[2,1] = e(se_effect_1)
matrix results_baseline[3,1] = e(effect_2)
matrix results_baseline[4,1] = e(se_effect_2)
matrix results_baseline[5,1] = e(effect_3)
matrix results_baseline[6,1] = e(se_effect_3)
matrix results_baseline[7,1] = e(placebo_1)
matrix results_baseline[8,1] = e(se_placebo_1)
matrix results_baseline[9,1] = e(placebo_2)
matrix results_baseline[10,1] = e(se_placebo_2)
matrix results_baseline[11,1] = e(N_effect_1)
matrix results_baseline[12,1] = `time_baseline'

********************************************************************************
* TEST 2: CONTROLS() with industry×time dummies
********************************************************************************

di ""
di "=============================================="
di "TEST 2: CONTROLS() with `N_industries' industry×time dummies"
di "=============================================="

timer clear 2
timer on 2

did_multiplegt_dyn Y group time D, effects(3) placebo(2) controls(`controls_list') cluster(group) graph_off

timer off 2
timer list 2
local time_controls = r(t2)

matrix results_controls = J(12,1,.)
matrix results_controls[1,1] = e(effect_1)
matrix results_controls[2,1] = e(se_effect_1)
matrix results_controls[3,1] = e(effect_2)
matrix results_controls[4,1] = e(se_effect_2)
matrix results_controls[5,1] = e(effect_3)
matrix results_controls[6,1] = e(se_effect_3)
matrix results_controls[7,1] = e(placebo_1)
matrix results_controls[8,1] = e(se_placebo_1)
matrix results_controls[9,1] = e(placebo_2)
matrix results_controls[10,1] = e(se_placebo_2)
matrix results_controls[11,1] = e(N_effect_1)
matrix results_controls[12,1] = `time_controls'

********************************************************************************
* TEST 3: HDFE(industry)
********************************************************************************

di ""
di "=============================================="
di "TEST 3: HDFE(industry) - `N_industries' categories"
di "=============================================="

timer clear 3
timer on 3

did_multiplegt_dyn Y group time D, effects(3) placebo(2) hdfe(industry) cluster(group) graph_off

timer off 3
timer list 3
local time_hdfe = r(t3)

matrix results_hdfe = J(12,1,.)
matrix results_hdfe[1,1] = e(effect_1)
matrix results_hdfe[2,1] = e(se_effect_1)
matrix results_hdfe[3,1] = e(effect_2)
matrix results_hdfe[4,1] = e(se_effect_2)
matrix results_hdfe[5,1] = e(effect_3)
matrix results_hdfe[6,1] = e(se_effect_3)
matrix results_hdfe[7,1] = e(placebo_1)
matrix results_hdfe[8,1] = e(se_placebo_1)
matrix results_hdfe[9,1] = e(placebo_2)
matrix results_hdfe[10,1] = e(se_placebo_2)
matrix results_hdfe[11,1] = e(N_effect_1)
matrix results_hdfe[12,1] = `time_hdfe'

********************************************************************************
* RESULTS COMPARISON - FULL PRECISION
********************************************************************************

di ""
di "=============================================="
di "RESULTS COMPARISON - FULL PRECISION"
di "=============================================="
di ""
di "Configuration:"
di "  Observations: " _N
di "  Groups: `N_groups'"
di "  Time periods: `T_periods'"
di "  FE categories: `N_industries'"
di "  True ATT: `true_att'"
di ""
di "=============================================="
di "POINT ESTIMATES (Full Precision)"
di "=============================================="
di ""
di "                    BASELINE              CONTROLS            HDFE"
di "-------------------------------------------------------------------"
di "Effect_1:     " %18.10f results_baseline[1,1] %18.10f results_controls[1,1] %18.10f results_hdfe[1,1]
di "Effect_2:     " %18.10f results_baseline[3,1] %18.10f results_controls[3,1] %18.10f results_hdfe[3,1]
di "Effect_3:     " %18.10f results_baseline[5,1] %18.10f results_controls[5,1] %18.10f results_hdfe[5,1]
di "Placebo_1:    " %18.10f results_baseline[7,1] %18.10f results_controls[7,1] %18.10f results_hdfe[7,1]
di "Placebo_2:    " %18.10f results_baseline[9,1] %18.10f results_controls[9,1] %18.10f results_hdfe[9,1]
di ""

di "=============================================="
di "STANDARD ERRORS (Full Precision)"
di "=============================================="
di ""
di "                    BASELINE              CONTROLS            HDFE"
di "-------------------------------------------------------------------"
di "SE_Effect_1:  " %18.10f results_baseline[2,1] %18.10f results_controls[2,1] %18.10f results_hdfe[2,1]
di "SE_Effect_2:  " %18.10f results_baseline[4,1] %18.10f results_controls[4,1] %18.10f results_hdfe[4,1]
di "SE_Effect_3:  " %18.10f results_baseline[6,1] %18.10f results_controls[6,1] %18.10f results_hdfe[6,1]
di "SE_Placebo_1: " %18.10f results_baseline[8,1] %18.10f results_controls[8,1] %18.10f results_hdfe[8,1]
di "SE_Placebo_2: " %18.10f results_baseline[10,1] %18.10f results_controls[10,1] %18.10f results_hdfe[10,1]
di ""

di "=============================================="
di "NUMERICAL DIFFERENCES: HDFE vs CONTROLS"
di "=============================================="
di ""
di "Effect_1 diff:    " %18.10f (results_hdfe[1,1] - results_controls[1,1])
di "Effect_2 diff:    " %18.10f (results_hdfe[3,1] - results_controls[3,1])
di "Effect_3 diff:    " %18.10f (results_hdfe[5,1] - results_controls[5,1])
di "Placebo_1 diff:   " %18.10f (results_hdfe[7,1] - results_controls[7,1])
di "Placebo_2 diff:   " %18.10f (results_hdfe[9,1] - results_controls[9,1])
di "SE_Effect_1 diff: " %18.10f (results_hdfe[2,1] - results_controls[2,1])
di ""

di "=============================================="
di "EXECUTION TIME COMPARISON"
di "=============================================="
di ""
di "Baseline:        " %8.2f results_baseline[12,1] " seconds"
di "Controls():      " %8.2f results_controls[12,1] " seconds"
di "HDFE():          " %8.2f results_hdfe[12,1] " seconds"
di ""
di "Speed improvement (Controls/HDFE): " %6.2f (results_controls[12,1]/results_hdfe[12,1]) "x faster"
di ""

di "=============================================="
di "SAMPLE SIZE"
di "=============================================="
di ""
di "N_effect_1: " %10.0f results_baseline[11,1]
di ""

di "=============================================="
di "TEST COMPLETED SUCCESSFULLY"
di "=============================================="
