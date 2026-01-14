********************************************************************************
* Test file for hdfe() option in did_multiplegt_dyn
* This file uses the LOCAL modified version of the ado file
********************************************************************************

clear all
set seed 12345
set more off

* IMPORTANT: Use the local modified version of the ado file
* First, drop any existing program definition
cap program drop did_multiplegt_dyn
cap program drop did_multiplegt_dyn_core_new

* Load the local version
adopath ++ "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata"

* Verify we're using the local version
which did_multiplegt_dyn

* Check if reghdfe is installed
cap which reghdfe
if _rc {
    ssc install reghdfe, replace
}

* Check if ftools is installed (required by reghdfe)
cap which ftools
if _rc {
    ssc install ftools, replace
}

********************************************************************************
* 1. Generate Synthetic Panel Data with Fixed Effects
********************************************************************************

* Parameters
local N_groups = 500        // Number of groups (e.g., firms)
local T_periods = 10        // Number of time periods
local N_industries = 50     // Number of industries (high-dimensional FE)
local N_regions = 20        // Number of regions

* Create panel structure
set obs `=`N_groups' * `T_periods''
gen group = ceil(_n / `T_periods')
bysort group: gen time = _n

* Assign time-invariant characteristics to groups
bysort group: gen industry = ceil(runiform() * `N_industries') if _n == 1
bysort group: replace industry = industry[1]

bysort group: gen region = ceil(runiform() * `N_regions') if _n == 1
bysort group: replace region = region[1]

* Generate treatment timing (staggered adoption)
bysort group: gen treat_time = . if _n == 1
bysort group: replace treat_time = ceil(runiform() * (`T_periods' + 3)) if _n == 1
bysort group: replace treat_time = treat_time[1]
replace treat_time = . if treat_time > `T_periods'

* Create treatment indicator
gen D = (time >= treat_time) & treat_time != .

* Generate fixed effects
* Industry fixed effects (time-invariant within industry)
gen industry_fe = 0
forvalues i = 1/`N_industries' {
    replace industry_fe = rnormal(0, 2) if industry == `i' & group == group[1]
}
bysort industry: replace industry_fe = industry_fe[1]

* Region fixed effects
gen region_fe = 0
forvalues r = 1/`N_regions' {
    replace region_fe = rnormal(0, 1.5) if region == `r' & group == group[1]
}
bysort region: replace region_fe = region_fe[1]

* Time fixed effects
gen time_fe = 0
forvalues t = 1/`T_periods' {
    replace time_fe = rnormal(0, 1) if time == `t' & _n == `t'
}
bysort time: replace time_fe = time_fe[1]

* Group fixed effects (unit-specific intercepts)
bysort group: gen group_fe = rnormal(0, 3) if _n == 1
bysort group: replace group_fe = group_fe[1]

* True treatment effect
local true_att = 2.5

* Generate outcome
gen Y = group_fe + time_fe + industry_fe + region_fe + `true_att' * D + rnormal(0, 1)

di ""
di "=========================================="
di "SYNTHETIC DATA SUMMARY"
di "=========================================="
di "Number of groups: `N_groups'"
di "Time periods: `T_periods'"
di "Number of industries: `N_industries'"
di "Number of regions: `N_regions'"
di "True ATT: `true_att'"
di ""

tab D
sum Y D industry region

********************************************************************************
* 2. Prepare dummy variables for controls() approach
********************************************************************************

* Create industry dummies interacted with time
tab industry, gen(ind_)

forvalues i = 1/`N_industries' {
    gen ind_time_`i' = ind_`i' * time
}

* Use subset for controls
local controls_subset ""
forvalues i = 1/20 {
    local controls_subset "`controls_subset' ind_time_`i'"
}

********************************************************************************
* 3. Test ORIGINAL approach: controls() with FE dummies
********************************************************************************

di ""
di "=========================================="
di "TEST 1: ORIGINAL APPROACH"
di "Using controls() with industry x time dummies"
di "=========================================="

timer clear 1
timer on 1

cap noi did_multiplegt_dyn Y group time D, effects(3) placebo(2) controls(`controls_subset') cluster(group) graph_off

timer off 1
timer list 1

* Store results
matrix results_controls = J(5,1,.)
if _rc == 0 {
    matrix results_controls[1,1] = e(effect_1)
    matrix results_controls[2,1] = e(se_effect_1)
    matrix results_controls[3,1] = e(placebo_1)
    matrix results_controls[4,1] = e(placebo_2)
    matrix results_controls[5,1] = 1  // success flag
}
else {
    matrix results_controls[5,1] = 0
}

********************************************************************************
* 4. Test NEW approach: hdfe() option
********************************************************************************

di ""
di "=========================================="
di "TEST 2: NEW APPROACH"
di "Using hdfe() with industry"
di "=========================================="

timer clear 2
timer on 2

cap noi did_multiplegt_dyn Y group time D, effects(3) placebo(2) hdfe(industry) cluster(group) graph_off

timer off 2
timer list 2

* Store results
matrix results_hdfe = J(5,1,.)
if _rc == 0 {
    matrix results_hdfe[1,1] = e(effect_1)
    matrix results_hdfe[2,1] = e(se_effect_1)
    matrix results_hdfe[3,1] = e(placebo_1)
    matrix results_hdfe[4,1] = e(placebo_2)
    matrix results_hdfe[5,1] = 1
}
else {
    di as error "HDFE approach failed with error code: " _rc
    matrix results_hdfe[5,1] = 0
}

********************************************************************************
* 5. Test with multiple HDFE variables
********************************************************************************

di ""
di "=========================================="
di "TEST 3: HDFE with multiple variables"
di "Using hdfe(industry region)"
di "=========================================="

timer clear 3
timer on 3

cap noi did_multiplegt_dyn Y group time D, effects(3) placebo(2) hdfe(industry region) cluster(group) graph_off

timer off 3
timer list 3

matrix results_hdfe2 = J(5,1,.)
if _rc == 0 {
    matrix results_hdfe2[1,1] = e(effect_1)
    matrix results_hdfe2[2,1] = e(se_effect_1)
    matrix results_hdfe2[3,1] = e(placebo_1)
    matrix results_hdfe2[4,1] = e(placebo_2)
    matrix results_hdfe2[5,1] = 1
}
else {
    di as error "HDFE approach with multiple vars failed: " _rc
    matrix results_hdfe2[5,1] = 0
}

********************************************************************************
* 6. Baseline comparison (no controls, no hdfe)
********************************************************************************

di ""
di "=========================================="
di "TEST 4: BASELINE (no controls, no hdfe)"
di "=========================================="

timer clear 4
timer on 4

cap noi did_multiplegt_dyn Y group time D, effects(3) placebo(2) cluster(group) graph_off

timer off 4
timer list 4

matrix results_baseline = J(5,1,.)
if _rc == 0 {
    matrix results_baseline[1,1] = e(effect_1)
    matrix results_baseline[2,1] = e(se_effect_1)
    matrix results_baseline[3,1] = e(placebo_1)
    matrix results_baseline[4,1] = e(placebo_2)
    matrix results_baseline[5,1] = 1
}

********************************************************************************
* 7. Summary Comparison
********************************************************************************

di ""
di "=========================================="
di "SUMMARY COMPARISON"
di "=========================================="
di ""
di "True ATT: `true_att'"
di ""
di "                      ATT        SE       Placebo1    Placebo2"
di "--------------------------------------------------------------"

if results_baseline[5,1] == 1 {
    di "Baseline       " %9.4f results_baseline[1,1] %9.4f results_baseline[2,1] %9.4f results_baseline[3,1] %9.4f results_baseline[4,1]
}

if results_controls[5,1] == 1 {
    di "Controls       " %9.4f results_controls[1,1] %9.4f results_controls[2,1] %9.4f results_controls[3,1] %9.4f results_controls[4,1]
}

if results_hdfe[5,1] == 1 {
    di "HDFE(ind)      " %9.4f results_hdfe[1,1] %9.4f results_hdfe[2,1] %9.4f results_hdfe[3,1] %9.4f results_hdfe[4,1]
}

if results_hdfe2[5,1] == 1 {
    di "HDFE(ind+reg)  " %9.4f results_hdfe2[1,1] %9.4f results_hdfe2[2,1] %9.4f results_hdfe2[3,1] %9.4f results_hdfe2[4,1]
}

di "--------------------------------------------------------------"
di ""
di "Timing comparison:"
timer list

di ""
di "=========================================="
di "TEST COMPLETED"
di "=========================================="
