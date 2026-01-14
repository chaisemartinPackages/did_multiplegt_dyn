********************************************************************************
* Test file for hdfe() option in did_multiplegt_dyn
* This file creates synthetic data with fixed effects and compares:
* 1. Original approach: FE as dummies interacted with time in controls()
* 2. New approach: hdfe() option using reghdfe
********************************************************************************

clear all
set seed 12345
set more off

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
* Industry and region are constant within groups (time-invariant)
bysort group: gen industry = ceil(runiform() * `N_industries') if _n == 1
bysort group: replace industry = industry[1]

bysort group: gen region = ceil(runiform() * `N_regions') if _n == 1
bysort group: replace region = region[1]

* Generate treatment timing (staggered adoption)
* Some groups never treated, others treated at different times
bysort group: gen treat_time = . if _n == 1
bysort group: replace treat_time = ceil(runiform() * (`T_periods' + 3)) if _n == 1
bysort group: replace treat_time = treat_time[1]
* Groups with treat_time > T_periods are never-treated
replace treat_time = . if treat_time > `T_periods'

* Create treatment indicator
gen D = (time >= treat_time) & treat_time != .

* Generate fixed effects (these are what we want to absorb)
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

* True treatment effect (constant for simplicity)
local true_att = 2.5

* Generate outcome with:
* Y = group_fe + time_fe + industry_fe + region_fe + treatment_effect + noise
gen Y = group_fe + time_fe + industry_fe + region_fe + `true_att' * D + rnormal(0, 1)

* Create time-varying covariate (for controls comparison)
gen X1 = rnormal(0, 1) + 0.5 * time

* Label variables
label var Y "Outcome"
label var group "Group ID"
label var time "Time period"
label var D "Treatment indicator"
label var industry "Industry code (time-invariant)"
label var region "Region code (time-invariant)"
label var treat_time "Treatment timing"

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

di ""
di "=========================================="
di "CREATING INDUSTRY x TIME DUMMIES"
di "=========================================="

* Create industry dummies interacted with time
* This is computationally expensive with many industries
tab industry, gen(ind_)

* Interact with time (as linear trend)
forvalues i = 1/`N_industries' {
    gen ind_time_`i' = ind_`i' * time
}

* Create region dummies interacted with time
tab region, gen(reg_)

forvalues r = 1/`N_regions' {
    gen reg_time_`r' = reg_`r' * time
}

di "Created `N_industries' industry x time interactions"
di "Created `N_regions' region x time interactions"

********************************************************************************
* 3. Save the dataset
********************************************************************************

save "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata/test_data_hdfe.dta", replace

********************************************************************************
* 4. Test ORIGINAL approach: controls() with FE dummies
********************************************************************************

di ""
di "=========================================="
di "TEST 1: ORIGINAL APPROACH"
di "Using controls() with industry x time dummies"
di "=========================================="

* Build the controls varlist (industry x time interactions)
local controls_list ""
forvalues i = 1/`N_industries' {
    local controls_list "`controls_list' ind_time_`i'"
}

* Note: Including all industry dummies would be very slow
* We'll use a subset for demonstration
local controls_subset ""
forvalues i = 1/20 {
    local controls_subset "`controls_subset' ind_time_`i'"
}

di ""
di "Running did_multiplegt_dyn with controls() approach..."
di "(Using subset of 20 industry x time interactions for speed)"
di ""

timer clear 1
timer on 1

cap noi did_multiplegt_dyn Y group time D, effects(3) placebo(2) controls(`controls_subset') cluster(group) graph_off

timer off 1
timer list 1

* Store results
if _rc == 0 {
    scalar att_controls = e(effect_1)
    scalar se_controls = e(se_effect_1)
    scalar placebo1_controls = e(placebo_1)
    scalar placebo2_controls = e(placebo_2)

    di ""
    di "RESULTS WITH CONTROLS() APPROACH:"
    di "Effect 1 (ATT): " %9.4f att_controls " (SE: " %9.4f se_controls ")"
    if placebo1_controls != . di "Placebo 1: " %9.4f placebo1_controls
    if placebo2_controls != . di "Placebo 2: " %9.4f placebo2_controls
}
else {
    di as error "Original approach failed with error code: " _rc
    scalar att_controls = .
    scalar se_controls = .
    scalar placebo1_controls = .
    scalar placebo2_controls = .
}

********************************************************************************
* 5. Test NEW approach: hdfe() option
********************************************************************************

di ""
di "=========================================="
di "TEST 2: NEW APPROACH"
di "Using hdfe() with industry and region"
di "=========================================="

di ""
di "Running did_multiplegt_dyn with hdfe() approach..."
di ""

timer clear 2
timer on 2

cap noi did_multiplegt_dyn Y group time D, effects(3) placebo(2) hdfe(industry) cluster(group) graph_off

timer off 2
timer list 2

* Store results
if _rc == 0 {
    scalar att_hdfe = e(effect_1)
    scalar se_hdfe = e(se_effect_1)
    scalar placebo1_hdfe = e(placebo_1)
    scalar placebo2_hdfe = e(placebo_2)

    di ""
    di "RESULTS WITH HDFE() APPROACH:"
    di "Effect 1 (ATT): " %9.4f att_hdfe " (SE: " %9.4f se_hdfe ")"
    if placebo1_hdfe != . di "Placebo 1: " %9.4f placebo1_hdfe
    if placebo2_hdfe != . di "Placebo 2: " %9.4f placebo2_hdfe
}
else {
    di as error "HDFE approach failed with error code: " _rc
    scalar att_hdfe = .
    scalar se_hdfe = .
    scalar placebo1_hdfe = .
    scalar placebo2_hdfe = .
}

********************************************************************************
* 6. Test with both industry and region in hdfe()
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

* Store results
if _rc == 0 {
    scalar att_hdfe2 = e(effect_1)
    scalar se_hdfe2 = e(se_effect_1)
    scalar placebo1_hdfe2 = e(placebo_1)
    scalar placebo2_hdfe2 = e(placebo_2)

    di ""
    di "RESULTS WITH HDFE(industry region) APPROACH:"
    di "Effect 1 (ATT): " %9.4f att_hdfe2 " (SE: " %9.4f se_hdfe2 ")"
    if placebo1_hdfe2 != . di "Placebo 1: " %9.4f placebo1_hdfe2
    if placebo2_hdfe2 != . di "Placebo 2: " %9.4f placebo2_hdfe2
}
else {
    di as error "HDFE approach with multiple vars failed with error code: " _rc
    scalar att_hdfe2 = .
    scalar se_hdfe2 = .
}

********************************************************************************
* 7. Baseline comparison (no controls, no hdfe)
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

* Store results
if _rc == 0 {
    scalar att_baseline = e(effect_1)
    scalar se_baseline = e(se_effect_1)
    scalar placebo1_baseline = e(placebo_1)
    scalar placebo2_baseline = e(placebo_2)

    di ""
    di "RESULTS BASELINE (no adjustments):"
    di "Effect 1 (ATT): " %9.4f att_baseline " (SE: " %9.4f se_baseline ")"
    if placebo1_baseline != . di "Placebo 1: " %9.4f placebo1_baseline
    if placebo2_baseline != . di "Placebo 2: " %9.4f placebo2_baseline
}

********************************************************************************
* 8. Summary Comparison
********************************************************************************

di ""
di "=========================================="
di "SUMMARY COMPARISON"
di "=========================================="
di ""
di "True ATT: `true_att'"
di ""
di _col(1) "Method" _col(25) "ATT" _col(35) "SE" _col(45) "Placebo 1" _col(60) "Placebo 2"
di "-----------------------------------------------------------------------"

if att_baseline != . {
    di _col(1) "Baseline" _col(22) %9.4f att_baseline _col(32) %9.4f se_baseline _col(45) %9.4f placebo1_baseline _col(60) %9.4f placebo2_baseline
}

if att_controls != . {
    di _col(1) "Controls (subset)" _col(22) %9.4f att_controls _col(32) %9.4f se_controls _col(45) %9.4f placebo1_controls _col(60) %9.4f placebo2_controls
}

if att_hdfe != . {
    di _col(1) "HDFE (industry)" _col(22) %9.4f att_hdfe _col(32) %9.4f se_hdfe _col(45) %9.4f placebo1_hdfe _col(60) %9.4f placebo2_hdfe
}

if att_hdfe2 != . {
    di _col(1) "HDFE (ind+reg)" _col(22) %9.4f att_hdfe2 _col(32) %9.4f se_hdfe2 _col(45) %9.4f placebo1_hdfe2 _col(60) %9.4f placebo2_hdfe2
}

di "-----------------------------------------------------------------------"
di ""
di "Timing comparison:"
timer list

di ""
di "=========================================="
di "TEST COMPLETED"
di "=========================================="

* Clean up
cap erase "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata/test_data_hdfe.dta"
