********************************************************************************
* Debug file for hdfe() option
* This will help identify why HDFE is not working
********************************************************************************

clear all
set seed 12345
set more off

* Load local version
cap program drop did_multiplegt_dyn
cap program drop did_multiplegt_dyn_core_new
adopath ++ "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata"
which did_multiplegt_dyn

* Check reghdfe
cap which reghdfe
if _rc {
    ssc install reghdfe, replace
}

********************************************************************************
* Create simple test data
********************************************************************************

local N_groups = 200
local T_periods = 8
local N_industries = 20

set obs `=`N_groups' * `T_periods''
gen group = ceil(_n / `T_periods')
bysort group: gen time = _n

* Industry - time invariant within groups
bysort group: gen industry = ceil(runiform() * `N_industries') if _n == 1
bysort group: replace industry = industry[1]

* Treatment timing
bysort group: gen treat_time = ceil(runiform() * (`T_periods' + 2)) if _n == 1
bysort group: replace treat_time = treat_time[1]
replace treat_time = . if treat_time > `T_periods'

gen D = (time >= treat_time) & treat_time != .

* Generate outcome with INDUSTRY-SPECIFIC TRENDS
* This is key: industry affects the TREND, not just the level
gen industry_trend = 0
forvalues i = 1/`N_industries' {
    * Each industry has a different trend coefficient
    replace industry_trend = (`i' - 10) * 0.3 if industry == `i'
}

* Group FE
bysort group: gen group_fe = rnormal(0, 2) if _n == 1
bysort group: replace group_fe = group_fe[1]

* True ATT
local true_att = 3.0

* Y = group_fe + industry_trend * time + ATT * D + noise
* The industry_trend * time creates differential trends by industry
gen Y = group_fe + industry_trend * time + `true_att' * D + rnormal(0, 1)

di ""
di "Data created with industry-specific trends"
di "True ATT: `true_att'"
di ""

* Check the industry trend effect
reg Y i.industry##c.time if D == 0
di "Industry trends are present in the data"

********************************************************************************
* Manual check: what does first-differenced data look like?
********************************************************************************

xtset group time
gen diff_Y = D.Y

di ""
di "=== CHECKING FIRST DIFFERENCES ==="
di ""

* Check if there's variation in diff_Y across industries
bysort industry: sum diff_Y
reg diff_Y i.industry if D == 0

di ""
di "If industry coefficients are significant, HDFE should matter"
di ""

********************************************************************************
* Test 1: Baseline (no HDFE)
********************************************************************************

di ""
di "=========================================="
di "TEST 1: BASELINE"
di "=========================================="

did_multiplegt_dyn Y group time D, effects(3) placebo(2) cluster(group) graph_off

********************************************************************************
* Test 2: HDFE with industry
********************************************************************************

di ""
di "=========================================="
di "TEST 2: HDFE(industry)"
di "=========================================="

did_multiplegt_dyn Y group time D, effects(3) placebo(2) hdfe(industry) cluster(group) graph_off

********************************************************************************
* Test 3: Manual HDFE check - what reghdfe actually does
********************************************************************************

di ""
di "=========================================="
di "TEST 3: MANUAL REGHDFE CHECK"
di "=========================================="

* Recreate what the ado should be doing
preserve

* Create first difference
xtset group time
cap drop diff_y_manual
gen diff_y_manual = D.Y

* Check original variation
sum diff_y_manual
local orig_sd = r(sd)
di "Original diff_Y SD: `orig_sd'"

* Run reghdfe to residualize
reghdfe diff_y_manual, absorb(i.industry i.time) residuals(resid_hdfe)

* Check residual variation
sum resid_hdfe
local resid_sd = r(sd)
di "Residualized diff_Y SD: `resid_sd'"

di ""
di "If SD changed significantly, HDFE is absorbing something"
di "Original SD: `orig_sd'"
di "Residual SD: `resid_sd'"
di "Difference: " %6.4f (`orig_sd' - `resid_sd')

restore

********************************************************************************
* Test 4: trends_nonparam for comparison
********************************************************************************

di ""
di "=========================================="
di "TEST 4: trends_nonparam(industry) for comparison"
di "=========================================="

did_multiplegt_dyn Y group time D, effects(3) placebo(2) trends_nonparam(industry) cluster(group) graph_off

di ""
di "=========================================="
di "DEBUG COMPLETE"
di "=========================================="
