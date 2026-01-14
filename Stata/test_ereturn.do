* Test what e() values are returned by did_multiplegt_dyn
clear all
set seed 12345

* Load local version
adopath ++ "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata"

* Create simple test data
set obs 500
gen group = ceil(_n / 10)
bysort group: gen time = _n
bysort group: gen industry = ceil(runiform() * 5) if _n == 1
bysort group: replace industry = industry[1]
bysort group: gen treat_time = ceil(runiform() * 12) if _n == 1
bysort group: replace treat_time = treat_time[1]
replace treat_time = . if treat_time > 10
gen D = (time >= treat_time) & treat_time != .
gen Y = rnormal(0, 1) + 3 * D

* Run command
did_multiplegt_dyn Y group time D, effects(2) placebo(1) hdfe(industry) cluster(group) graph_off

* Check what e() returns
di "=== E() SCALARS ==="
ereturn list

di ""
di "=== SPECIFIC VALUES ==="
di "e(effect_1) = " e(effect_1)
di "e(se_effect_1) = " e(se_effect_1)
di "e(N_effect_1) = " e(N_effect_1)

* Check if effects are in a matrix
di ""
di "=== MATRICES ==="
cap matrix list e(results)
cap matrix list e(effects)
