********************************************************************************
* Monte Carlo Simulation: Numerical Precision Test
* controls(FE×time dummies) vs hdfe()
* 500 simulations, 1000 observations each, 5 FE categories
********************************************************************************

clear all
set seed 98765
set more off

* Load local version
cap program drop did_multiplegt_dyn
cap program drop did_multiplegt_dyn_core_new
adopath ++ "/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/Stata"
which did_multiplegt_dyn

********************************************************************************
* Simulation Parameters
********************************************************************************

local N_sims = 500
local N_obs = 1000
local N_industries = 5
local T_periods = 10
local N_groups = ceil(`N_obs' / `T_periods')
local true_att = 3.0

di ""
di "=============================================="
di "MONTE CARLO SIMULATION PARAMETERS"
di "=============================================="
di "Number of simulations: `N_sims'"
di "Observations per simulation: `N_obs'"
di "Groups: `N_groups'"
di "Time periods: `T_periods'"
di "FE categories (industries): `N_industries'"
di "True ATT: `true_att'"
di "=============================================="
di ""

********************************************************************************
* Create matrix to store results
********************************************************************************

* Store differences: hdfe - controls for each metric
matrix results_diff = J(`N_sims', 10, .)
matrix colnames results_diff = diff_eff1 diff_eff2 diff_eff3 diff_plac1 diff_plac2 ///
                               diff_se_eff1 diff_se_eff2 diff_se_eff3 diff_se_plac1 diff_se_plac2

* Also store the actual estimates for reference
matrix results_controls = J(`N_sims', 10, .)
matrix results_hdfe = J(`N_sims', 10, .)

********************************************************************************
* Run Simulations
********************************************************************************

timer clear 1
timer on 1

forvalues sim = 1/`N_sims' {

    * Progress indicator
    if mod(`sim', 50) == 0 | `sim' == 1 {
        di "Simulation `sim' of `N_sims'..."
    }

    quietly {

        *-----------------------------------------------------------------------
        * Generate synthetic data
        *-----------------------------------------------------------------------

        clear
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
        gen industry_trend = (industry - `N_industries'/2) * 0.5

        * Group FE
        bysort group: gen group_fe = rnormal(0, 2) if _n == 1
        bysort group: replace group_fe = group_fe[1]

        * Y = group_fe + industry_trend * time + ATT * D + noise
        gen Y = group_fe + industry_trend * time + `true_att' * D + rnormal(0, 1)

        *-----------------------------------------------------------------------
        * Create industry×time dummies for controls() approach
        *-----------------------------------------------------------------------

        tab industry, gen(ind_)

        forvalues i = 1/`N_industries' {
            gen ind_time_`i' = ind_`i' * time
        }

        local controls_list ""
        forvalues i = 1/`N_industries' {
            local controls_list "`controls_list' ind_time_`i'"
        }

        *-----------------------------------------------------------------------
        * Run controls() version
        *-----------------------------------------------------------------------

        cap noi did_multiplegt_dyn Y group time D, effects(3) placebo(2) ///
            controls(`controls_list') cluster(group) graph_off

        if _rc == 0 {
            matrix results_controls[`sim', 1] = e(Effect_1)
            matrix results_controls[`sim', 2] = e(Effect_2)
            matrix results_controls[`sim', 3] = e(Effect_3)
            matrix results_controls[`sim', 4] = e(Placebo_1)
            matrix results_controls[`sim', 5] = e(Placebo_2)
            matrix results_controls[`sim', 6] = e(se_effect_1)
            matrix results_controls[`sim', 7] = e(se_effect_2)
            matrix results_controls[`sim', 8] = e(se_effect_3)
            matrix results_controls[`sim', 9] = e(se_placebo_1)
            matrix results_controls[`sim', 10] = e(se_placebo_2)
        }

        *-----------------------------------------------------------------------
        * Run hdfe() version
        *-----------------------------------------------------------------------

        cap noi did_multiplegt_dyn Y group time D, effects(3) placebo(2) ///
            hdfe(industry) cluster(group) graph_off

        if _rc == 0 {
            matrix results_hdfe[`sim', 1] = e(Effect_1)
            matrix results_hdfe[`sim', 2] = e(Effect_2)
            matrix results_hdfe[`sim', 3] = e(Effect_3)
            matrix results_hdfe[`sim', 4] = e(Placebo_1)
            matrix results_hdfe[`sim', 5] = e(Placebo_2)
            matrix results_hdfe[`sim', 6] = e(se_effect_1)
            matrix results_hdfe[`sim', 7] = e(se_effect_2)
            matrix results_hdfe[`sim', 8] = e(se_effect_3)
            matrix results_hdfe[`sim', 9] = e(se_placebo_1)
            matrix results_hdfe[`sim', 10] = e(se_placebo_2)
        }

        *-----------------------------------------------------------------------
        * Compute differences (hdfe - controls)
        *-----------------------------------------------------------------------

        forvalues j = 1/10 {
            matrix results_diff[`sim', `j'] = results_hdfe[`sim', `j'] - results_controls[`sim', `j']
        }

    }

}

timer off 1
timer list 1

di ""
di "=============================================="
di "SIMULATION COMPLETE"
di "=============================================="
di ""

********************************************************************************
* Convert matrices to dataset for analysis
********************************************************************************

clear
svmat results_diff, names(col)

label var diff_eff1 "Difference in Effect_1 (HDFE - Controls)"
label var diff_eff2 "Difference in Effect_2 (HDFE - Controls)"
label var diff_eff3 "Difference in Effect_3 (HDFE - Controls)"
label var diff_plac1 "Difference in Placebo_1 (HDFE - Controls)"
label var diff_plac2 "Difference in Placebo_2 (HDFE - Controls)"
label var diff_se_eff1 "Difference in SE_Effect_1 (HDFE - Controls)"
label var diff_se_eff2 "Difference in SE_Effect_2 (HDFE - Controls)"
label var diff_se_eff3 "Difference in SE_Effect_3 (HDFE - Controls)"
label var diff_se_plac1 "Difference in SE_Placebo_1 (HDFE - Controls)"
label var diff_se_plac2 "Difference in SE_Placebo_2 (HDFE - Controls)"

gen sim_id = _n

********************************************************************************
* Summary Statistics of Differences
********************************************************************************

di ""
di "=============================================="
di "SUMMARY STATISTICS OF DIFFERENCES (HDFE - Controls)"
di "=============================================="
di ""
di "If methods are numerically equivalent, all differences should be ~0"
di ""

di "--- POINT ESTIMATES ---"
sum diff_eff1 diff_eff2 diff_eff3 diff_plac1 diff_plac2, detail

di ""
di "--- STANDARD ERRORS ---"
sum diff_se_eff1 diff_se_eff2 diff_se_eff3 diff_se_plac1 diff_se_plac2, detail

********************************************************************************
* Detailed Statistics Table
********************************************************************************

di ""
di "=============================================="
di "DETAILED STATISTICS TABLE"
di "=============================================="
di ""
di "Variable          |    Mean         SD          Min         Max      |N|"
di "------------------|--------------------------------------------------|--|"

foreach var in diff_eff1 diff_eff2 diff_eff3 diff_plac1 diff_plac2 {
    sum `var', meanonly
    local n = r(N)
    sum `var'
    di "`var'" _col(20) %12.8f r(mean) %12.8f r(sd) %12.8f r(min) %12.8f r(max) _col(73) "`n'"
}

di ""
foreach var in diff_se_eff1 diff_se_eff2 diff_se_eff3 diff_se_plac1 diff_se_plac2 {
    sum `var', meanonly
    local n = r(N)
    sum `var'
    di "`var'" _col(20) %12.8f r(mean) %12.8f r(sd) %12.8f r(min) %12.8f r(max) _col(73) "`n'"
}

********************************************************************************
* Generate Histograms
********************************************************************************

di ""
di "=============================================="
di "GENERATING HISTOGRAMS"
di "=============================================="

set scheme s2color

* --- EFFECTS ---
histogram diff_eff1, title("Effect_1: HDFE - Controls") ///
    xtitle("Difference") ytitle("Frequency") ///
    note("N=`N_sims' sims, `N_obs' obs, `N_industries' FE categories") ///
    color(navy%70) xline(0, lcolor(red)) name(hist_eff1, replace)
graph export "histogram_diff_effect1.png", replace width(1200)

histogram diff_eff2, title("Effect_2: HDFE - Controls") ///
    xtitle("Difference") ytitle("Frequency") ///
    note("N=`N_sims' sims, `N_obs' obs, `N_industries' FE categories") ///
    color(navy%70) xline(0, lcolor(red)) name(hist_eff2, replace)
graph export "histogram_diff_effect2.png", replace width(1200)

histogram diff_eff3, title("Effect_3: HDFE - Controls") ///
    xtitle("Difference") ytitle("Frequency") ///
    note("N=`N_sims' sims, `N_obs' obs, `N_industries' FE categories") ///
    color(navy%70) xline(0, lcolor(red)) name(hist_eff3, replace)
graph export "histogram_diff_effect3.png", replace width(1200)

* --- PLACEBOS ---
histogram diff_plac1, title("Placebo_1: HDFE - Controls") ///
    xtitle("Difference") ytitle("Frequency") ///
    note("N=`N_sims' sims, `N_obs' obs, `N_industries' FE categories") ///
    color(forest_green%70) xline(0, lcolor(red)) name(hist_plac1, replace)
graph export "histogram_diff_placebo1.png", replace width(1200)

histogram diff_plac2, title("Placebo_2: HDFE - Controls") ///
    xtitle("Difference") ytitle("Frequency") ///
    note("N=`N_sims' sims, `N_obs' obs, `N_industries' FE categories") ///
    color(forest_green%70) xline(0, lcolor(red)) name(hist_plac2, replace)
graph export "histogram_diff_placebo2.png", replace width(1200)

* --- SE EFFECTS ---
histogram diff_se_eff1, title("SE(Effect_1): HDFE - Controls") ///
    xtitle("Difference") ytitle("Frequency") ///
    note("N=`N_sims' sims, `N_obs' obs, `N_industries' FE categories") ///
    color(maroon%70) xline(0, lcolor(red)) name(hist_se_eff1, replace)
graph export "histogram_diff_se_effect1.png", replace width(1200)

histogram diff_se_eff2, title("SE(Effect_2): HDFE - Controls") ///
    xtitle("Difference") ytitle("Frequency") ///
    note("N=`N_sims' sims, `N_obs' obs, `N_industries' FE categories") ///
    color(maroon%70) xline(0, lcolor(red)) name(hist_se_eff2, replace)
graph export "histogram_diff_se_effect2.png", replace width(1200)

histogram diff_se_eff3, title("SE(Effect_3): HDFE - Controls") ///
    xtitle("Difference") ytitle("Frequency") ///
    note("N=`N_sims' sims, `N_obs' obs, `N_industries' FE categories") ///
    color(maroon%70) xline(0, lcolor(red)) name(hist_se_eff3, replace)
graph export "histogram_diff_se_effect3.png", replace width(1200)

* --- SE PLACEBOS ---
histogram diff_se_plac1, title("SE(Placebo_1): HDFE - Controls") ///
    xtitle("Difference") ytitle("Frequency") ///
    note("N=`N_sims' sims, `N_obs' obs, `N_industries' FE categories") ///
    color(orange%70) xline(0, lcolor(red)) name(hist_se_plac1, replace)
graph export "histogram_diff_se_plac1.png", replace width(1200)

histogram diff_se_plac2, title("SE(Placebo_2): HDFE - Controls") ///
    xtitle("Difference") ytitle("Frequency") ///
    note("N=`N_sims' sims, `N_obs' obs, `N_industries' FE categories") ///
    color(orange%70) xline(0, lcolor(red)) name(hist_se_plac2, replace)
graph export "histogram_diff_se_plac2.png", replace width(1200)

********************************************************************************
* Combined Panel Graphs
********************************************************************************

graph combine hist_eff1 hist_eff2 hist_eff3, rows(1) cols(3) ///
    title("Point Estimates: HDFE vs Controls") ///
    note("Red line = 0 (perfect equivalence)") ///
    name(combined_effects, replace)
graph export "histogram_combined_effects.png", replace width(1800)

graph combine hist_plac1 hist_plac2, rows(1) cols(2) ///
    title("Placebos: HDFE vs Controls") ///
    note("Red line = 0 (perfect equivalence)") ///
    name(combined_placebos, replace)
graph export "histogram_combined_placebos.png", replace width(1400)

graph combine hist_se_eff1 hist_se_eff2 hist_se_eff3, rows(1) cols(3) ///
    title("SE of Effects: HDFE vs Controls") ///
    note("Red line = 0 (perfect equivalence)") ///
    name(combined_se_effects, replace)
graph export "histogram_combined_se_effects.png", replace width(1800)

graph combine hist_se_plac1 hist_se_plac2, rows(1) cols(2) ///
    title("SE of Placebos: HDFE vs Controls") ///
    note("Red line = 0 (perfect equivalence)") ///
    name(combined_se_placebos, replace)
graph export "histogram_combined_se_placebos.png", replace width(1400)

* Master graphs
graph combine combined_effects combined_placebos, rows(2) cols(1) ///
    title("Monte Carlo: Point Estimate Differences") ///
    subtitle("HDFE() vs Controls() - `N_sims' simulations") ///
    name(master_point, replace)
graph export "histogram_master_point_estimates.png", replace width(1800) height(1200)

graph combine combined_se_effects combined_se_placebos, rows(2) cols(1) ///
    title("Monte Carlo: Standard Error Differences") ///
    subtitle("HDFE() vs Controls() - `N_sims' simulations") ///
    name(master_se, replace)
graph export "histogram_master_se.png", replace width(1800) height(1200)

********************************************************************************
* Save Results
********************************************************************************

save "monte_carlo_results.dta", replace

di ""
di "=============================================="
di "FINAL SUMMARY"
di "=============================================="
di ""
di "Simulations: `N_sims' | Obs: `N_obs' | FE categories: `N_industries'"
di ""

di "Mean absolute differences (should be ~0 if equivalent):"
foreach var in diff_eff1 diff_eff2 diff_eff3 diff_plac1 diff_plac2 {
    qui sum `var'
    di "  `var': mean=" %14.12f r(mean) "  max_abs=" %14.12f max(abs(r(min)), abs(r(max)))
}

di ""
di "SE differences:"
foreach var in diff_se_eff1 diff_se_eff2 diff_se_eff3 diff_se_plac1 diff_se_plac2 {
    qui sum `var'
    di "  `var': mean=" %14.12f r(mean) "  max_abs=" %14.12f max(abs(r(min)), abs(r(max)))
}

di ""
di "=============================================="
di "MONTE CARLO COMPLETE"
di "=============================================="
