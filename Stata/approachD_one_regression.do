*! ===========================================================================
*! approachD_one_regression.do
*! ---------------------------------------------------------------------------
*! Same estimand as approachC (and as `did_multiplegt_dyn Y g t D, effects(2)`),
*! but recovered from a SINGLE reghdfe call -- no loop over (g,l), no hand
*! aggregation. The trick is stacked DID with two ingredients:
*!
*!   (i)  cell-by-unit and cell-by-time fixed effects, so each cell's within
*!        variation is a clean 2x2 DID.
*!   (ii) observation weights w = 1 for the switcher cohort and w = 1/K_c for
*!        each control in cell c. This makes every cell's weighted SSR of the
*!        demeaned regressor equal (= 1/4), so the stacked-OLS aggregator
*!        across cells collapses to the simple cross-switcher average that
*!        dCDH uses:  delta_l = (1/N_l) * sum_g S_g * DID_{g,l}.
*!
*! Without those weights, cells with more controls would get larger FWL
*! weight and the coefficient would NOT match did_multiplegt_dyn.
*!
*! Requires: reghdfe, ftools, did_multiplegt_dyn (already used by approachC).
*! ===========================================================================

clear all
set more off

* ===========================================================================
* PART A -- identical panel to approachC (same Y, same D path)
* ===========================================================================
clear
set obs 20
gen int g = ceil(_n / 4)
gen int t = _n - (g - 1) * 4

gen byte D = .
replace D = 0 if g==1
replace D = 1 if g==1 & t>=2
replace D = 0 if g==2
replace D = 1 if g==2 & t>=3
replace D = 0 if g==3
replace D = 1 if g==4 & t==1
replace D = 0 if g==4 & t>=2
replace D = 1 if g==5

gen double alpha = .
replace alpha =  1.3 if g==1
replace alpha = -2.1 if g==2
replace alpha =  0.7 if g==3
replace alpha =  3.4 if g==4
replace alpha = -1.8 if g==5
gen double Y = alpha + 1.5*t + 0.37*g - 0.21*t + 0.05*g*t
replace Y = Y + 4 if g==1 & t==2
replace Y = Y + 7 if g==1 & t==3
replace Y = Y + 4 if g==2 & t==3
replace Y = Y - 4 if g==4 & t==2
replace Y = Y - 7 if g==4 & t==3
drop alpha
tempfile panel
save `panel'

* ===========================================================================
* PART B -- derive design (D1g, Fg, Sg, Tg) from D, same as approachC
* ===========================================================================
qui sum t
local Tmax = r(max)

bysort g (t): gen byte D1g = D[1]
gen byte chg = (D != D1g)
bysort g (t): egen Fg = min(cond(chg==1, t, .))
replace Fg = `Tmax' + 1 if missing(Fg)
gen double _DatF = D if t==Fg
bysort g (t): egen DatFg = max(_DatF)
gen Sg = sign(DatFg - D1g)
bysort D1g: egen maxFstr = max(Fg)
gen Tg = maxFstr - 1
drop _DatF chg

forvalues gg = 1/5 {
    qui sum D1g if g==`gg', meanonly
    local D1_`gg' = r(mean)
    qui sum Fg  if g==`gg', meanonly
    local F_`gg' = r(mean)
    qui sum Tg  if g==`gg', meanonly
    local T_`gg' = r(mean)
    qui sum Sg  if g==`gg', meanonly
    local S_`gg' = r(mean)
}
local switchers ""
forvalues gg = 1/5 {
    if `F_`gg'' <= `Tmax' local switchers "`switchers' `gg'"
}

* ===========================================================================
* PART C -- BUILD THE STACKED DATASET
*   one row per (cell, original obs); cells are estimable (g_switcher, l)
*   pairs. Each cell carries its switcher cohort and valid controls observed
*   at t in {pb, pe}.
* ===========================================================================
local L = 2

tempfile stack
clear
set obs 0
gen int    cell_g = .
gen int    cell_l = .
gen int    g      = .
gen int    t      = .
gen double Y      = .
gen byte   treat  = .
gen byte   post   = .
gen byte   S      = .
gen double Kc     = .
save `stack', replace

foreach g of local switchers {
    forvalues l = 1/`L' {
        local pb = `F_`g'' - 1
        local pe = `F_`g'' - 1 + `l'
        if (`pe' <= `T_`g'') {

            * cohort-specific controls: same stratum, not yet switched at pe
            local ctrlcond ""
            local Kc = 0
            forvalues u = 1/5 {
                if (`u' != `g') & (`D1_`u'' == `D1_`g'') & (`F_`u'' > `pe') {
                    if "`ctrlcond'" == "" local ctrlcond "g==`u'"
                    else                   local ctrlcond "`ctrlcond' | g==`u'"
                    local Kc = `Kc' + 1
                }
            }

            preserve
                use `panel', clear
                keep if (g==`g') | (`ctrlcond')
                keep if t==`pb' | t==`pe'
                gen int    cell_g = `g'
                gen int    cell_l = `l'
                gen byte   treat  = (g==`g')
                gen byte   post   = (t==`pe')
                gen byte   S      = `S_`g''
                gen double Kc     = `Kc'
                keep cell_g cell_l g t Y treat post S Kc
                append using `stack'
                save `stack', replace
            restore
        }
    }
}

* ===========================================================================
* PART D -- ONE regression on the stacked data
* ===========================================================================
use `stack', clear

egen long cell_id      = group(cell_g cell_l)
egen long cell_unit_id = group(cell_id g)
egen long cell_time_id = group(cell_id t)

* equal-cell-weight trick: switcher gets w=1, each control gets w=1/K_c.
gen double w = cond(treat==1, 1, 1/Kc)

* signed treatment * post: flips sign for down-switchers (S = -1)
gen double SxTxP = S * treat * post

di _newline as txt "==== stacked-DID one-regression recovery of delta_l ===="
reghdfe Y c.SxTxP#i.cell_l [pw=w], absorb(cell_unit_id cell_time_id)

di _newline as txt "Read off:"
di as res "  coef on  1.cell_l#c.SxTxP  =  delta_1"
di as res "  coef on  2.cell_l#c.SxTxP  =  delta_2"

* ===========================================================================
* PART E -- sanity check against the official command
* ===========================================================================
use `panel', clear
di _newline as txt "==== did_multiplegt_dyn (target) ===="
did_multiplegt_dyn Y g t D, effects(2)

di _newline as txt "Expected: Effect_1 = 3.9750, Effect_2 = 4.6000."
di "Done."
