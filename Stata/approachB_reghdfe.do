*! ===========================================================================
*! approachB_reghdfe.do
*! ---------------------------------------------------------------------------
*! GOAL: get the de Chaisemartin-D'Haultfoeuille DID_{g,l} TWO ways on the SAME
*! data and show they coincide:
*!     (1) DATA OPERATIONS  -- the formula by hand (long-diffs and means)
*!     (2) reghdfe          -- a 2x2 two-way-FE regression on arranged data
*!
*! WHY THEY MATCH: DID_{g,l} is exactly a 2x2 difference-in-differences --
*! treated group g vs its cohort-specific controls, comparing the event period
*! (F_g-1+l) against the baseline (F_g-1). Regressing Y on D = 1{i=g}*1{t=event}
*! with unit and period fixed effects returns that 2x2 DID as the coefficient.
*! "Arranging the data" = subsetting to the right units and the right 2 periods.
*!
*! Requires reghdfe (and ftools):
*!     ssc install ftools, replace
*!     ssc install reghdfe, replace
*! Run:  do approachB_reghdfe.do
*! ===========================================================================

clear all
set more off
cap which reghdfe
if _rc {
    di as error "reghdfe not found. Run: ssc install ftools, replace ; ssc install reghdfe, replace"
    exit 198
}

* ---------------------------------------------------------------------------
* DESIGN (locals). Never-switchers get F = 9999 (a finite stand-in for +inf,
* so the comparison  F_u > event  is unambiguous). Only switchers carry S/T.
* ---------------------------------------------------------------------------
local L = 2
local switchers 1 2 4              // groups whose treatment switches at F_g
* period-one stratum D_{g,1}:
local D1_1 0
local D1_2 0
local D1_3 0
local D1_4 1
local D1_5 1
* first-switch period F_g  (9999 = never):
local F_1 2
local F_2 3
local F_3 9999
local F_4 2
local F_5 9999
* sign S_g  (+1 up, -1 down) and last-status-quo date T_g (switchers only):
local S_1  1
local S_2  1
local S_4 -1
local T_1 3
local T_2 3
local T_4 3

* ===========================================================================
* PART A -- build the panel (DETERMINISTIC, so the numbers are reproducible).
*   Y_{g,t} = alpha_g + 1.5 t + (0.37 g - 0.21 t + 0.05 g t) + injected effect
*   The 0.05*g*t piece does NOT difference out across cohort-specific controls,
*   so the DID_{g,l} are non-round -- a convincing test of the equivalence.
* ===========================================================================
clear
set obs 20
gen int g = ceil(_n / 4)
gen int t = _n - (g - 1) * 4
gen double alpha = .
replace alpha =  1.3 if g == 1
replace alpha = -2.1 if g == 2
replace alpha =  0.7 if g == 3
replace alpha =  3.4 if g == 4
replace alpha = -1.8 if g == 5
gen double Y = alpha + 1.5*t + 0.37*g - 0.21*t + 0.05*g*t
* inject signed effects at each switcher's event date (S_g * true effect):
replace Y = Y + 4 if g == 1 & t == 2      // g1 (S=+1), l=1, event t=2
replace Y = Y + 7 if g == 1 & t == 3      // g1 (S=+1), l=2, event t=3
replace Y = Y + 4 if g == 2 & t == 3      // g2 (S=+1), l=1, event t=3
replace Y = Y - 4 if g == 4 & t == 2      // g4 (S=-1), l=1, event t=2
replace Y = Y - 7 if g == 4 & t == 3      // g4 (S=-1), l=2, event t=3
drop alpha
tempfile panel
save `panel'

* ===========================================================================
* PART B -- for every estimable (g,l): compute DID_{g,l} by DATA OPS, and by
* reghdfe on the arranged 2x2 subset; store both. Also build a STACKED dataset
* and a crosswalk for PART C.
* ===========================================================================
tempfile resfile stack
postfile myh int g int l double did_ops double did_reg using "`resfile'", replace

local e = 0                        // experiment counter for the stack
local built = 0

foreach g of local switchers {
    forvalues l = 1/`L' {
        local pb = `F_`g'' - 1          // baseline period
        local pe = `F_`g'' - 1 + `l'    // event period
        * estimable iff event period does not exceed T_g
        if (`pe' <= `T_`g'') {

            * --- cohort-specific control set: same stratum, still status-quo ---
            local ctrlcond ""
            forvalues u = 1/5 {
                if (`u' != `g') & (`D1_`u'' == `D1_`g'') & (`F_`u'' > `pe') {
                    if "`ctrlcond'" == "" local ctrlcond "g==`u'"
                    else                   local ctrlcond "`ctrlcond' | g==`u'"
                }
            }

            * ============ (1) DATA OPERATIONS ============
            * treated long-difference  Y_{g,pe} - Y_{g,pb}
            qui sum Y if g==`g' & t==`pe', meanonly
            local tp = r(mean)
            qui sum Y if g==`g' & t==`pb', meanonly
            local tb = r(mean)
            * control mean long-difference  mean_u(Y_{u,pe}) - mean_u(Y_{u,pb})
            qui sum Y if (`ctrlcond') & t==`pe', meanonly
            local cp = r(mean)
            qui sum Y if (`ctrlcond') & t==`pb', meanonly
            local cb = r(mean)
            local did_ops = (`tp' - `tb') - (`cp' - `cb')

            * ============ (2) reghdfe on the arranged 2x2 ============
            preserve
                keep if g==`g' | (`ctrlcond')      // treated + its controls
                keep if t==`pb' | t==`pe'          // baseline + event period
                gen byte treat = (g==`g')
                gen byte post  = (t==`pe')
                gen byte D     = treat*post        // the DID treatment dummy
                qui reghdfe Y D, absorb(g t)       // unit + period fixed effects
                local did_reg = _b[D]
            restore

            post myh (`g') (`l') (`did_ops') (`did_reg')

            * ---- accumulate the STACKED dataset + crosswalk (for PART C) ----
            local ++e
            local e_g`e'   = `g'
            local e_l`e'   = `l'
            local e_S`e'   = `S_`g''
            local e_did`e' = `did_ops'
            preserve
                keep if g==`g' | (`ctrlcond')
                keep if t==`pb' | t==`pe'
                gen int  e     = `e'
                gen byte treat = (g==`g')
                gen byte post  = (t==`pe')
                gen byte D     = treat*post
                if `built'==0 {
                    save `stack', replace
                    local built 1
                }
                else {
                    append using `stack'
                    save `stack', replace
                }
            restore
        }
    }
}
postclose myh

* counts N_l (number of estimable groups at each horizon)
forvalues l = 1/`L' {
    local N`l' = 0
    forvalues ee = 1/`e' {
        if `e_l`ee'' == `l' local N`l' = `N`l'' + 1
    }
}

* ---- show the side-by-side comparison of DID_{g,l} ----
use "`resfile'", clear
gen double diff = did_ops - did_reg
format did_ops did_reg diff %9.4f
di _newline as txt "==== DID_{g,l}: data operations vs reghdfe (2x2) ===="
list g l did_ops did_reg diff, noobs abbrev(12) sepby(g)



* ===========================================================================
* PART C (bonus) -- ONE stacked reghdfe recovers ALL DID_{g,l} at once, as the
* experiment-specific slopes on D, with experiment-specific unit and period FE.
* Then lincom forms delta_l = (1/N_l) sum_g S_g DID_{g,l} (matching matrix A's
* rows) and reports a clustered SE. Compare to the data-ops delta_l.
* ===========================================================================
use `stack', clear
di _newline as txt "==== Stacked regression: one reghdfe for all DID_{g,l} ===="
reghdfe Y ibn.e#c.D, absorb(g#e t#e) vce(cluster g)

* per-experiment coefficient (= DID_{g,l}) vs data-ops value
di _newline as txt "experiment  (g,l)   stacked_b      did_ops"
forvalues ee = 1/`e' {
    local b = _b[`ee'.e#c.D]
    di as res "   e=`ee'    " %1.0f `e_g`ee'' "," %1.0f `e_l`ee'' ///
        "    " %9.4f `b' "   " %9.4f `e_did`ee''
}

* delta_l from the regression via lincom (signed, N_l-normalized)
di _newline as txt "==== delta_l from the stacked regression (lincom) ===="
forvalues l = 1/`L' {
    local expr ""
    forvalues ee = 1/`e' {
        if `e_l`ee'' == `l' {
            local s = `e_S`ee''
            if "`expr'" == "" local expr "(`s')*_b[`ee'.e#c.D]"
            else               local expr "`expr' + (`s')*_b[`ee'.e#c.D]"
        }
    }
    di as txt "delta_`l':"
    lincom (`expr') / `N`l''
}

* delta_l by pure data operations, for comparison
di _newline as txt "==== delta_l by data operations ===="
forvalues l = 1/`L' {
    local num = 0
    forvalues ee = 1/`e' {
        if `e_l`ee'' == `l' local num = `num' + (`e_S`ee'') * (`e_did`ee'')
    }
    di as res "  delta_`l' (data ops) = " %9.4f (`num' / `N`l'')
}

di _newline "Done."
