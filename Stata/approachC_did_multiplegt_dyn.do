*! ===========================================================================
*! approachC_did_multiplegt_dyn.do
*! ---------------------------------------------------------------------------
*! Reconcile THREE ways of getting the dCDH event-study effects on ONE dataset:
*!     (1) DATA OPERATIONS  -- the DID_{g,l} formula by hand, aggregated to d_l
*!     (2) reghdfe          -- the same DID_{g,l} as 2x2 two-way-FE regressions
*!     (3) did_multiplegt_dyn Y g t D, effects(2)   -- the official command
*!
*! KEY FIX vs the earlier files: we now DERIVE the design (D_{g,1}, F_g, S_g,
*! T_g) FROM the treatment path D, exactly as did_multiplegt_dyn does, instead
*! of hardcoding it. The earlier hand value T_g=3 was wrong: each stratum
*! contains a never-switcher (a valid control through the last period), so
*! T_g = T = 4. That makes (g=2, l=2) estimable, raising N_2 from 2 to 3 and
*! giving delta_2 = 4.6 (matching the command), not 6.95.
*!
*! Requires: reghdfe, ftools, and did_multiplegt_dyn:
*!     ssc install ftools, replace
*!     ssc install reghdfe, replace
*!     ssc install did_multiplegt_dyn, replace
*! ===========================================================================

clear all
set more off

* ===========================================================================
* PART A -- panel with the actual TREATMENT PATH D and the SAME deterministic Y
* as before. Treatment trajectories:
*   g1: 0 1 1 1 (up @2)   g2: 0 0 1 1 (up @3)   g3: 0 0 0 0 (never)
*   g4: 1 0 0 0 (down @2) g5: 1 1 1 1 (never)
* ===========================================================================
clear
set obs 20
gen int g = ceil(_n / 4)
gen int t = _n - (g - 1) * 4

* treatment path D_{g,t}
gen byte D = .
replace D = 0 if g==1
replace D = 1 if g==1 & t>=2
replace D = 0 if g==2
replace D = 1 if g==2 & t>=3
replace D = 0 if g==3
replace D = 1 if g==4 & t==1
replace D = 0 if g==4 & t>=2
replace D = 1 if g==5

* outcome (identical construction to the earlier files)
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
* PART B -- DERIVE the design from D (no hardcoding -> no T_g bug).
*   D1g  = D in period 1
*   Fg   = first period D differs from D1g (never -> T+1)
*   Sg   = sign of the jump at Fg (+1 up, -1 down; missing for never)
*   Tg   = (max over same-stratum Fg) - 1   [never-switchers carry F=T+1]
* ===========================================================================
qui sum t
local Tmax = r(max)

bysort g (t): gen byte D1g = D[1]                 // period-one treatment
gen byte chg = (D != D1g)
bysort g (t): egen Fg = min(cond(chg==1, t, .))   // first switch period
replace Fg = `Tmax' + 1 if missing(Fg)            // never-switchers -> T+1
gen double _DatF = D if t==Fg
bysort g (t): egen DatFg = max(_DatF)             // treatment value at Fg
gen Sg = sign(DatFg - D1g)                        // +1 / -1 / . (never)
bysort D1g: egen maxFstr = max(Fg)
gen Tg = maxFstr - 1                              // last period a control exists
drop _DatF chg

* pull the derived quantities into locals
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
di _newline as txt "Derived design (from D):"
foreach gg of local switchers {
    di as res "  g=`gg':  D1=`D1_`gg''  F=`F_`gg''  S=`S_`gg''  T=`T_`gg''"
}
di as txt "  switchers: `switchers'    (note T=`Tmax' for all -> (2,2) IS estimable)"

* ===========================================================================
* PART C -- DID_{g,l} by data ops AND by reghdfe; then aggregate to delta_l.
* ===========================================================================
tempfile resfile
postfile myh int g int l double did_ops double did_reg using "`resfile'", replace

* accumulators for delta_l: signed sum and count, per horizon
local L = 2
forvalues l = 1/`L' {
    local num`l' = 0
    local N`l' = 0
}

foreach g of local switchers {
    forvalues l = 1/`L' {
        local pb = `F_`g'' - 1
        local pe = `F_`g'' - 1 + `l'
        if (`pe' <= `T_`g'') {                        // estimable

            * cohort-specific controls: same stratum, not yet switched at pe
            local ctrlcond ""
            forvalues u = 1/5 {
                if (`u' != `g') & (`D1_`u'' == `D1_`g'') & (`F_`u'' > `pe') {
                    if "`ctrlcond'" == "" local ctrlcond "g==`u'"
                    else                   local ctrlcond "`ctrlcond' | g==`u'"
                }
            }

            * (1) data operations
            qui sum Y if g==`g' & t==`pe', meanonly
            local tp = r(mean)
            qui sum Y if g==`g' & t==`pb', meanonly
            local tb = r(mean)
            qui sum Y if (`ctrlcond') & t==`pe', meanonly
            local cp = r(mean)
            qui sum Y if (`ctrlcond') & t==`pb', meanonly
            local cb = r(mean)
            local did_ops = (`tp' - `tb') - (`cp' - `cb')

            * (2) reghdfe on the arranged 2x2
            preserve
                keep if g==`g' | (`ctrlcond')
                keep if t==`pb' | t==`pe'
                gen byte treat = (g==`g')
                gen byte post  = (t==`pe')
                gen byte Dd    = treat*post
                qui reghdfe Y Dd, absorb(g t)
                local did_reg = _b[Dd]
            restore

            post myh (`g') (`l') (`did_ops') (`did_reg')

            * accumulate signed average for delta_l
            local num`l' = `num`l'' + (`S_`g'') * (`did_ops')
            local N`l'   = `N`l'' + 1
        }
    }
}
postclose myh

* show DID_{g,l}: data ops vs reghdfe
use "`resfile'", clear
gen double diff = did_ops - did_reg
format did_ops did_reg diff %9.4f
di _newline as txt "==== DID_{g,l}: data operations vs reghdfe ===="
list g l did_ops did_reg diff, noobs abbrev(12) sepby(g)

* our delta_l (data ops)
di _newline as txt "==== delta_l (our data operations) ===="
forvalues l = 1/`L' {
    di as res "  delta_`l' = " %9.4f (`num`l'' / `N`l'') "   (N_`l' = `N`l'')"
}

* ===========================================================================
* PART D -- the official command on the ORIGINAL data.
*   Effect_1, Effect_2 in its table should equal our delta_1, delta_2.
* ===========================================================================
use `panel', clear
di _newline as txt "==== did_multiplegt_dyn Y g t D, effects(2) ===="
did_multiplegt_dyn Y g t D, effects(2)

* reveal stored results so you can read Effect_1 / Effect_2 programmatically
di _newline as txt "---- stored results (for the side-by-side) ----"
ereturn list
cap matrix list e(b)

di _newline "Expected: Effect_1 = 3.9750, Effect_2 = 4.6000 (matching our delta_l)."
di "Done."
