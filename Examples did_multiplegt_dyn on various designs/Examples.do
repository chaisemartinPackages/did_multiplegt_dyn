*===================================================
* PROGRAM:	DiD Examples for Various Designs       *
* PURPOSE:  did_multiplegt github examples         *
* DATE:	    May 19, 2025                           *
*===================================================

cd "/Users/renee/Desktop/Examples"

*=====================================================================
* Examples for One-shot Treatment: Deryugina (2017), Galagher (2014) *
*=====================================================================
******************************
****** Deryugina (2017) ******
******************************
use "data/Deryugina_2017.dta", clear
* Replication of Equation 1 & Figure 2 Panel A						
* Outcome of interest: log Per capita transfer from government (log_curr_trans_ind_gov_pc)
* Treatment: Hurricane 
* (g,t) panel: county, year

/* Original Paper: 
* Event Study Specification 
* ten years indicator of treatment, county and year fixed effects, hurricane occurrence outside the time interval of interest, year indicators interacted 1969 characteristics: land area; whether the county is coastal; population (in logs); the shares of the population that are black, under 20, or 65 and over; the employment rate; and per capita wages (in logs).
* Note the outcomes presented in the original paper apply spatially clustered standard errors, but the author also provides the following: if you do not have ols_spatial_HAC, you can obtain the paper's point estimates with standard errors clustered by county by running the following code (currently commented out)
*/ 

* Replication:
reghdfe log_curr_trans_ind_gov_pc hurricane hurr_*, absorb(i.county_fips i.year#i.coastal i.year#c.land_area1970 i.year#c.log_pop1969 i.year#c.frac_young1969 i.year#c.frac_old1969 i.year#c.frac_black1969 i.year#c.log_wage_pc1969 i.year#c.emp_rate1969) vce(cluster county_fips)

* Check one-shot treatment design
egen total_hurricane = total(hurricane), by(county_fips)
tab total_hurricane
/*total_hurri |
       cane |      Freq.     Percent        Cum.
------------+-----------------------------------
          0 |     33,975       68.36       68.36
          1 |     15,723       31.64      100.00
------------+-----------------------------------
      Total |     49,698      100.00
*/

* Time invariant controls interacted with time
local vars coastal land_area1970 log_pop1969 frac_young1969 frac_old1969 frac_black1969 log_wage_pc1969 emp_rate1969
foreach v of local vars {
    gen `v'_year = `v' * year
}

* did_multiplegt_dyn
* without controls
did_multiplegt_dyn log_curr_trans_ind_gov_pc county_fips year hurricane, effects(11) placebo(11) cluster(county_fips)
graph save "output/Deryugina_2017_g1.gph", replace
* with controls
did_multiplegt_dyn log_curr_trans_ind_gov_pc county_fips year hurricane, effects(11) placebo(11) controls(coastal_year land_area1970_year log_pop1969_year frac_young1969_year frac_old1969_year frac_black1969_year log_wage_pc1969_year emp_rate1969_year) cluster(county_fips)		   
graph save "output/Deryugina_2017_g2.gph", replace

graph combine "output/Deryugina_2017_g1.gph" "output/Deryugina_2017_g2.gph", ///
    col(2) xsize(10) ysize(5) ///
    saving("output/Deryugina_combined.gph", replace)
	
*****************************
****** Gallagher (2014) ******
*****************************
* Replication of Equation 1 & Figure 2				
* Outcome of interest: Food insurance take-up, ln_takeup_proxy1_dmn
* Treatment: Flood hit
* (g,t) panel: community (id_num), year (year2)
use "data/Gallagher_2014.dta", clear
*panel_1990_2007.dta
/* Original Paper: 
* Event Study Specification 
* community fixed effects, and state-by-year fixed effects
*/ 
* Replication: 
areg ln_takeup_proxy1_dmn hityear_m11_m17_dmn hityear_m10_dmn-hityear_m2_dmn hityear_dmn hityear_p1_dmn-hityear_p10_dmn hityear_p11_p17_dmn if panel9007 == 1, absorb(state_year) cluster(state)

* Check one-shot treatment design
order id_num year2 hityear hit90 hit91 hit92 hit93 hit94 hit95 hit96 hit97 hit98 hit99 hit00 hit01 hit02 hit03 hit04 hit05 hit06 hit07
egen total_flood = total(hityear), by(id_num)
tab total_flood

/*total_flood |      Freq.     Percent        Cum.
------------+-----------------------------------
          0 |     70,686       36.22       36.22
          1 |     55,656       28.52       64.74
          2 |     41,004       21.01       85.76
          3 |     17,856        9.15       94.91
          4 |      6,930        3.55       98.46
          5 |      2,106        1.08       99.54
          6 |        558        0.29       99.82
          7 |        126        0.06       99.89
          8 |         54        0.03       99.92
          9 |         90        0.05       99.96
         10 |         54        0.03       99.99
         11 |         18        0.01      100.00
------------+-----------------------------------
      Total |    195,138      100.00
*/

* did_multiplegt_dyn
encode st_fips, gen(state_id)
did_multiplegt_dyn ln_takeup_proxy1_dmn id_num year2 hityear if panel9007 == 1, effects(11) placebo(11) cluster(state_id) trends_nonparam(state_id) 
* trends_nonparam()
////////////////////////NOT SURE ABOUT THE STATE-YEAR-FIXED EFFECTS/////////////////////////////

*========================================================================================
* Examples for Several Treatments: Bradford et al. (2022), Rico-Straffon et al.  (2023) *
*========================================================================================
**********************************
****** Bradford et al. (2022) ****
**********************************
use "data/Bradford_et_al_2022.dta", clear

* rml_panel_data.dta
* Outcome: marijuana use in the past year
* Original Paper: 

* Check if two treatment in strict order, medical laws changes prior to recreational laws change
// keep only the age group 5 in the data. There are other age groups
keep if age_category == 5
gen Y = real(substr(year, 1, 4))
tabulate rm, generate(RM)
replace rm = RM2
drop RM1
drop RM2
tabulate mm, generate(MM)
replace mm = MM2
drop MM1
drop MM2
gen mm_rm_timing = mm-rm
sum mm_rm_timing // Medical Marijuana Legalization always before Recreational Marijuana Legalization

encode state, gen(state_encode)

* Rerun with did_multiplegt_dyn
// Estimate the Medical Marijuana treatment effect
preserve
*keep if rm ==0
did_multiplegt_dyn ln_mj_use_365 state_encode Y mm if rm==0, placebo(3) effects(3) 
restore
graph save "output/Bradford_et_al_2022_mm.gph", replace
// Estimate the Recreational Marijuana Law treatment effect
preserve
* Restrict sample to (g,t) s.t. D_g,t^1==1. Track the timing of F_g^1.
egen fg1 = min(cond(mm == 1, Y, .)), by(state_encode)
did_multiplegt_dyn ln_mj_use_365 state_encode Y rm if mm ==1, placebo(3) effects(3) trends_nonparam(fg1)
restore
graph save "output/Bradford_et_al_2022_rm.gph", replace

graph combine "output/Bradford_et_al_2022_mm.gph" "output/Bradford_et_al_2022_rm.gph", ///
    col(2) xsize(10) ysize(5) ///
    saving("output/Bradford_combined.gph", replace)

******************************************
****** Rico-Straffon et al.  (2023) ****** 
******************************************
* peru_concession_panel_for_fsc.dta
use "data/Rico-Straffon_et al_2023.dta", clear

* Replication of Equation 2 & Table 1 Column 3 						
* Outcome of interest:  Mapbiomas forest loss rate (%)
* Treatment: Logging concessions + eco-certifications of logging concessions (FSC)
* (g,t) panel: conecssion, year 

/* Original Paper: 
* TWFE Specification + DID_L (de Chaisemartin and D'Haultfœuille (2021; 2022b))
* Column 3 presents the FSC coefficient of running a TWFE regression in equation (2). In this case, FSC's TWFE estimator is a
weighted sum of two terms. The first term is a weighted sum of FSC's ATEs in each concession and year, while the second term is a weighted sum of the
concession effects (de Chaisemartin and D'Haultfœuille, 2022a). Column 3 shows that 4202/5937 ATEs in the second term receive negative weights
summing − 1.1. Thus, the TWFE estimator of the FSC effect is heavily contaminated by the effect of the concession treatment. Column 4 presents the
average DIDL estimator of the FSC effect, which is robust to heterogeneous effects and contamination biases. The joint significance test of the placebo
estimators rejects the null hypothesis (p-value =0.22). We clustered standard errors at concession level for FSC's impacts (Columns 3 and 4). 
*/ 

* did_multiplegt_dyn
* Check if two treatment in strict order, concession prior to FSC.
gen conc_fsc_timing = concession-FSC
tab conc_fsc_timing 
* Estimate Concession treatment effect
preserve 
did_multiplegt_dyn mapbiomasloss_p conc_id year concession if FSC==0, placebo(3) effects(3)
graph save "output/Rico-Straffon_et_al_2023_concession.gph", replace
restore 
* Estimate FSC treatment effect
preserve 
* Restrict sample to (g,t) s.t. D_g,t^=1. Track the timing of F_g^1.
egen fg1 = min(cond(concession == 1, year, .)), by(conc_id)
did_multiplegt_dyn mapbiomasloss_p conc_id year FSC if concession==1, placebo(3) effects(3) trends_nonparam(fg1)
restore
graph save "output/Rico-Straffon_et_al_2023_FSC.gph", replace

graph combine "output/Rico-Straffon_et_al_2023_concession.gph" "output/Rico-Straffon_et_al_2023_FSC.gph", ///
    col(2) xsize(10) ysize(5) ///
    saving("output/Rico-Straffon_combined.gph", replace)
	
*=========================================================================
* Examples for Treatment continuous dist in period 1: East et al. (2023) *
*=========================================================================
********************************
****** East et al. (2023) ****** 
********************************
use "data/East_et_al.dta", clear

* Replication of Equation 2 & Figure 7 (First generation outcome) & Section V. C. 
* Outcome of interest: First generation's birth weight 
* Treatment:  The percent of women of ­ child-bearing age eligible for Medicaid in the event of a pregnancy. 
* (g,t) panel: state, year

/* Original Paper: 
* Event study Specification
* Coefficient estimates are reported in percentage points. Estimated for infants born in ­ 1975–1988. ­ Preperiod trend is estimated 
and removed from all observations for each state prior to the event study estimation. For treated states, this is estimated using 
all ­ preperiod years for each state. For control states, we use the period ­ 1975–1981 to
estimate this trend. Regressions are weighted by birth cohort size and include state of birth and year of birth fixed
effects and controls for ­ state-year variables (unemployment rate, personal income per capita, maximum welfare
benefit for a family of four, indicators for state parental consent and notification laws and state Medicaid restric-
tions for abortion, and demographic controls for each state and year). Standard errors are clustered by infant's state
of birth.
*/ 

* Replication: 
global stcontrols stmarried stblack stother sthsdrop ///
sthsgrad stsomecoll pop0_4 pop5_17 pop18_24 pop25_44 ///
pop45_64 urate incpc maxafdc abortconsent abortmedr
global indvars evn6 evn5 evn4 evn3 evn2 omitted ev0 ev1 ev2 ev3 ev4
global graphvars evn5 evn4 evn3 evn2 omitted ev0 ev1 ev2 ev3
foreach x of var lbw {
reghdfe `x'_detrend81 $indvars $stcontrols [aweight=births], absorb(plborn dob_y_p) cluster(plborn)
estimates store baseline
coefplot (baseline,   color(black) ciopts(lcolor(black))), keep($graphvars) vertical omitted xline(5.5)  yline(0) xtitle("Year Relative to Expansion") ytitle("Coefficient Estimate")
*graph export `x'_detrend81_diffdiffcg.pdf, as(pdf) replace
}
*graph save "output/East_et_al_2022_replication.gph", replace

tab dob_y_p
* balanced panel, year 1975 to year 1988 
xtset plborn dob_y_p
gen d_simeli=D.simeligprenat
tab d_simeli
gen d_prenat=D.prenatelig 
tab d_prenat

*newsimeli is generated with the following rules with simeligprenat: average of 75 to 79 if never treated, average of 75 to 79 for year of birth before switch-in, actual simulated eligibility rate for switch-in year and afterwards. 

*did_multiplegt_dyn: 
*did_multiplegt_dyn lbw_detrend81 plborn dob_y_p newsimeli, effects(5) placebo(5) continuous(1) controls($stcontrols) weight(births)
*graph save "output/East_et_al_2022_g2.gph", replace
did_multiplegt_dyn lbw_detrend81 plborn dob_y_p newsimeli, effects(5) placebo(5) continuous(1) bootstrap(50,1) controls($stcontrols) weight(births) effects_equal("all")

did_multiplegt_dyn lbw_detrend81 plborn dob_y_p newsimeli, effects(5) placebo(5) continuous(1) bootstrap(50,1) controls($stcontrols) weight(births) normalized normalized_weights effects_equal("all")

graph save "output/East_et_al_2022_g3.gph", replace

graph combine "output/East_et_al_2022_g2.gph" "output/East_et_al_2022_g3.gph", ///
    col(2) xsize(10) ysize(5) ///
    saving("output/East_et_al_combined.gph", replace)

*=====================================================================
* Examples for On-and-off Discrete Treatment: Gentzkow et al.（2011） *
*=====================================================================
********************************
**** Gentzkow et al.（2011） ****
********************************
ssc desc cc_xd_didtextbook
net get cc_xd_didtextbook
dir *.zip
unzipfile cc_xd_didtextbook.zip
cd "cc_xd_didtextbook_2025_6_11/Data sets/Gentzkow et al 2011/"
use gentzkowetal_didtextbook.dta, clear

* Outcome of interest: voting turnouts
* Treatment:  Entries and exits of US daily newspapers.
* (g,t) panel: county, year
sort cnty90 year
xtset cnty90 year
bysort cnty90 (year): gen first_numdailies = numdailies[1]
tab first_numdailies

gen d_numdailies = .
bysort cnty90 (year): replace d_numdailies = numdailies - numdailies[_n-1]
tab d_numdailies

* Baseline regression .
did_multiplegt_dyn prestout cnty90 year numdailies, effects(4) placebo(4) 
* Understand the composition of the effect estimate with respect to treatment paths. 
did_multiplegt_dyn prestout cnty90 year numdailies, effects(2) design(0.8,console)
* Compute normalized actual - versus - status - quo event - study estimates. 
did_multiplegt_dyn prestout cnty90 year numdailies, effects(4) placebo(4) normalized effects_equal("all")
*Test that the first lagged treatment has no effect , and that treatment effects are constant over time .
did_multiplegt_dyn prestout cnty90 year numdailies if year <= first_change | same_treat_after_first_change==1, effects(2) effects_equal ("all") same_switchers graph_off


*================================================================
* Subsetting all of the datasets, Keeping variables needed only *
*================================================================
// PART 1
use "data/Deryugina_2017.dta", clear
keep log_curr_trans_ind_gov_pc hurricane hurr_* ///
     county_fips year ///
     coastal land_area1970 log_pop1969 frac_young1969 ///
     frac_old1969 frac_black1969 log_wage_pc1969 emp_rate1969
save "output/deryugina_2017.dta", replace
* source: Replication package on icspr, Final_dataset.dta

// PART 2
use "data/Bradford_et_al_2022.dta", clear
keep age_category year rm mm state ln_mj_use_365
keep if age_category == 5
gen new_year = real(substr(year, 1, 4))
drop year
rename new_year year
encode state, gen(state_encode)
drop state
rename state_encode state
save "output/hollingsworth_et_al_2022.dta", replace
* source: Replication package on icspr, rml_panel_data.dta

// PART 3
use "data/East_et_al_2022.dta", replace
keep plborn dob_y_p  newsimeli simeligprenat prenatelig lbw_detrend81 births stmarried stblack stother sthsdrop sthsgrad stsomecoll pop0_4 pop5_17 pop18_24 pop25_44 pop45_64 urate incpc maxafdc abortconsent abortmedr  
save "output/east_et_al_2023.dta", replace
* source: Replication package on icspr, data from NBER, self-defined newsimeli to ensure continous distribution yet with stayers  



