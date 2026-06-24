
********************************************************************************
*                                 PROGRAM 2                                    *
********************************************************************************
*Main goal: compute the variables U_Gg_plus_XX, U_Gg_minus_XX, U_Gg_var_plus_XX, and U_Gg_var_minus_XX,
*which are necessary to compute the DID_\ell estimators and their variance.

capture program drop did_multiplegt_dyn_core_new

program did_multiplegt_dyn_core_new, eclass
	version 12.0
	syntax varlist(min=4 max=4 numeric) [, effects(integer 1) placebo(integer 0) switchers(string) only_never_switchers controls(varlist numeric) trends_nonparam(varlist numeric) weight(varlist numeric) dont_drop_larger_lower NORMALIZED cluster(varlist numeric) graphoptions(string) SAVe_results(string) graph_off same_switchers same_switchers_pl drop_if_d_miss_before_first_switch trends_lin ci_level(integer 95) by(varlist numeric max=1) predict_het(string) design(string) date_first_switch(string) NORMALIZED_weights CONTinuous(string) switchers_core(string) less_conservative_se]
	qui{
		
////////// 1. Scalars initialization

///// Initializing the number of effects and placebos to estimate, depending on whether we consider switchers in or switchers out. 
///// Initializing also a scalar for whether the estimation is for switchers in or out.
	
if "`switchers_core'"=="in"{
scalar l_u_a_XX=min(`effects', L_u_XX)

if `placebo'!=0{
	scalar l_placebo_u_a_XX=min(`placebo', L_placebo_u_XX)
	}

	scalar increase_XX=1
}

if "`switchers_core'"=="out"{
scalar l_u_a_XX=min(`effects', L_a_XX)

if `placebo'!=0{
	scalar l_placebo_u_a_XX=min(`placebo', L_placebo_a_XX)
	}

scalar increase_XX=0
}

///// Initializing values of baseline treatment
levelsof d_sq_int_XX, local(levels_d_sq_XX)

////////// 2. Data preparation steps to generate variables necessary for computation of event-study effects	

///// Drop some variables we only create once
capture drop num_g_paths_0_XX
capture drop cohort_fullpath_0_XX
capture drop path_0_XX
capture drop d_fg0_XX 

///// Loop over the number of dynamic effects_equal

forvalue i=1/`=l_u_a_XX'{

///// Capture drop of variables created below

capture drop distance_to_switch_`i'_XX
capture drop never_change_d_`i'_XX
capture drop N`=increase_XX'_t_`i'_XX
capture drop N`=increase_XX'_t_`i'_XX_w
capture drop N`=increase_XX'_t_`i'_g_XX
capture drop N_gt_control_`i'_XX
capture drop diff_y_`i'_XX
capture drop dummy_U_Gg`i'_XX
capture drop U_Gg`i'_temp_XX
capture drop U_Gg`i'_XX
capture drop count`i'_core_XX
capture drop U_Gg`i'_temp_var_XX
capture drop U_Gg`i'_var_XX

capture drop never_change_d_`i'_wXX
capture drop distance_to_switch_`i'_wXX

capture drop d_fg`i'_XX
capture drop path_`i'_XX
capture drop num_g_paths_`i'_XX
capture drop cohort_fullpath_`i'_XX

capture drop dof_cohort_`i'_s_t_XX
capture drop dof_cohort_`i'_ns_t_XX
capture drop dof_cohort_`i'_s0_t_XX
capture drop dof_cohort_`i'_s1_t_XX
capture drop dof_cohort_`i'_s2_t_XX
capture drop count_cohort_`i'_s_t_XX
capture drop count_cohort_`i'_ns_t_XX
capture drop count_cohort_`i'_s0_t_XX
capture drop count_cohort_`i'_s1_t_XX
capture drop count_cohort_`i'_s2_t_XX
capture drop total_cohort_`i'_s_t_XX
capture drop total_cohort_`i'_ns_t_XX
capture drop total_cohort_`i'_s0_t_XX
capture drop total_cohort_`i'_s1_t_XX
capture drop total_cohort_`i'_s2_t_XX
capture drop mean_cohort_`i'_s_t_XX
capture drop mean_cohort_`i'_ns_t_XX
capture drop mean_cohort_`i'_s0_t_XX
capture drop mean_cohort_`i'_s1_t_XX
capture drop mean_cohort_`i'_s2_t_XX

capture drop E_hat_gt_`i'_XX
capture drop DOF_gt_`i'_XX

///// Creating long difference of outcome

xtset group_XX time_XX
bys group_XX : gen diff_y_`i'_XX = outcome_XX - L`i'.outcome_XX

///// Creating treatment paths if less_conservative_se option specified

if "`less_conservative_se'" != ""{

// Creating a time-invariant, group-level variable, containing g's treatment at F_g-1+\ell
gen d_fg_XX_temp=treatment_XX if time_XX==F_g_XX + `=-1+`i''
bys group_XX: gegen d_fg`i'_XX=mean(d_fg_XX_temp)

// This variable might be missing, for groups whose treatment never changes, and
// for groups not observed \ell periods after treatment change. Then, we impute
// their treatment at F_g-1+\ell-1. Inconsequential, just to avoid missing values.
// We also need to initialize a variable d_fg0_XX, when \ell=1.
if `i'==1{
gen d_fg0_XX=d_sq_XX
gegen path_0_XX=group(d_fg0_XX F_g_XX)
}
replace d_fg`i'_XX=d_fg`=`i'-1'_XX if d_fg`i'_XX==.
gegen path_`i'_XX = group(path_`=`i'-1'_XX d_fg`i'_XX)
drop d_fg_XX_temp 

// For each group, define a variable counting how many groups belong to the same cohort,
// with cohorts defined as d_fg0_XX F_g_XX, as well as the full path.
if `i'==1{
	bysort path_0_XX: gegen num_g_paths_0_XX=nunique(group_XX)
}
bysort path_`i'_XX: gegen num_g_paths_`i'_XX=nunique(group_XX)

// For each group, generate a dummy for whether at least two groups in their cohort.
gen cohort_fullpath_`i'_XX=(num_g_paths_`i'_XX>1)
if `i'==1{
gen cohort_fullpath_0_XX=(num_g_paths_0_XX>1) 
}
}

///// Identifying the control (g,t)s in the estimation of dynamic effect i 
bys group_XX: gen never_change_d_`i'_XX=(F_g_XX>time_XX) if diff_y_`i'_XX!=.

if "`only_never_switchers'" != "" {
	replace never_change_d_`i'_XX = 0 if F_g_XX > time_XX & F_g_XX < T_max_XX + 1 & diff_y_`i'_XX != .
}

///// Creating N^g_t:
///// number of control groups for g at t
gen never_change_d_`i'_wXX = never_change_d_`i'_XX*N_gt_XX
bys time_XX d_sq_XX `trends_nonparam': gegen N_gt_control_`i'_XX=total(never_change_d_`i'_wXX)

///// Creating binary variable indicating whether g is \ell periods away from switch

// If the same_switchers option is specified:

if ("`same_switchers'"!=""){

capture drop relevant_y_missing_XX
capture drop cum_fillin_XX
capture drop fillin_g_XX 
capture drop fillin_g_pl_XX
capture drop dum_fillin_temp_XX 
capture drop still_switcher_`i'_XX  

sort group_XX time_XX

// Modif. Diego 31-03-2025: long-difference and other variables to tailor distance_to_switch to the groups for which the all the effects can be computed
cap drop N_g_control_check_XX
gen N_g_control_check_XX = 0

forv j = 1/`=`effects'' {
	cap drop diff_y_last_XX
	cap drop never_change_d_last_XX
	cap drop never_change_d_last_wXX
	cap drop N_gt_control_last_XX
	cap drop N_g_control_last_temp_XX
	cap drop N_g_control_last_m_XX
	cap drop diff_y_relev_temp_XX
	cap drop diff_y_relev_XX

	xtset group_XX time_XX
	bys group_XX : gen diff_y_last_XX = outcome_XX - L`j'.outcome_XX
	bys group_XX: gen never_change_d_last_XX=(F_g_XX>time_XX) if diff_y_last_XX!=.
	if "`only_never_switchers'" != "" {
		replace never_change_d_last_XX = 0 if F_g_XX > time_XX & F_g_XX < T_max_XX + 1 & diff_y_last_XX != .
	}
	gen never_change_d_last_wXX = never_change_d_last_XX*N_gt_XX
	bys time_XX d_sq_XX `trends_nonparam': gegen N_gt_control_last_XX=total(never_change_d_last_wXX)

	gen N_g_control_last_temp_XX = N_gt_control_last_XX if time_XX == F_g_XX - 1 + `j'
	bys group_XX: egen N_g_control_last_m_XX = mean(N_g_control_last_temp_XX)

	gen diff_y_relev_temp_XX = diff_y_last_XX if time_XX == F_g_XX - 1 + `j'
	bys group_XX: egen diff_y_relev_XX = mean(diff_y_relev_temp_XX)

	replace N_g_control_check_XX = N_g_control_check_XX + (N_g_control_last_m_XX > 0 & diff_y_relev_XX != .)
}

* If the same_switchers_pl option is specified:

if ("`same_switchers_pl'"!=""){
cap drop N_g_control_check_pl_XX
gen N_g_control_check_pl_XX = 0
forv j = 1/`=`placebo'' {
	cap drop diff_y_last_XX
	cap drop never_change_d_last_XX
	cap drop never_change_d_last_wXX
	cap drop N_gt_control_last_XX
	cap drop N_g_control_last_temp_XX
	cap drop N_g_control_last_m_XX
	cap drop diff_y_relev_temp_XX
	cap drop diff_y_relev_XX

	xtset group_XX time_XX
	bys group_XX : gen diff_y_last_XX = outcome_XX - F`j'.outcome_XX
	bys group_XX: gen never_change_d_last_XX=(F_g_XX>time_XX) if diff_y_last_XX!=.
	if "`only_never_switchers'" != "" {
		replace never_change_d_last_XX = 0 if F_g_XX > time_XX & F_g_XX < T_max_XX + 1 & diff_y_last_XX != .
	}
	gen never_change_d_last_wXX = never_change_d_last_XX*N_gt_XX
	bys time_XX d_sq_XX `trends_nonparam': gegen N_gt_control_last_XX=total(never_change_d_last_wXX)

	gen N_g_control_last_temp_XX = N_gt_control_last_XX if time_XX == F_g_XX - 1 - `j'
	bys group_XX: egen N_g_control_last_m_XX = mean(N_g_control_last_temp_XX)

	gen diff_y_relev_temp_XX = diff_y_last_XX if time_XX == F_g_XX - 1 - `j'
	bys group_XX: egen diff_y_relev_XX = mean(diff_y_relev_temp_XX)

	replace N_g_control_check_pl_XX = N_g_control_check_pl_XX + (N_g_control_last_m_XX > 0 & diff_y_relev_XX != .)
}
* Generate a variable tagging the switchers that should be dropped
* Is the case if at least one of the placebos or effects we try to estimate is missing:
gen relevant_y_missing_XX=(outcome_XX==.&time_XX>=F_g_XX-1-`=`placebo''&time_XX<=F_g_XX-1+`=`effects'') 
* Or if some of the controls are missing:
if "`controls'" != ""{
replace relevant_y_missing_XX=1 if fd_X_all_non_missing_XX==0&time_XX>=F_g_XX-1-`=`placebo''&time_XX<=F_g_XX-1+`=`effects''
}

// Modif. Diego 19-04-25: make the same adjustments as below in case same_switcher_pl is specified
//bys group_XX: gen cum_fillin_XX = sum(relevant_y_missing_XX)
//gen dum_fillin_temp_XX = (cum_fillin_XX==0&time_XX==F_g_XX-1+`=`effects'')
//bys group_XX: gegen fillin_g_XX = total(dum_fillin_temp_XX)

//gen dum_fillin_temp_pl_XX = (cum_fillin_XX==0&time_XX==F_g_XX-1-`=`placebo'')
//bys group_XX: gegen fillin_g_pl_XX = total(dum_fillin_temp_pl_XX)
gen fillin_g_pl_XX = (N_g_control_check_pl_XX == `placebo')

capture drop dum_fillin_temp_XX
capture drop dum_fillin_temp_pl_XX

* tag switchers who have no missings from F_g_XX-1-`=`placebo'' to F_g_XX-1+`=`effects''
gen still_switcher_`i'_XX = (F_g_XX-1+`=`effects''<=T_g_XX & N_g_control_check_XX == `effects')  	

gen distance_to_switch_`i'_XX=(still_switcher_`i'_XX&time_XX==F_g_XX-1+`i'&`i'<=L_g_XX&S_g_XX==increase_XX&N_gt_control_`i'_XX>0&N_gt_control_`i'_XX!=.) if diff_y_`i'_XX!=. 
}

* If the same_switchers_pl option is not specified:
	
if ("`same_switchers_pl'"==""){
* Generate a variable tagging the switchers that should be dropped
* Is the case if at least one of the effects we try to estimate is missing:
gen relevant_y_missing_XX=(outcome_XX==.&time_XX>=F_g_XX-1&time_XX<=F_g_XX-1+`=`effects'') 
* Or if some of the controls are missing:
if "`controls'" != ""{
replace relevant_y_missing_XX=1 if fd_X_all_non_missing_XX==0&time_XX>=F_g_XX&time_XX<=F_g_XX-1+`=`effects''
}

// Modif. Diego 31-03-2025: we do not need this block if we use the changes above
//bys group_XX: gen cum_fillin_XX = sum(relevant_y_missing_XX)
//gen dum_fillin_temp_XX = (cum_fillin_XX==0&time_XX==F_g_XX-1+`=`effects'')
//bys group_XX: gegen fillin_g_XX = total(dum_fillin_temp_XX)

* tag switchers who have no missings from F_g_XX-1 to F_g_XX-1+`=`effects''
// Modif. Diego 31-03-2025: add check for control groups with the same baseline treat to still_switchers for all estimated periods

gen still_switcher_`i'_XX = (F_g_XX-1+`=`effects''<=T_g_XX & N_g_control_check_XX == `effects') 	
gen distance_to_switch_`i'_XX=(still_switcher_`i'_XX&time_XX==F_g_XX-1+`i'&`i'<=L_g_XX&S_g_XX==increase_XX&N_gt_control_`i'_XX>0&N_gt_control_`i'_XX!=.) if diff_y_`i'_XX!=.  

}
}

// If the same_switchers option is not specified

else{
gen distance_to_switch_`i'_XX=(time_XX==F_g_XX-1+`i'&`i'<=L_g_XX&S_g_XX==increase_XX&N_gt_control_`i'_XX>0&N_gt_control_`i'_XX!=.) if diff_y_`i'_XX!=. 
}

///// Creating a variable counting the number of groups \ell periods away from switch at t

gen distance_to_switch_`i'_wXX = distance_to_switch_`i'_XX*N_gt_XX
bys time_XX: gegen N`=increase_XX'_t_`i'_XX=total(distance_to_switch_`i'_wXX)

///// Computing N1_l,N0_l.

// Initializing the N1_`i'_XX/N0_`i'_XX scalar at 0. 
scalar N`=increase_XX'_`i'_XX =0

// Loop over t incrementing the scalar
forvalue t=`=t_min_XX'/`=T_max_XX'{
	sum N`=increase_XX'_t_`i'_XX if time_XX==`t'
	if !missing(r(mean)) scalar N`=increase_XX'_`i'_XX = N`=increase_XX'_`i'_XX + r(mean)
}

// Modif Felix weight
if "`weight'"!=""{

	bys time_XX: gegen N`=increase_XX'_t_`i'_XX_w=total(distance_to_switch_`i'_XX)
	
	scalar N`=increase_XX'_`i'_XX_weight =0

// Loop over t incrementing the scalar
forvalue t=`=t_min_XX'/`=T_max_XX'{
	sum N`=increase_XX'_t_`i'_XX_w if time_XX==`t'
	if !missing(r(mean)) scalar N`=increase_XX'_`i'_XX_weight = N`=increase_XX'_`i'_XX_weight + r(mean)
}
}

///// Creating N^1_{t,\ell,g}/N^0_{t,\ell,g}: 
///// Variable counting number of groups \ell periods away from switch at t, 
///// and with same D_{g,1} and trends_nonparam.

bys time_XX d_sq_XX `trends_nonparam': gegen N`=increase_XX'_t_`i'_g_XX=total(distance_to_switch_`i'_wXX)

///// Creating all the adjustment terms to compute estimators with controls, and their variances
 
if "`controls'" != ""{
	
// Initialize intermediate Variable needed later	
capture drop part2_switch`=increase_XX'_`i'_XX
gen part2_switch`=increase_XX'_`i'_XX=0

// generation of the T_d variable = max_{g:D_g,1 = d} F_g - 1: 
// last period when treatment effects can still be estimated for groups with baseline treatment equal to d

capture drop T_d_XX
gegen T_d_XX = max(F_g_XX), by(d_sq_int_XX)
replace T_d_XX = T_d_XX - 1

// Computing the long differences of the control variables (X_g_t - X_g_t-l)

local count_controls=0

foreach var of varlist `controls'{
		
local count_controls=`count_controls'+1

capture drop diff_X`count_controls'_`i'_XX

xtset group_XX time_XX
gen diff_X`count_controls'_`i'_XX=`var' - L`i'.`var'

// Computing N_g_t * (X_g_t - X_g_t-l)

capture drop diff_X`count_controls'_`i'_N_XX
gen diff_X`count_controls'_`i'_N_XX = N_gt_XX * diff_X`count_controls'_`i'_XX 

foreach l of local levels_d_sq_XX { // index l corresponds to d in the paper 

// intermediate variable to count the number of groups within each not yet switched cohort

// Changes Diego 27-03-25
//capture drop dummy_XX
//gen dummy_XX=(F_g_XX>time_XX & d_sq_int_XX == `l')

// Computing coordinates of vectors m_+_(g,d,l) and m_-_(g,d,l)

*Creating variable inside the summation across t in m^+_{g,d,\ell}/m^-_{g,d,\ell}
capture drop m`=increase_XX'_g_`count_controls'_`l'_`i'_XX 
gen m`=increase_XX'_g_`count_controls'_`l'_`i'_XX = (`i' <= T_g_XX-2 & d_sq_int_XX == `l')* (G_XX / N`=increase_XX'_`i'_XX) * ([distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * diff_X`count_controls'_`i'_N_XX)
*Summing that variable across t, and leaving one non missing observation per g	
capture drop m`=increase_XX'_`l'_`count_controls'_`i'_XX
bys group_XX: gegen m`=increase_XX'_`l'_`count_controls'_`i'_XX=total(m`=increase_XX'_g_`count_controls'_`l'_`i'_XX)
bys group_XX: replace m`=increase_XX'_`l'_`count_controls'_`i'_XX = . if _n != 1

// Computing coordinates of vectors M^+_{d,\ell} and M^-_{d,\ell}

capture drop M`=increase_XX'_`l'_`count_controls'_`i'_XX
egen M`=increase_XX'_`l'_`count_controls'_`i'_XX = total(m`=increase_XX'_`l'_`count_controls'_`i'_XX)
replace M`=increase_XX'_`l'_`count_controls'_`i'_XX = (1/G_XX)*M`=increase_XX'_`l'_`count_controls'_`i'_XX

// number of groups within each not yet switched cohort
capture drop E_hat_denom_`count_controls'_`l'_XX
//// CHANGE BELOW - tks + Changes Diego 27-03-25: replace dummy_XX_2 with dummy_XX

// Modif Clément 20/6/2025:
* Counting number of groups for DOF adjustment
{
cap drop dummy_XX
gen dummy_XX = 0
replace dummy_XX = (F_g_XX>time_XX & d_sq_int_XX == `l') if diff_y_XX < .
	if "`cluster'" == "" {
		bys time_XX : egen E_hat_denom_`count_controls'_`l'_XX = total(dummy_XX) if d_sq_int_XX == `l'
	}
	else {
		cap drop cluster_temp_XX
		gen cluster_temp_XX= `cluster' if dummy_XX == 1
		bys time_XX : gegen E_hat_denom_`count_controls'_`l'_XX = nunique(cluster_temp_XX) if !missing(cluster_temp_XX)
	}
}
// End Modif Clément 20/6/2025

// Add the indicator for at least two groups in the cohort to E_y_hat_gt_`l'_XX (demeaning is possible)
capture drop E_y_hat_gt_`l'_XX
gen E_y_hat_gt_`l'_XX=E_y_hat_gt_int_`l'_XX*(E_hat_denom_`count_controls'_`l'_XX>=2)


// Computing the summation from t=2 to F_g-1 that appears in the last term 
// of U^{+,var,X}_{g,\ell} and U^{-,var,X}_{g,\ell}, defined in the companion paper.
 
capture drop in_sum_temp_`count_controls'_`l'_XX
capture drop N_c_`l'_temp_XX
capture drop N_c_`l'_XX
//// CHANGE BELOW - tks, Changes Diego 27-03-25: Multiply by N_gt_XX
gen N_c_`l'_temp_XX = N_gt_XX * (d_sq_int_XX == `l' & time_XX >= 2 & time_XX <= T_d_XX & time_XX < F_g_XX & diff_y_XX < .)
////
egen N_c_`l'_XX = total(N_c_`l'_temp_XX)
/// Changes Diego 16-06-25: Adjust demeaning when E_hat_denom == 1
cap drop in_sum_temp_adj_`count_controls'_`l'_XX
gen in_sum_temp_adj_`count_controls'_`l'_XX = 0 if E_y_hat_gt_`l'_XX != .
replace in_sum_temp_adj_`count_controls'_`l'_XX = (sqrt((E_hat_denom_`count_controls'_`l'_XX)/(E_hat_denom_`count_controls'_`l'_XX - 1))-1) if E_y_hat_gt_`l'_XX != . & E_hat_denom_`count_controls'_`l'_XX > 1
gen in_sum_temp_`count_controls'_`l'_XX = (prod_X`count_controls'_Ngt_XX*(1+(E_hat_denom_`count_controls'_`l'_XX>=2)*in_sum_temp_adj_`count_controls'_`l'_XX)*(diff_y_XX-E_y_hat_gt_`l'_XX)*(time_XX>=2 & time_XX<=F_g_XX-1)) / N_c_`l'_XX


capture drop in_sum_`count_controls'_`l'_XX
capture drop in_sum_temp_adj_`count_controls'_`l'_XX
bys group_XX: gegen in_sum_`count_controls'_`l'_XX = total(in_sum_temp_`count_controls'_`l'_XX) 

}


// Residualize the outcome difference wrt control differences:
// Yg,t − Yg,t−ℓ − (Xg,t − Xg,t−ℓ)*θ_{Dg,1}

foreach l of local levels_d_sq_XX { 
	if (scalar(useful_res_`l'_XX)>1){
		replace diff_y_`i'_XX = diff_y_`i'_XX - coefs_sq_`l'_XX[`=`count_controls'',1]*diff_X`count_controls'_`i'_XX if d_sq_int_XX==`l' 
		* N.B. : in the above line, we do not add "&diff_X`count_controls'_`i'_XX!=." because we want to exclude from the estimation any first/long-difference for which the covariates are missing.

		// Initialize intermediate Variable needed later
		capture drop in_brackets_`l'_`count_controls'_XX	
		gen in_brackets_`l'_`count_controls'_XX=0

	}
}

}
}

///// Computing the variables used for the demeaning of outcome's long difference diff_y_`i',
///// which we will use to create the U_g^{var} variables, which are used to compute 
///// the estimators' variances. 

// Generate long-differences of outcomes time N_{g,t}, and dummy for (g,t) such that
// diff_y_`i' non missing and N_gt non missing.

capture drop diff_y_`i'_N_gt_XX
gen diff_y_`i'_N_gt_XX=N_gt_XX*diff_y_`i'_XX

capture drop dof_ns_`i'_XX
gen dof_ns_`i'_XX = (N_gt_XX != 0&diff_y_`i'_XX!=.&never_change_d_`i'_XX==1&N`=increase_XX'_t_`i'_XX>0&N`=increase_XX'_t_`i'_XX!=.)
capture drop dof_s_`i'_XX
gen dof_s_`i'_XX = (N_gt_XX != 0&distance_to_switch_`i'_XX==1)

// For never switchers, demeaning wrt to cohorts defined by D_{g,1}, `trends_nonparam' 
//(\mathcal{D}_k in companion paper)

// Modif Clément: we need to add by time
* Mean's denominator
bys d_sq_XX `trends_nonparam' time_XX : gegen count_cohort_`i'_ns_t_XX=total(N_gt_XX) if dof_ns_`i'_XX == 1

* Mean's numerator
bys d_sq_XX `trends_nonparam' time_XX: gegen total_cohort_`i'_ns_t_XX=total(diff_y_`i'_N_gt_XX) if dof_ns_`i'_XX == 1

* Mean 
gen mean_cohort_`i'_ns_t_XX=total_cohort_`i'_ns_t_XX/count_cohort_`i'_ns_t_XX

* Counting number of groups for DOF adjustment
{
	if "`cluster'" == "" {
		bys d_sq_XX `trends_nonparam' time_XX: gegen dof_cohort_`i'_ns_t_XX=total(dof_ns_`i'_XX) if dof_ns_`i'_XX == 1
	}
	else {
		cap drop cluster_dof_`i'_ns_XX
		gen cluster_dof_`i'_ns_XX = `cluster' if dof_ns_`i'_XX == 1
		bys d_sq_XX `trends_nonparam' time_XX: gegen dof_cohort_`i'_ns_t_XX = nunique(cluster_dof_`i'_ns_XX) if !missing(cluster_dof_`i'_ns_XX)
	}
}
// End modif Clément

// For switchers, if option less_conservative_se not specified, demeaning wrt to 
// cohorts defined by D_{g,1}, F_g, D_{g,F_g}, `trends_nonparam' (\mathcal{C}_k in companion paper).

if "`less_conservative_se'" == ""{
	
* Mean's denominator
bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam' : gegen count_cohort_`i'_s_t_XX=total(N_gt_XX) if dof_s_`i'_XX==1

* Mean's numerator
bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam' : gegen total_cohort_`i'_s_t_XX=total(diff_y_`i'_N_gt_XX) if dof_s_`i'_XX==1 

* Mean 
gen mean_cohort_`i'_s_t_XX=total_cohort_`i'_s_t_XX/count_cohort_`i'_s_t_XX	
		
* Counting number of groups for DOF adjustment
// Changes Diego 16-06-25: adjust DoF with clusters
{
	if "`cluster'" == "" {
		bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam' : gegen dof_cohort_`i'_s_t_XX=total(dof_s_`i'_XX) if dof_s_`i'_XX==1
	}
	else {
		cap drop cluster_dof_`i'_s_XX
		gen cluster_dof_`i'_s_XX = `cluster' if dof_s_`i'_XX == 1
		bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam': gegen dof_cohort_`i'_s_t_XX = nunique(cluster_dof_`i'_s_XX) if !missing(cluster_dof_`i'_s_XX)
	}
}
}
//save test_data, replace


// For switchers, if option less_conservative_se specified, demeaning wrt to cohorts 
// defined by D_{g,1} F_g, D_{g,F_g},..., D_{g,F_g+\ell}, if that cohort has at least two groups,
// if not: demeaning wrt to cohorts 
// defined by D_{g,1} F_g, D_{g,F_g}, if that cohort has at least two groups,
// if not: demeaning wrt D_{g,1} F_g.

if "`less_conservative_se'" != ""{
		
* Mean's denominator	
* by D_{g,1}, F_g, `trends_nonparam':
bys path_0_XX `trends_nonparam' : gegen count_cohort_`i'_s0_t_XX=total(N_gt_XX) if dof_s_`i'_XX==1
* by D_{g,1}, F_g, D_{g,F_g}, `trends_nonparam':
bys path_1_XX `trends_nonparam' : gegen count_cohort_`i'_s1_t_XX=total(N_gt_XX) if dof_s_`i'_XX==1
* by D_{g,1}, F_g, D_{g,F_g},..., D_{g,F_g+\ell}, `trends_nonparam':
bys path_`i'_XX `trends_nonparam' : gegen count_cohort_`i'_s2_t_XX=total(N_gt_XX) if dof_s_`i'_XX==1
	
* Mean's numerator
* by D_{g,1}, F_g, `trends_nonparam':
bys path_0_XX `trends_nonparam' : gegen total_cohort_`i'_s0_t_XX=total(diff_y_`i'_N_gt_XX) if dof_s_`i'_XX==1
* by D_{g,1}, F_g, D_{g,F_g}, `trends_nonparam':
bys path_1_XX `trends_nonparam' : gegen total_cohort_`i'_s1_t_XX=total(diff_y_`i'_N_gt_XX) if dof_s_`i'_XX==1	
* by D_{g,1}, F_g, D_{g,F_g},..., D_{g,F_g+\ell}, `trends_nonparam':
bys path_`i'_XX `trends_nonparam' : gegen total_cohort_`i'_s2_t_XX=total(diff_y_`i'_N_gt_XX) if dof_s_`i'_XX==1

* Counting number of groups for DOF adjustment
// Changes Diego 16-06-25: adjust DoF with clusters
{
	if "`cluster'" == "" {
		* by D_{g,1}, F_g, `trends_nonparam':
		bys path_0_XX `trends_nonparam' : gegen dof_cohort_`i'_s0_t_XX=total(dof_s_`i'_XX) if dof_s_`i'_XX==1
		* by D_{g,1}, F_g, D_{g,F_g}, `trends_nonparam':
		bys path_1_XX `trends_nonparam' : gegen dof_cohort_`i'_s1_t_XX=total(dof_s_`i'_XX) if dof_s_`i'_XX==1
		* by D_{g,1}, F_g, D_{g,F_g},..., D_{g,F_g+\ell}, `trends_nonparam':
		bys path_`i'_XX `trends_nonparam' : gegen dof_cohort_`i'_s2_t_XX=total(dof_s_`i'_XX) if dof_s_`i'_XX==1
	}
	else {
		cap drop cluster_dof_`i'_s_XX
		gen cluster_dof_`i'_s_XX = `cluster' if dof_s_`i'_XX == 1
		* by D_{g,1}, F_g, `trends_nonparam':
		bys path_0_XX `trends_nonparam' : gegen dof_cohort_`i'_s0_t_XX = nunique(cluster_dof_`i'_s_XX) if !missing(cluster_dof_`i'_s_XX)
		* by D_{g,1}, F_g, D_{g,F_g}, `trends_nonparam':
		bys path_1_XX `trends_nonparam' : gegen dof_cohort_`i'_s1_t_XX = nunique(cluster_dof_`i'_s_XX) if !missing(cluster_dof_`i'_s_XX)
		* by D_{g,1}, F_g, D_{g,F_g},..., D_{g,F_g+\ell}, `trends_nonparam':
		bys path_`i'_XX `trends_nonparam' : gegen dof_cohort_`i'_s2_t_XX = nunique(cluster_dof_`i'_s_XX) if !missing(cluster_dof_`i'_s_XX)

	}
}

*Attributing the right number of groups depending on which cohort will be used in demeaning:
gen dof_cohort_`i'_s_t_XX=dof_cohort_`i'_s2_t_XX if dof_cohort_`i'_s2_t_XX>=2
replace dof_cohort_`i'_s_t_XX=dof_cohort_`i'_s1_t_XX if dof_cohort_`i'_s2_t_XX<2&dof_cohort_`i'_s1_t_XX>=2 
replace dof_cohort_`i'_s_t_XX=dof_cohort_`i'_s0_t_XX if dof_cohort_`i'_s2_t_XX<2&dof_cohort_`i'_s1_t_XX<2 

* Mean
gen mean_cohort_`i'_s_t_XX=total_cohort_`i'_s2_t_XX/count_cohort_`i'_s2_t_XX if dof_cohort_`i'_s2_t_XX>=2
replace mean_cohort_`i'_s_t_XX=total_cohort_`i'_s1_t_XX/count_cohort_`i'_s1_t_XX if dof_cohort_`i'_s2_t_XX<2&dof_cohort_`i'_s1_t_XX>=2 
replace mean_cohort_`i'_s_t_XX=total_cohort_`i'_s0_t_XX/count_cohort_`i'_s0_t_XX if dof_cohort_`i'_s2_t_XX<2&dof_cohort_`i'_s1_t_XX<2

}

// Modif Clément 30/6/2025:
// If a switcher is the only one in their cohort or if a not-yet-switcher is the only one in their cohort, we demean wrt union of switchers and not-yet switchers, provided switchers and not-yet-switchers do not all come from the same cluster.

cap drop dof_ns_s_`i'_XX
cap drop count_cohort_`i'_ns_s_t_XX
cap drop total_cohort_`i'_ns_s_t_XX
cap drop mean_cohort_`i'_ns_s_t_XX
cap drop dof_cohort_`i'_ns_s_t_XX

gen dof_ns_s_`i'_XX=(dof_s_`i'_XX==1|dof_ns_`i'_XX ==1)

* Mean's denominator
bys d_sq_XX `trends_nonparam' time_XX : gegen count_cohort_`i'_ns_s_t_XX=total(N_gt_XX) if dof_ns_s_`i'_XX==1

* Mean's numerator
bys d_sq_XX `trends_nonparam' time_XX: gegen total_cohort_`i'_ns_s_t_XX=total(diff_y_`i'_N_gt_XX) if dof_ns_s_`i'_XX==1

* Mean 
gen mean_cohort_`i'_ns_s_t_XX=total_cohort_`i'_ns_s_t_XX/count_cohort_`i'_ns_s_t_XX

* Counting number of groups for DOF adjustment
{
	if "`cluster'" == "" {
		bys d_sq_XX `trends_nonparam' time_XX: gegen dof_cohort_`i'_ns_s_t_XX=total(dof_ns_s_`i'_XX) if dof_ns_s_`i'_XX==1
	}
	else {
		cap drop cluster_dof_`i'_ns_s_XX
		gen cluster_dof_`i'_ns_s_XX = `cluster' if dof_ns_s_`i'_XX==1
		bys d_sq_XX `trends_nonparam' time_XX: gegen dof_cohort_`i'_ns_s_t_XX = nunique(cluster_dof_`i'_ns_s_XX) if !missing(cluster_dof_`i'_ns_s_XX)
	}
}	

///// From those parts, generate variables for the demeaning and the DOF adjustment 
// E_hat_(g,t), defined from parts depending on the cohort definition
gen E_hat_gt_`i'_XX=0 if (time_XX<F_g_XX|F_g_XX-1+`i'==time_XX) 
replace E_hat_gt_`i'_XX=mean_cohort_`i'_ns_t_XX if (time_XX<F_g_XX&dof_cohort_`i'_ns_t_XX>=2)
replace E_hat_gt_`i'_XX=mean_cohort_`i'_s_t_XX if (F_g_XX-1+`i'==time_XX&dof_cohort_`i'_s_t_XX>=2)
replace E_hat_gt_`i'_XX=mean_cohort_`i'_ns_s_t_XX if dof_cohort_`i'_ns_s_t_XX>=2&((F_g_XX-1+`i'==time_XX&dof_cohort_`i'_s_t_XX==1)|(time_XX<F_g_XX&dof_cohort_`i'_ns_t_XX==1))

// DOF_(g,t) for DOF adjustement, defined from parts depending on the cohort definition 
// Diego - 02-03-25: when there is only 1 switcher, dof_cohort_`i'_s_t_XX = 1, hence the denominator in the expression below is 0
// The fraction is undefined and Stata puts it to missing
gen DOF_gt_`i'_XX=1 if (time_XX<F_g_XX|F_g_XX-1+`i'==time_XX)
replace DOF_gt_`i'_XX= sqrt(dof_cohort_`i'_s_t_XX/(dof_cohort_`i'_s_t_XX-1)) if (F_g_XX-1+`i'==time_XX & dof_cohort_`i'_s_t_XX > 1)
replace DOF_gt_`i'_XX=  sqrt(dof_cohort_`i'_ns_t_XX/(dof_cohort_`i'_ns_t_XX-1)) if (time_XX<F_g_XX & dof_cohort_`i'_ns_t_XX > 1)
replace DOF_gt_`i'_XX=  sqrt(dof_cohort_`i'_ns_s_t_XX/(dof_cohort_`i'_ns_s_t_XX-1)) if dof_cohort_`i'_ns_s_t_XX>=2&((F_g_XX-1+`i'==time_XX&dof_cohort_`i'_s_t_XX==1)|(time_XX<F_g_XX&dof_cohort_`i'_ns_t_XX==1))

// End modif Clément 30/6/2025.

////////// 3. Computing U_Gg_\ell variables

///// If the dynamic effect can be estimated (as there are switchers), we compute the U_Gg_\ell variables etc.

if (N`=increase_XX'_`i'_XX!=0){

// Creating a dummy variable indicating whether l<=T_g_XX-1

gen dummy_U_Gg`i'_XX = (`i'<=T_g_XX-1) 

// Computing U_+_(G,g,l)

gen U_Gg`i'_temp_XX = dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * N_gt_XX* [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] 

replace U_Gg`i'_temp_XX = U_Gg`i'_temp_XX* diff_y_`i'_XX 

bysort group_XX : gegen U_Gg`i'_XX=total(U_Gg`i'_temp_XX)

replace U_Gg`i'_XX = U_Gg`i'_XX*first_obs_by_gp_XX

// Counting the number of cells for which we can estimate U_Gg`i'_temp_XX, 
// to compute the "N" displayed by command

gen count`i'_core_XX=0
replace count`i'_core_XX=N_gt_XX if (U_Gg`i'_temp_XX!=.&U_Gg`i'_temp_XX!=0|(U_Gg`i'_temp_XX==0&diff_y_`i'_XX==0&(distance_to_switch_`i'_XX!=0|(N`=increase_XX'_t_`i'_g_XX!=0&never_change_d_`i'_XX!=0))))
 
// Computing U_(+,var)_(G,g,l)

gen U_Gg`i'_temp_var_XX = 0

// Final computation 

replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * N_gt_XX * DOF_gt_`i'_XX *(diff_y_`i'_XX-E_hat_gt_`i'_XX)

///// Adding the additional part of U^(+,var,X)_{G,g,l}/U^(-,var,X)_{G,g,l} when controls are included:
///// sum across values of baseline treatment d of M^+_(d,l)* a term in brackets in companion paper. 

if "`controls'"!=""{

// Loop over values of d_sq_int_XX:sum across values of baseline treatment  
levelsof d_sq_int_XX, local(levels_d_sq_XX)

foreach l of local levels_d_sq_XX {	
	// FELIX: Add this condition here, same for the placebo case!	
	if (scalar(useful_res_`l'_XX)>1){
	
capture drop combined`=increase_XX'_temp_`l'_`i'_XX	
gen combined`=increase_XX'_temp_`l'_`i'_XX=0
		
// Loop over controls: inner product of M^+_(d,l)* term in brackets in companion paper.
forvalues j=1/`count_controls'{

// Loop over k: computation of term in brackets in companion paper.
forvalues k=1/`count_controls'{		
	
// Computation of all cross products between elements of jth line of Den^{-1}_d and term multiplying it
capture drop in_brackets_`l'_`j'_`k'_temp_XX
//// CHANGE BELOW - tks + Changes Diego 27-03-25: remove N_c_`l'_XX from formula below, since we have adjusted inv_Denom above
gen in_brackets_`l'_`j'_`k'_temp_XX = inv_Denom_`l'_XX[`j',`k'] * in_sum_`k'_`l'_XX * (d_sq_int_XX == `l' & F_g_XX>=3) 
////

// Summing over k, to have jth coordinate of vector Den^{-1}_d*...

replace in_brackets_`l'_`j'_XX=in_brackets_`l'_`j'_XX+in_brackets_`l'_`j'_`k'_temp_XX
} // end loop over k

// Withdrawing theta_d
replace in_brackets_`l'_`j'_XX=in_brackets_`l'_`j'_XX - coefs_sq_`l'_XX[`j',1]

// Computation of all cross products between coordinates of M^+_(d,l) and those of term in brackets
capture drop combined`=increase_XX'_temp_`l'_`j'_`i'_XX
gen combined`=increase_XX'_temp_`l'_`j'_`i'_XX=M`=increase_XX'_`l'_`j'_`i'_XX*in_brackets_`l'_`j'_XX

// Summing over j
replace combined`=increase_XX'_temp_`l'_`i'_XX=combined`=increase_XX'_temp_`l'_`i'_XX+combined`=increase_XX'_temp_`l'_`j'_`i'_XX
} // end loop over j

// Final sum over the status quo treatment (outer sum over d in the formula)
// Modif. Diego: dropped  if d_sq_int_XX==`l'
replace part2_switch`=increase_XX'_`i'_XX=part2_switch`=increase_XX'_`i'_XX+combined`=increase_XX'_temp_`l'_`i'_XX
} // Modif FELIX: condition "useful residual" 
} // end loop over l

}	

// Summing the U^(var)_{G,g,l}s over time periods for each group
bys group_XX: gegen U_Gg`i'_var_XX=total(U_Gg`i'_temp_var_XX)

// Modeif Felix 21.03.2025 (adjust order of adding additional term and summing across t)
if "`controls'" != "" {
	// Making the adjustement to U^(+,var)_{G,g,l} when controls are included
	if `=increase_XX'==1{
	replace U_Gg`i'_var_XX=U_Gg`i'_var_XX - part2_switch1_`i'_XX 
	}

	if `=increase_XX'==0{
	//// CHANGE BELOW - tks (switching plus to minus in the sum since the final sign flip won't be until later)
	replace U_Gg`i'_var_XX=U_Gg`i'_var_XX - part2_switch0_`i'_XX 
	}
}

}

////////// 4. Computing adjustements for the normalized and trends_lin options 

///// Compute \delta^D for normalized option

if "`normalized'"!=""{
	
	capture drop sum_treat_until_`i'_XX
	capture drop delta_D_`i'_cum_temp_XX
	
	if `continuous'==0{
		bys group_XX : gegen sum_treat_until_`i'_XX = total(treatment_XX - d_sq_XX) if time_XX>=F_g_XX&time_XX<=F_g_XX-1+`i'&S_g_XX==increase_XX
	}
	// Redefine this with original treatment if continuous is defined (treatment was binarized and staggerized)
	else if `continuous'>0{
		bys group_XX : gegen sum_treat_until_`i'_XX = total(treatment_XX_orig - d_sq_XX_orig) if time_XX>=F_g_XX&time_XX<=F_g_XX-1+`i'&S_g_XX==increase_XX
	}
gen delta_D_`i'_cum_temp_XX = N_gt_XX/N`=increase_XX'_`i'_XX*[sum_treat_until_`i'_XX* S_g_XX + (1-S_g_XX)*(-sum_treat_until_`i'_XX)] if distance_to_switch_`i'_XX==1 

	sum  delta_D_`i'_cum_temp_XX
	scalar delta_norm_`i'_XX = r(sum) 

}

} // End of the loop over l_u_a_XX.
// At this point we have all the U_Gg_`i'_XX we need. Thus, we can sum them now for the trends_lin option.

///// trends_lin option

// As trends_lin relies on summing up the l effects we can only estimate effect l when we can estimate all prior l-1 effects, we verify if this condition is met.

scalar Ntrendslin=1
forvalue i=1/`=l_u_a_XX' {
scalar Ntrendslin=min(Ntrendslin , N`=increase_XX'_`i'_XX )
}

// Compute the U_Gg_\ell for trends_lin

if "`trends_lin'"!="" & Ntrendslin != 0  {
	
	capture drop U_Gg`=l_u_a_XX'_TL
	capture drop U_Gg`=l_u_a_XX'_var_TL
	
	// Initializing at 0
	gen U_Gg`=l_u_a_XX'_TL = 0
	gen U_Gg`=l_u_a_XX'_var_TL = 0

	// summing up the U_Gg's up to the l-th (current) effect
	forvalue i=1/`=l_u_a_XX'{
		replace U_Gg`=l_u_a_XX'_TL = U_Gg`=l_u_a_XX'_TL + U_Gg`i'_XX 		
		replace U_Gg`=l_u_a_XX'_var_TL =  U_Gg`=l_u_a_XX'_var_TL + U_Gg`i'_var_XX 
	}
	
	// replacing the U_Gg's with the one adjusted for group specific linear trends
	replace U_Gg`=l_u_a_XX'_XX = U_Gg`=l_u_a_XX'_TL
	replace U_Gg`=l_u_a_XX'_var_XX = U_Gg`=l_u_a_XX'_var_TL
		
}

////////// 5. Data preparation steps to generate variables necessary for computation of placebo effects	

if `placebo'!=0{
	if `=l_placebo_u_a_XX'>=1{
	
forvalue i=1/`=l_placebo_u_a_XX'{

///// Capture drop of variables created below
capture drop diff_y_pl_`i'_XX
capture drop U_Gg_pl_`i'_temp_XX
capture drop U_Gg_placebo_`i'_XX
capture drop U_Gg_pl_`i'_temp_var_XX
capture drop U_Gg_pl_`i'_var_XX
capture drop dist_to_switch_pl_`i'_XX
capture drop never_change_d_pl_`i'_XX
capture drop N`=increase_XX'_t_placebo_`i'_XX
capture drop N`=increase_XX'_t_pl_`i'_XX_w
capture drop N`=increase_XX'_t_placebo_`i'_g_XX
capture drop N_gt_control_placebo_`i'_XX
capture drop dummy_U_Gg_pl_`i'_XX 
capture drop never_change_d_pl_`i'_wXX
capture drop dist_to_switch_pl_`i'_wXX

capture drop dof_cohort_pl_`i'_s_t_XX
capture drop dof_cohort_pl_`i'_ns_t_XX
capture drop dof_cohort_pl_`i'_s0_t_XX
capture drop dof_cohort_pl_`i'_s1_t_XX
capture drop dof_cohort_pl_`i'_s2_t_XX
capture drop count_cohort_pl_`i'_s_t_XX
capture drop count_cohort_pl_`i'_ns_t_XX
capture drop count_cohort_pl_`i'_s0_t_XX
capture drop count_cohort_pl_`i'_s1_t_XX
capture drop count_cohort_pl_`i'_s2_t_XX
capture drop total_cohort_pl_`i'_s_t_XX
capture drop total_cohort_pl_`i'_ns_t_XX
capture drop total_cohort_pl_`i'_s0_t_XX
capture drop total_cohort_pl_`i'_s1_t_XX
capture drop total_cohort_pl_`i'_s2_t_XX
capture drop mean_cohort_pl_`i'_s_t_XX
capture drop mean_cohort_pl_`i'_ns_t_XX
capture drop mean_cohort_pl_`i'_s0_t_XX
capture drop mean_cohort_pl_`i'_s1_t_XX
capture drop mean_cohort_pl_`i'_s2_t_XX

capture drop E_hat_gt_pl_`i'_XX
capture drop DOF_gt_pl_`i'_XX

//The steps to compute the placebos are:
// 1. to place the corresponding outcome (y_{F_g-1} - y_{F_g - l - 1})) values in the same row of that (y_{F_g + l -1} - y_{F_g - 1}) of the symmetric DID_l. 
// 2. The other variables, such as N_gt, N0_l or N1_l, remain unchanged, except that we have to check if diff_y_placebo ( = y_{F_g - 2l -2}- y_{F_g - l -1}) exists. 

///// Computing the long differences for the placebos
xtset group_XX time_XX
bys group_XX : gen diff_y_pl_`i'_XX = L`=2*`i''.outcome_XX - L`i'.outcome_XX

///// Identifying the control (g,t)s in the estimation of placebo i 
bys group_XX: gen never_change_d_pl_`i'_XX=never_change_d_`i'_XX*(diff_y_pl_`i'_XX!=.) 
gen never_change_d_pl_`i'_wXX=never_change_d_pl_`i'_XX*N_gt_XX

///// number of control groups for g at t
bys time_XX d_sq_XX `trends_nonparam': gegen N_gt_control_placebo_`i'_XX=total(never_change_d_pl_`i'_wXX) 

///// Creating binary variable indicating whether g is \ell periods away from switch & (diff_y_pl_`i'_XX!=.) is well defined -> based on the already defined distance_to_switch_`i'_XX variable 
*Modif Clément 26/09/2025:
*gen dist_to_switch_pl_`i'_XX=distance_to_switch_`i'_XX*(diff_y_pl_`i'_XX!=.)
gen dist_to_switch_pl_`i'_XX=distance_to_switch_`i'_XX*(diff_y_pl_`i'_XX!=.)*(N_gt_control_placebo_`i'_XX>0&N_gt_control_placebo_`i'_XX!=.)
 
if "`same_switchers_pl'"!=""{
	replace dist_to_switch_pl_`i'_XX=dist_to_switch_pl_`i'_XX*fillin_g_pl_XX
}

///// Creating a variable counting the number of groups \ell periods away from switch at t
gen dist_to_switch_pl_`i'_wXX= dist_to_switch_pl_`i'_XX*N_gt_XX
bys time_XX: gegen N`=increase_XX'_t_placebo_`i'_XX=total(dist_to_switch_pl_`i'_wXX)

///// Computing N^1_\ell/N^0_\ell. for the placebos

// Initializing the N1_`i'_XX/N0_`i'_XX scalar at 0. 
scalar N`=increase_XX'_placebo_`i'_XX=0

// Loop over t incrementing the scalar
forvalue t=`=t_min_XX'/`=T_max_XX'{
	sum N`=increase_XX'_t_placebo_`i'_XX if time_XX==`t'
	if !missing(r(mean)) scalar N`=increase_XX'_placebo_`i'_XX = N`=increase_XX'_placebo_`i'_XX + r(mean) // Modif Diego
}

// Modif Felix weight
if "`weight'"!=""{

	bys time_XX: gegen N`=increase_XX'_t_pl_`i'_XX_w=total(dist_to_switch_pl_`i'_XX)
	
	scalar N`=increase_XX'_`i'_pl_XX_weight=0

// Loop over t incrementing the scalar
forvalue t=`=t_min_XX'/`=T_max_XX'{
	sum N`=increase_XX'_t_pl_`i'_XX_w if time_XX==`t'
	if !missing(r(mean)) scalar N`=increase_XX'_`i'_pl_XX_weight = N`=increase_XX'_`i'_pl_XX_weight + r(mean) // Modif Diego
}
}

///// Creating N^1_{t,\ell,g}/N^0_{t,\ell,g} for the placebos: 
///// Variable counting number of groups \ell periods away from switch at t, 
///// and with same D_{g,1} and trends_nonparam.

bys time_XX d_sq_XX `trends_nonparam': gegen N`=increase_XX'_t_placebo_`i'_g_XX=total(dist_to_switch_pl_`i'_wXX)

///// Creating all the adjustment terms to compute estimators with controls, and their variances 
if "`controls'" != ""{

// Initialize intermediate Variable needed later
capture drop part2_pl_switch`=increase_XX'_`i'_XX
gen part2_pl_switch`=increase_XX'_`i'_XX=0

local count_controls=0

// Computing the long differences of the control variables (X_g_t - X_g_t-l)

foreach var of varlist `controls'{

	local count_controls=`count_controls'+1

	capture drop diff_X_`count_controls'_placebo_`i'_XX

	xtset group_XX time_XX
	gen diff_X_`count_controls'_placebo_`i'_XX = L`=2*`i''.`var' - L`i'.`var'

	// Computing N_g_t * (X_g_t - X_g_t-l)

	capture drop diff_X`count_controls'_pl_`i'_N_XX
	gen diff_X`count_controls'_pl_`i'_N_XX = N_gt_XX * diff_X_`count_controls'_placebo_`i'_XX 

	foreach l of local levels_d_sq_XX { // index l corresponds to d in the paper 

		capture drop m`=increase_XX'_pl_g_`count_controls'_`l'_`i'_XX 
		gen m`=increase_XX'_pl_g_`count_controls'_`l'_`i'_XX = (`i' <= T_g_XX-2 & d_sq_int_XX == `l')* (G_XX / N`=increase_XX'_placebo_`i'_XX) * ([dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * diff_X`count_controls'_pl_`i'_N_XX)
		capture drop m_pl`=increase_XX'_`l'_`count_controls'_`i'_XX
		bys group_XX: gegen m_pl`=increase_XX'_`l'_`count_controls'_`i'_XX=total(m`=increase_XX'_pl_g_`count_controls'_`l'_`i'_XX)
		bys group_XX: replace m_pl`=increase_XX'_`l'_`count_controls'_`i'_XX = . if _n != 1
		capture drop M_pl`=increase_XX'_`l'_`count_controls'_`i'_XX 
		egen M_pl`=increase_XX'_`l'_`count_controls'_`i'_XX = total(m_pl`=increase_XX'_`l'_`count_controls'_`i'_XX)
		replace M_pl`=increase_XX'_`l'_`count_controls'_`i'_XX = (1/G_XX)*M_pl`=increase_XX'_`l'_`count_controls'_`i'_XX
	}

	foreach l of local levels_d_sq_XX {
		if (scalar(useful_res_`l'_XX)>1){ 
			replace diff_y_pl_`i'_XX = diff_y_pl_`i'_XX - coefs_sq_`l'_XX[`=`count_controls'',1]*diff_X_`count_controls'_placebo_`i'_XX if d_sq_int_XX==`l' 

			capture drop in_brackets_pl_`l'_`count_controls'_XX	
			gen in_brackets_pl_`l'_`count_controls'_XX=0
							
		}		
	}
}
}


///// Computing the variables used for the demeaning of outcome's placebo long difference diff_y_`i',
///// which we will use to create the U_g^{var} variables, which are used to compute 
///// the placebos' variances. 

// Generate placebo long-differences of outcomes time N_{g,t}, and dummy for (g,t) such that
// diff_y_`i' non missing and N_gt non missing.

capture drop diff_y_pl_`i'_N_gt_XX
gen diff_y_pl_`i'_N_gt_XX=N_gt_XX*diff_y_pl_`i'_XX
capture drop dof_ns_pl_`i'_XX
gen dof_ns_pl_`i'_XX = N_gt_XX!=0&diff_y_pl_`i'_XX!=.&never_change_d_pl_`i'_XX==1&N`=increase_XX'_t_placebo_`i'_XX>0&N`=increase_XX'_t_placebo_`i'_XX!=.
capture drop dof_s_pl_`i'_XX
gen dof_s_pl_`i'_XX = N_gt_XX!=0&dist_to_switch_pl_`i'_XX==1

// For never switchers, demeaning wrt to cohorts defined by D_{g,1}, `trends_nonparam' 
//(\mathcal{D}_k in companion paper)

// Modif Clément: we need to add by time
* Mean's denominator
bys d_sq_XX `trends_nonparam' time_XX : gegen count_cohort_pl_`i'_ns_t_XX=total(N_gt_XX) if dof_ns_pl_`i'_XX == 1

* Mean's numerator
bys d_sq_XX `trends_nonparam' time_XX : gegen total_cohort_pl_`i'_ns_t_XX=total(diff_y_pl_`i'_N_gt_XX) if dof_ns_pl_`i'_XX == 1

* Mean 
gen mean_cohort_pl_`i'_ns_t_XX=total_cohort_pl_`i'_ns_t_XX/count_cohort_pl_`i'_ns_t_XX

* Counting number of groups for DOF adjustment
{
	if "`cluster'" == "" {
		bys d_sq_XX `trends_nonparam' time_XX : gegen dof_cohort_pl_`i'_ns_t_XX=total(dof_ns_pl_`i'_XX) if dof_ns_pl_`i'_XX == 1
	}
	else {
		cap drop cluster_dof_pl_`i'_ns_XX
		gen cluster_dof_pl_`i'_ns_XX = `cluster' if dof_ns_pl_`i'_XX == 1
		bys d_sq_XX `trends_nonparam' time_XX: gegen dof_cohort_pl_`i'_ns_t_XX = nunique(cluster_dof_pl_`i'_ns_XX) if !missing(cluster_dof_pl_`i'_ns_XX)
	}
}
// End modif Clément

// For switchers, for placebos we no longer need to distinguish depending on whether the option
// less_conservative_se specified or not, we always demean wrt to 
// cohorts defined by D_{g,1}, F_g, D_{g,F_g}, `trends_nonparam' (\mathcal{C}_k in companion paper).

* Mean's denominator
bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam' : gegen count_cohort_pl_`i'_s_t_XX=total(N_gt_XX) if dof_s_pl_`i'_XX==1

* Mean's numerator
bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam' : gegen total_cohort_pl_`i'_s_t_XX=total(diff_y_pl_`i'_N_gt_XX) if dof_s_pl_`i'_XX==1 

* Mean 
gen mean_cohort_pl_`i'_s_t_XX=total_cohort_pl_`i'_s_t_XX/count_cohort_pl_`i'_s_t_XX	
		
* Counting number of groups for DOF adjustment
{
	if "`cluster'" == "" {
		bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam' : gegen dof_cohort_pl_`i'_s_t_XX=total(dof_s_pl_`i'_XX) if dof_s_pl_`i'_XX==1
	}
	else {
		cap drop cluster_dof_pl_`i'_s_XX
		gen cluster_dof_pl_`i'_s_XX = `cluster' if dof_s_pl_`i'_XX == 1
		bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam': gegen dof_cohort_pl_`i'_s_t_XX = nunique(cluster_dof_pl_`i'_s_XX) if !missing(cluster_dof_pl_`i'_s_XX)
	}
}

// Modif Clément 19/6/2025:
// If a switcher is the only one in their cohort or if a not-yet-switcher is the only one in their cohort, we demean wrt union of switchers and not-yet switchers.

cap drop dof_ns_s_pl_`i'_XX
cap drop count_cohort_pl_`i'_ns_s_t_XX
cap drop total_cohort_pl_`i'_ns_s_t_XX
cap drop mean_cohort_pl_`i'_ns_s_t_XX
cap drop dof_cohort_pl_`i'_ns_s_t_XX

gen dof_ns_s_pl_`i'_XX=(dof_s_pl_`i'_XX==1|dof_ns_pl_`i'_XX ==1)

* Mean's denominator
bys d_sq_XX `trends_nonparam' time_XX : gegen count_cohort_pl_`i'_ns_s_t_XX=total(N_gt_XX) if dof_ns_s_pl_`i'_XX==1

* Mean's numerator
bys d_sq_XX `trends_nonparam' time_XX: gegen total_cohort_pl_`i'_ns_s_t_XX=total(diff_y_pl_`i'_N_gt_XX) if dof_ns_s_pl_`i'_XX==1

* Mean 
gen mean_cohort_pl_`i'_ns_s_t_XX=total_cohort_pl_`i'_ns_s_t_XX/count_cohort_pl_`i'_ns_s_t_XX

* Counting number of groups for DOF adjustment
{
	if "`cluster'" == "" {
		bys d_sq_XX `trends_nonparam' time_XX: gegen dof_cohort_pl_`i'_ns_s_t_XX=total(dof_ns_s_pl_`i'_XX) if dof_ns_s_pl_`i'_XX==1
	}
	else {
		cap drop cluster_dof_pl_`i'_ns_s_XX
		gen cluster_dof_pl_`i'_ns_s_XX = `cluster' if dof_ns_s_pl_`i'_XX==1
		bys d_sq_XX `trends_nonparam' time_XX: gegen dof_cohort_pl_`i'_ns_s_t_XX = nunique(cluster_dof_pl_`i'_ns_s_XX) if !missing(cluster_dof_pl_`i'_ns_s_XX)
	}
}	

///// From those parts, generate variables for the demeaning and the DOF adjustment 
// E_hat_(g,t), defined from parts depending on the cohort definition 
gen E_hat_gt_pl_`i'_XX=0 if (time_XX<F_g_XX|F_g_XX-1+`i'==time_XX)
replace E_hat_gt_pl_`i'_XX=mean_cohort_pl_`i'_ns_t_XX if (dof_cohort_pl_`i'_ns_t_XX>=2&time_XX<F_g_XX)
replace E_hat_gt_pl_`i'_XX=mean_cohort_pl_`i'_s_t_XX if (dof_cohort_pl_`i'_s_t_XX>=2&F_g_XX-1+`i'==time_XX)
replace E_hat_gt_pl_`i'_XX=mean_cohort_pl_`i'_ns_s_t_XX if dof_cohort_pl_`i'_ns_s_t_XX>=2&(F_g_XX-1+`i'==time_XX&dof_cohort_pl_`i'_s_t_XX==1)|(time_XX<F_g_XX&dof_cohort_pl_`i'_ns_t_XX==1)

// DOF_(g,t) for DOF adjustement, defined from parts depending on the cohort definition 
// Diego - 02-03-25: when there is only 1 switcher, dof_cohort_`i'_s_t_XX = 1, hence the denominator in the expression below is 0
// The fraction is undefined and Stata puts it to missing
gen DOF_gt_pl_`i'_XX=1 if (time_XX<F_g_XX|F_g_XX-1+`i'==time_XX)
replace DOF_gt_pl_`i'_XX=  sqrt(dof_cohort_pl_`i'_s_t_XX/(dof_cohort_pl_`i'_s_t_XX-1)) if (F_g_XX-1+`i'==time_XX & dof_cohort_pl_`i'_s_t_XX > 1)
replace DOF_gt_pl_`i'_XX=  sqrt(dof_cohort_pl_`i'_ns_t_XX/(dof_cohort_pl_`i'_ns_t_XX-1)) if (time_XX<F_g_XX & dof_cohort_pl_`i'_ns_t_XX > 1)
replace DOF_gt_pl_`i'_XX=  sqrt(dof_cohort_pl_`i'_ns_s_t_XX/(dof_cohort_pl_`i'_ns_s_t_XX-1)) if dof_cohort_pl_`i'_ns_s_t_XX>=2&(F_g_XX-1+`i'==time_XX&dof_cohort_pl_`i'_s_t_XX==1)|(time_XX<F_g_XX&dof_cohort_pl_`i'_ns_t_XX==1)

// End modif Clément 19/6/2025.

////////// 6. Computing U_Gg_\ell variables for the placebos (similar to part 4, less commented)

gen dummy_U_Gg_pl_`i'_XX = (`i'<=T_g_XX-1)

if (N`=increase_XX'_placebo_`i'_XX!=0){

gen U_Gg_pl_`i'_temp_XX = dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * (time_XX>=`=`i'+1'&time_XX<=T_g_XX)* N_gt_XX * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] *diff_y_pl_`i'_XX 

bysort group_XX : gegen U_Gg_placebo_`i'_XX=total(U_Gg_pl_`i'_temp_XX)

replace U_Gg_placebo_`i'_XX= U_Gg_placebo_`i'_XX*first_obs_by_gp_XX

capture drop count`i'_pl_core_XX
gen count`i'_pl_core_XX=0
replace count`i'_pl_core_XX= N_gt_XX if (U_Gg_pl_`i'_temp_XX!=.&U_Gg_pl_`i'_temp_XX!=0|(U_Gg_pl_`i'_temp_XX==0&diff_y_pl_`i'_XX==0&(dist_to_switch_pl_`i'_XX!=0|(N`=increase_XX'_t_placebo_`i'_g_XX!=0&never_change_d_pl_`i'_XX!=0))))

gen U_Gg_pl_`i'_temp_var_XX =0

replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * N_gt_XX * DOF_gt_pl_`i'_XX *(diff_y_pl_`i'_XX-E_hat_gt_pl_`i'_XX)

if "`controls'"!=""{

levelsof d_sq_int_XX, local(levels_d_sq_XX)
foreach l of local levels_d_sq_XX {	
if (scalar(useful_res_`l'_XX)>1){ // MODIF FELIX
	
capture drop combined_pl`=increase_XX'_temp_`l'_`i'_XX	
gen combined_pl`=increase_XX'_temp_`l'_`i'_XX=0
		
forvalues j=1/`count_controls'{
forvalues k=1/`count_controls'{		
	
capture drop in_brackets_pl_`l'_`j'_`k'_temp_XX
//// CHANGE BELOW - tks + Changes Diego 27-03-25: remove N_c_`l'_XX from formula below, since we have adjusted inv_Denom above
gen in_brackets_pl_`l'_`j'_`k'_temp_XX = inv_Denom_`l'_XX[`j',`k'] * in_sum_`k'_`l'_XX * (d_sq_int_XX == `l' & F_g_XX>=3)
////
replace in_brackets_pl_`l'_`j'_XX=in_brackets_pl_`l'_`j'_XX+in_brackets_pl_`l'_`j'_`k'_temp_XX
}

replace in_brackets_pl_`l'_`j'_XX=in_brackets_pl_`l'_`j'_XX - coefs_sq_`l'_XX[`j',1]

capture drop combined_pl`=increase_XX'_temp_`l'_`j'_`i'_XX 
gen combined_pl`=increase_XX'_temp_`l'_`j'_`i'_XX=M_pl`=increase_XX'_`l'_`j'_`i'_XX*in_brackets_pl_`l'_`j'_XX

replace combined_pl`=increase_XX'_temp_`l'_`i'_XX=combined_pl`=increase_XX'_temp_`l'_`i'_XX+combined_pl`=increase_XX'_temp_`l'_`j'_`i'_XX
} 

replace part2_pl_switch`=increase_XX'_`i'_XX=part2_pl_switch`=increase_XX'_`i'_XX+combined_pl`=increase_XX'_temp_`l'_`i'_XX
} // MODIF FELIX
} // end loop oveer l

}	

bys group_XX: gegen U_Gg_pl_`i'_var_XX=total(U_Gg_pl_`i'_temp_var_XX)

// Modeif Felix 21.03.2025 (adjust order of adding additional term and summing across t)
if "`controls'"!=""{
	if `=increase_XX'==1{
	replace U_Gg_pl_`i'_var_XX=U_Gg_pl_`i'_var_XX - part2_pl_switch1_`i'_XX 
	}

	if `=increase_XX'==0{
	//// CHANGE BELOW - tks (change plus to minus since final sign flip won't be until later)
	replace U_Gg_pl_`i'_var_XX=U_Gg_pl_`i'_var_XX - part2_pl_switch0_`i'_XX 	
	}
}

}

////////// 7. Computing adjustements for the normalized and trends_lin options for placebos
////////// (similar to part 4, not commented) 

if "`normalized'"!=""{
	
	capture drop sum_treat_until_`i'_pl_XX
	capture drop delta_D_pl_`i'_cum_temp_XX

	if `continuous'==0{
	bys group_XX : gegen sum_treat_until_`i'_pl_XX = total(treatment_XX - d_sq_XX) if time_XX>=F_g_XX&time_XX<=F_g_XX-1+`i'& S_g_XX==increase_XX
	} // FELIX 18.03.2025 (delete the diff_y_pl_`i'_XX!=. condition)

	else if `continuous'>0{
		bys group_XX : gegen sum_treat_until_`i'_pl_XX = total(treatment_XX_orig - d_sq_XX_orig) if time_XX>=F_g_XX&time_XX<=F_g_XX-1+`i'&S_g_XX==increase_XX
	} // FELIX 18.03.2025 (delete the diff_y_pl_`i'_XX!=. condition)

gen delta_D_pl_`i'_cum_temp_XX = N_gt_XX/N`=increase_XX'_placebo_`i'_XX*[sum_treat_until_`i'_pl_XX* S_g_XX + (1-S_g_XX)*(-sum_treat_until_`i'_pl_XX)] if dist_to_switch_pl_`i'_XX==1

	sum  delta_D_pl_`i'_cum_temp_XX
	scalar delta_norm_pl_`i'_XX = r(sum)
}



} // End of the loop over l_placebo_u_a_XX.
}

scalar Ntrendslin_pl=1
forvalue i=1/`=l_placebo_u_a_XX' {
scalar Ntrendslin_pl=min(Ntrendslin_pl , N`=increase_XX'_placebo_`i'_XX )
}


if "`trends_lin'"!="" & Ntrendslin_pl != 0 {
	
	capture drop U_Gg_pl_`=l_placebo_u_a_XX'_TL
	capture drop U_Gg_pl_`=l_placebo_u_a_XX'_var_TL
	
	gen U_Gg_pl_`=l_placebo_u_a_XX'_TL = 0
	gen U_Gg_pl_`=l_placebo_u_a_XX'_var_TL = 0

	forvalue i=1/`=l_placebo_u_a_XX'{
		replace U_Gg_pl_`=l_placebo_u_a_XX'_TL = U_Gg_pl_`=l_placebo_u_a_XX'_TL + U_Gg_placebo_`i'_XX
		replace U_Gg_pl_`=l_placebo_u_a_XX'_var_TL = U_Gg_pl_`=l_placebo_u_a_XX'_var_TL + U_Gg_pl_`i'_var_XX
	
	}
	
	replace U_Gg_placebo_`=l_placebo_u_a_XX'_XX = U_Gg_pl_`=l_placebo_u_a_XX'_TL
	replace U_Gg_pl_`=l_placebo_u_a_XX'_var_XX =U_Gg_pl_`=l_placebo_u_a_XX'_var_TL

}
} //End of condition checking if at least one placebo requested.


////////// 8. Computing Average Total Effect estimator

///// Average total effect not estimated if trends_lin option requested:
if "`trends_lin'"==""{

///// Computing the sum of the N1_`i'_XX for the weights w. 
scalar sum_N`=increase_XX'_l_XX = 0
forvalue i=1/`=l_u_a_XX'{
	scalar sum_N`=increase_XX'_l_XX = sum_N`=increase_XX'_l_XX + N`=increase_XX'_`i'_XX
}	

///// Dropping/Initializing needed variables
capture drop U_Gg_XX
capture drop U_Gg_num_XX
capture drop U_Gg_den_XX
capture drop U_Gg_num_var_XX
capture drop  U_Gg_var_XX
gen U_Gg_num_XX=0
gen U_Gg_den_XX=0
gen U_Gg_num_var_XX=0


forvalue i=1/`=l_u_a_XX'{
	
if (N`=increase_XX'_`i'_XX!=0){
capture drop delta_D_`i'_temp_XX
capture drop delta_D_`i'_XX
		 
///// Computing the weights w_+,l
scalar w_`i'_XX = N`=increase_XX'_`i'_XX / sum_N`=increase_XX'_l_XX
		
///// Computing the delta_D_+,l", which enter in the denominator of \hat{\delta}_+/\hat{\delta}_-.

if `continuous'==0{
gen delta_D_`i'_temp_XX = N_gt_XX/N`=increase_XX'_`i'_XX*[(treatment_XX-d_sq_XX)* S_g_XX + (1-S_g_XX)*(d_sq_XX-treatment_XX)] if distance_to_switch_`i'_XX==1
}
// Redefine this with original treatment if continuous option specified (treatment was binarized and staggerized)
else if `continuous'>0{
gen delta_D_`i'_temp_XX = N_gt_XX/N`=increase_XX'_`i'_XX*[(treatment_XX_orig-d_sq_XX_orig)* S_g_XX + (1-S_g_XX)*(d_sq_XX_orig-treatment_XX_orig)] if distance_to_switch_`i'_XX==1	
}	

// For the average cumulative effect we want to keep delta_D_`i'_temp_XX and rescale it by N`=increase_XX'_`i'_XX/N_gt_XX to get an unweightened version of [D_{g,Fg-1+l} - D_{g,1}]	
// Maybe call it delta_D_g_`i'_XX instead of delta_D_`i'_temp_XX and then keep it?	
capture drop delta_D_g_`i'_XX
replace delta_D_`i'_temp_XX=0 if delta_D_`i'_temp_XX==.
gegen delta_D_`i'_XX = total(delta_D_`i'_temp_XX)
gen delta_D_g_`i'_XX=delta_D_`i'_temp_XX * (N`=increase_XX'_`i'_XX/N_gt_XX)
drop delta_D_`i'_temp_XX	
	
///// Computing the numerator of U^+_{G,g}: summing up the U_{G,g,l}s, after weighting them
bys group_XX : replace U_Gg_num_XX = U_Gg_num_XX + w_`i'_XX * U_Gg`i'_XX

///// Computing the numerator of the "alternative" U_{G,g}s for the variance : summing up the U_{G,g,l}_vars, after weighting them
bys group_XX : replace U_Gg_num_var_XX = U_Gg_num_var_XX + w_`i'_XX * U_Gg`i'_var_XX

///// Computing the denominator of U^+_{G,g}: summing up the delta^D_{+,l}s, after weighting them
bys group_XX : replace U_Gg_den_XX = U_Gg_den_XX + w_`i'_XX * delta_D_`i'_XX

}

}
		
///// Computing the U_+_(G,g)s".
bys group_XX : gen U_Gg_XX = U_Gg_num_XX/U_Gg_den_XX

///// Computing the U^+_{G,g}_vars.
bys group_XX : gen U_Gg_var_XX = U_Gg_num_var_XX/U_Gg_den_XX 

}
}
 // end of the quietly condition
	
end


