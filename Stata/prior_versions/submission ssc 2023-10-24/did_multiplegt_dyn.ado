*﻿* did_multiplegt_dyn : did_multiplegt, robust to dynamic effects and with asymptotic variances
** This version : September 19th, 2023
** Modified and commented by Clément.

********************************************************************************
*                                 PROGRAM 1                                    *
********************************************************************************

capture program drop did_multiplegt_dyn

program did_multiplegt_dyn, eclass
	version 12.0
	syntax varlist(min=4 max=4 numeric) [if] [in] [, effects(integer 1) placebo(integer 0) switchers(string) controls(varlist numeric) trends_nonparam(varlist numeric) weight(varlist numeric) drop_larger_lower NORMalized cluster(varlist numeric) graphoptions(string) SAVe_results(string) graph_off same_switchers effects_equal  drop_if_d_miss_before_first_switch ]
	
// Modif Felix:	
**** Check if gtools is installed, if not install it 
cap which gtools
if _rc ssc install gtools	

// Note from the gegen helpfile:
/*
The functions listed below have been compiled and hence will run very quickly. Functions not listed here hash the data and then call egen with by(varlist) set to the hash, which is often faster than calling egen directly, but not always.  Natively supported functions should always be faster, however.

-> So even if it is not optimized for every function this should not be a problem
*/
	
****The path of the initial dataset
local dataset_name_XX `c(filename)'
//>>
preserve


	qui{
				
capture drop outcome_XX
capture drop group_XX
capture drop time_XX
capture drop treatment_XX
capture drop d_sq_XX
capture drop diff_from_sq_XX
capture drop ever_change_d_XX
capture drop never_change_d_XX
capture drop temp_F_g_XX
capture drop F_g_XX
capture drop t_min_XX
capture drop T_max_XX
capture drop N_gt_XX
capture drop T_g_XX
capture drop min_F_g_XX
capture drop controls_time_XX
capture drop _fillin
capture drop treat_not_missing_XX
capture drop time_treat_not_miss_XX
capture drop time_d_nonmiss_XX
capture drop min_time_d_nonmiss_XX
capture drop max_treat_not_miss_XX
capture drop d_missing_startinghere_XX
capture drop d_missing_untilhere_XX
capture drop min_treat_not_miss_XX
capture drop first_switch_unknown_date_t_XX
capture drop first_switch_unknown_date_XX
capture drop d_F_g_temp_XX
capture drop d_F_g_XX
capture drop S_g_XX
capture drop L_g_XX
capture drop U_Gg_plus_XX
capture drop U_Gg_minus_XX
capture drop U_Gg_global_XX
capture drop w_plus_XX
capture drop first_obs_by_gp_XX
capture drop first_obs_by_clust_XX
capture drop avg_post_switch_treat_XX
capture drop count_time_post_switch_XX_temp
capture drop count_time_post_switch_XX
capture drop diff_from_sq_XX
capture drop ever_strict_increase_XX
capture drop ever_strict_decrease_XX
capture drop U_Gg_var_plus_XX
capture drop U_Gg_var_minus_XX
capture drop U_Gg_var_global_XX
capture drop U_Gg_var_global_2_XX
capture drop weight_XX
capture drop delta_XX


ereturn clear

tokenize `varlist'

///// Patching the cluster option: by default, command clusters at group level. 
///// If user specified clustering according to group, replace it by empty. 
if "`cluster'"=="`2'" local cluster ""

///// Selecting the sample 

//dropping observations with missing group or time
drop if `2'==.|`3'==.
//dropping observations with missing controls
	if "`controls'" !=""{
	foreach var of varlist `controls'{
	drop if `var'==.
	}
	}
//dropping observations not included in the if condition
	if "`if'" !=""{
	keep `if'
	}
// dropping observations not included in the in condition
	if "`in'" !=""{
	keep `in'
	} 
// dropping observations with a missing clustering variable
	if "`cluster'" !=""{
	drop if `cluster'==.
	}
//dropping groups with always missing treatment or outcome
capture drop mean_d_XX
capture drop mean_y_XX
bys `2': gegen mean_d_XX=mean(`4')
bys `2': gegen mean_y_XX=mean(`1')
drop if mean_d_XX==.|mean_y_XX==.
drop mean_d_XX mean_y_XX
	
///// Checking wether data has to be collapsed, because at a more disaggregated level than group*time. 
capture drop counter_XX
capture drop counter_temp_XX
gen counter_temp_XX=1
bys `2' `3' : gegen counter_XX=count(counter_temp_XX)
sum counter_XX
scalar aggregated_data=0
if r(max)==1{
scalar aggregated_data=1
}

///// Creating the weight variable.
if("`weight'"==""){
gen weight_XX = 1
}
else{
gen weight_XX = `weight'
}

replace weight_XX=0 if weight_XX==.

 ///// Collapsing the data when necessary

if scalar(aggregated_data)==0{
	
	replace weight_XX=0 if `1'==.
	
	if "`1'"!="`4'"{
				
		collapse (mean) `1' `4' `controls' `trends_nonparam' `cluster' /*`recat_treatment'*/ (count) weight_XX [pw=weight_XX], by(`2' `3')
	}
	
	if "`1'"=="`4'"{
				
		collapse (mean) `1' `controls' `trends_nonparam' `cluster' /*`recat_treatment'*/ (count) weight_XX [pw=weight_XX], by(`2' `3')
	}	
}
 
////// Creating variables that will be helpful to handle imbalanced panels//
gegen time_XX=group(`3')

gen time_d_nonmiss_XX=time_XX if `4'!=.
bys `2': gegen min_time_d_nonmiss_XX=min(time_d_nonmiss_XX)
capture drop max_time_d_nonmiss_XX
bys `2': gegen max_time_d_nonmiss_XX=max(time_d_nonmiss_XX)
capture drop time_y_nonmiss_XX
capture drop min_time_y_nonmiss_XX
gen time_y_nonmiss_XX=time_XX if `1'!=.
bys `2': gegen min_time_y_nonmiss_XX=min(time_y_nonmiss_XX)
capture drop time_d_miss_XX
capture drop min_time_d_miss_aft_ynm_XX
gen time_d_miss_XX=time_XX if `4'==.&time_XX>=min_time_y_nonmiss_XX
bys `2': gegen min_time_d_miss_aft_ynm_XX=min(time_d_miss_XX)
drop time_d_nonmiss_XX time_y_nonmiss_XX time_d_miss_XX

////// Creating baseline treatment: D_{g,1} in paper, redefined to account for imbalanced panels: g's treatment at first period where g's treatment non missing//

gen d_sq_XX_temp=`4' if time_XX==min_time_d_nonmiss_XX
bys `2': gegen d_sq_XX=mean(d_sq_XX_temp)
drop d_sq_XX_temp 

// Modif Doulo: Create a new variable containing the integer levels of d_sq_XX, in the case of non-integer averages of d_sq_XX.
capture drop d_sq_int_XX
gegen d_sq_int_XX = group(d_sq_XX)


////// If the option drop_larger_lower was specified, drop (g,t) cells such that at t, g has experienced both a strictly lower and a strictly higher treatment than its baseline treatment. //

gen diff_from_sq_XX=(`4'-d_sq_XX)

if "`drop_larger_lower'"!=""{		
	sort `2' `3'
	gen ever_strict_increase_XX=(diff_from_sq_XX>0&`4'!=.) 
	replace ever_strict_increase_XX=1 if ever_strict_increase_XX[_n-1]==1&`2'==`2'[_n-1]

	gen ever_strict_decrease_XX=(diff_from_sq_XX<0&`4'!=.) 
	replace ever_strict_decrease_XX=1 if ever_strict_decrease_XX[_n-1]==1&`2'==`2'[_n-1]

	drop if ever_strict_increase_XX==1 & ever_strict_decrease_XX==1
	
	drop ever_strict_increase_XX ever_strict_decrease_XX
		
}

////// Sort data //
sort `2' `3'

////// Creating Y D G variables //
gen outcome_XX=`1'
gegen group_XX=group(`2')
gen treatment_XX=`4'

////// Counting number of groups //
sum group_XX
scalar G_XX=r(max)

////// Ever changed treatment //
gen ever_change_d_XX=(abs(diff_from_sq_XX)>0&treatment_XX!=.)
replace ever_change_d_XX=1 if ever_change_d_XX[_n-1]==1&group_XX==group_XX[_n-1]

////// Date of the first treatment change //
sort group_XX time_XX
gen temp_F_g_XX=time_XX if ever_change_d_XX==1&ever_change_d_XX[_n-1]==0 
replace temp_F_g_XX=0 if temp_F_g_XX==.
bys group_XX: gegen F_g_XX=max(temp_F_g_XX)
drop temp_F_g_XX

///// Dropping the values of the baseline treatment such that no variance in F_g within those values.//
capture drop var_F_g_XX
bys d_sq_XX `trends_nonparam': gegen var_F_g_XX=sd(F_g_XX)
drop if var_F_g_XX==0
drop var_F_g_XX

count
if r(N)==0{
	di as error "No treatment effect can be estimated."
	di as error "This is because Assumption 1 in"
	di as error "de Chaisemartin & D'Haultfoeuille (2023)"
	di as error "is not satisfied in the data used for"
	di as error "estimation, given the options requested."

	exit
}

// Modif Mélitine : this to be put higher, so that T_max_XX is realistic
////// For each value of the baseline treatment, we drop time periods such that we do not have any control with that baseline treatment after that period //
//// Watch out: this means the panel is no longer balanced, though it is balanced within values of the baseline treatment //
gen never_change_d_XX=1-ever_change_d_XX
bys time_XX d_sq_XX `trends_nonparam': gegen controls_time_XX=max(never_change_d_XX)
drop if controls_time_XX==0

////// Computing t_min_XX, T_max_XX, and replacing F_g_XX by last period plus one for those that never change treatment //
sum time_XX
scalar t_min_XX=r(min)
scalar T_max_XX=r(max)
replace F_g_XX=T_max_XX+1 if F_g_XX==0 // defining F_g as T+1 for never-treated units


////// Dealing with missing treatments: most conservative option ///

*Let FMD_g denote the first date when g's treatment is missing while y has been not missing at least once, so that we know for sure that g already exists. If that date is before the first period when g's treatment changes, we do not know when g's treatment has changed for the first time. Then, a conservative option is to drop all of g's outcomes starting at FMD_g.

if "`drop_if_d_miss_before_first_switch'"!=""{
replace outcome_XX=. if min_time_d_miss_aft_ynm_XX<F_g_XX&time_XX>=min_time_d_miss_aft_ynm_XX
}

////// Dealing with missing treatments: most liberal option ///

*Let FD_g and LD_g respectively denote the first and last period where a group's treatment is non missing. Let FY_g denote the first period where a group's outcome is non missing.

*For groups that experience at least one treatment change, let LDBF_g denote the last date before F_g where g's treatment is non missing. We have FD_g<=LDBF_g<F_g<=LD_g, and we will deal with missing treatments depending on when they occur with respect to those four dates.

capture drop last_obs_D_bef_switch_t_XX
capture drop last_obs_D_bef_switch_XX
gen last_obs_D_bef_switch_t_XX=time_XX if time_XX<F_g_XX&treatment_XX!=.
bys group_XX: gegen last_obs_D_bef_switch_XX=max(last_obs_D_bef_switch_t_XX)

*For groups that do not experience a treatment change, we just have FD_g<=LD_g, and we will deal with missing treatments depending on when they occur with respect to those two dates.

*For t<FD_g, by default we are going to consider that g joins the panel at FD_g: any non-missing outcome before FD_g replaced as missing, but all non-missing outcomes after FD_g are kept. For groups such that FY_g<FD_g, this is a "liberal" convention: those groups exist before FD_g, so one could argue that their status quo treatment is missing and they should be dropped from the analysis. We give the user the option to do that, with drop_if_d_miss_before_first_switch option
  
replace outcome_XX=. if time_XX<min_time_d_nonmiss_XX 

*For groups that experience a treatment change, if D_gt missing at FD_g<t<LDBF_g, we replace their missing treatment by their status-quo treatment. Again, this is a liberal convention, so we give the user the option to not use those observations, with drop_if_d_miss_before_first_switch option.

replace treatment_XX=d_sq_XX if F_g_XX<T_max_XX+1&treatment_XX==.&time_XX<last_obs_D_bef_switch_XX&time_XX>min_time_d_nonmiss_XX

*For groups that experience a treatment change, if D_gt missing at LDBF_g<t<F_g (equivalent to LDBF_g<F_g-1), we cannot know the exact date when their treatment has changed, even in a binary and staggered design. Therefore, we set their outcomes at missing starting at LDBF_g+1. We also redefine their F_g as T+1 because they are effectively control groups. We also define the trunc_control_XX as LDBF_g+1 for them, because they can only be used as controls till that date.

replace outcome_XX=. if F_g_XX<T_max_XX+1&time_XX>last_obs_D_bef_switch_XX&last_obs_D_bef_switch_XX<F_g_XX-1
capture drop trunc_control_XX
gen trunc_control_XX=last_obs_D_bef_switch_XX+1 if F_g_XX<T_max_XX+1&last_obs_D_bef_switch_XX<F_g_XX-1
replace F_g_XX=T_max_XX+1 if F_g_XX<T_max_XX+1&last_obs_D_bef_switch_XX<F_g_XX-1

*For groups that experience a treatment change, if D_gt missing at F_g<t, we replace their missing treatment by D(g,F_g). This is again a liberal convention, but it is innocuous for the reduced-form parameters DID_l, so we do not give the user the option to overrule it (Note that overruling it could make the same_switchers option fail)

gen d_F_g_temp_XX = treatment_XX if time_XX==F_g_XX
bys group_XX: gegen d_F_g_XX = mean(d_F_g_temp_XX)
replace treatment_XX=d_F_g_XX if F_g_XX<T_max_XX+1&treatment_XX==.&time_XX>F_g_XX&last_obs_D_bef_switch_XX==F_g_XX-1

*For groups that do not experience a treatment change, if D_gt missing at FD_g<t<LD_g, we replace their missing treatment by D_g1. This is again a liberal convention, so we give the user the option to not use those observations, wi1th drop_if_d_miss_before_first_switch option.

replace treatment_XX=d_sq_XX if F_g_XX==T_max_XX+1&treatment_XX==.&time_XX>min_time_d_nonmiss_XX&time_XX<max_time_d_nonmiss_XX

*For groups that do not experience a treatment change, we replace all their outcomes by missing at t>LD_g. Even in a binary and staggered design, we cannot infer their treatment at t>LD_g.

replace outcome_XX=. if F_g_XX==T_max_XX+1&time_XX>max_time_d_nonmiss_XX
replace trunc_control_XX=max_time_d_nonmiss_XX+1 if F_g_XX==T_max_XX+1

////// Balancing the panel //
// Modification proposed by Diego, if there is already a variable name _merge in the user's data set 
cap rename _merge _merge_og

fillin group_XX time_XX
capture drop d_sq_XX_new
bys group_XX: gegen d_sq_XX_new=mean(d_sq_XX)
drop d_sq_XX
rename d_sq_XX_new d_sq_XX


////// Defining N_gt, the weight of each cell (g,t) //
gen N_gt_XX=1
replace N_gt_XX = N_gt_XX*weight_XX  //
replace N_gt_XX=0 if outcome_XX==.|treatment_XX==.

////// Determining T_g, last period where there is still a group with the same treatment as g's in period 1 and whose treatment has not changed since the start of the panel. Definition adapted from the paper, to account for imbalanced panel //

capture drop F_g_trunc_XX
gen F_g_trunc_XX=F_g_XX
replace F_g_trunc_XX=min(F_g_XX,trunc_control_XX) if trunc_control_XX!=.
bys d_sq_XX `trends_nonparam': gegen T_g_XX = max(F_g_trunc_XX) 

replace T_g_XX = T_g_XX-1

////// Defining S_g, an indicator variable for groups whose average post switch treatment value is larger than their initial value of treatment. They will be considered switchers in. If S_g==0, that means the group is a switcher out. For never-switchers, S_g will be undefined. //
bys group_XX : gegen avg_post_switch_treat_XX_temp = total(treatment_XX) if time_XX>=F_g_XX&time_XX<=T_g_XX
gen count_time_post_switch_XX_temp=(treatment_XX!=.) if time_XX>=F_g_XX&time_XX<=T_g_XX


bysort group_XX : gegen count_time_post_switch_XX=total(count_time_post_switch_XX_temp)

replace avg_post_switch_treat_XX_temp =avg_post_switch_treat_XX_temp/count_time_post_switch_XX

bys group_XX: gegen avg_post_switch_treat_XX=mean(avg_post_switch_treat_XX_temp)
drop avg_post_switch_treat_XX_temp

// When a group is a switching group, but its average post-treatment treatment value is exactly equal to its baseline treatment, we cannnot classify it as a swicher in or a switcher out, but it is not a control either. As such, we drop it from the estimation. Those groups are referred to as no-first-stage-switchers.
drop if avg_post_switch_treat_XX==d_sq_XX&F_g_XX!=T_g_XX+1

gen S_g_XX=(avg_post_switch_treat_XX>d_sq_XX) if F_g_XX!=T_max_XX+1
// Here, S_g_XX==. for never-switchers groups only.

///// Creating the variable L_g_XX = T_g_XX - F_g_XX, so that we can compute L_u or L_a afterwards. //
gen L_g_XX = T_g_XX - F_g_XX +1

///// Creating the equivalent variable L_g_placebo_XX for the placebos. //
if `placebo'!=0{
capture drop L_g_placebo_XX
gen L_g_placebo_XX = min(L_g_XX, F_g_XX-2) if F_g_XX>=3

}

///// Flagging first observation of each group_XX //
sort group_XX time_XX
bysort group_XX : gen first_obs_by_gp_XX = (_n==1)

// Flagging first observation of each cluster, and checking clustering variable weakly coarser than group variable. //
if "`cluster'"!=""{
	
	capture drop cluster_var_g_XX
	
	bysort `cluster' : gen first_obs_by_clust_XX = (_n==1)
	
	//Error message for clustering: non-nested case
	bysort group_XX: gegen cluster_var_g_XX = sd(`cluster')
	sum cluster_var_g_XX
	scalar max_cluster_var_XX = `r(max)'
	if (scalar(max_cluster_var_XX) >0){
	di as error "The group variable should be nested within the clustering variable."
	exit
	}
	
}

////// Declaring data as panel //
xtset group_XX time_XX

///// Creating first-differenced outcome and treatment//
capture drop diff_y_XX
capture drop diff_d_XX
gen diff_y_XX = d.outcome_XX
g diff_d_XX = d.treatment_XX

///// Dealing with controls : computing first differences of the controls and necessary variables for the residualization of the outcome variable //

if "`controls'" !=""{
local count_controls=0

capture drop fd_X_all_non_missing_XX
gen fd_X_all_non_missing_XX=1

foreach var of varlist `controls'{

local count_controls=`count_controls'+1

capture drop diff_X`count_controls'_XX
capture drop avg_diff_X`count_controls'_XX
capture drop resid_X`count_controls'_time_FE_XX

// Computing the first differences of the control variables
xtset group_XX time_XX
gen diff_X`count_controls'_XX=D.`var'
replace fd_X_all_non_missing_XX=0 if diff_X`count_controls'_XX==.
}

local count_controls=0
local mycontrols_XX ""
local prod_controls_y ""

foreach var of varlist `controls'{
local count_controls=`count_controls'+1

///// First step of the residualization : computing \Delta X_{.,t}: average of controls' first difference at time t, among groups whose treatment has not changed yet.

capture drop sum_weights_control_XX //So as to consider the weighted regressions. Do not move the capture drop outside the loop

bys time_XX d_sq_XX `trends_nonparam' : gegen sum_weights_control_XX = total(N_gt_XX) if ever_change_d_XX==0&diff_y_XX!=.&fd_X_all_non_missing_XX==1
bys time_XX d_sq_XX `trends_nonparam' : gegen avg_diff_X`count_controls'_XX = total(N_gt_XX*diff_X`count_controls'_XX) if ever_change_d_XX==0&diff_y_XX!=.&fd_X_all_non_missing_XX==1

bys time_XX d_sq_XX `trends_nonparam' : replace avg_diff_X`count_controls'_XX = avg_diff_X`count_controls'_XX/sum_weights_control_XX

// Computing \Delta\Dot{X}_{g,t}, the difference between the first differences of covariates and the average of their first-difference, which gives us the residuals of a regression of covariates on time fixed effects. 
//Multiply by sqrt(N_gt_XX) to replicate weighted regression

gen resid_X`count_controls'_time_FE_XX = sqrt(N_gt_XX)*(diff_X`count_controls'_XX - avg_diff_X`count_controls'_XX)

replace resid_X`count_controls'_time_FE_XX=0 if resid_X`count_controls'_time_FE_XX==.

// Storing the obtained residuals for the computation of theta_d
local mycontrols_XX "`mycontrols_XX' resid_X`count_controls'_time_FE_XX"
//

// Generating the product between \Delta\Dot{X}_{g,t} and \Delta Y_{g,t}
//Multiply by sqrt(N_gt_XX) to replicate weighted regression
capture drop prod_X`count_controls'_diff_y_temp_XX
capture drop prod_X`count_controls'_diff_y_XX
capture drop diff_y_wXX
gen diff_y_wXX = sqrt(N_gt_XX)*diff_y_XX

gen prod_X`count_controls'_diff_y_temp_XX = resid_X`count_controls'_time_FE_XX*diff_y_wXX if time_XX>=2&time_XX<F_g_XX

replace prod_X`count_controls'_diff_y_temp_XX = 0 if prod_X`count_controls'_diff_y_temp_XX ==.

// Computing the sum for each group to obtain the term \sum_{t=2}^{F_g-1}*N_{g,t}*\Delta \Dot{X}_{g,t}* \Delta Y_{g,t}
bys group_XX: gegen prod_X`count_controls'_diff_y_XX = total(prod_X`count_controls'_diff_y_temp_XX)

}

// Creating a local storing the status quos for which Denom_d is not defined 
local store_singular_XX ""

local store_noresidualization_XX ""

////// Storing the different possible values of the status quos //
//levelsof d_sq_XX, local(levels_d_sq_XX)

// Modif Doulo: use the new variable (d_sq_int_XX) containing integers levels of d_sq_XX, and replace all d_sq_XX==`l' by d_sq_int_XX==`l' 

levelsof d_sq_int_XX, local(levels_d_sq_XX)


foreach l of local levels_d_sq_XX {
	tempfile data_XX
	save "`data_XX'.dta", replace
// A baseline treatment is relevant iff it is taken by at least two groups with different values of F_g_XX
// and non-missing diff_y_XX, otherwise we do not need to perform the residualization for this specific baseline treatment.

scalar store_singular_`l'_XX = 0 
//tab F_g_XX if d_sq_XX==`l' 
tab F_g_XX if d_sq_int_XX==`l'  

scalar useful_res_`l'_XX = `r(r)'
	if (scalar(useful_res_`l'_XX)>1){

// Isolate the observations used in the computation of theta_d

keep if ever_change_d_XX==0&diff_y_XX!=.&fd_X_all_non_missing_XX==1&d_sq_int_XX==`l'


// Using the matrix accum function, to regress the first difference of outcome on the first differences of covariates. We will obtain the vectors of coefficients \theta_d s, where d indexes values of the baseline treatment.

capture matrix accum overall_XX = diff_y_wXX `mycontrols_XX'
scalar rc_XX=_rc

if scalar(rc_XX)!=0{
//local store_singular_XX = "`store_singular_XX' `l'" //Moved to line 515
scalar store_singular_`l'_XX = 1
local store_noresidualization_XX "`store_noresidualization_XX' `l'"
scalar useful_res_`l'_XX=1
}

else {
	
// Isolate the parts of the matrix we need
matrix didmgt_XX = overall_XX[2..`=`count_controls'+1',2..`=`count_controls'+1']
matrix didmgt_Xy = overall_XX[2..`=`count_controls'+1',1]

// Computing the vectors of coefficients \theta_d s for each value of the baseline treatment
matrix coefs_sq_`l'_XX = invsym(didmgt_XX)*didmgt_Xy

		// Computing the matrix Denom^{-1}
			//Check first if the matrix is invertible, invsym() inverts a matrix even if it is singular
			capture drop scalar det_XX
			scalar det_XX = det(didmgt_XX)

		if (abs(scalar(det_XX))<=10^(-16)){ 
					//local store_singular_XX = "`store_singular_XX' `l'" ////Moved to line 515
					scalar store_singular_`l'_XX = 1
					scalar drop det_XX
		}

		// if the matrix is invertible (i.e. det_XX!=0), compute Denom.
			matrix inv_Denom_`l'_XX = invsym(didmgt_XX)*G_XX
}

//>>
//restore
use "`data_XX'.dta", clear

}

}
//Modif Doulo: Fill up store_singular_XX, with correct values of statu quo and not the levels
levelsof d_sq_XX, local(levels_d_sq_bis_XX)
scalar index_sing_XX = 0
foreach l of local levels_d_sq_bis_XX {
scalar index_sing_XX = scalar(index_sing_XX)+1
if(scalar(store_singular_`=index_sing_XX'_XX) == 1){
local store_singular_XX = "`store_singular_XX' `l'"
}
}


//Display errors if one of the Denoms is not defined
if ("`store_singular_XX'"!=""){
	di as error "Some control variables are not taken into account for groups with baseline treatment equal to: [`store_singular_XX']"
	di as error "This may occur in the following situations:"
	di as error "1. For groups with those values of the baseline treatment,"
	di as error "the regression of the outcome first difference on the controls' first differences "
	di as error "and time fixed effects has fewer observations than variables."
	di as error "Note that for each value of the baseline treatment,"
	di as error "those regressions are estimated among (g,t)s such that g has not changed treatment yet at t."
	di as error "2. For groups with those values of the baseline treatment, "
	di as error "two or more of your control variables are perfectly collinear "
	di as error "in the sample where the regression is run, for instance because those control variables do not vary over time."
}

// Values of baseline treatment such that residualization could not be performed at all are dropped.
foreach l of local store_noresidualization_XX {
drop if d_sq_int_XX==`l'
}

} // end of the if "`controls'" !="" condition

///// Computing T_u/T_a, and L_u/L_a, to compare them to the number of effects requested //
// Initializing L_u_XX and L_a_XX with default values, otherwise some error messages are produced when they enter if conditions.
scalar L_u_XX=.
scalar L_a_XX=.

// For switchers in
if "`switchers'"==""|"`switchers'"=="in"{ 
sum L_g_XX if S_g_XX==1
scalar L_u_XX=r(max)
*For placebos
if `placebo'!=0{
sum L_g_placebo_XX if S_g_XX==1
scalar L_placebo_u_XX=r(max)
}
}

// For switchers out
if "`switchers'"==""|"`switchers'"=="out"{
sum L_g_XX if S_g_XX==0
scalar L_a_XX=r(max)
*For placebos
if `placebo'!=0{
sum L_g_placebo_XX if S_g_XX==0
scalar L_placebo_a_XX=r(max)
}
}

///// Error message if Assumption 1 in dCDH (2023) fails and no effect can be computed //

if ("`switchers'"=="in"&(L_u_XX==.|L_u_XX==0))|("`switchers'"=="out"&(L_a_XX==.|L_a_XX==0))|("`switchers'"==""&(L_u_XX==.|L_u_XX==0)&(L_a_XX==.|L_a_XX==0)){
		di as error "No treatment effect can be estimated."
		di as error "This is because Assumption 1 in"
		di as error "de Chaisemartin & D'Haultfoeuille (2023)"
		di as error "is not satisfied in the data used for"
		di as error "estimation, given the options requested."
		
		exit
	}

// Computing the final number of dynamic effects to be estimated. If the number asked by the user was too large, display an error message.
if "`switchers'"==""{
scalar l_XX=max(L_a_XX, L_u_XX)
scalar l_XX=min(l_XX, `effects')

if `placebo'!=0{

scalar l_placebo_XX=max(L_placebo_a_XX, L_placebo_u_XX)
scalar l_placebo_XX=min(l_placebo_XX, `placebo')
// The number of placebos computed cannot be greater than the number of effects computed:
scalar l_placebo_XX=min(l_placebo_XX, `effects')

}
else{
	scalar l_placebo_XX=0
}

}

if "`switchers'"=="in"{
	scalar l_XX=min(`effects', L_u_XX)
	if `placebo'!=0{
		scalar l_placebo_XX=min(`placebo', L_placebo_u_XX)
		// The number of placebos computed cannot be greater than the number of effects computed:
		scalar l_placebo_XX=min(l_placebo_XX, `effects') 
		}

else{
	scalar l_placebo_XX=0
}

}

if "`switchers'"=="out"{
	scalar l_XX=min(`effects', L_a_XX)
		if `placebo'!=0{
			scalar l_placebo_XX=min(`placebo', L_placebo_a_XX)
			// The number of placebos computed cannot be greater than the number of effects computed:
			scalar l_placebo_XX=min(l_placebo_XX, `effects')
		}

else{
	scalar l_placebo_XX=0
}

}


if l_XX<`effects'{
	di as error "The number of effects requested is too large."
	di as error "The number of effects which can be estimated is at most " l_XX "."
	di as error "The command will therefore try to estimate " l_XX " effect(s)."
}

if `placebo'!=0{
	if l_placebo_XX<`placebo'&`effects'>=`placebo'{
		di as error "The number of placebos which can be estimated is at most " l_placebo_XX "."
		di as error "The command will therefore try to estimate " l_placebo_XX " placebo(s)."
	}

	if `effects'<`placebo'{
		di as error "The number of placebo requested cannot be larger than the number of effects requested."
		di as error "The command cannot compute more than " l_placebo_XX " placebo(s)."
	}
}

///// Generating default values for the variables which will be aggregated. If there is a subsequent estimation, variables will be set equal to their estimated values instead. If there is no subsequent estimation, they will remain equal to zero but will not be missing. //
forvalue i=1/`=l_XX'{
	
	capture drop U_Gg`i'_plus_XX 
	capture drop U_Gg`i'_minus_XX
	capture drop count`i'_plus_XX
	capture drop count`i'_minus_XX
	capture drop U_Gg`i'_global_XX
	capture drop count`i'_global_XX
	capture drop N_effect_`i'_XX
	capture drop DID_`i'_XX
	capture drop U_Gg_var_`i'_in_XX
	capture drop U_Gg_var_`i'_out_XX
	capture drop U_Gg_var_glob_`i'_XX
	capture drop U_Gg_var_glob_eff`i'_sqrd_XX
	
	gen U_Gg`i'_plus_XX = 0
	gen U_Gg`i'_minus_XX = 0
	gen count`i'_plus_XX = 0
	gen count`i'_minus_XX = 0
	scalar N1_`i'_XX=0
	scalar N0_`i'_XX=0
	gen U_Gg_var_`i'_in_XX=0
	gen U_Gg_var_`i'_out_XX=0
	scalar sum_for_var_in_XX=0
	scalar sum_for_var_out_XX=0
	
	if "`normalized'"!=""{
		scalar delta_D_`i'_in_XX = 0
		scalar delta_D_`i'_out_XX = 0
	}
}

if l_placebo_XX!=0{

	forvalue i=1/`=l_placebo_XX'{	
	
	capture drop U_Gg_pl_`i'_plus_XX  
	capture drop U_Gg_pl_`i'_global_XX
	capture drop U_Gg_pl_`i'_minus_XX 
	capture drop DID_placebo_`i'_XX
	capture drop N_placebo_`i'_XX
	capture drop U_Gg_var_pl_`i'_in_XX
	capture drop U_Gg_var_pl_`i'_out_XX
	capture drop U_Gg_var_glob_pl_`i'_XX
	capture drop U_Gg_var_glob_pl_`i'_2_XX
	capture drop count`i'_pl_plus_XX
	capture drop  count`i'_pl_minus_XX
	capture drop count`i'_pl_global_XX
	
	gen U_Gg_pl_`i'_plus_XX   = 0
	gen U_Gg_pl_`i'_minus_XX   = 0
	scalar N1_placebo_`i'_XX=0
	scalar N0_placebo_`i'_XX=0
	gen U_Gg_var_pl_`i'_in_XX=0
	gen U_Gg_var_pl_`i'_out_XX=0
	scalar sum_for_var_placebo_in_XX=0
	scalar sum_for_var_placebo_out_XX=0
	gen count`i'_pl_plus_XX = 0
	gen count`i'_pl_minus_XX = 0

	if "`normalized'"!=""{
		scalar delta_D_pl_`i'_in_XX = 0
		scalar delta_D_pl_`i'_out_XX = 0
	}
	
}
}

gen U_Gg_plus_XX = 0
gen U_Gg_minus_XX = 0
scalar U_Gg_den_plus_XX = 0
scalar U_Gg_den_minus_XX = 0
scalar sum_N1_l_XX = 0
scalar sum_N0_l_XX = 0
gen U_Gg_var_plus_XX = 0
gen U_Gg_var_minus_XX = 0

////// Perform here the estimation, i.e. call the program did_multiplegt_dyn_core. //

// For switchers in
if ("`switchers'"==""|"`switchers'"=="in"){
if L_u_XX!=.&L_u_XX!=0{
	did_multiplegt_dyn_core outcome_XX group_XX time_XX treatment_XX, effects(`=l_XX') placebo(`=l_placebo_XX') switchers_core(in) controls(`controls') trends_nonparam(`trends_nonparam') `normalized' `same_switchers' `effects_equal'
	
// Store the results
forvalue i=1/`=l_u_a_XX'{
/////////// NB: in the case of unbalanced panels, it can happen that the U_Gg`i'_XX are not computed by program 2 (for example when y is missing). Consequently, for the command not to display an error message and continue running, we need to verify the variable is created, which is conditional on  N1_`i'_XX!=0.
if N1_`i'_XX!=0{
		replace U_Gg`i'_plus_XX = U_Gg`i'_XX
		replace count`i'_plus_XX= count`i'_core_XX
		replace U_Gg_var_`i'_in_XX=U_Gg`i'_var_XX
		
		if "`normalized'"!=""{
			scalar delta_D_`i'_in_XX = delta_norm_`i'_XX
		}
	}

}

if l_placebo_XX!=0{

	forvalue i=1/`=l_placebo_u_a_XX'{	
		if N1_placebo_`i'_XX!=0{
				replace U_Gg_pl_`i'_plus_XX  = U_Gg_placebo_`i'_XX
				replace count`i'_pl_plus_XX= count`i'_pl_core_XX
				replace U_Gg_var_pl_`i'_in_XX=U_Gg_pl_`i'_var_XX
				
				if "`normalized'"!=""{
					scalar delta_D_pl_`i'_in_XX = delta_norm_pl_`i'_XX
				}

		}
		
	}	
}
	
	if sum_N1_l_XX!=0{
	replace U_Gg_plus_XX = U_Gg_XX
	scalar U_Gg_den_plus_XX=U_Gg_den_XX
	replace U_Gg_var_plus_XX = U_Gg_var_XX
	}

}
	
} // end of the loop for switchers in

// For switchers out
if ("`switchers'"==""|"`switchers'"=="out"){
if L_a_XX!=.&L_a_XX!=0{
did_multiplegt_dyn_core outcome_XX group_XX time_XX treatment_XX, effects(`=l_XX') placebo(`=l_placebo_XX') switchers_core(out) controls(`controls')  trends_nonparam(`trends_nonparam') `normalized' `same_switchers' `effects_equal'

	
// Store the results
forvalue i=1/`=l_u_a_XX'{
if N0_`i'_XX!=0{
		replace U_Gg`i'_minus_XX = - U_Gg`i'_XX
		
		replace count`i'_minus_XX= count`i'_core_XX
		replace U_Gg_var_`i'_out_XX=U_Gg`i'_var_XX
		
		if "`normalized'"!=""{
			scalar delta_D_`i'_out_XX = delta_norm_`i'_XX
		}
		
	}
	
}


if l_placebo_XX!=0{

	forvalue i=1/`=l_placebo_u_a_XX'{
	if N0_placebo_`i'_XX!=0{
			replace U_Gg_pl_`i'_minus_XX  = -U_Gg_placebo_`i'_XX
			replace count`i'_pl_minus_XX= count`i'_pl_core_XX
			replace U_Gg_var_pl_`i'_out_XX=U_Gg_pl_`i'_var_XX
			
			if "`normalized'"!=""{
					scalar delta_D_pl_`i'_out_XX = delta_norm_pl_`i'_XX
	}

	}
	}		
}

	if sum_N0_l_XX!=0{
	replace U_Gg_minus_XX = - U_Gg_XX
	scalar U_Gg_den_minus_XX=U_Gg_den_XX
	replace U_Gg_var_minus_XX = U_Gg_var_XX
	}

}
	
} // end of the loop for switchers out 

///// Aggregating the results for switchers in and out. //

****************************DID_l*********************************

//Creation of the matrix which stores all the estimators (DID_l, DID_pl, delta, etc.), their sd and the CIs

matrix mat_res_XX = J(l_XX+l_placebo_XX+1,7,.) 

// Computing first the "global" U_Ggl s
forvalue i=1/`=l_XX'{
	
gen U_Gg`i'_global_XX = (N1_`i'_XX/(N1_`i'_XX+N0_`i'_XX))*U_Gg`i'_plus_XX +(N0_`i'_XX/(N1_`i'_XX+N0_`i'_XX))*U_Gg`i'_minus_XX
replace U_Gg`i'_global_XX=. if first_obs_by_gp_XX==0

gen count`i'_global_XX = max(count`i'_plus_XX, count`i'_minus_XX)

	if "`normalized'"!=""{
	scalar delta_D_`i'_global_XX = (N1_`i'_XX/(N1_`i'_XX+N0_`i'_XX))*delta_D_`i'_in_XX+(N0_`i'_XX/(N1_`i'_XX+N0_`i'_XX))*delta_D_`i'_out_XX
	}

// Storing the results and the number of switchers into the matrix mat_res_XX
local rownames "`rownames' Effect_`i'"

* Number of switchers
scalar N_switchers_effect_`i'_XX=N1_`i'_XX+N0_`i'_XX
matrix mat_res_XX[`i',6]=N_switchers_effect_`i'_XX
matrix mat_res_XX[`i',7]=`i'
ereturn scalar N_switchers_effect_`i' = N_switchers_effect_`i'_XX
* Number of observations used in the estimation
gegen N_effect_`i'_XX = total(count`i'_global_XX) 
scalar N_effect_`i'_XX = N_effect_`i'_XX
ereturn scalar N_effect_`i' = N_effect_`i'_XX
matrix mat_res_XX[`i',5]=N_effect_`i'_XX

if N_switchers_effect_`i'_XX==0|N_effect_`i'_XX==0{
	di as error "Effect_"`i' " cannot be estimated."
	di as error "There is no switcher or no control"
	di as error "for this effect."
}

* DID_l
gegen DID_`i'_XX = total(U_Gg`i'_global_XX) 
replace  DID_`i'_XX = DID_`i'_XX/G_XX
scalar DID_`i'_XX = DID_`i'_XX

if "`normalized'"!=""{
	scalar DID_`i'_XX = scalar(DID_`i'_XX)/scalar(delta_D_`i'_global_XX)
}

if ("`switchers'"==""&N1_`i'_XX==0&N0_`i'_XX==0)|("`switchers'"=="out"&N0_`i'_XX==0)|("`switchers'"=="in"&N1_`i'_XX==0){
	scalar DID_`i'_XX=.
}

ereturn scalar Effect_`i' = scalar(DID_`i'_XX)
matrix mat_res_XX[`i',1]= scalar(DID_`i'_XX) 

}

****************************Average Total Effect*********************************

// Computing the weight w_+.

if "`switchers'"==""{
	scalar w_plus_XX = U_Gg_den_plus_XX*sum_N1_l_XX/(U_Gg_den_plus_XX*sum_N1_l_XX+U_Gg_den_minus_XX*sum_N0_l_XX)
}

if "`switchers'"=="out"{
	scalar w_plus_XX=0
}

if "`switchers'"=="in"{
	scalar w_plus_XX=1
}

// Aggregating to obtain the "global" U_Gg s
gen U_Gg_global_XX = w_plus_XX*U_Gg_plus_XX +(1-w_plus_XX)*U_Gg_minus_XX
replace U_Gg_global_XX=. if first_obs_by_gp_XX==0

gegen delta_XX = total(U_Gg_global_XX)
replace delta_XX = delta_XX/G_XX
scalar delta_XX = delta_XX
ereturn scalar Av_tot_effect = scalar(delta_XX)

// Completing the results matrix
* Storing the results
matrix mat_res_XX[l_XX+1,1]= scalar(delta_XX)
local rownames "`rownames' Av_tot_eff"
* Number of switchers
scalar N_switchers_effect_XX=0
forvalue i=1/`=l_XX'{
	scalar N_switchers_effect_XX = N_switchers_effect_XX + N_switchers_effect_`i'_XX
}
matrix mat_res_XX[l_XX+1,6]=N_switchers_effect_XX
matrix mat_res_XX[l_XX+1,7]=0
ereturn scalar N_switchers_effect_average = N_switchers_effect_XX
 
* Number of observations used in the estimation

capture drop count_global_XX
capture drop N_effect_XX
gen count_global_XX=0
forvalue i=1/`=l_XX'{
	replace count_global_XX= max(count_global_XX, count`i'_global_XX)
}
sum count_global_XX
scalar N_effect_XX=r(sum)
matrix mat_res_XX[l_XX+1,5]=scalar(N_effect_XX)
ereturn scalar N_avg_total_effect = scalar(N_effect_XX)

****************************Placebos****************************
// Computing first the "global" U_Ggl s for the placebos

	
if l_placebo_XX!=0{

forvalue i=1/`=l_placebo_XX'{	

gen U_Gg_pl_`i'_global_XX = (N1_placebo_`i'_XX/(N1_placebo_`i'_XX+N0_placebo_`i'_XX))*U_Gg_pl_`i'_plus_XX  +(N0_placebo_`i'_XX/(N1_placebo_`i'_XX+N0_placebo_`i'_XX))*U_Gg_pl_`i'_minus_XX 
replace U_Gg_pl_`i'_global_XX=. if first_obs_by_gp_XX==0

gen count`i'_pl_global_XX=max(count`i'_pl_plus_XX, count`i'_pl_minus_XX)

if "`normalized'"!=""{
	scalar delta_D_pl_`i'_global_XX = (N1_placebo_`i'_XX/(N1_placebo_`i'_XX+N0_placebo_`i'_XX))*delta_D_pl_`i'_in_XX+(N0_placebo_`i'_XX/(N1_placebo_`i'_XX+N0_placebo_`i'_XX))*delta_D_pl_`i'_out_XX
		}

// Summing them to obtain the DID_pl

gegen DID_placebo_`i'_XX = total(U_Gg_pl_`i'_global_XX) 
replace  DID_placebo_`i'_XX= DID_placebo_`i'_XX/G_XX
scalar DID_placebo_`i'_XX = DID_placebo_`i'_XX


if "`normalized'"!=""{
	scalar DID_placebo_`i'_XX = scalar(DID_placebo_`i'_XX)/scalar(delta_D_pl_`i'_global_XX)
}

if ("`switchers'"==""&N1_placebo_`i'_XX==0&N0_placebo_`i'_XX==0)|("`switchers'"=="out"&N0_placebo_`i'_XX==0)|("`switchers'"=="in"&N1_placebo_`i'_XX==0){
	scalar DID_placebo_`i'_XX=.
}

ereturn scalar Placebo_`i' = scalar(DID_placebo_`i'_XX)

// Storing the results and the number of switchers into the matrix mat_res_XX
* Storing the results 
matrix mat_res_XX[`=l_XX'+ 1 + `i',1]=scalar(DID_placebo_`i'_XX)
local rownames "`rownames' Placebo_`i'"
* Number of switchers
scalar N_switchers_placebo_`i'_XX=N1_placebo_`i'_XX+N0_placebo_`i'_XX
matrix mat_res_XX[`=l_XX'+ 1 + `i',6]=N_switchers_placebo_`i'_XX
matrix mat_res_XX[`=l_XX'+ 1 + `i',7]= `=-`i''
ereturn scalar N_switchers_placebo_`i' = N_switchers_placebo_`i'_XX
* Number of observations used in the estimation
gegen N_placebo_`i'_XX = total(count`i'_pl_global_XX)
scalar N_placebo_`i'_XX = N_placebo_`i'_XX
ereturn scalar N_placebo_`i' = N_placebo_`i'_XX
matrix mat_res_XX[`=l_XX' + 1 + `i',5]=N_placebo_`i'_XX

if N_switchers_placebo_`i'_XX==0|N_placebo_`i'_XX==0{
	di as error "Placebo_"`i' " cannot be estimated."
	di as error "There is no switcher or no control"
	di as error "for this placebo."
}
}
}
//End of filling in the placebos matrix results

*****  Estimating the asymptotic variances *******************

///// Estimating \hat{\sigma}^2_l and the confidence intervals //
forvalue i=1/`=l_XX'{
if ("`switchers'"==""&(N1_`i'_XX!=0|N0_`i'_XX!=0))|("`switchers'"=="out"&N0_`i'_XX!=0)|("`switchers'"=="in"&N1_`i'_XX!=0){

// Estimating \hat{\sigma}^2_l

gen U_Gg_var_glob_`i'_XX = U_Gg_var_`i'_in_XX * (scalar(N1_`i'_XX)/(scalar(N1_`i'_XX)+scalar(N0_`i'_XX))) + U_Gg_var_`i'_out_XX* (scalar(N0_`i'_XX)/(scalar(N1_`i'_XX)+scalar(N0_`i'_XX)))


if "`cluster'"==""{
	gen U_Gg_var_glob_eff`i'_sqrd_XX = U_Gg_var_glob_`i'_XX^2*first_obs_by_gp_XX

	sum U_Gg_var_glob_eff`i'_sqrd_XX
	scalar sum_for_var_`i'_XX=r(sum)/G_XX^2
}


if "`cluster'"!=""{
		
	capture drop clust_U_Gg_var_glob_`i'_XX
	capture drop clust_U_Gg_var_glob_`i'_2_XX
	
	replace U_Gg_var_glob_`i'_XX = U_Gg_var_glob_`i'_XX*first_obs_by_gp_XX
	
	bys `cluster' : gegen clust_U_Gg_var_glob_`i'_XX = total(U_Gg_var_glob_`i'_XX)
	
	gen clust_U_Gg_var_glob_`i'_2_XX=clust_U_Gg_var_glob_`i'_XX^2*first_obs_by_clust_XX
	
	sum clust_U_Gg_var_glob_`i'_2_XX
	scalar sum_for_var_`i'_XX=r(sum)/G_XX^2
	
	replace U_Gg_var_glob_`i'_XX = clust_U_Gg_var_glob_`i'_XX
	
}

	scalar se_`i'_XX = sqrt(sum_for_var_`i'_XX)


if "`normalized'"!=""{
	scalar se_`i'_XX = scalar(se_`i'_XX)/scalar(delta_D_`i'_global_XX)
}

matrix mat_res_XX[`i',2]=scalar(se_`i'_XX)
ereturn scalar se_effect_`i'=scalar(se_`i'_XX)
	
// Lower bound of the 95% confidence interval
scalar LB_CI_`i'_XX = scalar(DID_`i'_XX) - 1.96*scalar(se_`i'_XX)
matrix mat_res_XX[`i',3]= scalar(LB_CI_`i'_XX)

// Upper bound of the 95% confidence interval
scalar UB_CI_`i'_XX = scalar(DID_`i'_XX) + 1.96*scalar(se_`i'_XX)
matrix mat_res_XX[`i',4]=scalar(UB_CI_`i'_XX)


}
}


****************************Placebos****************************
///// Estimating \hat{\sigma}^2_pl and the confidence intervals //


if l_placebo_XX!=0{

//Take into account the case where same_switchers leads to N0 or N1 = 0 for placebos
forvalue i=1/`=l_placebo_XX'{
if ("`switchers'"==""&(N1_placebo_`i'_XX!=0|N0_placebo_`i'_XX!=0))|("`switchers'"=="out"&N0_placebo_`i'_XX!=0)|("`switchers'"=="in"&N1_placebo_`i'_XX!=0){ 

// Estimating \hat{\sigma}^2_l

gen U_Gg_var_glob_pl_`i'_XX = U_Gg_var_pl_`i'_in_XX * (N1_placebo_`i'_XX/(N1_placebo_`i'_XX+N0_placebo_`i'_XX)) + U_Gg_var_pl_`i'_out_XX* (N0_placebo_`i'_XX/(N1_placebo_`i'_XX+N0_placebo_`i'_XX))

if "`cluster'"==""{

gen U_Gg_var_glob_pl_`i'_2_XX = U_Gg_var_glob_pl_`i'_XX^2*first_obs_by_gp_XX

sum U_Gg_var_glob_pl_`i'_2_XX
scalar sum_for_var_placebo_`i'_XX=r(sum)/G_XX^2
}

if "`cluster'"!=""{
	capture drop clust_U_Gg_var_glob_pl_`i'_XX
	capture drop clust_U_Gg_var_glob_pl_`i'_2_XX

	
	replace U_Gg_var_glob_pl_`i'_XX = U_Gg_var_glob_pl_`i'_XX*first_obs_by_gp_XX
	
	bys `cluster' : gegen clust_U_Gg_var_glob_pl_`i'_XX= total(U_Gg_var_glob_pl_`i'_XX)
	
	gen clust_U_Gg_var_glob_pl_`i'_2_XX=clust_U_Gg_var_glob_pl_`i'_XX^2*first_obs_by_clust_XX
	
	sum clust_U_Gg_var_glob_pl_`i'_2_XX
	scalar sum_for_var_placebo_`i'_XX=r(sum)/G_XX^2
	
	replace U_Gg_var_glob_pl_`i'_XX = clust_U_Gg_var_glob_pl_`i'_XX  
}


scalar se_placebo_`i'_XX = sqrt(sum_for_var_placebo_`i'_XX)


if "`normalized'"!=""{
	scalar se_placebo_`i'_XX = scalar(se_placebo_`i'_XX)/scalar(delta_D_pl_`i'_global_XX)
}

matrix mat_res_XX[`=l_XX' + 1 + `i',2]=scalar(se_placebo_`i'_XX)
ereturn scalar se_placebo_`i'=scalar(se_placebo_`i'_XX)
	
// Lower bound of the 95% confidence interval
scalar LB_CI_placebo_`i'_XX = scalar(DID_placebo_`i'_XX) - 1.96*scalar(se_placebo_`i'_XX)
matrix mat_res_XX[`=l_XX'+ 1 + `i', 3]= scalar(LB_CI_placebo_`i'_XX)

// Upper bound of the 95% confidence interval
scalar UB_CI_placebo_`i'_XX = scalar(DID_placebo_`i'_XX) + 1.96*scalar(se_placebo_`i'_XX)
matrix mat_res_XX[`=l_XX'+ 1 + `i',4]=scalar(UB_CI_placebo_`i'_XX)


}
}
}


****************************Average Effect****************************
///// Estimating \hat{\sigma}^2 and the confidence interval for the average effect //
if ("`switchers'"==""&(sum_N1_l_XX!=0|sum_N0_l_XX!=0))|("`switchers'"=="out"&sum_N0_l_XX!=0)|("`switchers'"=="in"&sum_N1_l_XX!=0){

// Estimating \hat{\sigma}^2

gen U_Gg_var_global_XX = w_plus_XX * U_Gg_var_plus_XX + (1 - w_plus_XX) * U_Gg_var_minus_XX

if "`cluster'"==""{
capture drop U_Gg_var_global_2_XX
gen U_Gg_var_global_2_XX = U_Gg_var_global_XX^2*first_obs_by_gp_XX

sum U_Gg_var_global_2_XX
scalar sum_for_var_XX=r(sum)/G_XX^2
}


if "`cluster'"!=""{
	capture drop clust_U_Gg_var_global_XX
	
	replace U_Gg_var_global_XX=U_Gg_var_global_XX*first_obs_by_gp_XX
	
	bys `cluster' : gegen clust_U_Gg_var_global_XX= total(U_Gg_var_global_XX)
	replace clust_U_Gg_var_global_XX=clust_U_Gg_var_global_XX^2*first_obs_by_clust_XX
	
	sum clust_U_Gg_var_global_XX
	scalar sum_for_var_XX=r(sum)/G_XX^2
}


scalar se_XX = sqrt(sum_for_var_XX)
matrix mat_res_XX[l_XX+1,2]=se_XX
ereturn scalar se_avg_total_effect = se_XX

// Lower bound of the 95% confidence interval
scalar LB_CI_XX = delta_XX - 1.96*se_XX
matrix mat_res_XX[l_XX+1,3]= LB_CI_XX

// Lower bound of the 95% confidence interval
scalar UB_CI_XX = delta_XX + 1.96*se_XX
matrix mat_res_XX[l_XX+1,4]= UB_CI_XX

}

 
 // Modif Felix 
*********************Matrix for event_plot******************************
matrix didmgt_b2=mat_res_XX[1...,1..1]

matrix didmgt_var2=J(`=l_placebo_XX'+`=l_XX'+1,1,0)
forvalue i=1/`=`=l_placebo_XX'+`=l_XX'+1'{
matrix didmgt_var2[`i',1]=mat_res_XX[`i',2]^2
}

*** add reference period
local rownames2 "`rownames' Effect_0"
matrix zero=J(1,1,0)
matrix didmgt_b2=didmgt_b2\zero
matrix didmgt_var2=didmgt_var2\zero

matrix rownames didmgt_var2= `rownames2'
// For compatibility with old do files:
matrix compatibility_variances = didmgt_var2

matrix rownames didmgt_b2 = `rownames2'
matrix rownames didmgt_var2 = `rownames2'

ereturn matrix estimates=didmgt_b2 
ereturn matrix didmgt_variances=didmgt_var2
ereturn matrix variances= compatibility_variances
ereturn local cmd "did_multiplegt_dyn"

****************************F-tests*************************************

// If the option cluster is specified, we have previously replaced U_Gg_var_glob_pl_`i'_XX by clust_U_Gg_var_glob_pl_`i'_XX, and U_Gg_var_glob_`i'_XX by clust_U_Gg_var_glob_`i'_XX. 
// Now, we must also replace first_obs_by_gp_XX by first_obs_by_clust_XX
if "`cluster'"!=""{

replace first_obs_by_gp_XX=first_obs_by_clust_XX

}

///// Performing a test to see whether all placebo effects are jointly equal to 0.

scalar all_Ns_pl_not_zero=.


if (l_placebo_XX!=0)&l_placebo_XX>1{
	
		scalar all_Ns_pl_not_zero=0
	
		forvalue i=1/`=l_placebo_XX'{
		if ("`switchers'"==""&(N1_placebo_`i'_XX!=0|N0_placebo_`i'_XX!=0))|("`switchers'"=="out"&N0_placebo_`i'_XX!=0)|("`switchers'"=="in"&N1_placebo_`i'_XX!=0){
			
			scalar all_Ns_pl_not_zero=all_Ns_pl_not_zero+1
			
		}
	}
	
	if all_Ns_pl_not_zero==l_placebo_XX{
	
	matrix didmgt_Placebo=J(l_placebo_XX,1,0)
	matrix didmgt_Var_Placebo=J(l_placebo_XX,l_placebo_XX,0)
	
	forvalue i=1/`=l_placebo_XX'{
		
		matrix didmgt_Placebo[`i',1]=scalar(DID_placebo_`i'_XX)
		matrix didmgt_Var_Placebo[`i',`i']= scalar(se_placebo_`i'_XX)^2
	
		if `i'<`=l_placebo_XX'{
		forvalue j=`=`i'+1'/`=l_placebo_XX'{
			
			capture drop U_Gg_var_pl_`i'_`j'_XX
			capture drop U_Gg_var_pl_`i'_`j'_2_XX
			
			if ("`normalized'"==""){
			gen U_Gg_var_pl_`i'_`j'_XX = U_Gg_var_glob_pl_`i'_XX + U_Gg_var_glob_pl_`j'_XX
		}
		
			if "`normalized'"!=""{
			gen U_Gg_var_pl_`i'_`j'_XX = U_Gg_var_glob_pl_`i'_XX/scalar(delta_D_pl_`i'_global_XX) + U_Gg_var_glob_pl_`j'_XX/scalar(delta_D_pl_`j'_global_XX)
}

			gen U_Gg_var_pl_`i'_`j'_2_XX = U_Gg_var_pl_`i'_`j'_XX^2*first_obs_by_gp_XX
			
			sum U_Gg_var_pl_`i'_`j'_2_XX
			scalar var_sum_pla_`i'_`j'_XX=r(sum)/G_XX^2
			
			scalar cov_pl_`i'_`j'_XX = (scalar(var_sum_pla_`i'_`j'_XX) - scalar(se_placebo_`i'_XX)^2 - scalar(se_placebo_`j'_XX)^2)/2
			
			matrix didmgt_Var_Placebo[`i',`j']= scalar(cov_pl_`i'_`j'_XX)
			matrix didmgt_Var_Placebo[`j',`i']= scalar(cov_pl_`i'_`j'_XX)

		}
	}
	
	}
	
	matrix didmgt_Var_Placebo_inv=invsym(didmgt_Var_Placebo)
	matrix didmgt_Placebo_t=didmgt_Placebo'
	matrix didmgt_chi2placebo=didmgt_Placebo_t*didmgt_Var_Placebo_inv*didmgt_Placebo
	scalar p_jointplacebo=1-chi2(l_placebo_XX,didmgt_chi2placebo[1,1])
	ereturn scalar p_jointplacebo=1-chi2(l_placebo_XX,didmgt_chi2placebo[1,1])

	}
	else{
		di as error "Some placebos could not be estimated."
		di as error "Therefore, the test of joint nullity of the placebos "
		di as error "could not be computed."
	}
	
}



///// Performing a test that all the DID_l effects are equal.

scalar all_Ns_not_zero=.

if ("`effects_equal'")!=""&l_XX>1{
	
	scalar all_Ns_not_zero=0
	
	forvalue i=1/`=l_XX'{
		if ("`switchers'"==""&(N1_`i'_XX!=0|N0_`i'_XX!=0))|("`switchers'"=="out"&N0_`i'_XX!=0)|("`switchers'"=="in"&N1_`i'_XX!=0){
			
			scalar all_Ns_not_zero=all_Ns_not_zero+1
			
		}
	}
	
	
	if all_Ns_not_zero==l_XX{
	
	// Creating a vector of all coefficient estimates
	matrix didmgt_Effects = mat_res_XX[1..l_XX, 1]
	
	// Creating a matrix where the variances and the covariances of the effects will be stored.
	matrix didmgt_Var_Effects=J(l_XX, l_XX, 0)
	
	matrix didmgt_identity = J(`=l_XX-1', l_XX,0)
	
	forvalue i=1/`=l_XX'{
				if ("`switchers'"==""&(N1_`i'_XX!=0|N0_`i'_XX!=0))|("`switchers'"=="out"&N0_`i'_XX!=0)|("`switchers'"=="in"&N1_`i'_XX!=0){

		
		// Storing the already computed variances.
		matrix didmgt_Var_Effects[`i',`i']= scalar(se_`i'_XX)^2
		
		
		if `i'<`=l_XX'{
			matrix didmgt_identity[`i',`i'] =1
		}
		
		
		// Computing and storing the covariances.
		if `i'<`=l_XX'{
		forvalue j=`=`i'+1'/`=l_XX'{
			
		capture drop U_Gg_var_`i'_`j'_XX
		capture drop U_Gg_var_`i'_`j'_2_XX
		
		
		if ("`normalized'"==""){
			gen U_Gg_var_`i'_`j'_XX = U_Gg_var_glob_`i'_XX+ U_Gg_var_glob_`j'_XX
		}
		
			if "`normalized'"!=""{
			gen U_Gg_var_`i'_`j'_XX = U_Gg_var_glob_`i'_XX/scalar(delta_D_`i'_global_XX) + U_Gg_var_glob_`j'_XX/scalar(delta_D_`j'_global_XX)
}
		

			gen U_Gg_var_`i'_`j'_2_XX = U_Gg_var_`i'_`j'_XX^2*first_obs_by_gp_XX

			sum U_Gg_var_`i'_`j'_2_XX
			scalar var_sum_`i'_`j'_XX=r(sum)/G_XX^2
			scalar cov_`i'_`j'_XX = (scalar(var_sum_`i'_`j'_XX) - scalar(se_`i'_XX)^2 - scalar(se_`j'_XX)^2)/2
			
			matrix didmgt_Var_Effects[`i',`j']= scalar(cov_`i'_`j'_XX)
			matrix didmgt_Var_Effects[`j',`i']= scalar(cov_`i'_`j'_XX)	
			
			
		}
		}
	}
	}
	
	
	// Computing the vector of recentered effects and its variance matrix
	matrix didmgt_D =didmgt_identity -J(`=l_XX-1',l_XX,(1/l_XX))
	matrix didmgt_test_effects = didmgt_D*didmgt_Effects
	matrix didmgt_test_var = didmgt_D*didmgt_Var_Effects*didmgt_D'
	matrix didmgt_test_var = (didmgt_test_var + didmgt_test_var')/2
		
	// Performing the test.	
	matrix didmgt_test_var_inv=invsym(didmgt_test_var)
	matrix didmgt_test_effects_t=didmgt_test_effects'
	matrix didmgt_chi2_equal_ef=didmgt_test_effects_t*didmgt_test_var_inv*didmgt_test_effects
	scalar p_equality_effects=1-chi2(`=l_XX-1',didmgt_chi2_equal_ef[1,1])
	ereturn scalar p_equality_effects=1-chi2(`=l_XX-1',didmgt_chi2_equal_ef[1,1])
	
	}
	else{
		di as error "Some effects could not be estimated."
		di as error "Therefore, the test of equality of the effects "
		di as error "could not be computed."
	}
	
	}


	} // end of the quietly condition
	
*****  Returning the results of the estimation. ******************************

matrix rownames mat_res_XX= `rownames'
matrix colnames mat_res_XX= "Estimate" "SE" "LB CI" "UB CI" "N" "Switchers" "time_to_treat"

/*
di " "
di "----------------------------------------------"
di "        Estimation of treatment effects       "
di "----------------------------------------------"
*/

display _newline
di as input "{hline 80}"
di as input _skip(13) "Estimation of treatment effects: Event-study effects"
di as input "{hline 80}"
noisily matlist mat_res_XX[1..l_XX, 1..6]
di as input "{hline 80}"
if l_XX>1&"`effects_equal'"!=""&all_Ns_not_zero==l_XX{
di as text "{it:Test of equality of the effects : p-value =} " scalar(p_equality_effects)
//}
}



/*
di " "
di "----------------------------------------------"
di "             Average total effect             "
di "----------------------------------------------"
*/

matrix mat_res_avg_XX=mat_res_XX[l_XX+1, 1..6]
matrix mat_res_avg_XX=(mat_res_avg_XX, .z )
matrix colnames mat_res_avg_XX= "Estimate" "SE" "LB CI" "UB CI" "N" "Switch" "x Periods"
display _newline
di as input "{hline 80}"
di as input _skip(4) "Estimation of treatment effects: Average total effect per treatment unit"
di as input "{hline 80}"
noisily matlist mat_res_avg_XX, nodotz
di as input "{hline 80}"


/*
di " "
di "----------------------------------------------"
di "         Testing the parallel trends          "
di "----------------------------------------------"
*/

if l_placebo_XX!=0{

display _newline
di as input "{hline 80}"
di as input _skip(10) "Testing the parallel trends and no anticipation assumptions"
di as input "{hline 80}"
matlist mat_res_XX[l_XX+2...,1..6]
di as input "{hline 80}"
if l_placebo_XX>1&all_Ns_pl_not_zero==l_placebo_XX{
di as text "{it:Test of joint nullity of the placebos : p-value =} " scalar(p_jointplacebo)

}
}


display _newline
di as text "The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899)."



*****  Producing a graph. ****************************************************
qui {
////Doulo: LIST all the variables name without using *?
//capture drop mat_res_XX*
capture drop mat_res_XX1 mat_res_XX2 mat_res_XX3 mat_res_XX4 mat_res_XX5 mat_res_XX6 mat_res_XX7
svmat mat_res_XX //rename the obtained variables ? yes!
////Doulo: use forvalue instead of foreach


forvalue index=1/6{
replace mat_res_XX`index'=0 if mat_res_XX7==0 // artificially replacing the result of average effect by 0
}

capture drop point_estimate
capture drop se_point_estimate
capture drop lb_CI_95
capture drop up_CI_95
capture drop N
capture drop N_switchers
capture drop time_to_treat

rename mat_res_XX1 point_estimate
rename mat_res_XX2 se_point_estimate
rename mat_res_XX3 lb_CI_95
rename mat_res_XX4 up_CI_95
rename mat_res_XX5 N
rename mat_res_XX6 N_switchers
rename mat_res_XX7 time_to_treat


sort time_to_treat

}

if ("`graph_off'"==""){
if "`graphoptions'"==""{
twoway (connected point_estimate time_to_treat, lpattern(solid)) (rcap up_CI_95 lb_CI_95 time_to_treat), xlabel(`=-l_placebo_XX'[1]`=l_XX') xtitle("Relative time to last period before treatment changes (t=0)", size(large)) title("DID, from last period before treatment changes (t=0) to t", size(large)) graphregion(color(white)) plotregion(color(white)) legend(off)
}
else{
global options="`graphoptions'"
twoway (connected point_estimate time_to_treat, lpattern(solid)) (rcap up_CI_95 lb_CI_95 time_to_treat),  $options
}
}

//Saving the results if requested 
if "`save_results'"!=""{
quietly{
keep if point_estimate!=.
}
keep time_to_treat point_estimate se_point_estimate lb_CI_95 up_CI_95 N N_switchers
label data "Stores did_multiplegt_dyn estimates' information: Type 'notes' for details." 

// notes : "INFORMATION ABOUT THE DATASET"

local command_options_XX effects(`effects') placebo(`placebo') switchers(`switchers') controls(`controls') trends_nonparam(`trends_nonparam') weight(`weight') `drop_larger_lower' `normalized' cluster(`cluster') `same_switchers' 
notes : "{bf:{ul:Date of Run:}} `c(current_date)' at `c(current_time)'"
notes : "{bf:{ul:Command Syntax:}} did_multiplegt `1' `2' `3' `4' `if' `in', `command_options_XX' "
*Use the local dataset_name_XX
notes : "{bf:{ul:Path of the Used Dataset:}} `dataset_name_XX' "
save "`save_results'", replace

}


//>>
restore

end


********************************************************************************
*                                 PROGRAM 2                                    *
********************************************************************************

capture program drop did_multiplegt_dyn_core

program did_multiplegt_dyn_core, eclass
	version 12.0
	syntax varlist(min=4 max=4 numeric) [, effects(integer 1) placebo(integer 0) switchers_core(string) controls(varlist numeric) trends_nonparam(varlist numeric) NORMalized same_switchers effects_equal]
	
	qui{
		
///// Initializing depending on whether we consider switchers in or switchers out: indicating the number of dynamic effects to estimate and the index for whether we compute N^1_l or N^0_l. //
	
if "`switchers_core'"=="in"{
scalar l_u_a_XX=min(`effects', L_u_XX)

if `placebo'!=0{
	if L_placebo_u_XX!=.&L_placebo_u_XX!=0{
	scalar l_placebo_u_a_XX=min(`placebo', L_placebo_u_XX)
	}
	}

	scalar increase_XX=1
}

if "`switchers_core'"=="out"{
scalar l_u_a_XX=min(`effects', L_a_XX)

if `placebo'!=0{
	if L_placebo_a_XX!=.&L_placebo_a_XX!=0{
	scalar l_placebo_u_a_XX=min(`placebo', L_placebo_a_XX)
	}
	}

scalar increase_XX=0
}

//levelsof d_sq_XX, local(levels_d_sq_XX)
levelsof d_sq_int_XX, local(levels_d_sq_XX)

*****  Estimating the DID_{+,l}s or the DID_{-,l}s *****************************

////// Loop to estimate the dynamic effects //


forvalue i=1/`=l_u_a_XX'{

capture drop distance_to_switch_`i'_XX
capture drop never_change_d_`i'_XX
capture drop N`=increase_XX'_t_`i'_XX
capture drop N`=increase_XX'_t_`i'_g_XX
capture drop N_gt_control_`i'_XX
capture drop diff_y_`i'_XX
capture drop diff_y_`i'_XX_temp
capture drop dummy_U_Gg`i'_XX
capture drop U_Gg`i'_temp_XX
capture drop U_Gg`i'_XX
capture drop count`i'_core_XX
capture drop mean_diff_y_`i'_nd_sq_t_XX
capture drop mean_diff_y_`i'_d_sq_t_XX
capture drop U_Gg`i'_temp_var_XX
capture drop U_Gg`i'_var_XX
capture drop U_Gg`i'_var_2_XX
capture drop count_var_`i'_ntreat_XX_temp
capture drop count_var_`i'_ntreat_XX
capture drop count_var_`i'_treat_XX_temp
capture drop count_var_`i'_treat_XX
capture drop avg_diff_y_`i'_tnp_XX
capture drop count_diff_y_`i'_nd_sq_t_XX
capture drop count_diff_y_`i'_d_sq_t_XX

capture drop never_change_d_`i'_wXX
capture drop distance_to_switch_`i'_wXX

////// Creating long difference of outcome //
xtset group_XX time_XX
bys group_XX : gen diff_y_`i'_XX = outcome_XX - L`i'.outcome_XX


///// Creating long differences of control variables //
if "`controls'" != ""{
		
local count_controls=0

// Computing the first differences of the control variables
foreach var of varlist `controls'{
		
local count_controls=`count_controls'+1

capture drop diff_X`count_controls'_`i'_XX

xtset group_XX time_XX
gen diff_X`count_controls'_`i'_XX=`var' - L`i'.`var'

foreach l of local levels_d_sq_XX { 
	if (scalar(useful_res_`l'_XX)>1){ 

replace diff_y_`i'_XX = diff_y_`i'_XX - coefs_sq_`l'_XX[`=`count_controls'',1]*diff_X`count_controls'_`i'_XX if d_sq_int_XX==`l' 

////// N.B. : in the above line, we do not add "&diff_X`count_controls'_`i'_XX!=." because we want to exclude from the estimation any first/long-difference for which the covariates are missing.
}
}

}

}

////// Identifying the control (g,t)s in the estimation of dynamic effect i //
bys group_XX: gen never_change_d_`i'_XX=(F_g_XX>time_XX) if diff_y_`i'_XX!=.

////// Counting the number of controls (g,t)s at each time period //

// N^g_t  //for trends_nonparam: we need controls defined as {D_{g',1} =  D_{g,1} and Sg′ =Sg}
gen never_change_d_`i'_wXX = never_change_d_`i'_XX*N_gt_XX
bys time_XX d_sq_XX `trends_nonparam': gegen N_gt_control_`i'_XX=total(never_change_d_`i'_wXX)

///// binary variable indicating whether group g is l periods away from switch //

//same_switchers options

if ("`same_switchers'"!=""){

capture drop relevant_y_missing_XX
capture drop cum_fillin_XX
capture drop fillin_g_XX 
capture drop dum_fillin_temp_XX 
capture drop still_switcher_`i'_XX  

sort group_XX time_XX

//Generate a variable tagging the switchers that should be dropped
gen relevant_y_missing_XX=(outcome_XX==.&time_XX>=F_g_XX-1&time_XX<=F_g_XX-1+`=`effects'')
if "`controls'" != ""{
replace relevant_y_missing_XX=1 if fd_X_all_non_missing_XX==0&time_XX>=F_g_XX&time_XX<=F_g_XX-1+`=`effects''
}

bys group_XX: gen cum_fillin_XX = sum(relevant_y_missing_XX)
gen dum_fillin_temp_XX = (cum_fillin_XX==0&time_XX==F_g_XX-1+`=`effects'')
bys group_XX: gegen fillin_g_XX = total(dum_fillin_temp_XX)
gen still_switcher_`i'_XX = (F_g_XX-1+`=`effects''<=T_g_XX&fillin_g_XX>0) 

gen distance_to_switch_`i'_XX=(still_switcher_`i'_XX&time_XX==F_g_XX-1+`i'&`i'<=L_g_XX&S_g_XX==increase_XX&N_gt_control_`i'_XX>0&N_gt_control_`i'_XX!=.) if diff_y_`i'_XX!=. 
}
else{
gen distance_to_switch_`i'_XX=(time_XX==F_g_XX-1+`i'&`i'<=L_g_XX&S_g_XX==increase_XX&N_gt_control_`i'_XX>0&N_gt_control_`i'_XX!=.) if diff_y_`i'_XX!=. 
}

///// Computing N^1_{t,l} or N^0_{t,l}. //
gen distance_to_switch_`i'_wXX = distance_to_switch_`i'_XX*N_gt_XX
bys time_XX: gegen N`=increase_XX'_t_`i'_XX=total(distance_to_switch_`i'_wXX)

///// Computing N^1_l or N^0_l. //
forvalue t=`=t_min_XX'/`=T_max_XX'{
	sum N`=increase_XX'_t_`i'_XX if time_XX==`t'
	scalar N`=increase_XX'_`i'_XX = N`=increase_XX'_`i'_XX + r(mean)
}

///// Computing N^0_{t,l,g} or N^1_{t,l,g}. //
bys time_XX d_sq_XX `trends_nonparam': gegen N`=increase_XX'_t_`i'_g_XX=total(distance_to_switch_`i'_wXX)

///// Computing the mean of first or long-differences of outcomes time N_{g,t} for non-treated and for treated separately - will be useful for the computation of the variance //

capture drop diff_y_`i'_N_gt_XX
gen diff_y_`i'_N_gt_XX=N_gt_XX*diff_y_`i'_XX
capture drop dof_diff_y_`i'_N_gt_XX
gen dof_diff_y_`i'_N_gt_XX=(N_gt_XX!=0&diff_y_`i'_XX!=.)

bys time_XX d_sq_XX `trends_nonparam' : gegen mean_diff_y_`i'_nd_sq_t_XX=mean(diff_y_`i'_N_gt_XX) if never_change_d_`i'_XX==1&N`=increase_XX'_t_`i'_XX>0&N`=increase_XX'_t_`i'_XX!=.

bys time_XX d_sq_XX `trends_nonparam' : gegen count_diff_y_`i'_nd_sq_t_XX=total(dof_diff_y_`i'_N_gt_XX) if diff_y_`i'_XX!=.&never_change_d_`i'_XX==1&N`=increase_XX'_t_`i'_XX>0&N`=increase_XX'_t_`i'_XX!=.

bys time_XX d_sq_XX `trends_nonparam' : gegen mean_diff_y_`i'_d_sq_t_XX=mean(diff_y_`i'_N_gt_XX) if distance_to_switch_`i'_XX==1 

bys time_XX d_sq_XX `trends_nonparam' : gegen count_diff_y_`i'_d_sq_t_XX=total(dof_diff_y_`i'_N_gt_XX) if distance_to_switch_`i'_XX==1

///// If the dynamic effect can be estimated (as there are switchers), we compute the U_Gg variables etc. //
if (N`=increase_XX'_`i'_XX!=0){

///// Computing the U^+_{G,g,l} // Make sure U^+_{G,g,l} is only defined for groups such that: [F_g_XX=t-l&increase_XX==l (==0 for switchers out) or F_g>t] otherwise set it ==.  

// Creating a dummy variable indicating whether l<=T_g_XX-1

gen dummy_U_Gg`i'_XX = (`i'<=T_g_XX-1) 

// Computing (g,t) cell U^+_{G,g,l}

gen U_Gg`i'_temp_XX = dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * N_gt_XX* [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] 

replace U_Gg`i'_temp_XX = U_Gg`i'_temp_XX* diff_y_`i'_XX 

bysort group_XX : gegen U_Gg`i'_XX=total(U_Gg`i'_temp_XX)

replace U_Gg`i'_XX = U_Gg`i'_XX*first_obs_by_gp_XX


// Counting the number of groups for which we can estimate U_Gg`i'_temp_XX - to help compute the "N" displayed by the command //

gen count`i'_core_XX=0
replace count`i'_core_XX=N_gt_XX if (U_Gg`i'_temp_XX!=.&U_Gg`i'_temp_XX!=0|(U_Gg`i'_temp_XX==0&diff_y_`i'_XX==0&(distance_to_switch_`i'_XX!=0|(N`=increase_XX'_t_`i'_g_XX!=0&never_change_d_`i'_XX!=0))))
 
///// Computing the "alternative" U_{G,g,l} which will be used for the computation of the variance only - these are like the above U_{G,g,l}s, except that the outcome differences are demeaned, and there is a DOF adjustment when possible//

gen U_Gg`i'_temp_var_XX =0

// For controls

replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * (diff_y_`i'_N_gt_XX-mean_diff_y_`i'_nd_sq_t_XX)*sqrt(count_diff_y_`i'_nd_sq_t_XX/(count_diff_y_`i'_nd_sq_t_XX-1)) if never_change_d_`i'_XX==1&count_diff_y_`i'_nd_sq_t_XX>1&count_diff_y_`i'_nd_sq_t_XX!=.

replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * diff_y_`i'_N_gt_XX  if never_change_d_`i'_XX==1&count_diff_y_`i'_nd_sq_t_XX==1

// For switchers
replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * (diff_y_`i'_N_gt_XX - mean_diff_y_`i'_d_sq_t_XX) * sqrt(count_diff_y_`i'_d_sq_t_XX/(count_diff_y_`i'_d_sq_t_XX-1)) if distance_to_switch_`i'_XX==1&count_diff_y_`i'_d_sq_t_XX>1&count_diff_y_`i'_d_sq_t_XX!=.

replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * diff_y_`i'_N_gt_XX if distance_to_switch_`i'_XX==1&count_diff_y_`i'_d_sq_t_XX==1

// Summing the U_{G,g,l}s over time periods for each group //
bys group_XX: gegen U_Gg`i'_var_XX=total(U_Gg`i'_temp_var_XX)

}

if "`normalized'"!=""{
	
	capture drop sum_treat_until_`i'_XX
	capture drop delta_D_`i'_cum_temp_XX
	
	bys group_XX : gegen sum_treat_until_`i'_XX = total(treatment_XX - d_sq_XX) if time_XX>=F_g_XX&time_XX<=F_g_XX-1+`i'&S_g_XX==increase_XX
	
gen delta_D_`i'_cum_temp_XX = N_gt_XX/N`=increase_XX'_`i'_XX*[sum_treat_until_`i'_XX* S_g_XX + (1-S_g_XX)*(-sum_treat_until_`i'_XX)] if distance_to_switch_`i'_XX==1 

	sum  delta_D_`i'_cum_temp_XX
	scalar delta_norm_`i'_XX = r(sum) 

}

}


****************************Placebos****************************
if `placebo'!=0{
	if `=l_placebo_u_a_XX'>=1{
forvalue i=1/`=l_placebo_u_a_XX'{

**Needed to compute placebos
capture drop diff_y_pl_`i'_XX
capture drop U_Gg_pl_`i'_temp_XX
capture drop U_Gg_placebo_`i'_XX
capture drop U_Gg_pl_`i'_temp_var_XX
capture drop U_Gg_pl_`i'_var_XX
capture drop mean_diff_y_pl_`i'_d_sq_t_XX
capture drop mean_diff_y_pl_`i'_nd_sq_t_XX
capture drop count_diff_y_pl_`i'_d_sq_t_XX
capture drop count_diff_y_pl_`i'_nd_sq_t_XX
capture drop dist_to_switch_pl_`i'_XX
capture drop never_change_d_pl_`i'_XX
capture drop N`=increase_XX'_t_placebo_`i'_XX
capture drop N`=increase_XX'_t_placebo_`i'_g_XX
capture drop N_gt_control_placebo_`i'_XX
capture drop dummy_U_Gg_pl_`i'_XX 
capture drop never_change_d_pl_`i'_wXX
capture drop dist_to_switch_pl_`i'_wXX
// Create long-differences for the placebos
xtset group_XX time_XX

//The main trick to computing the placebo point estimates is:
// 1. to place the corresponding outcome (y_{F_g-1} - y_{F_g - l - 1})) values in the same row of that (y_{F_g + l -1} - y_{F_g - 1}) of the symmetric DID_l. 
// 2. The other variables, such as N_gt, N0_l or N1_l, remain unchanged, except that we have to check if diff_y_placebo ( = y_{F_g - 2l -2}- y_{F_g - l -1}) exists. 
// 3. If y_{F_g - l -1} does not exist for a specific group, that group is excluded from the calculation, hence, for example, one always has #Switchers for DID_l>= #Switchers for DID_pl.

//Computing the diff_y_placebos
bys group_XX : gen diff_y_pl_`i'_XX = L`=2*`i''.outcome_XX - L`i'.outcome_XX

///// Creating long differences of control variables //
if "`controls'" != ""{

local count_controls=0
// Computing the first differences of the control variables
foreach var of varlist `controls'{

local count_controls=`count_controls'+1

capture drop diff_X_`count_controls'_placebo_`i'_XX

xtset group_XX time_XX

gen diff_X_`count_controls'_placebo_`i'_XX = L`=2*`i''.`var' - L`i'.`var'

foreach l of local levels_d_sq_XX {
	if (scalar(useful_res_`l'_XX)>1){ 

		replace diff_y_pl_`i'_XX = diff_y_pl_`i'_XX - coefs_sq_`l'_XX[`=`count_controls'',1]*diff_X_`count_controls'_placebo_`i'_XX if d_sq_int_XX==`l' 

	}
}

}

}

////// Identifying the control (g,t)s in the estimation of placebo i //

bys group_XX: gen never_change_d_pl_`i'_XX=never_change_d_`i'_XX*(diff_y_pl_`i'_XX!=.) 
gen never_change_d_pl_`i'_wXX=never_change_d_pl_`i'_XX*N_gt_XX

bys time_XX d_sq_XX `trends_nonparam': gegen N_gt_control_placebo_`i'_XX=total(never_change_d_pl_`i'_wXX) 

///// Binary variable indicating whether group g is l periods away from switch & (diff_y_pl_`i'_XX!=.) is well defined. //

gen dist_to_switch_pl_`i'_XX=distance_to_switch_`i'_XX*(diff_y_pl_`i'_XX!=.)
gen dist_to_switch_pl_`i'_wXX= dist_to_switch_pl_`i'_XX*N_gt_XX
///// Computing N^1_{t,l} or N^0_{t,l}. 

bys time_XX: gegen N`=increase_XX'_t_placebo_`i'_XX=total(dist_to_switch_pl_`i'_wXX)

///// Computing N^1_l or N^0_l. //

forvalue t=`=t_min_XX'/`=T_max_XX'{
	sum N`=increase_XX'_t_placebo_`i'_XX if time_XX==`t'
	scalar N`=increase_XX'_placebo_`i'_XX = N`=increase_XX'_placebo_`i'_XX + r(mean)
}


///// Computing N^0_{t,l,g} or N^1_{t,l,g}. 

bys time_XX d_sq_XX `trends_nonparam': gegen N`=increase_XX'_t_placebo_`i'_g_XX=total(dist_to_switch_pl_`i'_wXX)

///// Computing the mean of first or long-differences of outcomes for non-treated and for treated separately - will be useful for the computation of the variance // + //Weighting the means

capture drop diff_y_pl_`i'_N_gt_XX
gen diff_y_pl_`i'_N_gt_XX=N_gt_XX*diff_y_pl_`i'_XX
capture drop dof_diff_y_pl_`i'_N_gt_XX
gen dof_diff_y_pl_`i'_N_gt_XX=(N_gt_XX!=0&diff_y_pl_`i'_XX!=.)

bys time_XX d_sq_XX `trends_nonparam' : gegen mean_diff_y_pl_`i'_nd_sq_t_XX=mean(diff_y_pl_`i'_N_gt_XX) if never_change_d_pl_`i'_XX==1&N`=increase_XX'_t_placebo_`i'_XX>0&N`=increase_XX'_t_placebo_`i'_XX!=.

bys time_XX d_sq_XX `trends_nonparam' : gegen count_diff_y_pl_`i'_nd_sq_t_XX=total(dof_diff_y_pl_`i'_N_gt_XX) if diff_y_pl_`i'_XX!=.&never_change_d_pl_`i'_XX==1&N`=increase_XX'_t_placebo_`i'_XX>0&N`=increase_XX'_t_placebo_`i'_XX!=.

bys time_XX d_sq_XX `trends_nonparam' : gegen mean_diff_y_pl_`i'_d_sq_t_XX=mean(diff_y_pl_`i'_N_gt_XX) if dist_to_switch_pl_`i'_XX==1 

bys time_XX d_sq_XX `trends_nonparam' : gegen count_diff_y_pl_`i'_d_sq_t_XX=total(dof_diff_y_pl_`i'_N_gt_XX) if dist_to_switch_pl_`i'_XX==1

///// If the Placebos effect could be estimated (as there are switchers), we compute the U_Gg. //
gen dummy_U_Gg_pl_`i'_XX = (`i'<=T_g_XX-1)


if (N`=increase_XX'_placebo_`i'_XX!=0){
	gen U_Gg_pl_`i'_temp_XX = dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * N_gt_XX * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX]*diff_y_pl_`i'_XX 

	bysort group_XX : gegen U_Gg_placebo_`i'_XX=total(U_Gg_pl_`i'_temp_XX)

	replace U_Gg_placebo_`i'_XX= U_Gg_placebo_`i'_XX*first_obs_by_gp_XX

	// Counting the number of groups for which we can estimate U_Gg`i'_temp_XX - to help compute the "N" displayed by the command //
	capture drop count`i'_pl_core_XX

	gen count`i'_pl_core_XX=0
	replace count`i'_pl_core_XX= N_gt_XX if (U_Gg_pl_`i'_temp_XX!=.&U_Gg_pl_`i'_temp_XX!=0|(U_Gg_pl_`i'_temp_XX==0&diff_y_pl_`i'_XX==0&(dist_to_switch_pl_`i'_XX!=0|(N`=increase_XX'_t_placebo_`i'_g_XX!=0&never_change_d_pl_`i'_XX!=0))))

	// Initializing the U_{G,g,l}_var at 0 //
	gen U_Gg_pl_`i'_temp_var_XX =0
	
	replace U_Gg_pl_`i'_temp_var_XX =  dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * (diff_y_pl_`i'_N_gt_XX -mean_diff_y_pl_`i'_nd_sq_t_XX)*sqrt(count_diff_y_pl_`i'_nd_sq_t_XX/(count_diff_y_pl_`i'_nd_sq_t_XX-1)) if never_change_d_pl_`i'_XX==1&count_diff_y_pl_`i'_nd_sq_t_XX>1&count_diff_y_pl_`i'_nd_sq_t_XX!=.

	replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * diff_y_pl_`i'_N_gt_XX  if never_change_d_pl_`i'_XX==1&count_diff_y_pl_`i'_nd_sq_t_XX==1
	
	replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * (diff_y_pl_`i'_N_gt_XX - mean_diff_y_pl_`i'_d_sq_t_XX) * sqrt(count_diff_y_pl_`i'_d_sq_t_XX/(count_diff_y_pl_`i'_d_sq_t_XX-1)) if dist_to_switch_pl_`i'_XX==1&count_diff_y_pl_`i'_d_sq_t_XX>1&count_diff_y_pl_`i'_d_sq_t_XX!=.

	replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * diff_y_pl_`i'_N_gt_XX if dist_to_switch_pl_`i'_XX==1&count_diff_y_pl_`i'_d_sq_t_XX==1

	// Summing the U_{G,g,l}s over time periods for each group //
	bys group_XX: gegen U_Gg_pl_`i'_var_XX=total(U_Gg_pl_`i'_temp_var_XX)

	//End of placebos computation 

}

///// Normalization of the placebos//

if "`normalized'"!=""{
	
	capture drop sum_treat_until_`i'_pl_XX
	capture drop delta_D_pl_`i'_cum_temp_XX

	bys group_XX : gegen sum_treat_until_`i'_pl_XX = total(treatment_XX - d_sq_XX) if time_XX>=F_g_XX&time_XX<=F_g_XX-1+`i'&diff_y_pl_`i'_XX!=.&S_g_XX==increase_XX

gen delta_D_pl_`i'_cum_temp_XX = N_gt_XX/N`=increase_XX'_placebo_`i'_XX*[sum_treat_until_`i'_pl_XX* S_g_XX + (1-S_g_XX)*(-sum_treat_until_`i'_pl_XX)] if dist_to_switch_pl_`i'_XX==1

	sum  delta_D_pl_`i'_cum_temp_XX
	scalar delta_norm_pl_`i'_XX = r(sum)
}

} //End of loop on placebos
}
} //End of condition checking if we can compute at least one placebo

*** For the estimation of \hat{\delta} ***

///// Computing the sum of the N1_`i'_XX for the weights w. //
scalar sum_N`=increase_XX'_l_XX = 0
forvalue i=1/`=l_u_a_XX'{
	scalar sum_N`=increase_XX'_l_XX = sum_N`=increase_XX'_l_XX + N`=increase_XX'_`i'_XX
}	


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
	
	 
////// Computing the weights w. //
scalar w_`i'_XX = N`=increase_XX'_`i'_XX / sum_N`=increase_XX'_l_XX
	
	
///// Computing the delta^D_{+,l}s which are necessary for the denominator of \hat{\delta}_+. //

gen delta_D_`i'_temp_XX = N_gt_XX/N`=increase_XX'_`i'_XX*[(treatment_XX-d_sq_XX)* S_g_XX + (1-S_g_XX)*(d_sq_XX-treatment_XX)] if distance_to_switch_`i'_XX==1

replace delta_D_`i'_temp_XX=0 if delta_D_`i'_temp_XX==.
gegen delta_D_`i'_XX = total(delta_D_`i'_temp_XX)
drop delta_D_`i'_temp_XX	

	
///// Computing the numerator of U^+_{G,g}: summing up the U_{G,g,l}s, after weighting them. //
bys group_XX : replace U_Gg_num_XX = U_Gg_num_XX + w_`i'_XX * U_Gg`i'_XX


///// Computing the numerator of the "alternative" U_{G,g}s for the variance : summing up the U_{G,g,l}_vars, after weighting them. //
bys group_XX : replace U_Gg_num_var_XX = U_Gg_num_var_XX + w_`i'_XX * U_Gg`i'_var_XX

///// Computing the denominator of U^+_{G,g}: summing up the delta^D_{+,l}s, after weighting them.  //
bys group_XX : replace U_Gg_den_XX = U_Gg_den_XX + w_`i'_XX * delta_D_`i'_XX

}

}
		
///// Computing the U^+_{G,g}s. //
bys group_XX : gen U_Gg_XX = U_Gg_num_XX/U_Gg_den_XX

///// Computing the U^+_{G,g}_vars. //
bys group_XX : gen U_Gg_var_XX = U_Gg_num_var_XX/U_Gg_den_XX 

	} // end of the quietly condition
	
end


