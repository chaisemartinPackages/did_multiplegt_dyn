*﻿* did_multiplegt_dyn : did_multiplegt, robust to dynamic effects and with asymptotic variances
** This version : September 5th, 2023
** Modified and commented by Clément.


** In this version, we want to add the trends_lin option.
** The trends_lin option has no argument.


********************************************************************************
*                                 PROGRAM 1                                    *
********************************************************************************

capture program drop did_multiplegt_dyn_new

program did_multiplegt_dyn_new, eclass
	version 12.0
	syntax varlist(min=4 max=4 numeric) [if] [in] [, effects(integer 1) placebo(integer 0) switchers(string) controls(varlist numeric) trends_nonparam(varlist numeric) weight(varlist numeric) dont_drop_larger_lower NORMalized cluster(varlist numeric) graphoptions(string) SAVe_results(string) graph_off same_switchers same_switchers_pl effects_equal  drop_if_d_miss_before_first_switch trends_lin ci_level(integer 95) by(varlist numeric max=1) predict_het(string) design(string) date_first_switch(string) normalized_weights(string) CONTinuous(string) save_sample less_conservative_se]
	
// Modif Felix:	
**** Check if gtools is installed, if not present link to insall
qui cap which gtools
if _rc{
	di ""
	di as error "You have not installed the gtools package which is used within the did_multiplegt_dyn command."
	di `"{stata "ssc install gtools": Click here to install gtools}"'
	
	exit
}	

**** Add a Warning that same_switchers_pl only works when same_switchers is specified.
if "`same_switchers_pl'"!=""&"`same_switchers'"==""{
	di ""
	di as error "The same_switchers_pl option only works if same_switchers is specified as well!"
}

/////Doulo: Continous option syntax check:
scalar yes_cont_XX = 0
if("`continuous'"!=""){
	**Create the dummy equal to 1 if treatment is continuous
	scalar yes_cont_XX=1
if strpos("`continuous'", "pol") != 0{	
	if strpos("`continuous'", ",") == 0 {
		di as error ""
		di as error "Syntax error in continuous option: "
		di as error "You must specify the degree of the polynom by: pol, #."
		exit
	}
	else{
	scalar degree_pol = substr(subinstr("`continuous'", " ", "", .),strlen("pol,")+1,.)
	di "degree = "scalar(degree_pol)
	if (scalar(degree_pol)==""){
		di as error ""
		di as error "Syntax error in continuous option:"
		di as error "You must specify the degree of the polynomial function."
	}
	}
}
}

	
****The path of the initial dataset
local dataset_name_XX `c(filename)'

if "`save_sample'" != "" {
	cap drop _did_sample*
}

//>>
preserve


	qui{
		
capture drop outcome_non_diff_XX				
capture drop outcome_XX
capture drop group_XX
capture drop time_XX
capture drop treatment_XX
capture drop d_sq_XX
capture drop d_fg_XX
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
	
****************************************************************************// Modif Felix (from the by option): Add loop and constraint data if by is specified
local levels_byvar_XX 1 // define local so we run the loop once if by is not specified

scalar by_XX  = 1
if "`by'" !=""{
	xtset `2' `3'
	xtsum `by'
	if `r(sd_w)'!=0{
		display as error "The variable specified in the option by"
		display as error "is time-varying, the command will then"
		display as error "ignore the option by."
	 //scalar by_XX = 1
	 //local by ""
	 exit
	}
}
if "`by'" !=""{
// save data to load again in the multiple repetitions of the estimation
tempfile by_data_XX
save "`by_data_XX'.dta", replace

// store local with number of estimations we need to run
levelsof `by', local(levels_byvar_XX)
local `by'_lab_XX: value label `by'
}

// Modif Felix: loop to run estimation multiple times in case of by
} // end of the quietly condition
foreach l of local levels_byvar_XX {
	qui{

if "`by'" !=""{	
*** Check if we have value labels	
ds `by', has(vallabel)
local `by'_has_vallabel_XX=r(varlist)
if "``by'_has_vallabel_XX'"!="."{
	local val_lab_int_XX "`: label ``by'_lab_XX' `l''"
} 
if "``by'_has_vallabel_XX'"=="."{
	local val_lab_int_XX "`l'"
}
		
// Modif Felix: do sample restriction for the by option -> I think it is saver to do it right at the beginning where we also specify the if !!!!
	
	keep if `by' == `l'
}		
****************************************************************************	
	
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

************************************************************
* Modif Felix: Decompose the predict_het string 
if "`predict_het'"!=""{

if strpos("`predict_het'", ",") == 0 {
		di as error ""
		di as error "Syntax error in predict_het option"
		di as error "Comma required"
		
		exit
	}

local het_vars_XX = strtrim(substr("`predict_het'", 1, strpos("`predict_het'", ",") - 1))
local het_nums_XX = strtrim(substr("`predict_het'", strpos("`predict_het'", ",") + 1, .))
	
	local pred_het `het_vars_XX'
	local het_effects `het_nums_XX'	

* Check if predict_het and normalized are specified aimulataneously
if "`normalized'"!=""{
	local pred_het ""
	local het_effects ""
	local predict_het ""
	di as error ""
	di as error "The options normalized and predict_het can not be"
	di as error "specified at the same time, therefore predict_het"
	di as error "will be ignored."
}
}	

* Check if only ime constant variables are specified in predict_het
if "`pred_het'" !=""{
	xtset `2' `3'
	
	local predict_het_bad=""
	local predict_het_good=""
	
	foreach var in `pred_het'{
	xtsum `var'
	if `r(sd_w)'!=0{
		local predict_het_bad="`predict_het_bad' `var'"
	}
	if `r(sd_w)'==0{
		local predict_het_good="`predict_het_good' `var'"
	}
	}
	
	if "`predict_het_bad'"!=""{
		display as error "The variable(s) (`predict_het_bad') specified in the option predict_het"
	display as error "is(are) time-varying, the command will therefore ignore them."
	}
	
	//local pred_het "`predict_het_good'"
}
************************************************************
	
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
				
		collapse (mean) `1' `4' `controls' `trends_nonparam' `cluster' `by' /*`recat_treatment'*/ (count) weight_XX [pw=weight_XX], by(`2' `3') // Modif Felix: added `by'
	}
	
	if "`1'"=="`4'"{
				
		collapse (mean) `1' `controls' `trends_nonparam' `cluster'`by' /*`recat_treatment'*/ (count) weight_XX [pw=weight_XX], by(`2' `3')
	}	
}
 
// Modif Felix: Defining all inputs already here 
*** Also the continous part starts here!
////// Creating Y D G variables //
gen outcome_XX=`1'
gegen group_XX=group(`2')
gegen time_XX=group(`3')
gen treatment_XX=`4'

xtset group_XX time_XX

gen time_d_nonmiss_XX=time_XX if treatment_XX!=.
bys group_XX: gegen min_time_d_nonmiss_XX=min(time_d_nonmiss_XX)
capture drop max_time_d_nonmiss_XX
bys group_XX: gegen max_time_d_nonmiss_XX=max(time_d_nonmiss_XX)
capture drop time_y_nonmiss_XX
capture drop min_time_y_nonmiss_XX
gen time_y_nonmiss_XX=time_XX if outcome_XX!=.
bys group_XX: gegen min_time_y_nonmiss_XX=min(time_y_nonmiss_XX)

capture drop time_d_miss_XX
capture drop min_time_d_miss_aft_ynm_XX
gen time_d_miss_XX=time_XX if treatment_XX==.&time_XX>=min_time_y_nonmiss_XX
bys group_XX: gegen min_time_d_miss_aft_ynm_XX=min(time_d_miss_XX)
drop time_d_nonmiss_XX time_y_nonmiss_XX time_d_miss_XX 

////// Creating baseline treatment: D_{g,1} in paper, redefined to account for imbalanced panels: g's treatment at first period where g's treatment non missing//

gen d_sq_XX_temp=treatment_XX if time_XX==min_time_d_nonmiss_XX
bys group_XX: gegen d_sq_XX=mean(d_sq_XX_temp)
drop d_sq_XX_temp 

////// If the option drop_larger_lower was specified, drop (g,t) cells such that at t, g has experienced both a strictly lower and a strictly higher treatment than its baseline treatment. //
// Modif Felix: Changed it to be default to drop larger/lower cells. Now the option is called dont_drop_larger_lower and it basically reverses the version jow it is used to be. If it is not specified we drop those observations and if it is specified we do not drop the observations

gen diff_from_sq_XX=(treatment_XX-d_sq_XX)

if "`dont_drop_larger_lower'"==""{		
	sort group_XX time_XX
	gen ever_strict_increase_XX=(diff_from_sq_XX>0&treatment_XX!=.) 
	replace ever_strict_increase_XX=1 if ever_strict_increase_XX[_n-1]==1&group_XX==group_XX[_n-1]

	gen ever_strict_decrease_XX=(diff_from_sq_XX<0&treatment_XX!=.) 
	replace ever_strict_decrease_XX=1 if ever_strict_decrease_XX[_n-1]==1&group_XX==group_XX[_n-1]

	drop if ever_strict_increase_XX==1 & ever_strict_decrease_XX==1
	
	drop ever_strict_increase_XX ever_strict_decrease_XX
		
}

////// Sort data //
sort group_XX time_XX

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

// Modif Felix: shift this part here so we can then redefine d_sq_XX for the steps in between where we need the non-continuous version so that we do not run into some errors later
if "`continuous'"!=""{
// 1) generating polynomials of the baseline treatement	
forvalues pol_level = 1/`=degree_pol'{
	scalar pol_level_XX = `pol_level'
	capture drop d_sq_`pol_level'_XX 
	gen d_sq_`pol_level'_XX = d_sq_XX^scalar(pol_level_XX)
}

// 6) redefine d_sq_XX such that it matches the new treatment -> always 0 because with the new treatement definition all periods with no treatement change are 0 and therefore the first period always has to be 0.

capture drop d_sq_XX_orig
gen d_sq_XX_orig=d_sq_XX

replace d_sq_XX=0 // check again if there could be some issues with unbalanced panels, same as above where we replace all _n=0

}


// Modif Doulo: Create a new variable containing the integer levels of d_sq_XX, in the case of non-integer averages of d_sq_XX.
capture drop d_sq_int_XX
gegen d_sq_int_XX = group(d_sq_XX)

///// Dropping the values of the baseline treatment such that no variance in F_g within those values.//
capture drop var_F_g_XX
bys d_sq_XX `trends_nonparam': gegen var_F_g_XX=sd(F_g_XX)
drop if var_F_g_XX==0
drop var_F_g_XX

count
if r(N)==0{
	di as error ""
	di as error "No treatment effect can be estimated."
	di as error "This is because Assumption 1 in"
	di as error "de Chaisemartin & D'Haultfoeuille (2023)"
	di as error "is not satisfied in the data used for"
	di as error "estimation, given the options requested."
	di as error "If this is caused by your baseline treatement"
	di as error "being continuous you can  try using the"
	di as error "option continuous() which allows for a"
	di as error "continous period-one treatement."

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


// Modif Felix: Add this variable so we still can use the non FD version in the het_effects when trends_lin is specified
if "`predict_het_good'"!=""{
gen outcome_non_diff_XX=outcome_XX
}

// Modif Felix : Add the additional restrictions for the trends_lin option and define the new outcome as the firs differences
// When the trends_lin option is specified, drop units for which F_g_XX==2
if "`trends_lin'"!=""{
	drop if F_g_XX==2
	xtset group_XX time_XX
	capture drop FD_outcome_XX
	gen FD_outcome_XX =outcome_XX-L.outcome_XX
	replace outcome_XX = FD_outcome_XX
	drop FD_outcome_XX
	if "`controls'" !=""{
	foreach var of varlist `controls'{
	capture drop FD_`var'_XX
	gen FD_`var'_XX =`var' -L.`var'
	replace `var' = FD_`var'_XX
	drop FD_`var'_XX
	}
	}
	
	drop if time_XX==1
	
	sum time_XX
	scalar t_min_XX=r(min)
	
	///// N.B: if problems with unbalanced panels, maybe the t_min_XX is wrong, should maybe be redefined as T_min_XX, as first time period where a group is observed (only Y, ony D observed ? Or both Y and D should be observed?)
}

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

// When a group is a switching group, but its average post-treatment treatment value is exactly equal to its baseline treatment, we cannnot classify it as a swicher in or a switcher out, but it is not a control either. As such, we drop it from the estimation. Those groups are referred to as no-first-stage-switchers. -> shifted into continuous condition
//drop if avg_post_switch_treat_XX==d_sq_XX&F_g_XX!=T_g_XX+1

//// MODIF: ISSUE WITH CONTINUOUS HERE AS WE REDEFINED d_sq_XX TO 0 SO WE WILL ONLY HAVE SWITCHERS IN (if treatment non negative)

if "`continuous'"==""{
	drop if avg_post_switch_treat_XX==d_sq_XX&F_g_XX!=T_g_XX+1
	gen S_g_XX=(avg_post_switch_treat_XX>d_sq_XX) if F_g_XX!=T_max_XX+1
}
else if "`continuous'"!=""{
	drop if avg_post_switch_treat_XX==d_sq_XX_orig&F_g_XX!=T_g_XX+1
	gen S_g_XX=(avg_post_switch_treat_XX>d_sq_XX_orig) if F_g_XX!=T_max_XX+1
}	
// Here, S_g_XX==. for never-switchers groups only.

// Modif Felix: define the Version where S_g is -1 for switchers out which we need for estimating ß_het
if "`predict_het_good'"!="" | "`continuous'"!=""{
gen S_g_het_XX=S_g_XX
replace S_g_het_XX=-1 if S_g_het_XX==0
if "`predict_het_good'"!=""{
scalar heterogeneity_scalar_XX=1
}
}


// New place continuous
// Modif Felix: Modifications for the continuous option //
if "`continuous'"!=""{
	
// 2) Redefining the treatement -> use F-g and S_g for that
capture drop treatment_XX_temp
gen treatment_XX_temp = (F_g_XX<=time_XX)*S_g_het_XX if S_g_het_XX!=.

capture drop treatment_XX_orig
gen treatment_XX_orig=treatment_XX

replace treatment_XX = treatment_XX_temp 

// 3) generating dummy variables for the time fixed effects
tab time_XX, gen(time_fe_XX_)
local max_time_XX=r(r)
// Check if we may run in error there with r(r) but we should not as we group the time variable

// 4) generate interactions for the controls_time_XX -> exclude one
foreach var of varlist time_fe_XX_2-time_fe_XX_`max_time_XX'{
	forvalues pol_level = 1/`=degree_pol'{
		gen `var'_bt`pol_level'_XX = `var'*d_sq_`pol_level'_XX
		
		// 5) redefine the controls 
		local controls "`controls' `var'_bt`pol_level'_XX"
	}
	capture drop time_fe_XX_`var'
}

// check if we include all controls
di as error ""
di as error "`controls'"
}


//////Modif: Creating treatment at F_g: D_{g,F_g} //
gen d_fg_XX_temp=treatment_XX if time_XX==F_g_XX
bys group_XX: gegen d_fg_XX=mean(d_fg_XX_temp)
replace d_fg_XX=d_sq_XX if d_fg_XX==. & F_g_XX==T_max_XX+1
drop d_fg_XX_temp 


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
	di as error ""
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

// Modif Felix: create a general version of N_gt*delta_X*delta_Y because we need it in multiple places
gen prod_X`count_controls'_diff_y_int_XX = resid_X`count_controls'_time_FE_XX*diff_y_wXX
replace prod_X`count_controls'_diff_y_int_XX = 0 if prod_X`count_controls'_diff_y_int_XX ==.
// call this one later to produce the following variables 

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

//>>restore
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

// Modif Mélitine : Adjust placebo time horizon
// If the trends_lin option was specified, and we dropped the first time period, L_u_XX should be decreased by 1.
if "`trends_lin'"!=""{
	scalar L_placebo_u_XX=L_placebo_u_XX-1
}


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


if "`trends_lin'"!=""{
	scalar L_placebo_a_XX=L_placebo_a_XX-1
}


}
}

///// Error message if Assumption 1 in dCDH (2023) fails and no effect can be computed //

if ("`switchers'"=="in"&(L_u_XX==.|L_u_XX==0))|("`switchers'"=="out"&(L_a_XX==.|L_a_XX==0))|("`switchers'"==""&(L_u_XX==.|L_u_XX==0)&(L_a_XX==.|L_a_XX==0)){
		di as error ""
		di as error "No treatment effect can be estimated."
		di as error "This is because Assumption 1 in"
		di as error "de Chaisemartin & D'Haultfoeuille (2023)"
		di as error "is not satisfied in the data used for"
		di as error "estimation, given the options requested."
		di as error "If this is caused by your baseline treatement"
		di as error "being continuous you can  try using the"
		di as error "option continuous() which allows for a"
		di as error "continous period-one treatement."
		
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
	di as error ""
	di as error "The number of effects requested is too large."
	di as error "The number of effects which can be estimated is at most " l_XX "."
	di as error "The command will therefore try to estimate " l_XX " effect(s)."
}

if `placebo'!=0{
	if l_placebo_XX<`placebo'&`effects'>=`placebo'{
		di as error ""
		di as error "The number of placebos which can be estimated is at most " l_placebo_XX "."
		di as error "The command will therefore try to estimate " l_placebo_XX " placebo(s)."
	}

	if `effects'<`placebo'{
		di as error ""
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
	
	// Modif Felix -> In the Version from 03.10.23 on I use this scalar which is initialized in the same loop where we call the core program
	scalar N0_`i'_XX_new=0
	scalar N1_`i'_XX_new=0
	
	
	// Felix : do we have to create the covariance default variables -> No
	
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
	
	// Modif Felix -> In the Version from 03.10.23 on I use this scalar which is initialized in the same loop where we call the core program
	scalar N0_placebo_`i'_XX_new=0
	scalar N1_placebo_`i'_XX_new=0
	
	// For Felix : should we initialize covariance variables here -> No, only need them in the core

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
	
if "`trends_lin'"==""{
	did_multiplegt_dyn_core_new outcome_XX group_XX time_XX treatment_XX, effects(`=l_XX') placebo(`=l_placebo_XX') switchers_core(in) controls(`controls') trends_nonparam(`trends_nonparam') `normalized' `same_switchers' `same_switchers_pl' `effects_equal' continuous(`continuous') `less_conservative_se'
}	
	
// Store the results
forvalue i=1/`=l_XX'{

// Modif Mélitine : Call the core program inside of the loop if trends_lin is specified	(both for switchers in and out)
// Note that if the option trends_lin was specified, same_switchers must also be specified - HOW TO SPECIFY BOTH TRENDS_LIN and SAME_SWITCHERS?
if "`trends_lin'"!=""{
	did_multiplegt_dyn_core_new outcome_XX group_XX time_XX treatment_XX, effects(`i') switchers_core(in) controls(`controls') trends_nonparam(`trends_nonparam') `normalized' same_switchers `effects_equal' trends_lin continuous(`continuous') `less_conservative_se'
	
	/////////// NB: Much easier to "aggregate" the U_Gg_`i'_XX across horizons in the core program while they are still swithcers in / switchers out specific, because then the calling is easy. If we wait after aggregation, then the calling of the core program will be harder.
}

/////////// NB: in the case of unbalanced panels, it can happen that the U_Gg`i'_XX are not computed by program 2 (for example when y is missing). Consequently, for the command not to display an error message and continue running, we need to verify the variable is created, which is conditional on  N1_`i'_XX!=0.
if N1_`i'_XX!=0{
		replace U_Gg`i'_plus_XX = U_Gg`i'_XX
		replace count`i'_plus_XX= count`i'_core_XX
		replace U_Gg_var_`i'_in_XX=U_Gg`i'_var_XX
		
		// Modif Felix 
		scalar N1_`i'_XX_new=N1_`i'_XX
		
		// For Felix : should you store switchers-in covariances here ?
		
		if "`normalized'"!=""{
			scalar delta_D_`i'_in_XX = delta_norm_`i'_XX
		}
	}

}




if l_placebo_XX!=0{

	forvalue i=1/`=l_placebo_XX'{
		
		if "`trends_lin'"!=""{
	did_multiplegt_dyn_core_new outcome_XX group_XX time_XX treatment_XX, effects(`i') placebo(`i') switchers_core(in) controls(`controls') trends_nonparam(`trends_nonparam') `normalized' same_switchers same_switchers_pl `effects_equal' trends_lin continuous(`continuous') `less_conservative_se'
}
		
		
		if N1_placebo_`i'_XX!=0{
				replace U_Gg_pl_`i'_plus_XX  = U_Gg_placebo_`i'_XX
				replace count`i'_pl_plus_XX= count`i'_pl_core_XX
				replace U_Gg_var_pl_`i'_in_XX=U_Gg_pl_`i'_var_XX
				
				// Modif Felix 
				scalar N1_placebo_`i'_XX_new=N1_placebo_`i'_XX
				
				// should we store the covariance variables here ?
				
				if "`normalized'"!=""{
					scalar delta_D_pl_`i'_in_XX = delta_norm_pl_`i'_XX
				}

		}
		
	}	
}

	if "`trends_lin'"==""{
	if sum_N1_l_XX!=0{
	replace U_Gg_plus_XX = U_Gg_XX
	scalar U_Gg_den_plus_XX=U_Gg_den_XX
	replace U_Gg_var_plus_XX = U_Gg_var_XX
	}
	}

}
	
} // end of the loop for switchers in

// For switchers out
if ("`switchers'"==""|"`switchers'"=="out"){
if L_a_XX!=.&L_a_XX!=0{
	
if "`trends_lin'"==""{	
did_multiplegt_dyn_core_new outcome_XX group_XX time_XX treatment_XX, effects(`=l_XX') placebo(`=l_placebo_XX') switchers_core(out) controls(`controls')  trends_nonparam(`trends_nonparam') `normalized' `same_switchers' `same_switchers_pl' `effects_equal' continuous(`continuous') `less_conservative_se'
}
	
// Store the results
forvalue i=1/`=l_XX'{
	
if "`trends_lin'"!=""{	
	did_multiplegt_dyn_core_new outcome_XX group_XX time_XX treatment_XX, effects(`i') switchers_core(out) controls(`controls') trends_nonparam(`trends_nonparam') `normalized' same_switchers `effects_equal' trends_lin continuous(`continuous') `less_conservative_se'
}
	
if N0_`i'_XX!=0{
		replace U_Gg`i'_minus_XX = - U_Gg`i'_XX
		replace count`i'_minus_XX= count`i'_core_XX
		replace U_Gg_var_`i'_out_XX=U_Gg`i'_var_XX
		
		// For Felix : should you store switchers-out covariances here ?

		// Modif Felix 
		scalar N0_`i'_XX_new=N0_`i'_XX
		
		if "`normalized'"!=""{
			scalar delta_D_`i'_out_XX = delta_norm_`i'_XX
		}
		
	}
	
}


if l_placebo_XX!=0{

	forvalue i=1/`=l_placebo_XX'{
		
if "`trends_lin'"!=""{	
	did_multiplegt_dyn_core_new outcome_XX group_XX time_XX treatment_XX, effects(`i') placebo(`i') switchers_core(out) controls(`controls') trends_nonparam(`trends_nonparam') `normalized' same_switchers same_switchers_pl `effects_equal' trends_lin continuous(`continuous') `less_conservative_se'
}
		
	if N0_placebo_`i'_XX!=0{
			replace U_Gg_pl_`i'_minus_XX  = -U_Gg_placebo_`i'_XX
			replace count`i'_pl_minus_XX= count`i'_pl_core_XX
			replace U_Gg_var_pl_`i'_out_XX=U_Gg_pl_`i'_var_XX
			
			// Modif Felix 
			scalar N0_placebo_`i'_XX_new=N0_placebo_`i'_XX
			
			if "`normalized'"!=""{
					scalar delta_D_pl_`i'_out_XX = delta_norm_pl_`i'_XX
	}

	}
	}		
}

	if "`trends_lin'"==""{
	if sum_N0_l_XX!=0{
	replace U_Gg_minus_XX = - U_Gg_XX
	scalar U_Gg_den_minus_XX=U_Gg_den_XX
	replace U_Gg_var_minus_XX = U_Gg_var_XX
	}
	}

}
	
} // end of the loop for switchers out 

	

///// Aggregating the results for switchers in and out. //

****************************DID_l*********************************

//Creation of the matrix which stores all the estimators (DID_l, DID_pl, delta, etc.), their sd and the CIs

matrix mat_res_XX = J(l_XX+l_placebo_XX+1,7,.) 

// Computing first the "global" U_Ggl s
forvalue i=1/`=l_XX'{
	
gen U_Gg`i'_global_XX = (N1_`i'_XX_new/(N1_`i'_XX_new+N0_`i'_XX_new))*U_Gg`i'_plus_XX +(N0_`i'_XX_new/(N1_`i'_XX_new+N0_`i'_XX_new))*U_Gg`i'_minus_XX
replace U_Gg`i'_global_XX=. if first_obs_by_gp_XX==0

gen count`i'_global_XX = max(count`i'_plus_XX, count`i'_minus_XX)

	if "`normalized'"!=""{
	scalar delta_D_`i'_global_XX = (N1_`i'_XX_new/(N1_`i'_XX_new+N0_`i'_XX_new))*delta_D_`i'_in_XX+(N0_`i'_XX_new/(N1_`i'_XX_new+N0_`i'_XX_new))*delta_D_`i'_out_XX
	}

// Storing the results and the number of switchers into the matrix mat_res_XX
local rownames "`rownames' Effect_`i'"

* Number of switchers
scalar N_switchers_effect_`i'_XX=N1_`i'_XX_new+N0_`i'_XX_new // Modif Felix
matrix mat_res_XX[`i',6]=N_switchers_effect_`i'_XX
matrix mat_res_XX[`i',7]=`i'
ereturn scalar N_switchers_effect_`i' = N_switchers_effect_`i'_XX
* Number of observations used in the estimation
gegen N_effect_`i'_XX = total(count`i'_global_XX) 
scalar N_effect_`i'_XX = N_effect_`i'_XX
ereturn scalar N_effect_`i' = N_effect_`i'_XX
matrix mat_res_XX[`i',5]=N_effect_`i'_XX

if N_switchers_effect_`i'_XX==0|N_effect_`i'_XX==0{
	di as error ""
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

if ("`switchers'"==""&N1_`i'_XX_new==0&N0_`i'_XX_new==0)|("`switchers'"=="out"&N0_`i'_XX_new==0)|("`switchers'"=="in"&N1_`i'_XX_new==0){
	scalar DID_`i'_XX=.
}

ereturn scalar Effect_`i' = scalar(DID_`i'_XX)
matrix mat_res_XX[`i',1]= scalar(DID_`i'_XX) 

}

****************************Average Total Effect*********************************


if "`trends_lin'"==""{ // Problem with N placebo is somwhere here!

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

}


// Modif Felix: add time_to_treat so graph will be correct
if "`trends_lin'"!=""{
	local rownames "`rownames' Av_tot_eff" 
	matrix mat_res_XX[l_XX+1,7]=0
}	
****************************Placebos****************************
// Computing first the "global" U_Ggl s for the placebos

	
if l_placebo_XX!=0{

forvalue i=1/`=l_placebo_XX'{	

gen U_Gg_pl_`i'_global_XX = (N1_placebo_`i'_XX_new/(N1_placebo_`i'_XX_new+N0_placebo_`i'_XX_new))*U_Gg_pl_`i'_plus_XX  +(N0_placebo_`i'_XX_new/(N1_placebo_`i'_XX_new+N0_placebo_`i'_XX_new))*U_Gg_pl_`i'_minus_XX 
replace U_Gg_pl_`i'_global_XX=. if first_obs_by_gp_XX==0

gen count`i'_pl_global_XX=max(count`i'_pl_plus_XX, count`i'_pl_minus_XX)

if "`normalized'"!=""{
	scalar delta_D_pl_`i'_global_XX = (N1_placebo_`i'_XX_new/(N1_placebo_`i'_XX_new+N0_placebo_`i'_XX_new))*delta_D_pl_`i'_in_XX+(N0_placebo_`i'_XX_new/(N1_placebo_`i'_XX_new+N0_placebo_`i'_XX_new))*delta_D_pl_`i'_out_XX
		}

// Summing them to obtain the DID_pl

gegen DID_placebo_`i'_XX = total(U_Gg_pl_`i'_global_XX) 
replace  DID_placebo_`i'_XX= DID_placebo_`i'_XX/G_XX
scalar DID_placebo_`i'_XX = DID_placebo_`i'_XX


if "`normalized'"!=""{
	scalar DID_placebo_`i'_XX = scalar(DID_placebo_`i'_XX)/scalar(delta_D_pl_`i'_global_XX)
}

if ("`switchers'"==""&N1_placebo_`i'_XX_new==0&N0_placebo_`i'_XX_new==0)|("`switchers'"=="out"&N0_placebo_`i'_XX_new==0)|("`switchers'"=="in"&N1_placebo_`i'_XX_new==0){
	scalar DID_placebo_`i'_XX=.
}

ereturn scalar Placebo_`i' = scalar(DID_placebo_`i'_XX)

// Storing the results and the number of switchers into the matrix mat_res_XX
* Storing the results 
matrix mat_res_XX[`=l_XX'+ 1 + `i',1]=scalar(DID_placebo_`i'_XX)

local rownames "`rownames' Placebo_`i'"
* Number of switchers
scalar N_switchers_placebo_`i'_XX=N1_placebo_`i'_XX_new+N0_placebo_`i'_XX_new
matrix mat_res_XX[`=l_XX'+ 1 + `i',6]=N_switchers_placebo_`i'_XX
matrix mat_res_XX[`=l_XX'+ 1 + `i',7]= `=-`i''
ereturn scalar N_switchers_placebo_`i' = N_switchers_placebo_`i'_XX
* Number of observations used in the estimation
gegen N_placebo_`i'_XX = total(count`i'_pl_global_XX)
scalar N_placebo_`i'_XX = N_placebo_`i'_XX
ereturn scalar N_placebo_`i' = N_placebo_`i'_XX
matrix mat_res_XX[`=l_XX' + 1 + `i',5]=N_placebo_`i'_XX

if N_switchers_placebo_`i'_XX==0|N_placebo_`i'_XX==0{
	di as error ""
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
if ("`switchers'"==""&(N1_`i'_XX_new!=0|N0_`i'_XX_new!=0))|("`switchers'"=="out"&N0_`i'_XX_new!=0)|("`switchers'"=="in"&N1_`i'_XX_new!=0){

// Estimating \hat{\sigma}^2_l

gen U_Gg_var_glob_`i'_XX = U_Gg_var_`i'_in_XX * (scalar(N1_`i'_XX_new)/(scalar(N1_`i'_XX_new)+scalar(N0_`i'_XX_new))) + U_Gg_var_`i'_out_XX* (scalar(N0_`i'_XX_new)/(scalar(N1_`i'_XX_new)+scalar(N0_`i'_XX_new)))

// For Felix : equivalent aggregation for the covariances ? Same weighting scheme as for the U_Gg_var_glob.

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

//ci_level
scalar ci_level = `ci_level'/100
scalar z_level = invnormal(scalar(ci_level) + (1-scalar(ci_level))/2)
	
// Lower bound of the 95% confidence interval
scalar LB_CI_`i'_XX = scalar(DID_`i'_XX) - scalar(z_level)*scalar(se_`i'_XX)
matrix mat_res_XX[`i',3]= scalar(LB_CI_`i'_XX)

// Upper bound of the 95% confidence interval
scalar UB_CI_`i'_XX = scalar(DID_`i'_XX) + scalar(z_level)*scalar(se_`i'_XX)
matrix mat_res_XX[`i',4]=scalar(UB_CI_`i'_XX)


}
}


****************************Placebos****************************
///// Estimating \hat{\sigma}^2_pl and the confidence intervals //


if l_placebo_XX!=0{

//Take into account the case where same_switchers leads to N0 or N1 = 0 for placebos
forvalue i=1/`=l_placebo_XX'{
if ("`switchers'"==""&(N1_placebo_`i'_XX_new!=0|N0_placebo_`i'_XX_new!=0))|("`switchers'"=="out"&N0_placebo_`i'_XX_new!=0)|("`switchers'"=="in"&N1_placebo_`i'_XX_new!=0){ 

// Estimating \hat{\sigma}^2_l

gen U_Gg_var_glob_pl_`i'_XX = U_Gg_var_pl_`i'_in_XX * (N1_placebo_`i'_XX_new/(N1_placebo_`i'_XX_new+N0_placebo_`i'_XX_new)) + U_Gg_var_pl_`i'_out_XX* (N0_placebo_`i'_XX_new/(N1_placebo_`i'_XX_new+N0_placebo_`i'_XX_new))

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
scalar LB_CI_placebo_`i'_XX = scalar(DID_placebo_`i'_XX) - scalar(z_level)*scalar(se_placebo_`i'_XX)
matrix mat_res_XX[`=l_XX'+ 1 + `i', 3]= scalar(LB_CI_placebo_`i'_XX)

// Upper bound of the 95% confidence interval
scalar UB_CI_placebo_`i'_XX = scalar(DID_placebo_`i'_XX) + scalar(z_level)*scalar(se_placebo_`i'_XX)
matrix mat_res_XX[`=l_XX'+ 1 + `i',4]=scalar(UB_CI_placebo_`i'_XX)


}
}
}


****************************Average Effect****************************
///// Estimating \hat{\sigma}^2 and the confidence interval for the average effect //

if "`trends_lin'"==""{

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
scalar LB_CI_XX = delta_XX - scalar(z_level)*se_XX
matrix mat_res_XX[l_XX+1,3]= LB_CI_XX

// Lower bound of the 95% confidence interval
scalar UB_CI_XX = delta_XX + scalar(z_level)*se_XX
matrix mat_res_XX[l_XX+1,4]= UB_CI_XX

}

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
		if ("`switchers'"==""&(N1_placebo_`i'_XX_new!=0|N0_placebo_`i'_XX_new!=0))|("`switchers'"=="out"&N0_placebo_`i'_XX_new!=0)|("`switchers'"=="in"&N1_placebo_`i'_XX_new!=0){
			
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
			
			// For Felix : here, you can find the covariances for the placebos. To be put in thr core program where indicated.
			
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
		di as error ""
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
		if ("`switchers'"==""&(N1_`i'_XX_new!=0|N0_`i'_XX_new!=0))|("`switchers'"=="out"&N0_`i'_XX_new!=0)|("`switchers'"=="in"&N1_`i'_XX_new!=0){
			
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
				if ("`switchers'"==""&(N1_`i'_XX_new!=0|N0_`i'_XX_new!=0))|("`switchers'"=="out"&N0_`i'_XX_new!=0)|("`switchers'"=="in"&N1_`i'_XX_new!=0){

		
		// Storing the already computed variances.
		matrix didmgt_Var_Effects[`i',`i']= scalar(se_`i'_XX)^2
		
		
		if `i'<`=l_XX'{
			matrix didmgt_identity[`i',`i'] =1
		}
		
		
		// Computing and storing the covariances.
		if `i'<`=l_XX'{
		forvalue j=`=`i'+1'/`=l_XX'{
			
			// For Felix : here are the covariances for the effects. To be displaced.
			
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
		di as error ""
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
if "`by'" !=""{	
di as input _skip(35) "By: `by' = `val_lab_int_XX'"
}	
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

if "`trends_lin'"==""{

matrix mat_res_avg_XX=mat_res_XX[l_XX+1, 1..6]
matrix mat_res_avg_XX=(mat_res_avg_XX, .z )
matrix colnames mat_res_avg_XX= "Estimate" "SE" "LB CI" "UB CI" "N" "Switch" "x Periods"
display _newline
di as input "{hline 80}"
di as input _skip(4) "Estimation of treatment effects: Average total effect per treatment unit"
if "`by'" !=""{	
di as input _skip(35) "By: `by' = `val_lab_int_XX'"
}	
di as input "{hline 80}"
noisily matlist mat_res_avg_XX, nodotz
di as input "{hline 80}"

}

if "`trends_lin'"!=""{
display _newline
di as input "{hline 80}"
di as input _skip(4) "When the option {it:trends_lin} is specified no average effects are reported"
di as input "{hline 80}"
}	

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
if "`by'" !=""{	
di as input _skip(35) "By: `by' = `val_lab_int_XX'"
}	
di as input "{hline 80}"
matlist mat_res_XX[l_XX+2...,1..6]
di as input "{hline 80}"
if l_placebo_XX>1&all_Ns_pl_not_zero==l_placebo_XX{
di as text "{it:Test of joint nullity of the placebos : p-value =} " scalar(p_jointplacebo)

}
}

/*
di " "
di "----------------------------------------------"
di "        Predicting effect heterogeneity       "
di "----------------------------------------------"
*/


if "`predict_het_good'"!=""{ 
qui{
	
sort group_XX time_XX	
	
//  Modif Felix: Now generate the correct long differences (that will have to be inside the loop over the l effects (in the code `i'))
* Isolate the two Y values in Yg,Fg−1+ℓ − Yg,Fg−1, they are fixed within group and ℓ

*** Define number of effects we want to calculate
if "`het_effects'"==""{
local all_effects_XX ""	
forvalues i=1/`=l_XX'{
local all_effects_XX "`all_effects_XX' `i'"
}
}

if "`het_effects'"!=""{ // allow to only show some effects 
local all_effects_XX "`het_effects'"
local test_eff_XX : subinstr local het_effects " " ",", all 
local max_test_eff_XX = max(`test_eff_XX')
	if `max_test_eff_XX'>`=l_XX'{
		di as error ""
		di as error "You specified some numbers in predict_het that exceed the number of effects possible to estimate!"
		di as error "Please specify only numbers that are smaller or equal to the number you request in effects()."
		exit
	}
}

foreach i in `all_effects_XX'{
	
capture drop Yg_Fg_min_1_XX
capture drop Yg_Fg_`i'_XX
capture drop prod_het_`i'_XX
capture drop interact_var_XX
capture drop group_interact_var_XX
capture drop d_sq_group_XX	

* Yg,Fg−1
gen Yg_Fg_min1_XX_temp=outcome_non_diff_XX if time_XX==F_g_XX-1
bys group_XX: gegen Yg_Fg_min_1_XX=mean(Yg_Fg_min1_XX_temp)
capture drop Yg_Fg_min1_XX_temp

* Yg,Fg−1+ℓ
gen Yg_Fg_`i'_XX_temp=outcome_non_diff_XX if time_XX==F_g_XX-1+`i'
bys group_XX: gegen Yg_Fg_`i'_XX=mean(Yg_Fg_`i'_XX_temp)
capture drop Yg_Fg_`i'_XX_temp

* Now we can generate Sg*(Yg,Fg−1+ℓ − Yg,Fg−1)
gen prod_het_`i'_XX=S_g_het_XX*(Yg_Fg_`i'_XX - Yg_Fg_min_1_XX)

* keep one observation by group to not artificially increase sample
bys group_XX: replace prod_het_`i'_XX = . if _n != 1

* F_g_XX#d_sq_XX#S_g_het_XX
gegen d_sq_group_XX=group(d_sq_XX)

*** Regression of interest 
reg prod_het_`i'_XX `predict_het_good' F_g_XX#d_sq_group_XX#S_g_XX if F_g_XX-1+`i'<=T_g_XX, level(`ci_level') vce(robust)

forvalues j=1/`:word count `predict_het_good''{
scalar beta_het_`j'_eff`i'_hat_XX=_b["`: word `j' of `predict_het_good''"]
scalar se_het_`j'_eff`i'_hat_XX=_se["`: word `j' of `predict_het_good''"]
scalar lb_het_`j'_eff`i'_hat_XX=_r_lb["`: word `j' of `predict_het_good''"]
scalar ub_het_`j'_eff`i'_hat_XX=_r_ub["`: word `j' of `predict_het_good''"]
scalar t_het_`j'_eff`i'_hat_XX=scalar(beta_het_`j'_eff`i'_hat_XX)/scalar(se_het_`j'_eff`i'_hat_XX)
}
scalar N_het_`i'_hat_XX=e(N)

// Add F-test on joint significance
test `predict_het_good'
scalar p_het_`i'_hat_XX=r(p)
scalar F_het_`i'_hat_XX=r(F)
}
}

**** Output Part ****	

foreach i in `all_effects_XX'{
matrix effect_het_`i'_XX=J(`:word count `predict_het_good'',6,.)

local effect_het_rownames_`i'_XX ""

forvalue j=1/`:word count `predict_het_good''{
matrix effect_het_`i'_XX[`j',1]=scalar(beta_het_`j'_eff`i'_hat_XX)
matrix effect_het_`i'_XX[`j',2]=scalar(se_het_`j'_eff`i'_hat_XX)
matrix effect_het_`i'_XX[`j',3]=scalar(t_het_`j'_eff`i'_hat_XX)
matrix effect_het_`i'_XX[`j',4]=scalar(lb_het_`j'_eff`i'_hat_XX)
matrix effect_het_`i'_XX[`j',5]=scalar(ub_het_`j'_eff`i'_hat_XX)
matrix effect_het_`i'_XX[`j',6]=scalar(N_het_`i'_hat_XX)

*local effect_het_rownames_`i'_XX "`effect_het_rownames_`i'_XX' "Effect_`i'_`: word `j' of `predict_het_good''""
local effect_het_rownames_`i'_XX "`effect_het_rownames_`i'_XX' "`: word `j' of `predict_het_good''""
}


matrix colnames effect_het_`i'_XX= "Estimate" "SE" "t" "LB CI" "UB CI" "N"
matrix rownames effect_het_`i'_XX=`effect_het_rownames_`i'_XX'

}	



display _newline
di as input "{hline 80}"
di as input _skip(25) "Predicting effect heterogeneity"
if "`by'" !=""{	
di as input _skip(35) "By: `by' = `val_lab_int_XX'"
}	
di as input "{hline 80}"

foreach i in `all_effects_XX'{
display _newline
di as input "{hline 80}"	
//di as input _skip(36) "Effect ℓ=`i'"	
di as input _skip(37) "Effect_`i'"
di as input "{hline 80}"
matlist effect_het_`i'_XX
di as input "{hline 80}"
di as text "{it:Test of joint nullity of the estimates : p-value =} " p_het_`i'_hat_XX
}
}

** Normalized weights **
if "`normalized_weights'" != "" {
	qui {
		if "`normalized'" == "" {
			noi di as err "normalized option required to compute normalized_weights"
			exit
		}

	local lag_d = strtrim("`normalized_weights'")

	if !inlist("`lag_d'", "by_k", "by_calendar") {
		noi di as err "First argument of normalized_weights incorrectly specified"
		noi di as err "normalized_weights requires by_k or by_calendar as arguments"
		exit
	}	

	mat define weight_mat = J(`=l_XX + 1', `=l_XX', .)
	local cols ""
	forv i = 1/`=l_XX'{	
		local m = "ℓ=`i'"
		local cols `cols' `m'
		forv k = 0/`=`i' - 1' {
			if "`lag_d'" == "by_k" {
				local row = `k' + 1
			}
			else {
				local row = `i' - `k' // Visualization by Fg - 1 + (l - k)
			}
			// Modif Felix: adjust for continuous option
			if "`continuous'"==""{
			gen delta_`i'_`k' = abs(treatment_XX - d_sq_XX) if time_XX == F_g_XX - 1 + `i' - `k' & F_g_XX - 1 + `i' <= T_g_XX
			}
			else if "`continuous'"!=""{
			gen delta_`i'_`k' = abs(treatment_XX_orig - d_sq_XX_orig) if time_XX == F_g_XX - 1 + `i' - `k' & F_g_XX - 1 + `i' <= T_g_XX
			}
			sum delta_`i'_`k'
			mat weight_mat[`row', `i'] = (`r(sum)' / scalar(delta_D_`i'_global_XX)) / scalar(N_switchers_effect_`i'_XX)
		}		
	}
	

	local rows ""
	if "`lag_d'" == "by_k" {
		forv j = 1/`=l_XX' {
				local r = "k=`=`j'-1'"
				local rows `rows' `r'
		}	
	}
	else {
		forv j = 1/`=l_XX' {
			if `j' == 1 {
				local r = "D_Fg"	
			}
			else {
				local r = "D_Fg+`=`j'-1'"
			}
			local rows `rows' `r'
		}
	}
	
	
	local rows `rows' "Total"
	
	matrix define mat_temp = J(`=l_XX', `=l_XX', 0)	
	forv i = 1/`=l_XX' {
		forv j = 1/`=l_XX' {
			if weight_mat[`i', `j'] != . {
				mat mat_temp[`i', `j'] = weight_mat[`i', `j']
			}
		}
	}	
	matrix define mat_total = mat_temp' * J(`=l_XX', 1, 1)
	forv i = 1/`=l_XX' {
		mat weight_mat[`=l_XX+1', `i'] = mat_total[`i', 1]
	}
	
	mat rown weight_mat = `rows'
	mat coln weight_mat = `cols'
	local by_opt_lag = subinstr("`lag_d'", "by_", "", .)
	}
	noi display _newline
	di as input "{hline 70}"
	di as input _skip(15) "Weights on treatment lags - by `by_opt_lag'"
	if "`by'" !=""{	
	di as input _skip(25) "By: `by' = `val_lab_int_XX'"
	}	
	di as input "{hline 70}"
	noi matlist weight_mat, format(%9.4fc) lines(rowtotal)
	noi display _newline		
}

*** Modif Diego: Design detection option

// Modif Diego (for design): Saving data for further options //
qui tempfile data_predesign
qui save `data_predesign', replace

if "`save_sample'" != "" {
	qui {
		keep if !missing(`2') & !missing(`3') // Drop _fillin output //
		keep `2' `3' S_g_XX
		if "`by'"==""{
		rename S_g_XX _did_sample
		replace _did_sample = -1 if _did_sample == 0
		replace _did_sample = 0 if _did_sample == .
		cap label drop switch_lab_XX 
		label define switch_lab_XX 0 "Never-switcher" 1 "Switcher-in" -1 "Switcher-out"
		label values _did_sample switch_lab_XX
		tempfile did_sample
		save `did_sample', replace
		}
		else if "`by'"!=""{
		rename S_g_XX _did_sample_`l'
		replace _did_sample_`l' = -1 if _did_sample_`l' == 0
		replace _did_sample_`l' = 0 if _did_sample_`l' == .
		cap label drop switch_lab_XX 
		label define switch_lab_XX 0 "Never-switcher" 1 "Switcher-in" -1 "Switcher-out"
		label values _did_sample_`l' switch_lab_XX
		tempfile did_sample_`l'
		save `did_sample_`l'', replace
		}
	}
}

// Syntax 										     : design(% over total, dataset path)
// Syntax for console only (with max treatment paths): design(, console)
if "`design'" != "" {	
	qui {
	
	use `data_predesign', clear
	
	if strpos("`design'", ",") == 0 {
		di as error ""
		di as error "Syntax error in design option"
		di as error "Comma required"
		
		exit
	}
	local des_p = strtrim(substr("`design'", 1, strpos("`design'", ",") - 1))
	local des_path = strtrim(substr("`design'", strpos("`design'", ",") + 1, .))
	local des_n = l_XX 
	
	if missing("`des_p'") {
		local des_p = 1
	}
	// New
	local des_per = `des_p' * 100
	local des_per : di %9.2fc `des_per'
	
	gen F_g_plus_n_XX = F_g_XX + `des_n' - 1
	keep if time_XX >= F_g_XX - 1 & time_XX <= F_g_plus_n_XX
	sort group_XX time_XX
	bys group_XX : gen time_l_XX = _n

	// Aggregate weights by group //
	keep group_XX time_l_XX weight_XX treatment_XX F_g_XX
	{
		if "`weight'" != "" {
			gegen g_weight_XX = sum(weight_XX), by(group_XX)
		}
		else {
			gen g_weight_XX = 1
		}
	}
	drop weight_XX
	greshape wide treatment_XX, i(group_XX g_weight_XX F_g_XX) j(time_l_XX)
	
	// Drop missing treatments //
	foreach v of varlist treatment_XX* {
		drop if `v' == .
	}
	
	* Remove non switchers *
	count if F_g_XX > T_max_XX
	local non_switchers = r(N)
	drop if F_g_XX > T_max_XX
	
	gen N_XX = 1
	sum g_weight_XX
	gen N_w_XX =  (g_weight_XX * N_XX) / r(sum)
	drop group_XX g_weight_XX
	gcollapse (sum) N_XX N_w_XX, by(treatment_XX*)
	order N* treat*	
	sum N_XX
	local tot_switch = r(sum)
	
	// Keep the observations amounting to p% of the detected treatment paths //
	gen neg_N_XX = -N_XX
	sort neg_N_XX treatment_XX*
	gen cum_sum_XX = sum(N_w_XX)
	gen in_table_XX = (cum_sum_XX <= `des_p')
	sort in_table_XX cum_sum_XX
	bys in_table_XX: gen id_XX = _n
	keep if (in_table_XX == 1) | (in_table_XX == 0 & id_XX == 1) // Keep the first observation exceeding the p%	
	{
		if `des_p' < 1 {
			sum cum_sum_XX if in_table == 0
			local last_p : di %9.2fc 100 * r(min)
		}
		else {
			local last_p : di %9.2fc 100
		}
	}
	sort neg_N_XX treatment_XX*
	
	drop neg_N_XX in_table_XX id_XX cum_sum_XX
	
	mat define desmat = J(`=_N', `= 2 + 1 + `des_n'', .)
	local colnames
	label var N_XX "#Groups"
	label var N_w_XX "%Groups"
	replace N_w_XX=100*N_w_XX // Modif Felix: change to %
	local j = 1
	foreach v of varlist _all {
		local m : var label `v'
		if strpos("`v'", "treatment_XX") > 0 {
			local k = substr("`v'", length("treatment_XX") + 1, .)
			local m = "ℓ=`=`k'-1'"
		}
		local colnames `colnames' `m'	
		forv i = 1/`=_N' {
			matrix desmat[`i', `j'] = `v'[`i']
		}
		local j = `j' + 1
	}
	local rownames
	forv i = 1/`=_N' {
		local rownames `rownames' "TreatPath`i'"
	}
	foreach v in col row {
		mat `v'n desmat = ``v'names'
	}

	if "`des_path'" == "console" {		
		noi display _newline
		di as input "{hline 80}"
		di as input _skip(10) "Detection of treatment paths - `des_n' periods after first switch"
		if "`by'" !=""{	
			di as input _skip(35) "By: `by' = `val_lab_int_XX'"
		}	
		di as input "{hline 80}"
		noi matlist desmat, format(%9.4gc)
		di as input "{hline 80}"
		// New
		noi di "{it: Treatment paths detected in at least `=strtrim("`des_per'")'% of the switching groups `tot_switch' for which `effects' effects could be estimated.}"
		noi di "{it: Total % = `=strtrim("`last_p'")' %}"
		
		// Reference line //
		local n_groups = N_XX[1]
		local d_start = treatment_XX1[1]
		if `des_n' == 1 {
			local d_vec : di %9.2gc treatment_XX2[1]
			local d_vec = strtrim("`d_vec'")
		}
		else {
			local d_vec = "vector ("
			forv j = 1/`des_n' {
				local k : di %9.2gc treatment_XX`=`j' + 1'[1]
				local k = strtrim("`k'")
				local d_vec = "`d_vec'`k', "
			}
			local d_vec = substr("`d_vec'", 1, length("`d_vec'")-2) + ")"
		}
		noi di as text "{it: Design interpretation (first row):}"
		noi di as text "{it: `n_groups' groups started with treatment `d_start' and then experienced treatment `d_vec'.}"		
		//noi di _newline
	}
	else if "`des_path'" != "console" & "`des_path'" !=  "" {
		if "`by'"==""{
		 putexcel set "`des_path'", replace sheet("Design", replace)
		 putexcel A1 = matrix(desmat, names)
		 noi di as input "Design exported to `des_path'"
		}
		else if "`by'"!=""{
			putexcel set "`des_path'", modify sheet("Design `by' = `val_lab_int_XX'", replace)
		 putexcel A1 = matrix(desmat, names)
		 noi di as input "Design exported to `des_path'"
	}  
	}
	}
}

*** Modif Diego: Date first switch option
if "`date_first_switch'" != "" {
	qui {

	use `data_predesign', clear
	//keep if L_g_XX >= l_XX //New
	
	local dfs_opt = strtrim(substr("`date_first_switch'", 1, strpos("`date_first_switch'", ",") - 1))
	local dfs_path = strtrim(substr("`date_first_switch'", strpos("`date_first_switch'", ",") + 1, .))
	
	if "`dfs_opt'" != "" & "`dfs_opt'" != "by_baseline_treat" {
		di as error ""
		di as error "Only option {bf:by_baseline_treat} allowed"
		exit
	}
	
	// Drop non switchers
	drop if F_g_XX == T_max_XX+1 | F_g_XX == .
	keep if time_XX == F_g_XX // Keep just 1 observation per group
	keep `2' `3' F_g_XX d_sq_XX
	
	
	if missing("`dfs_opt'") {
		gcollapse (count) `2', by(`3')
		tostring `3', replace
		sum `2'
		local tot_s = r(sum)
		gen share_XX = `2'/`tot_s' * 100
		
		mkmat `2' share_XX, matrix(dfs) rownames(`3')
		mat colnames dfs = #Groups %Groups
		if "`dfs_path'" == "console" {		
			noi display _newline
			di as input "{hline 48}"
			di as input _skip(10) "Switching dates"
			if "`by'" !=""{	
				di as input _skip(10) "By: `by' = `val_lab_int_XX'"
			}	
			di as input "{hline 48}"
			noi matlist dfs, format(%9.4gc) rowtitle("Any status quo treat.") twidth(25) // New text
			di as input "{hline 48}"		
			ereturn scalar switch_dates = `=_N' // New return instead of text
			//noi di _newline
		}
		else if "`dfs_path'" != "console" & "`dfs_path'" !=  "" {
			if "`by'"==""{
				putexcel set "`dfs_path'", replace sheet("Switching dates", replace)
				putexcel A1 = matrix(dfs, names)
				noi di as input "Switching dates exported to `dfs_path'"
				//noi di _newline
			}
			else if "`by'"!=""{
				putexcel set "`dfs_path'", modify sheet("Switching dates `by' = `val_lab_int_XX'", replace)
				putexcel A1 = matrix(dfs, names)
				noi di as input "Switching dates exported to `dfs_path'"
				//noi di _newline	
			}
		}		
	}
	
	else {
		gcollapse (count) `2', by(`3' d_sq_XX)
		tostring `3', replace
		gegen tot_share_XX = sum(`2'), by(d_sq_XX)
		gen share_XX = `2'/tot_share_XX * 100
		sort d_sq_XX `3'
		
		noi display _newline
		if "`dfs_path'" == "console" {
			di as input "{hline 48}"
			di as input _skip(10) "Switching dates"
			di as input "{hline 48}"
		}
		
		tostring d_sq_XX, replace
		levelsof d_sq_XX, local(levels_d_sq_XX)
		local i = 1
		foreach m of local levels_d_sq_XX {
			mkmat `2' share_XX if d_sq_XX == "`m'", matrix(dfs`i') rownames(`3')
			mat colnames dfs`i' = #Groups %Groups
			sum `2'
			local tot_s = r(sum)
			if "`dfs_path'" == "console" {		
				noi matlist dfs`i', format(%9.4gc) rowtitle("Status quo treat. = `m'") twidth(25) // New
				di as input "{hline 48}"
				sum `2' if d_sq_XX == "`m'"
				local tot_s = r(sum)
				count if d_sq_XX == "`m'"
				ereturn scalar switch_dates_Dsq_`l' = `r(N)' // New
			}
			else if "`dfs_path'" != "console" & "`dfs_path'" !=  "" {
				if `i' == 1 {
					local rep_opt "replace"
				}
				else {
					local rep_opt "modify"
				}
				
				if "`by'"==""{
					putexcel set "`dfs_path'", `rep_opt' sheet("Switch. Dates Base treat `m'", replace)
					putexcel A1 = matrix(dfs`i', names)
					noi di _newline
				}
				else if "`by'"!=""{
					putexcel set "`dfs_path'", modify sheet("Dates Base treat=`m' & `by' = `val_lab_int_XX'", replace) //`by' = `val_lab_int_XX'
					putexcel A1 = matrix(dfs`i', names)
					noi di _newline
				}
			}	
			local i = `i' + 1			
		}
		if "`dfs_path'" != "console" & "`dfs_path'" != "" {
			noi di as input "Switching dates exported to `dfs_path'"
		}
	}
	//noi di _newline
	}
}

use `data_predesign', clear

// Modif Felix: Reloade the data in case we deleted something
if "`by'" !=""{	
use "`by_data_XX'.dta", clear
local rownames ""

// Save the matrices to generate multiple graphs
matrix mat_res_XX_by_`l'=mat_res_XX
}

} // Modif Felix: End of the byvar loop



*****  Producing a graph. ****************************************************
qui {

// Create own color pattern order -> simplify if CI's same color
local col1 "midblue"
local col2 "midblue" //navy
local col3 "red"  
local col4 "red" //maroon    
local col5 "green" 
local col6 "green"    
local col7 "magenta" 
local col8 "magenta"   
local col9 "gold"   
local col10 "gold"    
local col11 "lime" 
local col12 "lime"      

local graph_input ""
local graph_options ""
local graph_options_int ""

if "`by'" ==""|scalar(by_XX)==0{	
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

local graph_input "(connected point_estimate time_to_treat, lpattern(solid)) (rcap up_CI_95 lb_CI_95 time_to_treat)"
local graph_options "legend(off)"

}

// Modif Felix: Prepare data for multiple graphs 
if "`by'" !=""{
		
*** Save tempfile because we will do sorting afterwards 
tempfile before_graph_XX
save "`before_graph_XX'.dta", replace		
		
local cnt_XX=1
local merge_count=1
foreach l of local levels_byvar_XX {
use "`before_graph_XX'.dta", clear		
	
capture drop mat_res_XX_by_`l'1 mat_res_XX_by_`l'2 mat_res_XX_by_`l'3 mat_res_XX_by_`l'4 mat_res_XX_by_`l'5 mat_res_XX_by_`l'6 mat_res_XX_by_`l'7
svmat mat_res_XX_by_`l' //rename the obtained variables ? yes!
////Doulo: use forvalue instead of foreach


forvalue index=1/6{
replace mat_res_XX_by_`l'`index'=0 if mat_res_XX_by_`l'7==0 // artificially replacing the result of average effect by 0
}

capture drop point_estimate_`l'
capture drop se_point_estimate_`l'
capture drop lb_CI_95_`l'
capture drop up_CI_95_`l'
capture drop N_`l'
capture drop N_switchers_`l'
capture drop time_to_treat

rename mat_res_XX_by_`l'1 point_estimate_`l'
rename mat_res_XX_by_`l'2 se_point_estimate_`l'
rename mat_res_XX_by_`l'3 lb_CI_95_`l'
rename mat_res_XX_by_`l'4 up_CI_95_`l'
rename mat_res_XX_by_`l'5 N_`l'
rename mat_res_XX_by_`l'6 N_switchers_`l'
rename mat_res_XX_by_`l'7 time_to_treat

local col_CI=`cnt_XX'+1

// Create a local to just spicify graph once
local graph_input "`graph_input' (connected point_estimate_`l' time_to_treat, lpattern(solid) lcolor(`col`cnt_XX'') mcolor(`col`cnt_XX'')) (rcap up_CI_95_`l' lb_CI_95_`l' time_to_treat, lcolor(`col`col_CI''))" 

*** Adapt if no value labels
ds `by', has(vallabel)
local `by'_has_vallabel_XX=r(varlist)
if "``by'_has_vallabel_XX'"!="."{
	local val_lab_int_XX "`: label ``by'_lab_XX' `l''"
} 
if "``by'_has_vallabel_XX'"=="."{
	local val_lab_int_XX "`l'"
}

local graph_options_int "`graph_options_int' `cnt_XX' "`by' = `val_lab_int_XX'""

local cnt_XX=`cnt_XX'+2 // +2 to skip the CI's

*** Save the specific dataset 
keep point_estimate_`l' se_point_estimate_`l' lb_CI_95_`l' up_CI_95_`l' N_`l' N_switchers_`l' time_to_treat

tempfile graph_`l'_XX
	//Doulo:
keep if time_to_treat!=. //to allow time_to_treat as id of the dataset and use merge 1:1 (safer)

save "`graph_`l'_XX'.dta", replace	

if `merge_count'==1{
	use "`graph_`l'_XX'.dta", clear
	sort time_to_treat
	tempfile graph_merge_XX
	save "`graph_merge_XX'.dta", replace
}

if `merge_count'>1{
	use "`graph_merge_XX'.dta", clear
	merge 1:1 time_to_treat using "`graph_`l'_XX'.dta"
	capture drop _merge
	tempfile graph_merge_XX
	save "`graph_merge_XX'.dta", replace
}	

local merge_count=`merge_count'+1

}	

local graph_options "legend(pos(6) order(`graph_options_int') rows(1))"

}	

}

if ("`graph_off'"==""){
if "`graphoptions'"==""{
twoway `graph_input', xlabel(`=-l_placebo_XX'[1]`=l_XX') xtitle("Relative time to last period before treatment changes (t=0)", size(large)) title("DID, from last period before treatment changes (t=0) to t", size(large)) graphregion(color(white)) plotregion(color(white)) `graph_options'
}
else{
global options="`graphoptions'"
twoway `graph_input', $options
}
}


//Saving the results if requested -> Modif Felix: maybe put in "if by" because the Varaibles have other names when by is specified
if "`save_results'"!=""{
quietly{
//keep if point_estimate!=.
keep if time_to_treat!=.
}

if "`by'"==""{
keep time_to_treat point_estimate se_point_estimate lb_CI_95 up_CI_95 N N_switchers
}
else if "`by'"!=""{
keep time_to_treat point_estimate* se_point_estimate* lb_CI_95* up_CI_95* N* N_switchers*
}

label data "Stores did_multiplegt_dyn estimates' information: Type 'notes' for details." 

// notes : "INFORMATION ABOUT THE DATASET"

local command_options_XX effects(`effects') placebo(`placebo') switchers(`switchers') controls(`controls') trends_nonparam(`trends_nonparam') weight(`weight') `dont_drop_larger_lower' `normalized' cluster(`cluster') `same_switchers' `trends_lin' by(`by') predict_het(`predict_het') ci_level(`ci_level') design(`design') date_first_switch(`date_first_switch') continuous(`continuous') `less_conservative_se'

notes : "{bf:{ul:Date of Run:}} `c(current_date)' at `c(current_time)'"
notes : "{bf:{ul:Command Syntax:}} did_multiplegt `1' `2' `3' `4' `if' `in', `command_options_XX' "
*Use the local dataset_name_XX
notes : "{bf:{ul:Path of the Used Dataset:}} `dataset_name_XX' "
save "`save_results'", replace

}

display _newline
di as text "The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899)."

//>>
restore

if "`save_sample'" != "" { // foreach l of local levels_byvar_XX {
	if "`by'"==""{
	qui merge m:1 `2' `3' using `did_sample', gen(_merge_XX)
	qui drop if _merge_XX == 2
	qui drop _merge_XX
	}
	else if "`by'"!=""{
	foreach l of local levels_byvar_XX {	
	qui merge m:1 `2' `3' using `did_sample_`l'', gen(_merge_XX)
	qui drop if _merge_XX == 2
	qui drop _merge_XX	
	}
	}
}
end

********************************************************************************
*                                 PROGRAM 2                                    *
********************************************************************************

capture program drop did_multiplegt_dyn_core_new

program did_multiplegt_dyn_core_new, eclass
	version 12.0
	syntax varlist(min=4 max=4 numeric) [, effects(integer 1) placebo(integer 0) switchers(string) controls(varlist numeric) trends_nonparam(varlist numeric) weight(varlist numeric) dont_drop_larger_lower NORMalized cluster(varlist numeric) graphoptions(string) SAVe_results(string) graph_off same_switchers same_switchers_pl effects_equal  drop_if_d_miss_before_first_switch trends_lin ci_level(integer 95) by(varlist numeric max=1) predict_het(string) design(string) date_first_switch(string) normalized_weights(string) CONTinuous(string) switchers_core(string) less_conservative_se]
	
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

// Drop some variables we only create once
capture drop num_g_paths_0_XX
capture drop cohort_fullpath_0_XX

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
capture drop U_Gg`i'_temp_var_XX
capture drop U_Gg`i'_var_XX
capture drop U_Gg`i'_var_2_XX
capture drop count_var_`i'_ntreat_XX_temp
capture drop count_var_`i'_ntreat_XX
capture drop count_var_`i'_treat_XX_temp
capture drop count_var_`i'_treat_XX
capture drop avg_diff_y_`i'_tnp_XX
capture drop dof_cohort_`i'_ns_t_XX
capture drop dof_cohort_`i'_s_t_XX
capture drop count_cohort_`i'_ns_t_XX
capture drop count_cohort_`i'_s_t_XX
capture drop total_cohort_`i'_ns_t_XX
capture drop total_cohort_`i'_s_t_XX
capture drop mean_cohort_`i'_ns_t_XX
capture drop mean_cohort_`i'_s_t_XX

capture drop never_change_d_`i'_wXX
capture drop distance_to_switch_`i'_wXX

// Modif: Doulo Paths
capture drop d_fg`i'_XX
capture drop path_0_XX 
capture drop path_`i'_XX
capture drop d_fg0_XX
capture drop num_g_paths_`i'_XX
capture drop cohort_fullpath_`i'_XX

capture drop dof_cohort_`i'_ns_t_XX
capture drop dof_cohort_`i'_s0_t_XX
capture drop dof_cohort_`i'_s1_t_XX
capture drop dof_cohort_`i'_s2_t_XX
capture drop count_cohort_`i'_ns_t_XX
capture drop count_cohort_`i'_s0_t_XX
capture drop count_cohort_`i'_s1_t_XX
capture drop count_cohort_`i'_s2_t_XX
capture drop total_cohort_`i'_ns_t_XX
capture drop total_cohort_`i'_s0_t_XX
capture drop total_cohort_`i'_s1_t_XX
capture drop total_cohort_`i'_s2_t_XX
capture drop mean_cohort_`i'_ns_t_XX
capture drop mean_cohort_`i'_s0_t_XX
capture drop mean_cohort_`i'_s1_t_XX
capture drop mean_cohort_`i'_s2_t_XX

////// Creating long difference of outcome //
xtset group_XX time_XX
bys group_XX : gen diff_y_`i'_XX = outcome_XX - L`i'.outcome_XX

// Modif: Doulo Paths
if "`less_conservative_se'" != ""{
gen d_fg0_XX=d_sq_XX

gen d_fg_XX_temp=treatment_XX if time_XX==F_g_XX + `=-1+`i''
bys group_XX: gegen d_fg`i'_XX=mean(d_fg_XX_temp)
replace d_fg`i'_XX=d_fg`=`i'-1'_XX if d_fg`i'_XX==. //& F_g_XX==T_max_XX+1
drop d_fg_XX_temp 

gegen path_0_XX=group(d_fg0_XX F_g_XX) // path_0_XX= group(d_fg0_XX F_g_XX) ?
gegen path_`i'_XX = group(path_`=`i'-1'_XX d_fg`i'_XX) if path_`=`i'-1'_XX!=.

// For each group_XX define a variable which tells us how many other groups have the same path i.e. are in the same cohort for each of the definitions
if `i'==1{ // generate for path_0 once
	bysort path_0_XX: gegen num_g_paths_0_XX=nunique(group_XX)
}
bysort path_`i'_XX: gegen num_g_paths_`i'_XX=nunique(group_XX)

// generate a dummy for each group determining wether they are in a chohrt of size>1 and then build the E_hat paths in different if conditions based on that (for each group but this does not need to be additionally specified as this dummy is at the group level)
bys group_XX: gen cohort_fullpath_`i'_XX=(num_g_paths_`i'_XX>1) // I think bys not needed
if `i'==1{ // generate for path_0 once
	bys group_XX: gen cohort_fullpath_0_XX=(num_g_paths_0_XX>1) 
}

**** Doublecheck that I did the numbering in this explaination correct, is path_0 actually path_1 or reverse? ****
// generate the same dummy for the reduced paths 
// 0) For version a) (D_{g,1},F_g, D_{g,F_g},...,D_{g,F_g-1+\ell}) -> if cohort_fullpath_`i'_XX==1
// 1) For version b) (D_{g,1},F_g, D_{g,F_g}) we can just rely on cohort_fullpath_1_XX and then set path to just D_{g,F_g} (D_{g,1} and F_g will always be in bysort and not a part of `paths' -> actually I think D_{g,1} is included in paths as path_0_XX so we only need to add F_g_XX) -> if cohort_fullpath_1_XX==1
// 2) If I redefine path_0_XX then I could just use path_0_XX here -> if cohort_fullpath_0_XX==1
// 3) In the last case I can just use all that have num_g_paths_0_XX==1, so they have a 0 in the dummy defined in 2) -> if cohort_fullpath_0_XX==0, can use the case without demeaning we should already have when there is only one group with this path which should be the case here
// Implement this starting by the most general case (3) and then replace if the more strict cases are met
	
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

if ("`same_switchers_pl'"!=""){ // Modif Felix: new option for trends_lin
//Generate a variable tagging the switchers that should be dropped
gen relevant_y_missing_XX=(outcome_XX==.&time_XX>=F_g_XX-1-`=`placebo''&time_XX<=F_g_XX-1+`=`effects'') // at least one of the effects we try to estimate is missing
if "`controls'" != ""{
replace relevant_y_missing_XX=1 if fd_X_all_non_missing_XX==0&time_XX>=F_g_XX-`=`placebo''&time_XX<=F_g_XX-1+`=`effects'' // Is a -1 missing here???????? at F_g_XX-`=`placebo'' ????
}

bys group_XX: gen cum_fillin_XX = sum(relevant_y_missing_XX)
gen dum_fillin_temp_XX = (cum_fillin_XX==0&time_XX==F_g_XX-1+`=`effects'')
bys group_XX: gegen fillin_g_XX = total(dum_fillin_temp_XX)


gen still_switcher_`i'_XX = (F_g_XX-1+`=`effects''<=T_g_XX&fillin_g_XX>0) // tag switchers who are observable for all periods up to the specified number of effects 	

gen distance_to_switch_`i'_XX=(still_switcher_`i'_XX&time_XX==F_g_XX-1+`i'&`i'<=L_g_XX&S_g_XX==increase_XX&N_gt_control_`i'_XX>0&N_gt_control_`i'_XX!=.) if diff_y_`i'_XX!=. 
}	
if ("`same_switchers_pl'"==""){
//Generate a variable tagging the switchers that should be dropped
gen relevant_y_missing_XX=(outcome_XX==.&time_XX>=F_g_XX-1&time_XX<=F_g_XX-1+`=`effects'') // at least one of the effects we try to estimate is missing
if "`controls'" != ""{
replace relevant_y_missing_XX=1 if fd_X_all_non_missing_XX==0&time_XX>=F_g_XX&time_XX<=F_g_XX-1+`=`effects''
}

bys group_XX: gen cum_fillin_XX = sum(relevant_y_missing_XX)
gen dum_fillin_temp_XX = (cum_fillin_XX==0&time_XX==F_g_XX-1+`=`effects'')
bys group_XX: gegen fillin_g_XX = total(dum_fillin_temp_XX)


gen still_switcher_`i'_XX = (F_g_XX-1+`=`effects''<=T_g_XX&fillin_g_XX>0) // tag switchers who are observable for all periods up to the specified number of effects 	

gen distance_to_switch_`i'_XX=(still_switcher_`i'_XX&time_XX==F_g_XX-1+`i'&`i'<=L_g_XX&S_g_XX==increase_XX&N_gt_control_`i'_XX>0&N_gt_control_`i'_XX!=.) if diff_y_`i'_XX!=.  
}

}
else{
gen distance_to_switch_`i'_XX=(time_XX==F_g_XX-1+`i'&`i'<=L_g_XX&S_g_XX==increase_XX&N_gt_control_`i'_XX>0&N_gt_control_`i'_XX!=.) if diff_y_`i'_XX!=. 
}

///// Computing N^1_{t,l} or N^0_{t,l}. //
gen distance_to_switch_`i'_wXX = distance_to_switch_`i'_XX*N_gt_XX
bys time_XX: gegen N`=increase_XX'_t_`i'_XX=total(distance_to_switch_`i'_wXX)

///// Computing N^1_l or N^0_l. //

// Initializing the N1_`i'_XX/N0_`i'_XX scalar at 0. 
scalar N`=increase_XX'_`i'_XX =0 // Modif Felix: see if we do not initialize it again after each repitition -> thats not what we want

forvalue t=`=t_min_XX'/`=T_max_XX'{
	sum N`=increase_XX'_t_`i'_XX if time_XX==`t'
	scalar N`=increase_XX'_`i'_XX = N`=increase_XX'_`i'_XX + r(mean)
}

///// Computing N^0_{t,l,g} or N^1_{t,l,g}. //
bys time_XX d_sq_XX `trends_nonparam': gegen N`=increase_XX'_t_`i'_g_XX=total(distance_to_switch_`i'_wXX)


///// Creating long differences of control variables //
if "`controls'" != ""{

*********************************** New variances controls (Modif Felix)	
capture drop part2_switch`=increase_XX'_`i'_XX
gen part2_switch`=increase_XX'_`i'_XX=0

// Number of control groups with same status-quo treatement that did not change treatement yet	
capture drop dummy_ctrl_U_Gg`i'_XX

** (Modif Diego) Deprecated since D_g_1 should equate the levels of sq treatment (after grouping) 
// Creating a dummy variable indicating whether l<=T_g_XX-2 & D_g_1=d 
//gen dummy_ctrl_U_Gg`i'_XX = (`i'<=T_g_XX-2 & d_sq_XX==treatment_XX) // replace d_sq_XX==treatment_XX by d_sq_XX==`l' and put it in the corresponding loop

** Modif. Diego: generation of the T_d variable = max_g:D_g,1 = d F_d - 1
capture drop T_d_XX
gegen T_d_XX = max(F_g_XX), by(d_sq_int_XX)
replace T_d_XX = T_d_XX - 1

**********************************************************	
		
local count_controls=0

// Computing the first differences of the control variables
foreach var of varlist `controls'{
		
local count_controls=`count_controls'+1

capture drop diff_X`count_controls'_`i'_XX

xtset group_XX time_XX
gen diff_X`count_controls'_`i'_XX=`var' - L`i'.`var'


*********************************** New variances controls (Modif Felix)
*** Note: I will have to adapt the placebos as well!!!
**# Variance Control Setup

* intermediate steps
*N_g_t * (X_g_t - X_g_t-l)
capture drop diff_X`count_controls'_`i'_N_XX
gen diff_X`count_controls'_`i'_N_XX = N_gt_XX * diff_X`count_controls'_`i'_XX 

foreach l of local levels_d_sq_XX { //l is d in the model //

capture drop dummy_XX
capture drop N_d_t_XX
gen dummy_XX=(F_g_XX>time_XX & d_sq_int_XX == `l')
gegen N_d_t_XX=total(dummy_XX) if d_sq_int_XX == `l' // We never use this !!!

* small m
/*
** Building blocks of m^*_gdl (Debug) **
gen block1_`=increase_XX'_g_`count_controls'_`l'_`i'_XX = (`i' <= T_g_XX-2 & d_sq_int_XX == `l')
gen block2_`=increase_XX'_g_`count_controls'_`l'_`i'_XX = (G_XX / N`=increase_XX'_`i'_XX)
gen block3_`=increase_XX'_g_`count_controls'_`l'_`i'_XX = [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX]
gen block4_`=increase_XX'_g_`count_controls'_`l'_`i'_XX = (time_XX>=`=`i'+1'&time_XX<=T_g_XX)
gen block5_`=increase_XX'_g_`count_controls'_`l'_`i'_XX = diff_X`count_controls'_`i'_N_XX
*/

capture drop m`=increase_XX'_g_`count_controls'_`l'_`i'_XX // new
gen m`=increase_XX'_g_`count_controls'_`l'_`i'_XX = (`i' <= T_g_XX-2 & d_sq_int_XX == `l')* (G_XX / N`=increase_XX'_`i'_XX) * ([distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * diff_X`count_controls'_`i'_N_XX)	

* capital M
capture drop m`=increase_XX'_`l'_`count_controls'_`i'_XX // new
bys group_XX: gegen m`=increase_XX'_`l'_`count_controls'_`i'_XX=total(m`=increase_XX'_g_`count_controls'_`l'_`i'_XX) // use mean and delete second part

bys group_XX: replace m`=increase_XX'_`l'_`count_controls'_`i'_XX = . if _n != 1

capture drop M`=increase_XX'_`l'_`count_controls'_`i'_XX // new
egen M`=increase_XX'_`l'_`count_controls'_`i'_XX = total(m`=increase_XX'_`l'_`count_controls'_`i'_XX)

* is this the correct bysort or should it be bysort time_XX d_sq_XX???
replace M`=increase_XX'_`l'_`count_controls'_`i'_XX = (1/G_XX)*M`=increase_XX'_`l'_`count_controls'_`i'_XX

// the E_hat defined above

capture drop E_hat_t`count_controls'_`l'_XX	
capture drop E_hat_denom_`count_controls'_`l'_XX
capture drop E_hat_t`count_controls'_`l'_temp_XX
bys time_XX : egen E_hat_denom_`count_controls'_`l'_XX = total(dummy_XX) if d_sq_int_XX == `l'

gen E_hat_t`count_controls'_`l'_temp_XX = (prod_X`count_controls'_diff_y_int_XX * dummy_XX) / E_hat_denom_`count_controls'_`l'_XX
bys d_sq_int_XX time_XX : egen E_hat_t`count_controls'_`l'_XX = total(E_hat_t`count_controls'_`l'_temp_XX)

//bys G: gegen E_hat_t`count_controls'_`l'_XX = mean(prod_X`count_controls'_diff_y_int_XX) if d_sq_int_XX == `l' & time_XX < F_g_XX

//bysort d_sq_XX time_XX: gegen E_hat_d_t`count_controls'_XX=mean(prod_X`count_controls'_diff_y_int_XX) if time_XX<F_g_XX & treatment_XX==d_sq_XX // do I also need to include the treatement level condition of the sum somehow (latter part of the if condition)???


// define part inside the sum (sum [t=2 -> F_g-1])	
// sum up the count conrol specific parts as those are Vector products
*capture drop in_sum_temp_`count_controls'_XX
*gen in_sum_temp_`count_controls'_XX=prod_X`count_controls'_diff_y_int_XX-(N_d_t_XX>=2)*sqrt((N_d_t_XX)/(N_d_t_XX-1))*E_hat_d_t`count_controls'_XX
capture drop in_sum_temp_`count_controls'_`l'_XX
//gen in_sum_temp_`count_controls'_`l'_XX=prod_X`count_controls'_diff_y_int_XX-(N_d_t_XX>=2)*sqrt((N_d_t_XX)/(N_d_t_XX-1))*E_hat_t`count_controls'_`l'_XX

capture drop N_c_`l'_temp_XX
capture drop N_c_`l'_XX
gen N_c_`l'_temp_XX = d_sq_int_XX == `l' & time_XX >= 2 & time_XX <= T_d_XX & time_XX < F_g_XX
egen N_c_`l'_XX = total(N_c_`l'_temp_XX)
gen in_sum_temp_`count_controls'_`l'_XX = (prod_X`count_controls'_diff_y_int_XX-(E_hat_denom_`count_controls'_`l'_XX>=2)*sqrt((E_hat_denom_`count_controls'_`l'_XX)/(E_hat_denom_`count_controls'_`l'_XX - 1))*E_hat_t`count_controls'_`l'_XX * (time_XX>=2 & time_XX<=F_g_XX-1)) / N_c_`l'_XX
	
// define the inner sum (sum [t=2 -> F_g-1])
capture drop in_sum_`count_controls'_`l'_XX

bys group_XX: gegen in_sum_`count_controls'_`l'_XX = total(in_sum_temp_`count_controls'_`l'_XX) 
}

**********************************************************


foreach l of local levels_d_sq_XX { 
	if (scalar(useful_res_`l'_XX)>1){ 

replace diff_y_`i'_XX = diff_y_`i'_XX - coefs_sq_`l'_XX[`=`count_controls'',1]*diff_X`count_controls'_`i'_XX if d_sq_int_XX==`l' 

////// N.B. : in the above line, we do not add "&diff_X`count_controls'_`i'_XX!=." because we want to exclude from the estimation any first/long-difference for which the covariates are missing.

*********************************** New variances controls (Modif Felix)
capture drop in_brackets_`l'_`count_controls'_XX	
gen in_brackets_`l'_`count_controls'_XX=0
**********************************************************

}
}

}

}


///// Computing the mean of first or long-differences of outcomes time N_{g,t} for non-treated and for treated separately - will be useful for the computation of the variance //

// Modif Felix: change names because now we will generate two things, the cohorts and the baseline cohorts

capture drop diff_y_`i'_N_gt_XX
gen diff_y_`i'_N_gt_XX=N_gt_XX*diff_y_`i'_XX
capture drop dof_diff_y_`i'_N_gt_XX
gen dof_diff_y_`i'_N_gt_XX=(N_gt_XX!=0&diff_y_`i'_XX!=.)

// Modif: Define number of groups in each cohort / baseline cohort (seperately for switchers and for controls)

// Modif: Define the expected outcome evolution in those cohorts for the demeaning 

*** Cohorts never switchers
* DOF
bys d_sq_XX `trends_nonparam' : gegen dof_cohort_`i'_ns_t_XX=total(dof_diff_y_`i'_N_gt_XX) if diff_y_`i'_XX!=.&never_change_d_`i'_XX==1&N`=increase_XX'_t_`i'_XX>0&N`=increase_XX'_t_`i'_XX!=.

* Denominator
bys d_sq_XX `trends_nonparam' : gegen count_cohort_`i'_ns_t_XX=total(N_gt_XX) if diff_y_`i'_XX!=.&never_change_d_`i'_XX==1&N`=increase_XX'_t_`i'_XX>0&N`=increase_XX'_t_`i'_XX!=.

* Numerator
bys d_sq_XX `trends_nonparam' : gegen total_cohort_`i'_ns_t_XX=total(diff_y_`i'_N_gt_XX) if never_change_d_`i'_XX==1&N`=increase_XX'_t_`i'_XX>0&N`=increase_XX'_t_`i'_XX!=.

* Estimator for the expectation (no need for bysort or any conditioning as the Numerator and denominator are generated along the same set of conditions)
gen mean_cohort_`i'_ns_t_XX=total_cohort_`i'_ns_t_XX/count_cohort_`i'_ns_t_XX


*** Cohorts switchers
if "`less_conservative_se'" == ""{
* DOF
bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam' : gegen dof_cohort_`i'_s_t_XX=total(dof_diff_y_`i'_N_gt_XX) if distance_to_switch_`i'_XX==1

* Denominator
bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam' : gegen count_cohort_`i'_s_t_XX=total(N_gt_XX) if distance_to_switch_`i'_XX==1

* Numerator
bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam' : gegen total_cohort_`i'_s_t_XX=total(diff_y_`i'_N_gt_XX) if distance_to_switch_`i'_XX==1 

* Estimator for the expectation (no need for bysort or any conditioning as the Numerator and denominator are generated along the same set of conditions)
gen mean_cohort_`i'_s_t_XX=total_cohort_`i'_s_t_XX/count_cohort_`i'_s_t_XX
}

if "`less_conservative_se'" != ""{
** DOF 
* switcher path0 (s0)
bys path_0_XX `trends_nonparam' : gegen dof_cohort_`i'_s0_t_XX=total(dof_diff_y_`i'_N_gt_XX) if distance_to_switch_`i'_XX==1&cohort_fullpath_0_XX==1&cohort_fullpath_1_XX==0
if `i'>1{ // otherwise equivalent to fullpath
* switcher path1 (s1)
bys path_1_XX `trends_nonparam' : gegen dof_cohort_`i'_s1_t_XX=total(dof_diff_y_`i'_N_gt_XX) if distance_to_switch_`i'_XX==1&cohort_fullpath_1_XX==1&cohort_fullpath_`i'_XX==0	
}
* switcher fullpath (s2)
bys path_`i'_XX `trends_nonparam' : gegen dof_cohort_`i'_s2_t_XX=total(dof_diff_y_`i'_N_gt_XX) if distance_to_switch_`i'_XX==1&cohort_fullpath_`i'_XX==1

** Denominator	
* switcher path0 (s0)
bys path_0_XX `trends_nonparam' : gegen count_cohort_`i'_s0_t_XX=total(N_gt_XX) if distance_to_switch_`i'_XX==1&cohort_fullpath_0_XX==1&cohort_fullpath_1_XX==0
if `i'>1{ // otherwise equivalent to fullpath
* switcher path1 (s1)
bys path_1_XX `trends_nonparam' : gegen count_cohort_`i'_s1_t_XX=total(N_gt_XX) if distance_to_switch_`i'_XX==1&cohort_fullpath_1_XX==1&cohort_fullpath_`i'_XX==0	
}
* switcher fullpath (s2)
bys path_`i'_XX `trends_nonparam' : gegen count_cohort_`i'_s2_t_XX=total(N_gt_XX) if distance_to_switch_`i'_XX==1&cohort_fullpath_`i'_XX==1	
	
** Numerator
* switcher path0 (s0)
bys path_0_XX `trends_nonparam' : gegen total_cohort_`i'_s0_t_XX=total(diff_y_`i'_N_gt_XX) if distance_to_switch_`i'_XX==1&cohort_fullpath_0_XX==1&cohort_fullpath_1_XX==0
if `i'>1{ // otherwise equivalent to fullpath
* switcher path1 (s1)
bys path_1_XX `trends_nonparam' : gegen total_cohort_`i'_s1_t_XX=total(diff_y_`i'_N_gt_XX) if distance_to_switch_`i'_XX==1&cohort_fullpath_1_XX==1&cohort_fullpath_`i'_XX==0	
}
* switcher fullpath (s2)
bys path_`i'_XX `trends_nonparam' : gegen total_cohort_`i'_s2_t_XX=total(diff_y_`i'_N_gt_XX) if distance_to_switch_`i'_XX==1&cohort_fullpath_`i'_XX==1

** Estimator for the expectation (no need for bysort or any conditioning as the Numerator and denominator are generated along the same set of conditions)
* switcher path0 (s0)
gen mean_cohort_`i'_s0_t_XX=total_cohort_`i'_s0_t_XX/count_cohort_`i'_s0_t_XX
* switcher path1 (s1)
if `i'>1{
gen mean_cohort_`i'_s1_t_XX=total_cohort_`i'_s1_t_XX/count_cohort_`i'_s1_t_XX
}
* switcher fullüpath (s2)
gen mean_cohort_`i'_s2_t_XX=total_cohort_`i'_s2_t_XX/count_cohort_`i'_s2_t_XX
	
}	

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
 
// Modif Felix: Adjusted this part to the new cohorts/cohorts 
///// Computing the "alternative" U_{G,g,l} which will be used for the computation of the variance only - these are like the above U_{G,g,l}s, except that the outcome differences are demeaned, and there is a DOF adjustment when possible//

gen U_Gg`i'_temp_var_XX = 0

// For controls
replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * diff_y_`i'_N_gt_XX  if never_change_d_`i'_XX==1&dof_cohort_`i'_ns_t_XX==1 // Do I need to change cohort to dof in the conditions as well? -> Yes, see indicator function

* For all t < Fg -> only need this case for controls
replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * N_gt_XX * (diff_y_`i'_XX - (never_change_d_`i'_XX * sqrt(dof_cohort_`i'_ns_t_XX/(dof_cohort_`i'_ns_t_XX-1)) * mean_cohort_`i'_ns_t_XX)) if never_change_d_`i'_XX==1&dof_cohort_`i'_ns_t_XX>1&dof_cohort_`i'_ns_t_XX!=.


// For switchers
if "`less_conservative_se'" == ""{
replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * diff_y_`i'_N_gt_XX if distance_to_switch_`i'_XX==1&dof_cohort_`i'_s_t_XX==1

* For t = Fg − 1 + ℓ -> only need this case for switchers
replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * N_gt_XX * (diff_y_`i'_XX - (distance_to_switch_`i'_XX * sqrt(dof_cohort_`i'_s_t_XX/(dof_cohort_`i'_s_t_XX-1)) * mean_cohort_`i'_s_t_XX)) if distance_to_switch_`i'_XX==1&dof_cohort_`i'_s_t_XX>1&dof_cohort_`i'_s_t_XX!=.
}

if "`less_conservative_se'" != ""{
* Case where we only have one group in a pathway, even with the least strict definition -> no demeaning
if `i'==1{
replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * diff_y_`i'_N_gt_XX if distance_to_switch_`i'_XX==1&(dof_cohort_`i'_s0_t_XX==.&dof_cohort_`i'_s2_t_XX==.)
}

if `i'>1{
replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * diff_y_`i'_N_gt_XX if distance_to_switch_`i'_XX==1&(dof_cohort_`i'_s0_t_XX==.&dof_cohort_`i'_s1_t_XX==.&dof_cohort_`i'_s2_t_XX==.)
}
// should I do conditions with cohort_fullpath or with dof_cohort missing/non-missing? -> If everything in the prior code is correct then all the dof should just be defined if there are at least 2 groups in each pathway given the corresponding definition

* Case where we have multiple groups per pathways when conditioning on D_{g,1},F_g
replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * N_gt_XX * (diff_y_`i'_XX - (distance_to_switch_`i'_XX * sqrt(dof_cohort_`i'_s0_t_XX/(dof_cohort_`i'_s0_t_XX-1)) * mean_cohort_`i'_s0_t_XX)) if distance_to_switch_`i'_XX==1&dof_cohort_`i'_s0_t_XX>1&dof_cohort_`i'_s0_t_XX!=.

* Case where we have multiple groups per pathways when conditioning on D_{g,1},F_g, D_{g,F_g}
if `i'>1{
replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * N_gt_XX * (diff_y_`i'_XX - (distance_to_switch_`i'_XX * sqrt(dof_cohort_`i'_s1_t_XX/(dof_cohort_`i'_s1_t_XX-1)) * mean_cohort_`i'_s1_t_XX)) if distance_to_switch_`i'_XX==1&dof_cohort_`i'_s1_t_XX>1&dof_cohort_`i'_s1_t_XX!=.
}

* Case where we have multiple groups per pathways only when conditioning on D_{g,1},F_g, D_{g,F_g},...,D_{g,F_g-1+\ell}
replace U_Gg`i'_temp_var_XX= dummy_U_Gg`i'_XX*(G_XX / N`=increase_XX'_`i'_XX) * [distance_to_switch_`i'_XX - (N`=increase_XX'_t_`i'_g_XX/N_gt_control_`i'_XX) * never_change_d_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * N_gt_XX * (diff_y_`i'_XX - (distance_to_switch_`i'_XX * sqrt(dof_cohort_`i'_s2_t_XX/(dof_cohort_`i'_s2_t_XX-1)) * mean_cohort_`i'_s2_t_XX)) if distance_to_switch_`i'_XX==1&dof_cohort_`i'_s2_t_XX>1&dof_cohort_`i'_s2_t_XX!=.
}

****************************************************************************
// Modif Felix: Add the additional part with controls 
* How to deal with the `l' and `count_controls' here, we have several variables which dpend on the status-quo level of d (denoted by l) and all the control(X) specific variables are numbered from 1 to k and are not a kx1 vector
if "`controls'"!=""{

levelsof d_sq_int_XX, local(levels_d_sq_XX)
foreach l of local levels_d_sq_XX {	
	
capture drop combined`=increase_XX'_temp_`l'_`i'_XX	
gen combined`=increase_XX'_temp_`l'_`i'_XX=0
		
// loops to do the summations
forvalues j=1/`count_controls'{
forvalues k=1/`count_controls'{		
	
// define the part in brackets -> this is the first matrix sum where we sum the inv_denom rows with the kx1 coulmn vectors!
capture drop in_brackets_`l'_`j'_`k'_temp_XX
gen in_brackets_`l'_`j'_`k'_temp_XX = inv_Denom_`l'_XX[`j',`k'] * in_sum_`k'_`l'_XX * (d_sq_int_XX == `l' & F_g_XX>=3) // G does not show up because it is part of inv_Denom_ and N does not show up because it is already included in the in_sum part (we use mean and not total there)
// step where we have the first sum over the k
*capture drop in_brackets_`l'_`j'_XX
replace in_brackets_`l'_`j'_XX=in_brackets_`l'_`j'_XX+in_brackets_`l'_`j'_`k'_temp_XX

} // end loop over k

replace in_brackets_`l'_`j'_XX=in_brackets_`l'_`j'_XX - coefs_sq_`l'_XX[`j',1]
// Now in this part put the summing with M and the rest which is indexed by j
capture drop combined`=increase_XX'_temp_`l'_`j'_`i'_XX // new
gen combined`=increase_XX'_temp_`l'_`j'_`i'_XX=M`=increase_XX'_`l'_`j'_`i'_XX*in_brackets_`l'_`j'_XX

// step where we have the second sum over the j
replace combined`=increase_XX'_temp_`l'_`i'_XX=combined`=increase_XX'_temp_`l'_`i'_XX+combined`=increase_XX'_temp_`l'_`j'_`i'_XX

} // end loop over j

// Final sum over the status quo treatment (outer sum in the formula)
replace part2_switch`=increase_XX'_`i'_XX=part2_switch`=increase_XX'_`i'_XX+combined`=increase_XX'_temp_`l'_`i'_XX if d_sq_int_XX==`l' // does the if condition here even matter because those are seperate variables depending on l

} // loop over l

/*
// define the second part which is added with controls as a whole
bysort group_XX: gegen part_2_`l'_XX=total(M`=increase_XX'_d_`i'_XX*in_brackets_XX) 
*/


// Making the adjustement
if `=increase_XX'==1{
replace U_Gg`i'_temp_var_XX=U_Gg`i'_temp_var_XX - part2_switch1_`i'_XX 
}

if `=increase_XX'==0{
replace U_Gg`i'_temp_var_XX=U_Gg`i'_temp_var_XX + part2_switch0_`i'_XX 
}
	
}	

// Summing the U_{G,g,l}s over time periods for each group //
bys group_XX: gegen U_Gg`i'_var_XX=total(U_Gg`i'_temp_var_XX)

}
****************************************************************************



if "`normalized'"!=""{
	
	capture drop sum_treat_until_`i'_XX
	capture drop delta_D_`i'_cum_temp_XX
	
	// Modif Felix: redefine this with original treatment if continuous is defined
	if "`continuous'"==""{
		bys group_XX : gegen sum_treat_until_`i'_XX = total(treatment_XX - d_sq_XX) if time_XX>=F_g_XX&time_XX<=F_g_XX-1+`i'&S_g_XX==increase_XX
	}
	else if "`continuous'"!=""{
		bys group_XX : gegen sum_treat_until_`i'_XX = total(treatment_XX_orig - d_sq_XX_orig) if time_XX>=F_g_XX&time_XX<=F_g_XX-1+`i'&S_g_XX==increase_XX
	}
	
gen delta_D_`i'_cum_temp_XX = N_gt_XX/N`=increase_XX'_`i'_XX*[sum_treat_until_`i'_XX* S_g_XX + (1-S_g_XX)*(-sum_treat_until_`i'_XX)] if distance_to_switch_`i'_XX==1 

	sum  delta_D_`i'_cum_temp_XX
	scalar delta_norm_`i'_XX = r(sum) 

}

} 
// End of the loop over l_u_a_XX, at this point we have all the U_Gg_`i'_XX that we need. Thus, we can sum them now for the trends_lin option.

// Modif Felix: Only include if we have all placebos
scalar Ntrendslin=1
forvalue i=1/`=l_u_a_XX' {
scalar Ntrendslin=min(Ntrendslin , N`=increase_XX'_`i'_XX )
}

// Modif Felix: Compute the variances adjusted for trends_lin
if "`trends_lin'"!="" & Ntrendslin != 0  {
	
// Modif Felix: put the computation of the covariances after this loop
// Try to put the covariance calculation here 
// U_Gg_var_glob_`i'_XX not defiend here, use underlying in/out specific variable

/* // Actually we do not need this part here anymore 
forvalue i=1/`=l_u_a_XX'{
if `i'<`=l_u_a_XX'{ 
		forvalue j=`=`i'+1'/`=l_u_a_XX'{

capture drop U_Gg_var_`i'_`j'_XX
capture drop U_Gg_var_`i'_`j'_2_XX
		
// Modif Felix: changed sum to product here	-> Issue, very large values or lots of 0s	
		if ("`normalized'"==""){
			gen U_Gg_var_`i'_`j'_XX = U_Gg`i'_var_XX * U_Gg`j'_var_XX // having a product here instead of a summation
		}
		
			if "`normalized'"!=""{
			gen U_Gg_var_`i'_`j'_XX = U_Gg`i'_var_XX/scalar(delta_norm_`i'_XX) * U_Gg`j'_var_XX/scalar(delta_norm_`j'_XX) // having a product here instead of a summation 
}

		}
}
}
*/	
	
	capture drop U_Gg`=l_u_a_XX'_TL
	capture drop U_Gg`=l_u_a_XX'_var_TL
	
	* Initializing at 0
	gen U_Gg`=l_u_a_XX'_TL = 0
	gen U_Gg`=l_u_a_XX'_var_TL = 0

	// Add some condition that l_u_a_XX>1? 
	
	
	forvalue i=1/`=l_u_a_XX'{
		replace U_Gg`=l_u_a_XX'_TL = U_Gg`=l_u_a_XX'_TL + U_Gg`i'_XX 		
		replace U_Gg`=l_u_a_XX'_var_TL =  U_Gg`=l_u_a_XX'_var_TL + U_Gg`i'_var_XX 
		/*
		if `i'<`=l_u_a_XX'{
		forvalue j=`=`i'+1'/`=l_u_a_XX'{			
		replace U_Gg`=l_u_a_XX'_var_TL = U_Gg`=l_u_a_XX'_var_TL + 2*U_Gg_var_`i'_`j'_XX
		
		}
	}
*/
	}
	
	replace U_Gg`=l_u_a_XX'_XX = U_Gg`=l_u_a_XX'_TL
	replace U_Gg`=l_u_a_XX'_var_XX = U_Gg`=l_u_a_XX'_var_TL
	
	// For Felix : Here, the computation of the U_Gg_var is wrong : we have to add the covariances as in the pdf web_appendix
	
	// the variances itself are already added up here, just covariances missing
	
	// can we just put the U_Gg_var_`i'_`j'_XX here as if they were the covariances and the U_Gg`i'_var_XX as if they were the variances and then compute U_Gg's for the variance of the sum following the formula ???
	// Var(sum(X1,X2)) = sum(Var(X1),Var(X2)) + 2*Cov(X1,X2)
	
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
capture drop dist_to_switch_pl_`i'_XX
capture drop never_change_d_pl_`i'_XX
capture drop N`=increase_XX'_t_placebo_`i'_XX
capture drop N`=increase_XX'_t_placebo_`i'_g_XX
capture drop N_gt_control_placebo_`i'_XX
capture drop dummy_U_Gg_pl_`i'_XX 
capture drop never_change_d_pl_`i'_wXX
capture drop dist_to_switch_pl_`i'_wXX
capture drop dof_cohort_pl_`i'_ns_t_XX
capture drop dof_cohort_pl_`i'_s_t_XX
capture drop count_cohort_pl_`i'_ns_t_XX
capture drop count_cohort_pl_`i'_s_t_XX
capture drop total_cohort_pl_`i'_ns_t_XX
capture drop total_cohort_pl_`i'_s_t_XX
capture drop mean_cohort_pl_`i'_ns_t_XX
capture drop mean_cohort_pl_`i'_s_t_XX

// Modif Felix: new variables control path 
capture drop dof_cohort_pl_`i'_ns_t_XX
capture drop dof_cohort_pl_`i'_s0_t_XX
capture drop dof_cohort_pl_`i'_s1_t_XX
capture drop dof_cohort_pl_`i'_s2_t_XX
capture drop count_cohort_pl_`i'_ns_t_XX
capture drop count_cohort_pl_`i'_s0_t_XX
capture drop count_cohort_pl_`i'_s1_t_XX
capture drop count_cohort_pl_`i'_s2_t_XX
capture drop total_cohort_pl_`i'_ns_t_XX
capture drop total_cohort_pl_`i'_s0_t_XX
capture drop total_cohort_pl_`i'_s1_t_XX
capture drop total_cohort_pl_`i'_s2_t_XX
capture drop mean_cohort_pl_`i'_ns_t_XX
capture drop mean_cohort_pl_`i'_s0_t_XX
capture drop mean_cohort_pl_`i'_s1_t_XX
capture drop mean_cohort_pl_`i'_s2_t_XX

// Create long-differences for the placebos
xtset group_XX time_XX

//The main trick to computing the placebo point estimates is:
// 1. to place the corresponding outcome (y_{F_g-1} - y_{F_g - l - 1})) values in the same row of that (y_{F_g + l -1} - y_{F_g - 1}) of the symmetric DID_l. 
// 2. The other variables, such as N_gt, N0_l or N1_l, remain unchanged, except that we have to check if diff_y_placebo ( = y_{F_g - 2l -2}- y_{F_g - l -1}) exists. 
// 3. If y_{F_g - l -1} does not exist for a specific group, that group is excluded from the calculation, hence, for example, one always has #Switchers for DID_l>= #Switchers for DID_pl.

//Computing the diff_y_placebos
bys group_XX : gen diff_y_pl_`i'_XX = L`=2*`i''.outcome_XX - L`i'.outcome_XX

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
scalar N`=increase_XX'_placebo_`i'_XX=0
forvalue t=`=t_min_XX'/`=T_max_XX'{
	sum N`=increase_XX'_t_placebo_`i'_XX if time_XX==`t'
	scalar N`=increase_XX'_placebo_`i'_XX = N`=increase_XX'_placebo_`i'_XX + r(mean)
}


///// Computing N^0_{t,l,g} or N^1_{t,l,g}. 

bys time_XX d_sq_XX `trends_nonparam': gegen N`=increase_XX'_t_placebo_`i'_g_XX=total(dist_to_switch_pl_`i'_wXX)


///// Creating long differences of control variables //
if "`controls'" != ""{

*********************************** New variances controls (Modif Felix)	
capture drop part2_pl_switch`=increase_XX'_`i'_XX
gen part2_pl_switch`=increase_XX'_`i'_XX=0

** Modif. Diego: generation of the T_d variable = max_g:D_g,1 = d F_d - 1
capture drop T_d_XX
gegen T_d_XX = max(F_g_XX), by(d_sq_int_XX)
replace T_d_XX = T_d_XX - 1

**********************************************************	

local count_controls=0
// Computing the first differences of the control variables
foreach var of varlist `controls'{

local count_controls=`count_controls'+1

capture drop diff_X_`count_controls'_placebo_`i'_XX

xtset group_XX time_XX

gen diff_X_`count_controls'_placebo_`i'_XX = L`=2*`i''.`var' - L`i'.`var'


*********************************** New variances controls (Modif Felix)
*** Note: I will have to adapt the placebos as well!!!
**# Variance Control Setup

* intermediate steps
*N_g_t * (X_g_t - X_g_t-l)
capture drop diff_X`count_controls'_pl_`i'_N_XX
gen diff_X`count_controls'_pl_`i'_N_XX = N_gt_XX * diff_X_`count_controls'_placebo_`i'_XX // should be correct N_gt, see diff_y_pl below

foreach l of local levels_d_sq_XX { //l is d in the model //

capture drop dummy_XX
capture drop N_d_t_XX
gen dummy_XX=(F_g_XX>time_XX & d_sq_int_XX == `l')
gegen N_d_t_XX=total(dummy_XX) if d_sq_int_XX == `l' // We never use this !!!

* small m
capture drop m`=increase_XX'_pl_g_`count_controls'_`l'_`i'_XX // new
gen m`=increase_XX'_pl_g_`count_controls'_`l'_`i'_XX = (`i' <= T_g_XX-2 & d_sq_int_XX == `l')* (G_XX / N`=increase_XX'_placebo_`i'_XX) * ([dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+1'&time_XX<=T_g_XX) * diff_X`count_controls'_pl_`i'_N_XX)

* capital M
capture drop m_pl`=increase_XX'_`l'_`count_controls'_`i'_XX // new
bys group_XX: gegen m_pl`=increase_XX'_`l'_`count_controls'_`i'_XX=total(m`=increase_XX'_pl_g_`count_controls'_`l'_`i'_XX) // use mean and delete second part

bys group_XX: replace m_pl`=increase_XX'_`l'_`count_controls'_`i'_XX = . if _n != 1

capture drop M_pl`=increase_XX'_`l'_`count_controls'_`i'_XX // new
egen M_pl`=increase_XX'_`l'_`count_controls'_`i'_XX = total(m_pl`=increase_XX'_`l'_`count_controls'_`i'_XX)

* is this the correct bysort or should it be bysort time_XX d_sq_XX???
replace M_pl`=increase_XX'_`l'_`count_controls'_`i'_XX = (1/G_XX)*M_pl`=increase_XX'_`l'_`count_controls'_`i'_XX

// the E_hat defined above -> is there some reason why it should differ for the placebos? I do not think so
* Any need to change prod_X`count_controls'_diff_y_int_XX??

capture drop E_hat_t`count_controls'_`l'_XX	
capture drop E_hat_denom_`count_controls'_`l'_XX
capture drop E_hat_t`count_controls'_`l'_temp_XX
bys time_XX : egen E_hat_denom_`count_controls'_`l'_XX = total(dummy_XX) if d_sq_int_XX == `l'

gen E_hat_t`count_controls'_`l'_temp_XX = (prod_X`count_controls'_diff_y_int_XX * dummy_XX) / E_hat_denom_`count_controls'_`l'_XX
bys d_sq_int_XX time_XX : egen E_hat_t`count_controls'_`l'_XX = total(E_hat_t`count_controls'_`l'_temp_XX)


// define part inside the sum (sum [t=2 -> F_g-1])	
// sum up the count conrol specific parts as those are Vector products
*capture drop in_sum_temp_`count_controls'_XX
*gen in_sum_temp_`count_controls'_XX=prod_X`count_controls'_diff_y_int_XX-(N_d_t_XX>=2)*sqrt((N_d_t_XX)/(N_d_t_XX-1))*E_hat_d_t`count_controls'_XX
capture drop in_sum_temp_pl_`count_controls'_`l'_XX
//gen in_sum_temp_`count_controls'_`l'_XX=prod_X`count_controls'_diff_y_int_XX-(N_d_t_XX>=2)*sqrt((N_d_t_XX)/(N_d_t_XX-1))*E_hat_t`count_controls'_`l'_XX

capture drop N_c_`l'_temp_XX
capture drop N_c_`l'_XX
gen N_c_`l'_temp_XX = d_sq_int_XX == `l' & time_XX >= 2 & time_XX <= T_d_XX & time_XX < F_g_XX
egen N_c_`l'_XX = total(N_c_`l'_temp_XX)
gen in_sum_temp_pl_`count_controls'_`l'_XX = (prod_X`count_controls'_diff_y_int_XX-(E_hat_denom_`count_controls'_`l'_XX>=2)*sqrt((E_hat_denom_`count_controls'_`l'_XX)/(E_hat_denom_`count_controls'_`l'_XX - 1))*E_hat_t`count_controls'_`l'_XX * (time_XX>=2 & time_XX<=F_g_XX-1)) / N_c_`l'_XX
	
// define the inner sum (sum [t=2 -> F_g-1])
capture drop in_sum_pl_`count_controls'_`l'_XX

bys group_XX: gegen in_sum_pl_`count_controls'_`l'_XX = total(in_sum_temp_pl_`count_controls'_`l'_XX) 
}

**********************************************************



foreach l of local levels_d_sq_XX {
	if (scalar(useful_res_`l'_XX)>1){ 

		replace diff_y_pl_`i'_XX = diff_y_pl_`i'_XX - coefs_sq_`l'_XX[`=`count_controls'',1]*diff_X_`count_controls'_placebo_`i'_XX if d_sq_int_XX==`l' 

*********************************** New variances controls (Modif Felix)
capture drop in_brackets_pl_`l'_`count_controls'_XX	
gen in_brackets_pl_`l'_`count_controls'_XX=0
**********************************************************		
		
}		
}

}

}


///// Computing the mean of first or long-differences of outcomes for non-treated and for treated separately - will be useful for the computation of the variance // + //Weighting the means

// Modif Felix: change names because now we will generate two things, the cohorts and the baseline cohorts

capture drop diff_y_pl_`i'_N_gt_XX
gen diff_y_pl_`i'_N_gt_XX=N_gt_XX*diff_y_pl_`i'_XX
capture drop dof_diff_y_pl_`i'_N_gt_XX
gen dof_diff_y_pl_`i'_N_gt_XX=(N_gt_XX!=0&diff_y_pl_`i'_XX!=.)

*** Cohorts never switchers
* DOF
bys d_sq_XX `trends_nonparam' : gegen dof_cohort_pl_`i'_ns_t_XX=total(dof_diff_y_pl_`i'_N_gt_XX) if diff_y_pl_`i'_XX!=.&never_change_d_pl_`i'_XX==1&N`=increase_XX'_t_placebo_`i'_XX>0&N`=increase_XX'_t_placebo_`i'_XX!=.

* Denominator
bys d_sq_XX `trends_nonparam' : gegen count_cohort_pl_`i'_ns_t_XX=total(N_gt_XX) if diff_y_pl_`i'_XX!=.&never_change_d_pl_`i'_XX==1&N`=increase_XX'_t_placebo_`i'_XX>0&N`=increase_XX'_t_placebo_`i'_XX!=.

* Numerator
bys d_sq_XX `trends_nonparam' : gegen total_cohort_pl_`i'_ns_t_XX=total(diff_y_pl_`i'_N_gt_XX) if never_change_d_pl_`i'_XX==1&N`=increase_XX'_t_placebo_`i'_XX>0&N`=increase_XX'_t_placebo_`i'_XX!=.

* Estimator for the expectation (no need for bysort or any conditioning as the Numerator and denominator are generated along the same set of conditions)
gen mean_cohort_pl_`i'_ns_t_XX=total_cohort_pl_`i'_ns_t_XX/count_cohort_pl_`i'_ns_t_XX


*** Cohorts switchers
if "`less_conservative_se'" == ""{
* DOF
bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam' : gegen dof_cohort_pl_`i'_s_t_XX=total(dof_diff_y_pl_`i'_N_gt_XX) if dist_to_switch_pl_`i'_XX==1 

* Denominator
bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam' : gegen count_cohort_pl_`i'_s_t_XX=total(N_gt_XX) if dist_to_switch_pl_`i'_XX==1 

* Numerator
bys d_sq_XX F_g_XX d_fg_XX `trends_nonparam' : gegen total_cohort_pl_`i'_s_t_XX=total(diff_y_pl_`i'_N_gt_XX) if dist_to_switch_pl_`i'_XX==1 

* Estimator for the expectation (no need for bysort or any conditioning as the Numerator and denominator are generated along the same set of conditions)
gen mean_cohort_pl_`i'_s_t_XX=total_cohort_pl_`i'_s_t_XX/count_cohort_pl_`i'_s_t_XX
}


if "`less_conservative_se'" != ""{
** DOF 
* switcher path0 (s0)
bys path_0_XX `trends_nonparam' : gegen dof_cohort_pl_`i'_s0_t_XX=total(dof_diff_y_pl_`i'_N_gt_XX) if dist_to_switch_pl_`i'_XX==1&cohort_fullpath_0_XX==1&cohort_fullpath_1_XX==0
if `i'>1{ // otherwise equivalent to fullpath
* switcher path1 (s1)
bys path_1_XX `trends_nonparam' : gegen dof_cohort_pl_`i'_s1_t_XX=total(dof_diff_y_pl_`i'_N_gt_XX) if dist_to_switch_pl_`i'_XX==1&cohort_fullpath_1_XX==1&cohort_fullpath_`i'_XX==0	
}
* switcher fullpath (s2)
bys path_`i'_XX `trends_nonparam' : gegen dof_cohort_pl_`i'_s2_t_XX=total(dof_diff_y_pl_`i'_N_gt_XX) if dist_to_switch_pl_`i'_XX==1&cohort_fullpath_`i'_XX==1

** Denominator	
* switcher path0 (s0)
bys path_0_XX `trends_nonparam' : gegen count_cohort_pl_`i'_s0_t_XX=total(N_gt_XX) if dist_to_switch_pl_`i'_XX==1&cohort_fullpath_0_XX==1&cohort_fullpath_1_XX==0
if `i'>1{ // otherwise equivalent to fullpath
* switcher path1 (s1)
bys path_1_XX `trends_nonparam' : gegen count_cohort_pl_`i'_s1_t_XX=total(N_gt_XX) if dist_to_switch_pl_`i'_XX==1&cohort_fullpath_1_XX==1&cohort_fullpath_`i'_XX==0	
}
* switcher fullpath (s2)
bys path_`i'_XX `trends_nonparam' : gegen count_cohort_pl_`i'_s2_t_XX=total(N_gt_XX) if dist_to_switch_pl_`i'_XX==1&cohort_fullpath_`i'_XX==1	
	
** Numerator
* switcher path0 (s0)
bys path_0_XX `trends_nonparam' : gegen total_cohort_pl_`i'_s0_t_XX=total(diff_y_pl_`i'_N_gt_XX) if dist_to_switch_pl_`i'_XX==1&cohort_fullpath_0_XX==1&cohort_fullpath_1_XX==0
if `i'>1{ // otherwise equivalent to fullpath
* switcher path1 (s1)
bys path_1_XX `trends_nonparam' : gegen total_cohort_pl_`i'_s1_t_XX=total(diff_y_pl_`i'_N_gt_XX) if dist_to_switch_pl_`i'_XX==1&cohort_fullpath_1_XX==1&cohort_fullpath_`i'_XX==0	
}
* switcher fullpath (s2)
bys path_`i'_XX `trends_nonparam' : gegen total_cohort_pl_`i'_s2_t_XX=total(diff_y_pl_`i'_N_gt_XX) if dist_to_switch_pl_`i'_XX==1&cohort_fullpath_`i'_XX==1

** Estimator for the expectation (no need for bysort or any conditioning as the Numerator and denominator are generated along the same set of conditions)
* switcher path0 (s0)
gen mean_cohort_pl_`i'_s0_t_XX=total_cohort_pl_`i'_s0_t_XX/count_cohort_pl_`i'_s0_t_XX
* switcher path1 (s1)
if `i'>1{
gen mean_cohort_pl_`i'_s1_t_XX=total_cohort_pl_`i'_s1_t_XX/count_cohort_pl_`i'_s1_t_XX
}
* switcher fullüpath (s2)
gen mean_cohort_pl_`i'_s2_t_XX=total_cohort_pl_`i'_s2_t_XX/count_cohort_pl_`i'_s2_t_XX
	
}


///// If the Placebos effect could be estimated (as there are switchers), we compute the U_Gg. //
gen dummy_U_Gg_pl_`i'_XX = (`i'<=T_g_XX-1)


*** This if condition seems to not be satisfied for i>1 if trends_lin is specified !!!!!! (in some specific cases)
if (N`=increase_XX'_placebo_`i'_XX!=0){
	gen U_Gg_pl_`i'_temp_XX = dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * N_gt_XX * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX]*diff_y_pl_`i'_XX 

	bysort group_XX : gegen U_Gg_placebo_`i'_XX=total(U_Gg_pl_`i'_temp_XX)

	replace U_Gg_placebo_`i'_XX= U_Gg_placebo_`i'_XX*first_obs_by_gp_XX

	// Counting the number of groups for which we can estimate U_Gg`i'_temp_XX - to help compute the "N" displayed by the command //
	capture drop count`i'_pl_core_XX

	gen count`i'_pl_core_XX=0
	replace count`i'_pl_core_XX= N_gt_XX if (U_Gg_pl_`i'_temp_XX!=.&U_Gg_pl_`i'_temp_XX!=0|(U_Gg_pl_`i'_temp_XX==0&diff_y_pl_`i'_XX==0&(dist_to_switch_pl_`i'_XX!=0|(N`=increase_XX'_t_placebo_`i'_g_XX!=0&never_change_d_pl_`i'_XX!=0))))


// Modif Felix: Adjusted this part to the new cohorts/cohorts 
///// Computing the "alternative" U_{G,g,l} which will be used for the computation of the variance only - these are like the above U_{G,g,l}s, except that the outcome differences are demeaned, and there is a DOF adjustment when possible//

// Initializing the U_{G,g,l}_var at 0 //
gen U_Gg_pl_`i'_temp_var_XX =0
	
*****
	
// For controls
replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * diff_y_pl_`i'_N_gt_XX  if never_change_d_pl_`i'_XX==1&dof_cohort_pl_`i'_ns_t_XX==1

* For all t < Fg -> only need this case for controls
replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * N_gt_XX * (diff_y_pl_`i'_XX - (never_change_d_pl_`i'_XX * sqrt(dof_cohort_pl_`i'_ns_t_XX/(dof_cohort_pl_`i'_ns_t_XX-1)) * mean_cohort_pl_`i'_ns_t_XX)) if never_change_d_pl_`i'_XX==1&dof_cohort_pl_`i'_ns_t_XX>1&dof_cohort_pl_`i'_ns_t_XX!=.


// For switchers
if "`less_conservative_se'" == ""{
replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * diff_y_pl_`i'_N_gt_XX if dist_to_switch_pl_`i'_XX==1&dof_cohort_pl_`i'_s_t_XX==1 

* For t = Fg − 1 + ℓ -> only need this case for switchers
replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * N_gt_XX * (diff_y_pl_`i'_XX - (dist_to_switch_pl_`i'_XX * sqrt(dof_cohort_pl_`i'_s_t_XX/(dof_cohort_pl_`i'_s_t_XX-1)) * mean_cohort_pl_`i'_s_t_XX)) if dist_to_switch_pl_`i'_XX==1&dof_cohort_pl_`i'_s_t_XX>1&dof_cohort_pl_`i'_s_t_XX!=.		
}

if "`less_conservative_se'" != ""{
* Case where we only have one group in a pathway, even with the least strict definition -> no demeaning
if `i'==1{
replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * diff_y_pl_`i'_N_gt_XX if dist_to_switch_pl_`i'_XX==1&(dof_cohort_pl_`i'_s0_t_XX==.&dof_cohort_pl_`i'_s2_t_XX==.)
}

if `i'>1{
replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * diff_y_pl_`i'_N_gt_XX if dist_to_switch_pl_`i'_XX==1&(dof_cohort_pl_`i'_s0_t_XX==.&dof_cohort_pl_`i'_s1_t_XX==.&dof_cohort_pl_`i'_s2_t_XX==.)
}
// should I do conditions with cohort_fullpath or with dof_cohort missing/non-missing? -> If everything in the prior code is correct then all the dof should just be defined if there are at least 2 groups in each pathway given the corresponding definition

* Case where we have multiple groups per pathways when conditioning on D_{g,1},F_g
replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * N_gt_XX * (diff_y_pl_`i'_XX - (dist_to_switch_pl_`i'_XX * sqrt(dof_cohort_pl_`i'_s0_t_XX/(dof_cohort_pl_`i'_s0_t_XX-1)) * mean_cohort_pl_`i'_s0_t_XX)) if dist_to_switch_pl_`i'_XX==1&dof_cohort_pl_`i'_s0_t_XX>1&dof_cohort_pl_`i'_s0_t_XX!=.

* Case where we have multiple groups per pathways when conditioning on D_{g,1},F_g, D_{g,F_g}
if `i'>1{
replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * N_gt_XX * (diff_y_pl_`i'_XX - (dist_to_switch_pl_`i'_XX * sqrt(dof_cohort_pl_`i'_s1_t_XX/(dof_cohort_pl_`i'_s1_t_XX-1)) * mean_cohort_pl_`i'_s1_t_XX)) if dist_to_switch_pl_`i'_XX==1&dof_cohort_pl_`i'_s1_t_XX>1&dof_cohort_pl_`i'_s1_t_XX!=.
}

* Case where we have multiple groups per pathways only when conditioning on D_{g,1},F_g, D_{g,F_g},...,D_{g,F_g-1+\ell}
replace U_Gg_pl_`i'_temp_var_XX= dummy_U_Gg_pl_`i'_XX*(G_XX / N`=increase_XX'_placebo_`i'_XX) * [dist_to_switch_pl_`i'_XX - (N`=increase_XX'_t_placebo_`i'_g_XX/N_gt_control_placebo_`i'_XX) * never_change_d_pl_`i'_XX] * (time_XX>=`=`i'+2'&time_XX<=T_g_XX) * N_gt_XX * (diff_y_pl_`i'_XX - (dist_to_switch_pl_`i'_XX * sqrt(dof_cohort_pl_`i'_s2_t_XX/(dof_cohort_pl_`i'_s2_t_XX-1)) * mean_cohort_pl_`i'_s2_t_XX)) if dist_to_switch_pl_`i'_XX==1&dof_cohort_pl_`i'_s2_t_XX>1&dof_cohort_pl_`i'_s2_t_XX!=.
}	
	
	
****************************************************************************


// Modif Felix: Add the additional part with controls 
if "`controls'"!=""{

levelsof d_sq_int_XX, local(levels_d_sq_XX)
foreach l of local levels_d_sq_XX {	
	
capture drop combined_pl`=increase_XX'_temp_`l'_`i'_XX	
gen combined_pl`=increase_XX'_temp_`l'_`i'_XX=0
		
// loops to do the summations
forvalues j=1/`count_controls'{
forvalues k=1/`count_controls'{		
	
// define the part in brackets -> this is the first matrix sum where we sum the inv_denom rows with the kx1 coulmn vectors!
capture drop in_brackets_pl_`l'_`j'_`k'_temp_XX
gen in_brackets_pl_`l'_`j'_`k'_temp_XX = inv_Denom_`l'_XX[`j',`k'] * in_sum_pl_`k'_`l'_XX * (d_sq_int_XX == `l' & F_g_XX>=3) // G does not show up because it is part of inv_Denom_ and N does not show up because it is already included in the in_sum part (we use mean and not total there)
// step where we have the first sum over the k
*capture drop in_brackets_`l'_`j'_XX
replace in_brackets_pl_`l'_`j'_XX=in_brackets_pl_`l'_`j'_XX+in_brackets_pl_`l'_`j'_`k'_temp_XX

} // end loop over k

replace in_brackets_pl_`l'_`j'_XX=in_brackets_pl_`l'_`j'_XX - coefs_sq_`l'_XX[`j',1]
// Now in this part put the summing with M and the rest which is indexed by j
capture drop combined_pl`=increase_XX'_temp_`l'_`j'_`i'_XX // new
gen combined_pl`=increase_XX'_temp_`l'_`j'_`i'_XX=M_pl`=increase_XX'_`l'_`j'_`i'_XX*in_brackets_pl_`l'_`j'_XX

// step where we have the second sum over the j
replace combined_pl`=increase_XX'_temp_`l'_`i'_XX=combined_pl`=increase_XX'_temp_`l'_`i'_XX+combined_pl`=increase_XX'_temp_`l'_`j'_`i'_XX

} // end loop over j

// Final sum over the status quo treatment (outer sum in the formula)
replace part2_pl_switch`=increase_XX'_`i'_XX=part2_pl_switch`=increase_XX'_`i'_XX+combined_pl`=increase_XX'_temp_`l'_`i'_XX if d_sq_int_XX==`l' // does the if condition here even matter because those are seperate variables depending on l

} // loop over l

// Making the adjustement
if `=increase_XX'==1{
replace U_Gg_pl_`i'_temp_var_XX=U_Gg_pl_`i'_temp_var_XX - part2_pl_switch1_`i'_XX 
}

if `=increase_XX'==0{
replace U_Gg_pl_`i'_temp_var_XX=U_Gg_pl_`i'_temp_var_XX + part2_pl_switch0_`i'_XX 
}
	
}	

// Summing the U_{G,g,l}s over time periods for each group //
bys group_XX: gegen U_Gg_pl_`i'_var_XX=total(U_Gg_pl_`i'_temp_var_XX)

}
****************************************************************************

	//End of placebos computation 



///// Normalization of the placebos//

if "`normalized'"!=""{
	
	capture drop sum_treat_until_`i'_pl_XX
	capture drop delta_D_pl_`i'_cum_temp_XX

	if "`continuous'"==""{
	bys group_XX : gegen sum_treat_until_`i'_pl_XX = total(treatment_XX - d_sq_XX) if time_XX>=F_g_XX&time_XX<=F_g_XX-1+`i'&diff_y_pl_`i'_XX!=.&S_g_XX==increase_XX
	}
		else if "`continuous'"!=""{
		bys group_XX : gegen sum_treat_until_`i'_pl_XX = total(treatment_XX_orig - d_sq_XX_orig) if time_XX>=F_g_XX&time_XX<=F_g_XX-1+`i'&S_g_XX==increase_XX&diff_y_pl_`i'_XX!=.
	}

gen delta_D_pl_`i'_cum_temp_XX = N_gt_XX/N`=increase_XX'_placebo_`i'_XX*[sum_treat_until_`i'_pl_XX* S_g_XX + (1-S_g_XX)*(-sum_treat_until_`i'_pl_XX)] if dist_to_switch_pl_`i'_XX==1

	sum  delta_D_pl_`i'_cum_temp_XX
	scalar delta_norm_pl_`i'_XX = r(sum)
}



} //End of loop on placebos
}

// Modif Felix: adapted the part analogous to the one above (without placebos)

scalar Ntrendslin_pl=1
forvalue i=1/`=l_placebo_u_a_XX' {
scalar Ntrendslin_pl=min(Ntrendslin_pl , N`=increase_XX'_placebo_`i'_XX )
}

*** Modif Felix: Solve Issue with U_Gg which are created but due to the condition not in the sum and then get called by the main program again (thats why we had placebo 3 even if placebo 2 was undefined). Idea for now, verify this agian!!! -> Drop all U_Gg after the effect we can not use anymore 
if "`trends_lin'"!="" & Ntrendslin_pl == 0 { 
	capture drop U_Gg_placebo_`=l_placebo_u_a_XX'_XX
	scalar N`=increase_XX'_placebo_`=l_placebo_u_a_XX'_XX=0
}	

if "`trends_lin'"!="" & Ntrendslin_pl != 0 {
	
	capture drop U_Gg_pl_`=l_placebo_u_a_XX'_TL
	capture drop U_Gg_pl_`=l_placebo_u_a_XX'_var_TL
	
	gen U_Gg_pl_`=l_placebo_u_a_XX'_TL = 0
	gen U_Gg_pl_`=l_placebo_u_a_XX'_var_TL = 0

	// Add some condition that l_u_a_XX>1?
	
	forvalue i=1/`=l_placebo_u_a_XX'{
		*if (N`=increase_XX'_placebo_`i'_XX!=0){ // add this for now but discuss again how to deal with the problem of having no observations and therefore U_Gg_placebo_`i' undefined so we get an error here when trying to add it -> maybe check for it some time earlier in the code!
		replace U_Gg_pl_`=l_placebo_u_a_XX'_TL = U_Gg_pl_`=l_placebo_u_a_XX'_TL + U_Gg_placebo_`i'_XX
		replace U_Gg_pl_`=l_placebo_u_a_XX'_var_TL = U_Gg_pl_`=l_placebo_u_a_XX'_var_TL + U_Gg_pl_`i'_var_XX
	*}
	}
	replace U_Gg_placebo_`=l_placebo_u_a_XX'_XX = U_Gg_pl_`=l_placebo_u_a_XX'_TL
		replace U_Gg_pl_`=l_placebo_u_a_XX'_var_XX =U_Gg_pl_`=l_placebo_u_a_XX'_var_TL
		
		// For Felix : as above for the effects, the U_Gg_vars have to be modified.

	}
	

}
//End of condition checking if we can compute at least one placebo

*** For the estimation of \hat{\delta} ***

///// Computing the sum of the N1_`i'_XX for the weights w. //

if "`trends_lin'"==""{

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

// Modif Felix: Adjust for the continuous option
if "`continuous'"==""{
gen delta_D_`i'_temp_XX = N_gt_XX/N`=increase_XX'_`i'_XX*[(treatment_XX-d_sq_XX)* S_g_XX + (1-S_g_XX)*(d_sq_XX-treatment_XX)] if distance_to_switch_`i'_XX==1
}
else if "`continuous'"!=""{
gen delta_D_`i'_temp_XX = N_gt_XX/N`=increase_XX'_`i'_XX*[(treatment_XX_orig-d_sq_XX_orig)* S_g_XX + (1-S_g_XX)*(d_sq_XX_orig-treatment_XX_orig)] if distance_to_switch_`i'_XX==1	
}	

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


// For Felix : here is where we compute the U_Gg_XX for the average effect. When they will be called in the main program, they will be switchers-in/switchers-out specific. Can we avoid this, and directly aggregate the U_Gg_global_XX of the main prohram to have the avg effect ? Do we still need the w_`i'_XX though ?

}
}
 // end of the quietly condition
	
end


