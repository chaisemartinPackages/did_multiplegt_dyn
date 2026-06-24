***** File to run standardized tests for "did_multiplegt_dyn" before updating versions*****

						*** Version: 22.05.25 ***


							*** Preamble ***
********************************************************************************
// Set the directory to load the datasets and to store the resulting log files




//David
* data path
global data_path="C:\Users\dearb\OneDrive - Universidad de los Andes\Really Credible\debug did_multiplegt_dyn\did_multiplegt_dyn_test\testing_data"

* save path
global save_path="C:\Users\dearb\OneDrive - Universidad de los Andes\Really Credible\debug did_multiplegt_dyn\did_multiplegt_dyn_test\test_dyn\log"

global gitupd = "C:\Users\dearb\OneDrive - Universidad de los Andes\Really Credible\debug did_multiplegt_dyn\did_multiplegt_dyn_test\GitHub_updates"
// Felix




// Felix
* data path
*global data_path="C:\Users\Felix\Desktop\DID_programs\did_multiplegt_dyn\testing_data"

* save path
*global save_path="C:\Users\Felix\Desktop\DID_programs\did_multiplegt_dyn\test_dyn\log"

// Diego
* data path

* save path


// Romain
* data path

* save path


// Clément
* data path

* save path


* Specify the used dataset 
local dataname = "wagepan"
local dataset = "`dataname'.dta"
* Specify version to be saved with the log file 
local version = "03_10_2025"

********************************************************************************

*** Initialize the version of the command one wants to test 
do "$gitupd/`version'\did_multiplegt_dyn.ado"


**# Bookmark #1
*** Load data and initialize all variables for the corresponding dataset ***
if "`dataset'"=="wagepan.dta"{
	
	use "$data_path/`dataset'",clear
	
	* Base variables
	local y = "lwage"
	local g = "nr"
	local t = "year"
	local d = "union"
	
	* Number of effects/placebos
	local n_eff = 5
	local n_pl = 2
	
	* Bootstrap options
	local b_reps = 20
	local b_seed = 1234
	
	* Controls
	local cont = "hours"
	
	* Nonparametric trends by
	local nonparam = "black"
	
	* Weighting variable 
	local wght = "educ"
	
	* Cluster variables
	local clust = "hisp"
	
	* by variable
	local by_var = "black"
	
	* predict_het variable 
	local het_var = "black"
	
	* effects to test in effects-equal
	local eff_eq="all"
}

if "`dataset'"=="favara_imbs.dta"{
	
	use "$data_path/`dataset'",clear
	
	* Base variables
	local y = ""
	local g = ""
	local t = ""
	local d = ""
	
	* Number of effects/placebos
	local n_eff = 5
	local n_pl = 2
	
	* Bootstrap options
	local b_reps = 20
	local b_seed = 1234
	
	* Controls
	local cont = ""
	
	* Nonparametric trends by
	local nonparam = ""
	
	* Weighting variable 
	local wght = ""
	
	* Cluster variables
	local clust = ""
}

if "`dataset'"=="gentzkowetal_didtextbook.dta"{
	
	use "$data_path/`dataset'",clear
	
	* Base variables
	local y = ""
	local g = ""
	local t = ""
	local d = ""

	* Number of effects/placebos
	local n_eff = 5
	local n_pl = 2
	
	* Bootstrap options
	local b_reps = 20
	local b_seed = 1234	
	
	* Controls
	local cont = ""
	
	* Nonparametric trends by
	local nonparam = ""
	
	* Weighting variable 
	local wght = ""
	
	* Cluster variables
	local clust = ""
}

********************************************************************************

							*** Run the actual tests ***

**# Bookmark #2
** Options without additional outputs
{
cap matrix drop store_results
cap matrix drop A
mat A=.

* Start saving log file 
log using "$save_path/log_main_`version'_`dataname'", replace

// Base program, no placebos or other optios (only graph_off to save time)
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') 
mat rownames A="Baseline"
mat store_results=A
mat store_results=A\e(estimates)[1..`n_eff'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

*****************************************************

// Placebos 
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl')
mat rownames A="Placebos"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// Normalized 
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') normalized
mat rownames A="Normalized"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// Controls 
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') controls(`cont')
mat rownames A="Controls"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// Trends_nonparam 
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') trends_nonparam(`nonparam')
mat rownames A="Trends_Nonparam"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// trends_lin
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') trends_lin
mat rownames A="Trends_Lin"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// Continuous (maybe exclude for now)
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') continuous(1)
mat rownames A="Continuous"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// weight
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') weight(`wght')
mat rownames A="Weight"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// cluster
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') cluster(`clust')
mat rownames A="Cluster"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// same_switchers
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') same_switchers
mat rownames A="Same_Switchers"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// same_switchers_placebo
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') same_switchers same_switchers_pl
mat rownames A="Same_SwitcherS_Placebo"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// switchers in 
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') switchers(in)
mat rownames A="Switchers_In"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// switchers out
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') switchers(out)
mat rownames A="Switchers_Out"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// only never switchers
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') only_never_switchers 
mat rownames A="Only_Never_Switchers"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// CI level
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') ci_level(90)
mat rownames A="CI_Level_90"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int


did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') ci_level(99)
mat rownames A="CI_Level_99"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// less conservative se's
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') less_conservative_se
mat rownames A="Less_Conservative_SE"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// bootstrap
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') bootstrap(`b_reps',`b_seed')
mat rownames A="Bootstrap"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// dont_drop_larger_lower
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') dont_drop_larger_lower
mat rownames A="Dont_Drop_Larger_Lower"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

*****************************************************

// effects-equal
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') effects_equal(`eff_eq')
mat rownames A="Effects_Equal"
mat store_results=store_results\A
mat store_results=store_results\e(estimates)[1..`n_eff'+`n_pl'+1,1]
mat store_results=store_results\e(variances)[1..`n_eff'+`n_pl'+1,1]

forvalues j=1/`n_eff'{
	matrix xx_int=e(N_effect_`j')
	matrix rownames xx_int = "N_effect_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_effect_`j')
	matrix rownames xx_int = "N_switchers_effect_`j'" 
	mat store_results=store_results\xx_int
}
// Do we not report p-value from joint test of significance effects as ereturn???

matrix xx_int=e(N_avg_total_effect)
matrix rownames xx_int = "N_avg_total_effect" 
mat store_results=store_results\xx_int

matrix xx_int=e(N_switchers_effect_average)
matrix rownames xx_int = "N_switchers_effect_average" 
mat store_results=store_results\xx_int

forvalues j=1/`n_pl'{
	
	matrix xx_int=e(N_placebo_`j')
	matrix rownames xx_int = "N_placebo_`j'" 
	mat store_results=store_results\xx_int
	
	matrix xx_int=e(N_switchers_placebo_`j')
	matrix rownames xx_int = "N_switchers_placebo_`j'" 	
	mat store_results=store_results\xx_int
}

matrix xx_int=e(p_jointeffects)
matrix rownames xx_int = "p_jointeffects" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_jointplacebo)
matrix rownames xx_int = "p_jointplacebo" 
mat store_results=store_results\xx_int

matrix xx_int=e(p_equality_effects)
matrix rownames xx_int = "p_equality_effects" 
mat store_results=store_results\xx_int

*****************************************************

* End main log 
log close
translate "$save_path/log_main_`version'_`dataname'.smcl" "$save_path/log_main_`version'_`dataname'.pdf" 
}

*** Translate results vector into dataset 
svmat2 store_results, rnames(Effect)
rename store_results results_`version'_`dataname'
preserve 
	keep Effect results_`version'_`dataname'
	order Effect results_`version'_`dataname'
	local dim (`= rowsof(store_results)')
	drop if _n>`dim'
	save "$save_path/datasets/data_`version'_`dataname'.dta" // FELIX: DELETE THE REPLACE AT THE END TO NOT OVERWRITE DATASET BY ACCIDENT This is just for testing purpose
restore 

// Working on storing the results, exclude the second part for now
/*
**# Bookmark #3
** Options with additional outputs
{
	
* Start saving log file 
log using "$save_path/log_add_`version'_`dataname'", replace	
	
// design
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') design(,console)	

// by 
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') by(`by_var')

// by_path
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') by_path(4)
// Error in this option!!!!!! -> also does not seem to match with results from design option	
	
// predict_het
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') predict_het(`het_var',all)

// date_first_switch
did_multiplegt_dyn `y' `g' `t' `d', graph_off effects(`n_eff') placebo(`n_pl') date_first_switch(by_baseline_treat,console)
	
* End add log 
log close
translate "$save_path/log_add_`version'_`dataname'.smcl" "$save_path/log_add_`version'_`dataname'.pdf" 	
}
*/

** What about other options with no direct outputs (graph options, save sample, etc)?












