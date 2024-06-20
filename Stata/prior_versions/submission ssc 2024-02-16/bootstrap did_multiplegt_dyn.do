//Install command and load example data
net get did_multiplegt_dyn
use favara_imbs_did_multiplegt_dyn.dta, clear

*** Wagepan to test for now
use "C:\Users\fe-kn\C DE CHAISEMARTIN Dropbox\RAs De Chaisemartin\RAs Really Credible DID-TWFE\did_multiplegt\test_dofiles_and_datasets\WAGEPAN.DTA", clear

// do "did_multiplegt_dyn.do"
// Initialize number of bootsrap repetitions (adapt to your desired number)
local b_reps=10
local n_effects=5
local n_placebos=2

//qui levelsof county
qui levelsof nr
local lvl_g=r(r)
di `lvl_g'

//Store the coefficients to bootstrap
	local coefs " "
	forvalues i=1/`n_effects'{
		local coefs = "`coefs' e(Effect_`i')"
	}
	forvalues i=1/`n_placebos'{
		local coefs = "`coefs' e(Placebo_`i')"
	}
//Run the bootstrap replications (adapt to your specification)
	set seed 1234
	//bootstrap "`coefs'", reps(`b_reps') cluster(county) size(`lvl_g'): did_multiplegt_dyn Dl_vloans_b county year inter_bra, effects(`n_effects') placebo(`n_placebos') continuous(1) graph_off
	bootstrap "`coefs'", reps(`b_reps') cluster(nr) size(`lvl_g'): did_multiplegt_dyn lwage nr year union, effects(`n_effects') placebo(`n_placebos') continuous(1) graph_off //by(black)

// Note: The treatment in this example dataset is not continuous, therefore you get a warning message when running this example.