* Check Monte Carlo results
use "monte_carlo_results.dta", clear
describe
sum diff_eff1 diff_eff2 diff_eff3 diff_plac1 diff_plac2
sum diff_se_eff1 diff_se_eff2 diff_se_eff3 diff_se_plac1 diff_se_plac2
list diff_eff1 diff_se_eff1 in 1/20
