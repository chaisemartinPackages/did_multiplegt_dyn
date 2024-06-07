clear
set seed 123
set obs 1500
gen G = mod(_n-1, 250) + 1
bys G: gen T = _n
sum T
gen D = uniform() > 0.5 & T >= r(max) - 1
gen X = uniform()
gen Y = uniform() * (1 + D) + 2*X

qui do "C:/Users/39380/C DE CHAISEMARTIN Dropbox/RAs De Chaisemartin/RAs Really Credible DID-TWFE/GitHub repo - local/did_multiplegt_dyn/Stata/did_multiplegt_dyn.ado"
est clear
did_multiplegt_dyn Y G T D, effects(1) placebo(1) graph_off
eststo

qui do "C:/Users/39380/C DE CHAISEMARTIN Dropbox/RAs De Chaisemartin/RAs Really Credible DID-TWFE/GitHub repo - local/did_multiplegt_dyn/Stata/mp_did_multiplegt_dyn.ado"
mp_did_multiplegt_dyn Y G T D, effects(2) placebo(4) controls(X)

eststo
esttab, order(Effect_* Placebo_* Av_*) se