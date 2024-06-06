clear
set seed 123
set obs 1000
gen G = mod(_n-1, 250) + 1
bys G: gen T = _n
gen D = uniform() > 0.5 & T == 4
gen Y = uniform() * (1 + D)

est clear
did_multiplegt_dyn Y G T D, effects(1) placebo(1) graph_off
eststo

qui do "mp_did_multiplegt_dyn.ado"
mp_did_multiplegt_dyn Y G T D, effects(1) placebo(2)
eststo
esttab, order(Effect_* Placebo_* Av_*)