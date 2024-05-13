clear
global dir0 : pwd
global dir1 = subinstr("$dir0", "\vignettes\src", "", .)

qui do "$dir1/did_multiplegt_dyn_13_05_24.do"
clear
set seed 123
scalar TT = 20
scalar GG = 1000
set obs `= TT * GG'
gen G = mod(_n-1,GG) + 1
gen T = floor((_n-1)/GG)
sort G T

gen D = uniform() > 0.5
gen X = uniform() * T
gen Y = uniform() * (1 + D + X)

est clear
did_multiplegt_dyn Y G T D, effects(5) graph_off effects_equal
estadd scalar p_joint = e(p_equality_effects) 
estadd local controls = "No"
est sto model_1

did_multiplegt_dyn Y G T D, effects(5) placebo(3) effects_equal controls(X) graph_off
estadd scalar p_joint = e(p_equality_effects) 
estadd scalar p_placebo = e(p_jointplacebo)
estadd local controls = "Yes"
est sto model_2

// Basic usage
esttab model_* using "$dir1/vignettes/assets/reg1.tex", replace booktabs se s(p_joint p_placebo controls) standalone

// Refined format
esttab model_* using "$dir1/vignettes/assets/reg1a.tex", replace booktabs se s(control p_joint p_placebo, label("Controls" "Joint Eq. Effects" "Joint Sig. Placebo"))  b(%9.5fc) coeflabels(Effect_1 "$\hat{\delta}_1$" Effect_2 "$\hat{\delta}_2$" Effect_3 "$\hat{\delta}_3$" Effect_4 "$\hat{\delta}_4$" Effect_5 "$\hat{\delta}_5$" Avg_Tot_Effect "$\hat{\delta}$" Placebo_1 "$\hat{\delta}_1^{pl}$" Placebo_2 "$\hat{\delta}_2^{pl}$" Placebo_3 "$\hat{\delta}_3^{pl}$") substitute(\_ _) mlabels(,none) collabels(,none) standalone

est clear

gen H1 = G/GG
gen H2 = mod(G, TT)

// Predict Het
did_multiplegt_dyn Y G T D, predict_het(H1 H2, all) graph_off effects(3)
est sto model_1

esttab model_* using "$dir1/vignettes/assets/reg2.tex", replace booktabs se noobs standalone