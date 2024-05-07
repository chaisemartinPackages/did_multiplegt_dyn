// A toy example //
/*
clear
set obs 9
gen Party = char(64 + ceil(_n/3))
egen Partyid = group(Party)
bys Party: gen Year = 2003 + (_n-1)*2
mat define Share = (0.4\0.35\0.25\0.4\0.45\0.50\0.15\0.1\0.05)
mat define LeadChangeYr = (2004\2004\2004\.\.\.\2005\2005\2005)
foreach v in Share LeadChangeYr {
    mat coln `v' = "`v'"
    svmat `v', n(col)
}
gen Treatment = Year >= LeadChangeYr

mat define res = J(4, 1, .)
did_multiplegt_dyn Share Partyid Year Treatment if inlist(Party, "A", "B"), effects(2) graph_off
mat res[2,1] = e(Effect_1)
mat res[4,1] = e(Effect_2)
did_multiplegt_dyn Share Partyid Year Treatment if inlist(Party, "C", "B"), effects(2) graph_off
mat res[1,1] = e(Effect_1)
mat res[3,1] = e(Effect_2)
svmat res

gen ell = _n
line res1 ell in 1/4, lc(black) || ///
scatter res1 ell in 1/4 if mod(_n, 2) == 0, mcolor(red) msize(1.5)  || ///
scatter res1 ell in 1/4 if mod(_n, 2) == 1, mcolor(blue) msize(1.5) || ///
, legend(order(2 "Party A's dynamic effects" 3 "Party B's dynamic effects") pos(6) col(2)) ///
ytitle("dynamic effects") xtitle("# periods after first switch") ylabel(-0.3(0.1)0.1, ) yline(0)
gr export "assets/vignette_1_Stata_fig1.jpg", replace
*/


// General case //
/*
clear
set seed 123
scalar TT = 20
scalar GG = 5
set obs `= TT * GG'
gen G = mod(_n-1,GG) + 1
gen T = floor((_n-1)/GG)
sort G T

gen D = 0
forv j=2/5 {
    replace D = 1 if G == `j' & T == `j' + 2
}
bys G: gen D_stag = sum(D)
gen Y = uniform() * (1 + 100*D_stag) if mod(T,4) == 0
drop D_stag*

browse

did_multiplegt_dyn Y G T D, effects(5) graph_off

bys G: gen D0 = D[1]
gen D_change = abs(D - D0) != 0
bys G: gen at_least_one_D_change = sum(D_change)

bys G: egen never_treated = max(at_least_one_D_change)
replace never_treated = 1 - never_treated
bys G: egen F_g_temp = min(T * D_change) if D_change != 0
bys G: egen F_g = mean(F_g_temp)
sum T
replace F_g = r(max) + 1 if missing(F_g)

gen subsample = mod(F_g, 4) + 1
replace subsample = 0 if at_least_one_D_change == 0

sort G T
browse

keep if !missing(Y)

// Naive run
did_multiplegt_dyn Y G T at_least_one_D_change, effects(2)

scalar effects = 2
mat define res = J(4 * effects, 6, .)
forv j=1/4 {
    did_multiplegt_dyn Y G T at_least_one_D_change ///
     if inlist(subsample, 0, `j'), effects(`=effects') graph_off
    forv i = 1/`=effects' {
        mat adj = mat_res_XX[`i',1..6]
        forv c =1/6 {
            mat res[`j'+(`i'-1)*4,`c'] = adj[1, `c']
        }
    }
}
mat li res
*/

// Bigger case for graph //

clear
set seed 123
scalar TT = 20
scalar GG = 1000
set obs `= TT * GG'
gen G = mod(_n-1,GG) + 1
gen T = floor((_n-1)/GG)
sort G T

gen D = 0
replace D = 1 if T == mod(G - 1, 5) + 3 & mod(G-1, 5) != 0
bys G: gen D_stag_temp = sum(D)
bys G: gen D_stag = sum(D_stag_temp)
gen Y = uniform() * (1 + D_stag)
replace Y = . if mod(T, 4) != 0
drop D_stag*

bys G: gen D0 = D[1]
gen D_change = abs(D - D0) != 0
bys G: gen at_least_one_D_change = sum(D_change)

bys G: egen never_treated = max(at_least_one_D_change)
replace never_treated = 1 - never_treated
bys G: egen F_g_temp = min(T * D_change) if D_change != 0
bys G: egen F_g = mean(F_g_temp)
sum T
replace F_g = r(max) + 1 if missing(F_g)

gen subsample = (4 - mod(F_g, 4)) * (mod(F_g, 4) != 0) + 1
replace subsample = 0 if at_least_one_D_change == 0

sort G T

keep if !missing(Y)

scalar effects = 2
mat define res = J(4 * effects, 6, .)
forv j=1/4 {
    did_multiplegt_dyn Y G T at_least_one_D_change ///
     if inlist(subsample, 0, `j'), effects(`=effects') graph_off
    forv i = 1/`=effects' {
        mat adj = mat_res_XX[`i',1..6]
        forv c =1/6 {
            mat res[`j'+(`i'-1)*4,`c'] = adj[1, `c']
        }
    }
}
mat li res

did_multiplegt_dyn Y G T at_least_one_D_change, graph_off

// Add the 0-period to the saved matrix
mat res = (0,0,0,0,0,0) \ res
svmat res
gen rel_time = _n-1 if !missing(res1)
tw rcap res3 res4 rel_time, lc(blue) || ///
connected res1 rel_time, mc(blue) lc(blue) || ///
, xtitle("# periods after first switch") yline(0) ytitle("dynamic effects") leg(off)
gr export "assets/vignette_1_Stata_fig2.jpg", replace
