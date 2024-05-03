// A toy example //
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
ytitle("dynamic effects") xtitle("# periods after first switch")
gr export "assets/vignette_1_Stata_fig1.jpg", replace

