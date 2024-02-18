//ssc install did_multiplegt_dyn


qui bcuse wagepan, clear
sort nr year
gen clus = mod(nr, 3)
gen over_gr = mod(nr, 2)

did_multiplegt_dyn lwage nr year union, effects(5) graph_off placebo(2) controls(married exper) predict_het(clus, ) switchers(in)

did_multiplegt_dyn lwage nr year union, effects(5) graph_off placebo(2) trends_nonparam(over_gr) weight(hours) dont_drop_larger_lower

did_multiplegt_dyn lwage nr year union, effects(5) graph_off placebo(2) normalized drop_if_d_miss_before_first_switch normalized_weights cluster(clus)

did_multiplegt_dyn lwage nr year union, effects(5) graph_off placebo(2) continuous(pol, 1) controls(hours) trends_lin

did_multiplegt_dyn lwage nr year union, effects(5) effects_equal placebo(1) graph_off design(1, console) date_first_switch(by_baseline_treat, console) same_switchers same_switchers_pl