library(DIDmultiplegtDYN)
library(wooldridge)

data("wagepan")
wagepan$over_gr <- as.numeric(wagepan$nr %% 2)
wagepan$clus <- as.numeric(wagepan$nr %% 3)

# Some brute force consistency checks for the package #
summary(did_multiplegt_dyn(df = wagepan, Y = "lwage", G = "nr", T = "year", D = "union", effects = 5, placebo = 2, graph_off = TRUE, controls = c("married", "exper"), predict_het = list("clus",-1), switchers = "in"))

summary(did_multiplegt_dyn(df = wagepan, Y = "lwage", G = "nr", T = "year", D = "union", effects = 5, placebo = 2, graph_off = TRUE, trends_nonparam = "over_gr", weight = "hours", dont_drop_larger_lower = TRUE))

summary(did_multiplegt_dyn(df = wagepan, Y = "lwage", G = "nr", T = "year", D = "union", effects = 5, placebo = 2, graph_off = TRUE, normalized = TRUE, normalized_weights = TRUE, drop_if_d_miss_before_first_switch = TRUE, cluster = "clus"))

summary(did_multiplegt_dyn(df = wagepan, Y = "lwage", G = "nr", T = "year", D = "union", effects = 5, placebo = 2, graph_off = TRUE, continuous = 1, trends_lin = TRUE, controls = "hours"))

summary(did_multiplegt_dyn(df = wagepan, Y = "lwage", G = "nr", T = "year", D = "union", effects = 5, graph_off = TRUE, placebo = 2, design = list(1, "console"), date_first_switch = list("by_baseline_treat", "console"), same_switchers = TRUE, same_switchers_pl = TRUE, effects_equal = TRUE))