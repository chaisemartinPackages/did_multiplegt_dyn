rm(list=ls())

## Dependencies - did_multiplegt_dyn
library(dplyr)
library(ggplot2)
library(matlib, include.only = c("Ginv"))
library(xlsx, include.only = c("write.xlsx"))
library(plm, include.only = c("pdata.frame", "make.pbalanced"))
library(data.table, include.only = c("shift", "setnames"))

library(wooldridge)
library(devtools)

dir0 <- "C:/Users/39380/C DE CHAISEMARTIN Dropbox/RAs De Chaisemartin/RAs Really Credible DID-TWFE/did_multiplegt_dyn/_R_did_multiplegt_dyn"
devtools::load_all(paste0(dir0,"/DIDmultiplegtDYN"))

data("wagepan")
wagepan$over_gr <- as.numeric(wagepan$nr %% 3)
did_multiplegt_dyn <- did_multiplegt_dyn(
    df = wagepan, 
    Y = "lwage", 
    G = "nr", 
    T = "year", 
    D = "union", 
    effects = 5, 
    placebo = 2, 
    effects_equal = TRUE, 
    save_sample = TRUE, 
    graph_off = TRUE, 
    normalized = TRUE, 
    normalized_weights = "by_k", 
    design = c(0.5, paste0(dir0, "/test_des.xlsx")), 
    date_first_switch = c("by_baseline_treat",paste0(dir0, "/test.xlsx")),
    by = "over_gr",
    same_switchers = TRUE
)
print(did_multiplegt_dyn)


#BUG TO SOLVE - CASE WITH NO OUT-SWITCHERS - SOLVED

#wagepan <- data.frame(matrix(0, nrow = 1000, ncol = 4))
#names(wagepan) <- c("Y", "D", "G", "T")
#wagepan$T <- sapply(1:dim(wagepan)[1], function(x) x %% 4 + 1)
#wagepan <- wagepan[order(wagepan$T), ]
#wagepan$G <- sapply(1:dim(wagepan)[1], function(x) x %% (dim(wagepan)[1] / 4) + 1)
#wagepan$D <- wagepan$G > dim(wagepan)[1]/8  &  wagepan$T >= 2

#wagepan$Y <- runif(dim(wagepan)[1]) * (1 + wagepan$D)
#wagepan <- wagepan[order(wagepan$G, wagepan$T), ]
#write.csv(wagepan, "test_data.csv")
#did_multiplegt_dyn <- did_multiplegt_dyn(wagepan, "Y", "G", "T", "D", effects = 2)
#print(did_multiplegt_dyn)

#wagepan <- wagepan[order(wagepan$nr, wagepan$year), ]
#wagepan$clus <- wagepan$nr %% 3
#wagepan$wr_clus <- round(runif(n = nrow(wagepan)), 2)
#wagepan$union <- 2 * wagepan$union
#res <- did_multiplegt_dyn(wagepan, "lwage", "nr", "year", "union", effects = 6, placebo = 4)
#res <- did_multiplegt_dyn(wagepan, "lwage", "nr", "year", "union", cluster = "clus", effects = 4, trends_nonparam = "exper", same_switchers = TRUE, controls = c("married", "hours"), effects_equal = TRUE, placebo = 3, same_switchers_pl = TRUE) 

#wagepan$year[2] <- 1980
#wagepan$married[2] <- NA
#wagepan$union[1] <- NA
#wagepan$union[50] <- NA
#wagepan$union[224] <- NA
#wagepan$union[1042] <- NA
#wagepan$union[33] <- NA

#did_multiplegt_dyn(wagepan, "lwage", "nr", "year", "union", effects = 5, placebo = 2, normalized = TRUE, date_first_switch = c("by_baseline_treat", "console"), normalized_weights = "by_calendar", trends_nonparam = "exper", effects_equal = TRUE, controls = c("married", "hours"), same_switchers = TRUE)

# Examples dont_drop_larger_lower
#res <- did_multiplegt_dyn(wagepan, "lwage", "nr", "year", "hours")
#res <- did_multiplegt_dyn(wagepan, "lwage", "nr", "year", "hours", dont_drop_larger_lower = TRUE)

#View(res)
# First union missing (discrepancy with Stata)
#wagepan <- wagepan[order(wagepan$nr, wagepan$year), ]
#wagepan <- wagepan %>% group_by(.data$nr) %>% mutate(first_obs = row_number() == 1)
#wagepan$union[wagepan$first_obs == 1] <- NA
#wagepan$union[1] <- NA
#wagepan$union[50] <- NA
#wagepan$union[224] <- NA
#wagepan$union[1042] <- NA
#wagepan$union[33] <- NA
#View(wagepan)

#wagepan <- subset(wagepan, (wagepan$year > 1983 & wagepan$nr < 50) | wagepan$nr>=50)
#results <- did_multiplegt_dyn(wagepan, "lwage", "nr", "year", "union", effects = 3)
#View(results)

#did_trends_nonparam <- function(varlist) {
    #trends_nonparam_XX <- ""
    #for (i in varlist) {
        #trends_nonparam_XX <- paste0(trends_nonparam_XX, ".data$", i, ", ")
    #}
    #trends_nonparam_XX <- substr(trends_nonparam_XX, 1, nchar(trends_nonparam_XX) - 2)
    #assign("trends_nonparam_XX", trends_nonparam_XX, envir = globalenv())
#}

#joint_trends <- function(df, var, trends_nonparam) {
    #joint_vars <- c(trends_nonparam, var)
    #print(joint_vars)
    #df$joint_trends_XX <- df[joint_vars]
    #df
#}

#trends_nonparam <- c("nr", "married")
#wagepan <- joint_trends(wagepan, "year", trends_nonparam)
#wagepan <- wagepan %>% 
    #group_by(.data$joint_trends_XX) %>% 
    #mutate(avg_union = mean(union, na.rm = TRUE))
#View(wagepan)

