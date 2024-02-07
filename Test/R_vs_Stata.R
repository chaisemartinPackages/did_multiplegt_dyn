rm(list = ls())
library(devtools)
#devtools::install_github('lbraglia/RStata')
library(RStata)

diff_print <- function(robj, sobj) {
    stata_effects <- subset(sobj, sobj$time_to_treat > 0)
    stata_effects <- stata_effects[order(stata_effects$time_to_treat), ]
    stata_effects$time_to_treat <- NULL
    print(round(abs(stata_effects - robj$results$Effects), 4))

    stata_placebo <- subset(sobj, sobj$time_to_treat < 0)
    stata_placebo$time_to_treat <- -stata_placebo$time_to_treat
    stata_placebo <- stata_placebo[order(stata_placebo$time_to_treat), ]
    stata_placebo$time_to_treat <- NULL
    print(round(abs(stata_placebo - robj$results$Placebo), 4))
}

setwd(gsub("Test", "DIDmultiplegtDYN", dirname(sys.frame(1)$ofile)))
devtools::load_all()
options("RStata.StataPath" = "\"C:\\Program Files\\Stata18\\StataMP-64\"")
options("RStata.StataVersion" = 18)

## Test 1: Example from Favara-Imbs
setwd(gsub("/DIDmultiplegtDYN", "", getwd()))
repo <- "chaisemartinPackages/ApplicationData/main" 
file <- "favara_imbs_did_multiplegt_dyn.dta"
url <- paste("https://raw.githubusercontent.com", repo, file, sep = "/")
data <-  haven::read_dta(url)
did_multiplegt_dyn <- did_multiplegt_dyn(
    df = data, Y = "Dl_vloans_b", G = "county", T = "year", D = "inter_bra",
    effects = 8, placebo = 3, graph_off = TRUE)
#print(did_multiplegt_dyn)
stata_mat <- stata("qui do did_multiplegt_dyn_working_24_01_24.do
       qui did_multiplegt_dyn_new Dl_vloans_b county year inter_bra, effects(8) placebo(3) save_results(res) graph_off
       qui use res.dta, clear", data.in = data, data.out = TRUE)                         
diff_print(did_multiplegt_dyn, stata_mat)
stata_mat <- NULL
stata("clear")

# Test 2: Example from WAGEPAN data
library(wooldridge)
data("wagepan")
setwd(gsub("/DIDmultiplegtDYN", "", getwd()))
wagepan$over_gr <- wagepan$nr %% 3
did_multiplegt_dyn <- did_multiplegt_dyn(
    df = wagepan, Y = "lwage", G = "nr", T = "year", D = "union",
    effects = 5, placebo = 2, graph_off = TRUE, cluster = "over_gr")
stata_mat <- stata("qui do did_multiplegt_dyn_working_24_01_24.do
       qui did_multiplegt_dyn_new lwage nr year union, effects(5) placebo(2) save_results(res) graph_off cluster(over_gr)
       qui use res.dta, clear", data.in = wagepan, data.out = TRUE)                         
diff_print(did_multiplegt_dyn, stata_mat)
stata_mat <- NULL
stata("clear")

library(wooldridge)
data("wagepan")
setwd(gsub("/DIDmultiplegtDYN", "", getwd()))
wagepan$over_gr <- wagepan$nr %% 3
did_multiplegt_dyn <- did_multiplegt_dyn(
    df = wagepan, Y = "lwage", G = "nr", T = "year", D = "union",
    effects = 5, placebo = 2, graph_off = TRUE, cluster = "over_gr")
stata_mat <- stata("qui do did_multiplegt_dyn_working_24_01_24.do
       qui did_multiplegt_dyn_new lwage nr year union, effects(5) placebo(2) save_results(res) graph_off cluster(over_gr)
       qui use res.dta, clear", data.in = wagepan, data.out = TRUE)                         
diff_print(did_multiplegt_dyn, stata_mat)
stata_mat <- NULL
stata("clear")

unlink("res.dta")






