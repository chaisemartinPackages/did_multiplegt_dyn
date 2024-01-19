library(devtools)
#devtools::install_github('lbraglia/RStata')
library(RStata)

setwd(gsub("Test", "DIDmultiplegtDYN", dirname(sys.frame(1)$ofile)))
devtools::load_all()

if (is.null(options("RStata.StataPath")[[1]])) {
    chooseStataBin()  # Locate your Stata executable from the popup window
}
options("RStata.StataVersion" = 18)

## Test 1: Example from Favara-Imbs
#setwd(gsub("/DIDmultiplegtDYN", "", getwd()))
#data =  haven::read_dta("favara_imbs_did_multiplegt_dyn.dta")
#did_multiplegt_dyn <- did_multiplegt_dyn(
    #df = data, Y = "Dl_vloans_b", G = "county", T = "year", D = "inter_bra",
    #effects = 8, placebo = 3, graph_off = TRUE)
##print(did_multiplegt_dyn)
#setwd(paste0(getwd(), "/Stata"))
#stata <- stata("qui do did_multiplegt_dyn_new_path_option.do
       #qui did_multiplegt_dyn Dl_vloans_b county year inter_bra, effects(8) placebo(3) save_results(res) graph_off
       #qui use res.dta, clear
       #"
#, data.in = data, data.out = TRUE)                         

#stata_effects <- subset(stata, stata$time_to_treat > 0)
#stata_effects <- stata_effects[order(stata_effects$time_to_treat), ]
#stata_effects$time_to_treat <- NULL
#cat("Mean abs. diff, Effects: ", mean(abs(unlist(stata_effects - did_multiplegt_dyn$results$Effects))), "\n")

# Test 2: Example from WAGEPAN data
library(wooldridge)
data("wagepan")
setwd(gsub("/DIDmultiplegtDYN", "", getwd()))
wagepan$over_gr <- wagepan$nr %% 3
did_multiplegt_dyn <- did_multiplegt_dyn(
    df = wagepan, Y = "lwage", G = "nr", T = "year", D = "union",
    effects = 5, placebo = 2, graph_off = TRUE, cluster = "over_gr")
#print(did_multiplegt_dyn)
setwd(paste0(getwd(), "/Stata"))
stata <- stata("qui do did_multiplegt_dyn_new_path_option.do
       qui did_multiplegt_dyn lwage nr year union, effects(5) placebo(2) save_results(res) graph_off cluster(over_gr)
       qui use res.dta, clear
       "
, data.in = wagepan, data.out = TRUE)                         

stata_effects <- subset(stata, stata$time_to_treat > 0)
stata_effects <- stata_effects[order(stata_effects$time_to_treat), ]
stata_effects$time_to_treat <- NULL
cat("Mean abs. diff, Effects: ", mean(abs(unlist(stata_effects - did_multiplegt_dyn$results$Effects))), "\n")
stata_placebo <- subset(stata, stata$time_to_treat < 0)
stata_placebo$time_to_treat <- -stata_placebo$time_to_treat
stata_placebo <- stata_placebo[order(stata_placebo$time_to_treat), ]
stata_placebo$time_to_treat <- NULL
cat("Mean abs. diff, Placebo: ", mean(abs(unlist(stata_placebo - did_multiplegt_dyn$results$Placebo))), "\n")

unlink("res.dta")






