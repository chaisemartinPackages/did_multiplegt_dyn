#install_github("chaisemartinPackages/did_multiplegt_dyn/R", force = TRUE) 
#install_github("chaisemartinPackages/did_multiplegt_dyn/did_multiplegt_dyn_all_pl/R", force = TRUE) 

library(devtools)

library(haven)
library(DIDmultiplegtDYN)
library(DIDmultiplegtDYNallPL)
data <- read_dta("https://raw.githubusercontent.com/chaisemartinPackages/did_multiplegt_dyn/main/vignettes/assets/vignette_3_data.dta")
did <- did_multiplegt_dyn_all_pl(
    df = data,
    outcome = "Y",
    group = "G",
    time = "T",
    treatment = "D",
    effects = 1,
    placebo = 2
)
print(did)
