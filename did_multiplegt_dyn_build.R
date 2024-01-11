library(devtools)
library(usethis)
library(roxygen2)

setwd(paste0(dirname(sys.frame(1)$ofile), "/R"))
document()
#devtools::check()