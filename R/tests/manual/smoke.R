# quick smoke test - small data, small bootstrap count, no callr to start
suppressPackageStartupMessages({
    library(DIDmultiplegtDYN)
    library(polars)
})

data(favara_imbs)
df <- as.data.frame(favara_imbs)
cat("nrow:", nrow(df), " ncol:", ncol(df), "\n")
cat("groups:", length(unique(df$county)), "\n")

# Force in-process to surface errors quickly
options(
    DID_BOOTSTRAP_METHOD = "subprocess",
    DID_BOOTSTRAP_SAMPLE_DIR = NULL,
    DID_BOOTSTRAP_LOAD_SAMPLES = FALSE
)

t0 <- Sys.time()
res <- did_multiplegt_dyn(
    df = df, outcome = "Dl_hpi", group = "county",
    time = "year", treatment = "inter_bra",
    effects = 3, placebo = 2, cluster = "state_n",
    graph_off = TRUE, bootstrap = list(5, 42)
)
cat(sprintf("elapsed: %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))
print(res$results$Effects)
print(res$results$ATE)
print(res$results$Placebos)
