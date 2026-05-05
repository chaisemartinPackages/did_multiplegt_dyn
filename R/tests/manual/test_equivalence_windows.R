# ============================================================================
# Test 1 — equivalence on favara_imbs (small, ~30s)
#
# Drives both bootstrap paths off the SAME unit-selection CSVs and verifies
# the resulting Effects/ATE/Placebos match to floating-point precision.
# ============================================================================
suppressPackageStartupMessages({
    library(DIDmultiplegtDYN)
    library(polars)
})

cat("R version:        ", R.version.string, "\n", sep = "")
cat("DIDmultiplegtDYN: ", as.character(packageVersion("DIDmultiplegtDYN")), "\n\n", sep = "")

data(favara_imbs)
df <- as.data.frame(favara_imbs)

run <- function(method, sample_dir, load) {
    options(
        DID_BOOTSTRAP_METHOD       = method,
        DID_BOOTSTRAP_SAMPLE_DIR   = sample_dir,
        DID_BOOTSTRAP_LOAD_SAMPLES = load
    )
    did_multiplegt_dyn(
        df = df, outcome = "Dl_hpi", group = "county",
        time = "year", treatment = "inter_bra",
        effects = 3, placebo = 2, cluster = "state_n",
        graph_off = TRUE, bootstrap = list(8, 999)
    )
}

sd <- tempfile("boot_eq_"); dir.create(sd)
cat("Sample dir: ", sd, "\n", sep = "")

cat("\n--- NEW (subprocess + weights, writes CSVs) ---\n")
new <- run("subprocess",       sd, FALSE)

cat("\n--- LEGACY (row replication, replays CSVs) ---\n")
old <- run("row_replication",  sd, TRUE)

cmp <- function(label, a, b, tol = 1e-8) {
    if (is.null(a) && is.null(b)) {
        cat(sprintf("[skip] %s\n", label)); return(invisible(NULL))
    }
    diff <- max(abs(as.matrix(a) - as.matrix(b)), na.rm = TRUE)
    cat(sprintf("[%s] %s: max abs diff = %.3e\n",
                if (diff < tol) "OK" else "FAIL", label, diff))
}

cat("\n--- COMPARISON ---\n")
cmp("Effects",  new$results$Effects,  old$results$Effects)
cmp("ATE",      new$results$ATE,      old$results$ATE)
cmp("Placebos", new$results$Placebos, old$results$Placebos)

cat("\nNEW Effects:\n");    print(new$results$Effects)
cat("\nLEGACY Effects:\n"); print(old$results$Effects)
cat("\nNEW ATE:\n");        print(new$results$ATE)
cat("\nLEGACY ATE:\n");     print(old$results$ATE)
cat("\nNEW Placebos:\n");   print(new$results$Placebos)
cat("\nLEGACY Placebos:\n"); print(old$results$Placebos)
