# ----------------------------------------------------------------------------
# Equivalence test: row-replication vs. weight-based bootstrap on identical
# bootstrap draws. The two paths are driven through the *same* CSVs of unit
# selections, so any mismatch in Effects/ATE/Placebos points to a bug in the
# new implementation.
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(DIDmultiplegtDYN)
    library(polars)
})

cat("R version:        ", R.version.string, "\n", sep = "")
cat("DIDmultiplegtDYN: ", as.character(packageVersion("DIDmultiplegtDYN")), "\n", sep = "")
cat("polars:           ", as.character(packageVersion("polars")), "\n\n", sep = "")

# Load east_et_al_2023 from the user's tutorial repo (cached locally to avoid
# repeated network hits).
cache_dir <- normalizePath(
    file.path("~/.cache/did_multiplegt_dyn"), winslash = "/", mustWork = FALSE
)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
data_url <- "https://raw.githubusercontent.com/Credible-Answers/did_multiplegt_dyn_tutorial/main/data/east_et_al_2023.RData"
data_file <- file.path(cache_dir, "east_et_al_2023.RData")
if (!file.exists(data_file)) {
    cat("Downloading east_et_al_2023 ...\n")
    utils::download.file(data_url, data_file, mode = "wb", quiet = TRUE)
}
load(data_file)

stcontrols <- c("stmarried", "stblack", "stother", "sthsdrop", "sthsgrad",
                "stsomecoll", "pop0_4", "pop5_17", "pop18_24", "pop25_44",
                "pop45_64", "urate", "incpc", "maxafdc", "abortconsent", "abortmedr")

run_args <- list(
    df            = df,
    outcome       = "lbw_detrend81",
    group         = "plborn",
    time          = "dob_y_p",
    treatment     = "newsimeli",
    effects       = 5,
    placebo       = 5,
    continuous    = 1,
    controls      = stcontrols,
    weight        = "births",
    effects_equal = "all",
    graph_off     = TRUE,
    bootstrap     = list(30, 12345)
)

sample_dir <- tempfile("boot_samples_")
cat("Sample dir: ", sample_dir, "\n", sep = "")

# ---- 1. New (subprocess + weights) path, writing samples ---------------------
cat("\n=== NEW path (subprocess + weights) =====================================\n")
options(
    DID_BOOTSTRAP_METHOD       = "subprocess",
    DID_BOOTSTRAP_SAMPLE_DIR   = sample_dir,
    DID_BOOTSTRAP_LOAD_SAMPLES = FALSE
)
t0 <- Sys.time()
res_new <- do.call(did_multiplegt_dyn, run_args)
cat(sprintf("NEW elapsed: %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))

# ---- 2. Legacy (row-replication) path, reading the same samples --------------
cat("\n=== LEGACY path (row replication, replaying CSVs) =======================\n")
options(
    DID_BOOTSTRAP_METHOD       = "row_replication",
    DID_BOOTSTRAP_SAMPLE_DIR   = sample_dir,
    DID_BOOTSTRAP_LOAD_SAMPLES = TRUE
)
t0 <- Sys.time()
res_old <- do.call(did_multiplegt_dyn, run_args)
cat(sprintf("LEGACY elapsed: %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))

# ---- 3. Compare ---------------------------------------------------------------
options(
    DID_BOOTSTRAP_METHOD       = "subprocess",
    DID_BOOTSTRAP_SAMPLE_DIR   = NULL,
    DID_BOOTSTRAP_LOAD_SAMPLES = FALSE
)

cmp <- function(label, a, b, tol = 1e-8) {
    if (is.null(a) && is.null(b)) {
        cat(sprintf("[skip] %s: both NULL\n", label)); return(invisible(NULL))
    }
    a <- as.matrix(a); b <- as.matrix(b)
    if (!all(dim(a) == dim(b))) {
        cat(sprintf("[FAIL] %s: dim mismatch %s vs %s\n",
                    label, paste(dim(a), collapse = "x"),
                    paste(dim(b), collapse = "x")))
        return(invisible(NULL))
    }
    diff <- max(abs(a - b), na.rm = TRUE)
    cat(sprintf("[%s] %s: max abs diff = %.3e\n",
                if (diff < tol) "OK" else "FAIL", label, diff))
}

cat("\n=== COMPARISON ==========================================================\n")
cmp("Effects",  res_new$results$Effects,  res_old$results$Effects)
cmp("ATE",      res_new$results$ATE,      res_old$results$ATE)
cmp("Placebos", res_new$results$Placebos, res_old$results$Placebos)

cat("\nNEW Effects:\n");   print(res_new$results$Effects)
cat("\nLEGACY Effects:\n"); print(res_old$results$Effects)
cat("\nNEW ATE:\n");       print(res_new$results$ATE)
cat("\nLEGACY ATE:\n");    print(res_old$results$ATE)
cat("\nNEW Placebos:\n");  print(res_new$results$Placebos)
cat("\nLEGACY Placebos:\n"); print(res_old$results$Placebos)
