# ----------------------------------------------------------------------------
# Memory test: favara_imbs scaled up 50x. The original (row-replication) path
# leaks polars buffers across bootstrap iterations because they live in
# Rust's allocator and gc()/rm() can't return those pages to the OS. The new
# (subprocess) path runs each engine call in callr::r() so memory is
# reclaimed when the worker exits.
#
# This script runs Test A (one bootstrap=30 run) and Test B (3 sequential
# runs without bootstrap) under both methods and prints RSS at each step.
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(DIDmultiplegtDYN)
    library(polars)
    library(data.table)
})

cat("R version:        ", R.version.string, "\n", sep = "")
cat("DIDmultiplegtDYN: ", as.character(packageVersion("DIDmultiplegtDYN")), "\n", sep = "")

data(favara_imbs)
fi_big <- rbindlist(lapply(seq_len(50), function(i) {
    tmp <- as.data.table(as.data.frame(favara_imbs))
    tmp[, county := county + (i - 1L) * 100000L]
    tmp
}))
fi_big <- as.data.frame(fi_big)
cat("scaled rows:    ", nrow(fi_big), "\n")
cat("scaled groups:  ", length(unique(fi_big$county)), "\n\n")

pid <- Sys.getpid()
get_rss_gb <- function() {
    # macOS: ps -o rss= -p $PID  -> RSS in KB
    out <- tryCatch(
        system(sprintf("ps -o rss= -p %d", pid), intern = TRUE),
        error = function(e) NA_character_
    )
    val <- suppressWarnings(as.numeric(trimws(out)))
    if (length(val) == 0 || is.na(val)) return(NA_real_)
    val[1] / (1024 * 1024)
}

run_one <- function(method_label, method_value) {
    cat("\n############################################################\n")
    cat("##  Method: ", method_label, "\n", sep = "")
    cat("############################################################\n")
    options(
        DID_BOOTSTRAP_METHOD       = method_value,
        DID_BOOTSTRAP_SAMPLE_DIR   = NULL,
        DID_BOOTSTRAP_LOAD_SAMPLES = FALSE
    )

    # --- Test A: bootstrap leak -------------------------------------------
    gc(full = TRUE, verbose = FALSE)
    cat(sprintf("[A] RSS before:                 %5.2f GB\n", get_rss_gb()))

    t0 <- Sys.time()
    res <- did_multiplegt_dyn(
        df = fi_big, outcome = "Dl_hpi", group = "county",
        time = "year", treatment = "inter_bra",
        effects = 3, placebo = 2, cluster = "state_n",
        graph_off = TRUE, bootstrap = list(30, 1)
    )
    el <- as.numeric(Sys.time() - t0, units = "secs")

    cat(sprintf("[A] elapsed:                    %5.1f s\n", el))
    cat(sprintf("[A] RSS after bootstrap=30:     %5.2f GB\n", get_rss_gb()))
    rm(res); gc(full = TRUE, verbose = FALSE)
    cat(sprintf("[A] RSS after gc:               %5.2f GB\n", get_rss_gb()))

    # --- Test B: between-run leak (no bootstrap) --------------------------
    for (i in 1:3) {
        res <- did_multiplegt_dyn(
            df = fi_big, outcome = "Dl_hpi", group = "county",
            time = "year", treatment = "inter_bra",
            effects = 3, placebo = 2, cluster = "state_n",
            graph_off = TRUE
        )
        cat(sprintf("[B] RSS after run %d:            %5.2f GB\n", i, get_rss_gb()))
        rm(res); gc(full = TRUE, verbose = FALSE)
        cat(sprintf("[B] RSS after gc:               %5.2f GB\n", get_rss_gb()))
    }
}

# Run new path first; then legacy.
run_one("subprocess (NEW)",      "subprocess")
run_one("row_replication (OLD)", "row_replication")
