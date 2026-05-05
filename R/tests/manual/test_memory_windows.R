# ============================================================================
# Test 3 — memory leak (Windows-specific RSS reader)
#
# favara_imbs scaled 50x by re-labeling counties. Two passes in the same R
# session: NEW (subprocess + weights) first, then LEGACY (row replication).
# Watch the RSS-after-gc number after the bootstrap; on Windows it should
# stay bounded for NEW and climb for LEGACY.
# ============================================================================
suppressPackageStartupMessages({
    library(DIDmultiplegtDYN)
    library(polars)
    library(data.table)
})

cat("R version:        ", R.version.string, "\n", sep = "")
cat("DIDmultiplegtDYN: ", as.character(packageVersion("DIDmultiplegtDYN")), "\n", sep = "")

# Scale favara_imbs 50x by relabeling counties.
data(favara_imbs)
fi_big <- rbindlist(lapply(seq_len(50), function(i) {
    tmp <- as.data.table(as.data.frame(favara_imbs))
    tmp[, county := county + (i - 1L) * 100000L]
    tmp
}))
fi_big <- as.data.frame(fi_big)
cat("scaled rows:    ", nrow(fi_big), "\n")
cat("scaled groups:  ", length(unique(fi_big$county)), "\n\n")

# Windows-compatible RSS reader. Uses PowerShell's Get-Process for the
# parent R process. Returns NA on non-Windows hosts so the script still
# runs (it just prints NA in the RSS column).
pid <- Sys.getpid()
get_rss_gb <- function() {
    out <- tryCatch(
        suppressWarnings(system(
            paste0('powershell -NoProfile -Command "(Get-Process -Id ', pid,
                   ').WorkingSet64 / 1GB"'),
            intern = TRUE
        )),
        error = function(e) NA_character_
    )
    val <- suppressWarnings(as.numeric(trimws(out)))
    if (length(val) == 0 || is.na(val[1])) return(NA_real_)
    val[1]
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

# NEW first (so its baseline is fresh), LEGACY second.
run_one("subprocess (NEW)",      "subprocess")
run_one("row_replication (OLD)", "row_replication")
