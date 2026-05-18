#!/usr/bin/env Rscript
# =============================================================================
# Single-method runner. Runs ONE spec with ONE bootstrap method in this R
# process and saves the result (or error) as an RDS file.
#
# Usage:
#   Rscript run_one_path.R <spec_name> <method> <sample_dir> <load_samples> <out_rds>
#
#   method        = "subprocess" | "row_replication"
#   sample_dir    = directory used by DID_BOOTSTRAP_SAMPLE_DIR
#   load_samples  = "TRUE" | "FALSE"   (sets DID_BOOTSTRAP_LOAD_SAMPLES)
#   out_rds       = path to write an RDS file with the captured result
#
# RDS payload:
#   list(
#     status   = "ok" | "error",
#     err_msg  = NA | error message,
#     Effects  = matrix or NULL,
#     ATE      = matrix or NULL,
#     Placebos = matrix or NULL,
#     elapsed  = seconds (numeric)
#   )
# =============================================================================

suppressPackageStartupMessages({
    library(DIDmultiplegtDYN)
    library(polars)
    library(callr)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) {
    stop("Usage: Rscript run_one_path.R <spec_name> <method> <sample_dir> <load_samples> <out_rds>")
}
spec_name    <- args[[1L]]
method       <- args[[2L]]
sample_dir   <- args[[3L]]
load_samples <- as.logical(args[[4L]])
out_rds      <- args[[5L]]

script_dir <- tryCatch({
    cargs <- commandArgs(trailingOnly = FALSE)
    f <- sub("--file=", "", cargs[grep("--file=", cargs)])
    if (length(f)) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
source(file.path(script_dir, "specs.R"))

if (is.null(SPECS_BY_NAME[[spec_name]])) {
    stop("Unknown spec '", spec_name, "'")
}
spec <- SPECS_BY_NAME[[spec_name]]

BOOT_REPS <- 4   # double on purpose for did_multiplegt_dyn() validator
BOOT_SEED <- 1234L

loader_fn <- get(spec$loader, mode = "function")
df <- loader_fn()
if (!is.null(spec$preprocess)) {
    pre_fn <- get(spec$preprocess, mode = "function")
    df <- pre_fn(df)
}

options(
    DID_BOOTSTRAP_METHOD       = method,
    DID_BOOTSTRAP_SAMPLE_DIR   = sample_dir,
    DID_BOOTSTRAP_LOAD_SAMPLES = load_samples
)

t0 <- Sys.time()
res <- tryCatch({
    full <- c(list(df = df, graph_off = TRUE,
                   bootstrap = list(BOOT_REPS, BOOT_SEED)),
              spec$args)
    do.call(did_multiplegt_dyn, full)
}, error = function(e) e)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")

extract_result <- function(res) {
    inner <- res$results
    if (is.null(inner)) inner <- res
    list(
        Effects  = if (!is.null(inner$Effects))  as.matrix(inner$Effects)  else NULL,
        ATE      = if (!is.null(inner$ATE))      as.matrix(inner$ATE)      else NULL,
        Placebos = if (!is.null(inner$Placebos)) as.matrix(inner$Placebos) else NULL
    )
}

if (inherits(res, "error")) {
    payload <- list(
        status   = "error",
        err_msg  = conditionMessage(res),
        Effects  = NULL, ATE = NULL, Placebos = NULL,
        elapsed  = elapsed
    )
} else {
    parts <- extract_result(res)
    payload <- list(
        status   = "ok",
        err_msg  = NA_character_,
        Effects  = parts$Effects,
        ATE      = parts$ATE,
        Placebos = parts$Placebos,
        elapsed  = elapsed
    )
}

dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(payload, out_rds)
cat(sprintf("  [%s] %s -> %s (%.1fs)\n",
            payload$status, method, basename(out_rds), elapsed))
if (payload$status == "error") {
    cat("    err: ", payload$err_msg, "\n", sep = "")
}
