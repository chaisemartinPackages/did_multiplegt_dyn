#!/usr/bin/env Rscript
# =============================================================================
# Bootstrap equivalence runner -- runs ONE spec by name, using TWO separate
# Rscript subprocesses (one per method) for maximum memory isolation.
#
# Usage:
#   Rscript run_one_spec.R <spec_name> [out_csv]
#
# Per-spec process layout:
#   1. mkdir shared CSV sample dir
#   2. Rscript run_one_path.R <spec> subprocess <sd> FALSE   new.rds
#         -> writes per-replicate unit selections to <sd>
#   3. Rscript run_one_path.R <spec> row_replication <sd> TRUE  old.rds
#         -> reads the same CSVs so the two paths see identical draws
#   4. Load both RDS files, compare Effects/ATE/Placebos, write a CSV summary.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
    stop("Usage: Rscript run_one_spec.R <spec_name> [out_csv]")
}
spec_name <- args[[1L]]
out_csv   <- if (length(args) >= 2L) args[[2L]] else NA_character_

script_dir <- tryCatch({
    cargs <- commandArgs(trailingOnly = FALSE)
    f <- sub("--file=", "", cargs[grep("--file=", cargs)])
    if (length(f)) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())

TOL <- 1e-8

cat(sprintf("\n=== %s ===\n", spec_name))

# ---- per-spec scratch dir --------------------------------------------------
sd       <- tempfile(paste0("boot_eq_", spec_name, "_"))
work_dir <- tempfile(paste0("boot_eq_results_", spec_name, "_"))
dir.create(sd, showWarnings = FALSE)
dir.create(work_dir, showWarnings = FALSE)
new_rds <- file.path(work_dir, "new.rds")
old_rds <- file.path(work_dir, "old.rds")

# ---- run a method in its own Rscript ---------------------------------------
runner_path <- file.path(script_dir, "run_one_path.R")
invoke_one <- function(method, load_samples, out_rds) {
    cmd <- c(
        "--no-save", "--no-restore", runner_path,
        spec_name, method, sd,
        if (load_samples) "TRUE" else "FALSE",
        out_rds
    )
    # Show child output live so the driver log is informative.
    system2("Rscript", cmd, stdout = "", stderr = "")
}

# Order matters: new path writes the CSVs, legacy path loads them.
invoke_one("subprocess",       load_samples = FALSE, out_rds = new_rds)
invoke_one("row_replication",  load_samples = TRUE,  out_rds = old_rds)

# ---- compare ---------------------------------------------------------------
if (!file.exists(new_rds) || !file.exists(old_rds)) {
    rows <- data.frame(
        spec     = spec_name, status = "ERROR",
        quantity = NA_character_, diff = NA_real_,
        t_new = NA_real_, t_old = NA_real_,
        stringsAsFactors = FALSE
    )
} else {
    new_payload <- readRDS(new_rds)
    old_payload <- readRDS(old_rds)

    max_abs_diff <- function(a, b) {
        if (is.null(a) && is.null(b)) return(0)
        if (is.null(a) || is.null(b)) return(NA_real_)
        if (!identical(dim(a), dim(b))) return(NA_real_)
        diff <- abs(as.numeric(a) - as.numeric(b))
        both_na <- is.na(as.numeric(a)) & is.na(as.numeric(b))
        diff[both_na] <- 0
        diff[is.na(diff)] <- Inf
        if (length(diff) == 0) 0 else max(diff)
    }

    cmp_one <- function(label, a, b) {
        d <- max_abs_diff(a, b)
        st <- if (is.na(d)) "ERROR" else if (d < TOL) "OK" else "FAIL"
        data.frame(spec = spec_name, status = st, quantity = label,
                   diff = d, stringsAsFactors = FALSE)
    }

    if (new_payload$status == "error" && old_payload$status == "error") {
        # Both paths errored. As long as the error message matches, this is a
        # consistent outcome -- "same exact result" in the sense the user asked
        # for. The misfit specs (data does not support the option combo) land
        # here naturally.
        same_err <- identical(trimws(new_payload$err_msg),
                              trimws(old_payload$err_msg))
        st <- if (same_err) "OK" else "FAIL"
        rows <- data.frame(
            spec = spec_name, status = st, quantity = "error",
            diff = if (same_err) 0 else NA_real_,
            stringsAsFactors = FALSE
        )
        cat(sprintf("  Both paths errored. messages %s\n",
                    if (same_err) "MATCH (consistent failure -> OK)"
                                 else "DIFFER (FAIL)"))
        if (!same_err) {
            cat("    new: ", new_payload$err_msg, "\n", sep = "")
            cat("    old: ", old_payload$err_msg, "\n", sep = "")
        }
    } else if (new_payload$status == "error" || old_payload$status == "error") {
        rows <- data.frame(
            spec = spec_name, status = "FAIL", quantity = "error",
            diff = NA_real_, stringsAsFactors = FALSE
        )
        cat(sprintf("  Path mismatch: new=%s old=%s\n",
                    new_payload$status, old_payload$status))
        if (new_payload$status == "error")
            cat("    new err: ", new_payload$err_msg, "\n", sep = "")
        if (old_payload$status == "error")
            cat("    old err: ", old_payload$err_msg, "\n", sep = "")
    } else {
        parts <- list(
            cmp_one("Effects",  new_payload$Effects,  old_payload$Effects),
            cmp_one("ATE",      new_payload$ATE,      old_payload$ATE),
            cmp_one("Placebos", new_payload$Placebos, old_payload$Placebos)
        )
        rows <- do.call(rbind, parts)
    }
    rows$t_new <- new_payload$elapsed
    rows$t_old <- old_payload$elapsed
}

# ---- report ----------------------------------------------------------------
for (i in seq_len(nrow(rows))) {
    cat(sprintf("  [%s] %s: max abs diff = %s\n",
                rows$status[i],
                ifelse(is.na(rows$quantity[i]), "<all>", rows$quantity[i]),
                if (is.na(rows$diff[i])) "NA" else sprintf("%.3e", rows$diff[i])))
}
cat(sprintf("  (NEW=%.1fs LEGACY=%.1fs)\n",
            rows$t_new[1], rows$t_old[1]))

if (!is.na(out_csv) && nzchar(out_csv)) {
    dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(rows, out_csv, row.names = FALSE)
}

unlink(sd, recursive = TRUE)
unlink(work_dir, recursive = TRUE)

bad <- sum(rows$status != "OK")
if (bad > 0) quit(status = 1) else quit(status = 0)
