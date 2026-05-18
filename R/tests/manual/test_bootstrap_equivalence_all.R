#!/usr/bin/env Rscript
# =============================================================================
# Bootstrap equivalence test: NEW (subprocess + weights) vs LEGACY (row replication)
# -----------------------------------------------------------------------------
# Drives every R example from the tutorial chapters (ch01-ch04) plus the option
# sweep translated from tests_stata.do through both bootstrap paths, sharing a
# CSV of unit selections so the two paths see identical draws. Asserts that
# Effects / ATE / Placebos match to floating-point tolerance.
#
# Run as:   Rscript tests/manual/test_bootstrap_equivalence_all.R
# Log:      written to tests/manual/logs/test_bootstrap_equivalence_all.log
# =============================================================================

suppressPackageStartupMessages({
    library(DIDmultiplegtDYN)
    library(polars)
    library(callr)
    library(data.table)
})

# ---- Paths -----------------------------------------------------------------
script_dir <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f)) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())

tutorial_dir <- normalizePath(file.path(script_dir, "..", "did_multiplegt_dyn_tutorial"),
                              mustWork = FALSE)
data_dir     <- file.path(tutorial_dir, "data")

logs_dir <- file.path(script_dir, "logs")
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Bootstrap settings ----------------------------------------------------
BOOT_REPS <- 8L
BOOT_SEED <- 1234L
TOL       <- 1e-8

# ---- Dataset loaders -------------------------------------------------------
load_deryugina <- function() {
    e <- new.env()
    load(file.path(data_dir, "deryugina_2017.RData"), envir = e)
    e$df
}
load_east <- function() {
    e <- new.env()
    load(file.path(data_dir, "east_et_al_2023.RData"), envir = e)
    e$df
}
load_gentzkow <- function() {
    e <- new.env()
    load(file.path(data_dir, "gentzkowetal_didtextbook.RData"), envir = e)
    e$df
}
load_hollingsworth <- function() {
    e <- new.env()
    load(file.path(data_dir, "hollingsworth_et_al_2022.RData"), envir = e)
    e$df
}

# ---- Pre-processing helpers (mirroring the tutorial chapter code) ----------
prep_deryugina_controls <- function(df) {
    vars <- c("coastal", "land_area1970", "log_pop1969", "frac_young1969",
              "frac_old1969", "frac_black1969", "log_wage_pc1969", "emp_rate1969")
    for (v in vars) df[[paste0(v, "_year")]] <- df[[v]] * df$year
    df
}

prep_gentzkow_t11 <- function(df) {
    df[df$year <= df$first_change | df$same_treat_after_first_change == 1, ]
}

prep_hollingsworth_t13 <- function(df) df[df$rm == 0, ]
prep_hollingsworth_t14 <- function(df) {
    df$fg1 <- ave(ifelse(df$mm == 1, df$year, NA), df$state,
                  FUN = function(x) suppressWarnings(min(x, na.rm = TRUE)))
    df$fg1[!is.finite(df$fg1)] <- NA
    df[df$mm == 1, ]
}

# ---- east controls -------------------------------------------------------
EAST_CONTROLS <- c("stmarried", "stblack", "stother", "sthsdrop", "sthsgrad",
                   "stsomecoll", "pop0_4", "pop5_17", "pop18_24", "pop25_44",
                   "pop45_64", "urate", "incpc", "maxafdc", "abortconsent",
                   "abortmedr")

# =============================================================================
#  Test specs.
#  Each is a list with: name, loader, preprocess (optional), args (NAMED list
#  of did_multiplegt_dyn() arguments, EXCLUDING df / bootstrap / graph_off).
# =============================================================================
SPECS <- list(

    # -------- Chapter 01: deryugina_2017 ----------------------------------
    list(name = "ch01_T2_deryugina_basic",
         loader = load_deryugina,
         args = list(
             outcome = "log_curr_trans_ind_gov_pc",
             group = "county_fips", time = "year", treatment = "hurricane",
             effects = 3, placebo = 3, cluster = "county_fips"
         )),

    list(name = "ch01_T3_deryugina_controls",
         loader = load_deryugina,
         preprocess = prep_deryugina_controls,
         args = list(
             outcome = "log_curr_trans_ind_gov_pc",
             group = "county_fips", time = "year", treatment = "hurricane",
             effects = 3, placebo = 3, cluster = "county_fips",
             controls = c("coastal_year", "land_area1970_year", "log_pop1969_year",
                          "frac_young1969_year", "frac_old1969_year",
                          "frac_black1969_year", "log_wage_pc1969_year",
                          "emp_rate1969_year")
         )),

    # -------- Chapter 02: east_et_al_2023 ---------------------------------
    list(name = "ch02_T5_east_continuous_weight",
         loader = load_east,
         args = list(
             outcome = "lbw_detrend81", group = "plborn",
             time = "dob_y_p", treatment = "newsimeli",
             effects = 5, placebo = 5, continuous = 1,
             controls = EAST_CONTROLS, weight = "births",
             effects_equal = "all"
         )),

    list(name = "ch02_T6_east_normalized",
         loader = load_east,
         args = list(
             outcome = "lbw_detrend81", group = "plborn",
             time = "dob_y_p", treatment = "newsimeli",
             effects = 5, placebo = 5, continuous = 1,
             controls = EAST_CONTROLS, weight = "births",
             normalized = TRUE, effects_equal = "all"
         )),

    # -------- Chapter 03: gentzkow ----------------------------------------
    list(name = "ch03_T8_gentzkow_effects_equal",
         loader = load_gentzkow,
         args = list(
             outcome = "prestout", group = "cnty90",
             time = "year", treatment = "numdailies",
             effects = 4, placebo = 4, effects_equal = "all"
         )),

    list(name = "ch03_T10_gentzkow_normalized",
         loader = load_gentzkow,
         args = list(
             outcome = "prestout", group = "cnty90",
             time = "year", treatment = "numdailies",
             effects = 4, placebo = 4,
             normalized = TRUE, effects_equal = "all"
         )),

    list(name = "ch03_T11_gentzkow_same_switchers",
         loader = load_gentzkow,
         preprocess = prep_gentzkow_t11,
         args = list(
             outcome = "prestout", group = "cnty90",
             time = "year", treatment = "numdailies",
             effects = 2, effects_equal = "all", same_switchers = TRUE
         )),

    # -------- Chapter 04: hollingsworth -----------------------------------
    list(name = "ch04_T13_marijuana_medical",
         loader = load_hollingsworth,
         preprocess = prep_hollingsworth_t13,
         args = list(
             outcome = "ln_mj_use_365", group = "state",
             time = "year", treatment = "mm",
             effects = 3, placebo = 3, cluster = "state"
         )),

    list(name = "ch04_T14_marijuana_recreational_trends_nonparam",
         loader = load_hollingsworth,
         preprocess = prep_hollingsworth_t14,
         args = list(
             outcome = "ln_mj_use_365", group = "state",
             time = "year", treatment = "rm",
             effects = 3, placebo = 3,
             trends_nonparam = "fg1", cluster = "state"
         )),

    # =====================================================================
    #  Option sweep translated from tests_stata.do
    #  Stata wagepan dataset is not available; we use the deryugina dataset
    #  with column roles substituted as below. The point is option coverage,
    #  not effect interpretation.
    #
    #    y     -> log_curr_trans_ind_gov_pc
    #    g     -> county_fips
    #    t     -> year
    #    d     -> hurricane            (binary)
    #    cont  -> log_pop1969          (continuous control)
    #    npar  -> coastal              (binary nonparam trend)
    #    wght  -> emp_rate1969         (positive)
    #    clust -> county_fips          (cluster on group; we don't have state)
    # =====================================================================
    list(name = "stata_baseline",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5)),

    list(name = "stata_placebos",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2)),

    list(name = "stata_normalized",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, normalized=TRUE)),

    list(name = "stata_controls",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, controls="log_pop1969")),

    list(name = "stata_trends_nonparam",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, trends_nonparam="coastal")),

    list(name = "stata_trends_lin",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, trends_lin=TRUE)),

    list(name = "stata_continuous",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, continuous=1)),

    list(name = "stata_weight",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, weight="emp_rate1969")),

    list(name = "stata_cluster",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, cluster="county_fips")),

    list(name = "stata_same_switchers",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, same_switchers=TRUE)),

    list(name = "stata_same_switchers_pl",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2,
                     same_switchers=TRUE, same_switchers_pl=TRUE)),

    list(name = "stata_switchers_in",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, switchers="in")),

    list(name = "stata_switchers_out",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, switchers="out")),

    list(name = "stata_only_never_switchers",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, only_never_switchers=TRUE)),

    list(name = "stata_ci_level_90",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, ci_level=90)),

    list(name = "stata_ci_level_99",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, ci_level=99)),

    list(name = "stata_less_conservative_se",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, less_conservative_se=TRUE)),

    list(name = "stata_dont_drop_larger_lower",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, dont_drop_larger_lower=TRUE)),

    list(name = "stata_effects_equal",
         loader = load_deryugina,
         args = list(outcome="log_curr_trans_ind_gov_pc", group="county_fips",
                     time="year", treatment="hurricane",
                     effects=5, placebo=2, effects_equal="all"))
)

# =============================================================================
#  Test driver
# =============================================================================

# Single side-by-side run.
run_one_method <- function(method, df, args, sample_dir, load_samples) {
    options(
        DID_BOOTSTRAP_METHOD       = method,
        DID_BOOTSTRAP_SAMPLE_DIR   = sample_dir,
        DID_BOOTSTRAP_LOAD_SAMPLES = load_samples
    )
    full <- c(list(df = df, graph_off = TRUE,
                   bootstrap = list(BOOT_REPS, BOOT_SEED)),
              args)
    do.call(did_multiplegt_dyn, full)
}

# Pull comparable parts out of the result object.
extract_results <- function(res) {
    inner <- res$results
    if (is.null(inner)) inner <- res
    list(
        Effects  = as.matrix(inner$Effects),
        ATE      = as.matrix(inner$ATE),
        Placebos = if (!is.null(inner$Placebos)) as.matrix(inner$Placebos) else NULL
    )
}

# Max abs diff between two matrices (NA-safe). Returns NA if either is NULL.
max_abs_diff <- function(a, b) {
    if (is.null(a) && is.null(b)) return(0)
    if (is.null(a) || is.null(b)) return(NA_real_)
    if (any(dim(a) != dim(b))) return(NA_real_)
    diff <- abs(as.numeric(a) - as.numeric(b))
    diff[is.na(diff)] <- 0  # treat aligned NAs as equal
    if (length(diff) == 0) 0 else max(diff)
}

cmp_one <- function(label, a, b) {
    d <- max_abs_diff(a, b)
    status <- if (is.na(d)) "ERROR" else if (d < TOL) "OK" else "FAIL"
    list(status = status, diff = d, label = label)
}

run_spec <- function(spec) {
    cat(sprintf("\n=== %s ===\n", spec$name))
    df <- spec$loader()
    if (!is.null(spec$preprocess)) df <- spec$preprocess(df)

    sd <- tempfile(paste0("boot_eq_", spec$name, "_"))
    dir.create(sd, showWarnings = FALSE)
    on.exit(unlink(sd, recursive = TRUE), add = TRUE)

    t0 <- Sys.time()
    new <- tryCatch(
        run_one_method("subprocess", df, spec$args, sd, FALSE),
        error = function(e) {
            cat(sprintf("  [error in subprocess path] %s\n", conditionMessage(e)))
            NULL
        }
    )
    el_new <- as.numeric(Sys.time() - t0, units = "secs")

    t0 <- Sys.time()
    old <- tryCatch(
        run_one_method("row_replication", df, spec$args, sd, TRUE),
        error = function(e) {
            cat(sprintf("  [error in row_replication path] %s\n", conditionMessage(e)))
            NULL
        }
    )
    el_old <- as.numeric(Sys.time() - t0, units = "secs")

    if (is.null(new) || is.null(old)) {
        return(data.frame(
            spec = spec$name, status = "ERROR", quantity = NA_character_,
            diff = NA_real_, t_new = el_new, t_old = el_old,
            stringsAsFactors = FALSE
        ))
    }

    r_new <- extract_results(new)
    r_old <- extract_results(old)

    results <- list(
        cmp_one("Effects",  r_new$Effects,  r_old$Effects),
        cmp_one("ATE",      r_new$ATE,      r_old$ATE),
        cmp_one("Placebos", r_new$Placebos, r_old$Placebos)
    )
    out <- do.call(rbind, lapply(results, function(r) {
        cat(sprintf("  [%s] %s: max abs diff = %.3e\n",
                    r$status, r$label,
                    if (is.na(r$diff)) -1 else r$diff))
        data.frame(spec = spec$name, status = r$status,
                   quantity = r$label, diff = r$diff,
                   t_new = el_new, t_old = el_old,
                   stringsAsFactors = FALSE)
    }))
    cat(sprintf("  (NEW=%.1fs LEGACY=%.1fs)\n", el_new, el_old))
    out
}

# =============================================================================
#  Run every spec
# =============================================================================
cat(sprintf("Running %d bootstrap equivalence specs (reps=%d seed=%d)\n",
            length(SPECS), BOOT_REPS, BOOT_SEED))
cat(sprintf("Tolerance: %g\n", TOL))

t_start <- Sys.time()
all_rows <- do.call(rbind, lapply(SPECS, function(s) {
    tryCatch(run_spec(s), error = function(e) {
        cat(sprintf("  [hard error in %s] %s\n", s$name, conditionMessage(e)))
        data.frame(spec = s$name, status = "ERROR",
                   quantity = NA_character_, diff = NA_real_,
                   t_new = NA_real_, t_old = NA_real_,
                   stringsAsFactors = FALSE)
    })
}))
t_total <- as.numeric(Sys.time() - t_start, units = "secs")

# ----- Summary --------------------------------------------------------------
cat("\n\n===================== SUMMARY =====================\n")
ok    <- sum(all_rows$status == "OK", na.rm = TRUE)
fail  <- sum(all_rows$status == "FAIL", na.rm = TRUE)
errs  <- sum(all_rows$status == "ERROR", na.rm = TRUE)
total <- nrow(all_rows)
cat(sprintf("Total comparisons: %d   OK: %d   FAIL: %d   ERROR: %d\n",
            total, ok, fail, errs))
cat(sprintf("Elapsed: %.1f s\n", t_total))

bad <- all_rows[all_rows$status != "OK", ]
if (nrow(bad) > 0) {
    cat("\nFailing rows:\n")
    print(bad, row.names = FALSE)
}

# Persist machine-readable result for follow-up.
out_csv <- file.path(logs_dir, "test_bootstrap_equivalence_all.csv")
utils::write.csv(all_rows, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_csv))

if (fail + errs > 0) {
    quit(status = 1)
} else {
    cat("\nAll OK.\n")
}
