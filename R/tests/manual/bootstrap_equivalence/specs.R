# =============================================================================
# Bootstrap-equivalence test specs.
# Defines `SPECS` -- the list of test cases driving the per-spec runner.
# This file is sourced by run_one_spec.R and run_all_specs.sh (via a list dump).
# =============================================================================

# ---- Paths -----------------------------------------------------------------
# Locate specs.R's own directory robustly: when source()'d we cannot rely on
# commandArgs("--file="), so try the sourced-file ofile attribute first.
.specs_dir <- local({
    # 1) source() exposes the path via sys.frames -- search for it.
    found <- NULL
    for (i in seq_along(sys.frames())) {
        of <- sys.frame(i)$ofile
        if (!is.null(of) && length(of) == 1 && nzchar(of)) {
            found <- of; break
        }
    }
    if (!is.null(found)) return(dirname(normalizePath(found)))
    # 2) Fall back to commandArgs (when invoked directly with Rscript).
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f)) dirname(normalizePath(f)) else getwd()
})

# tests/manual/bootstrap_equivalence/specs.R -> ../../did_multiplegt_dyn_tutorial/data
TUTORIAL_DIR <- normalizePath(
    file.path(.specs_dir, "..", "..", "did_multiplegt_dyn_tutorial"),
    mustWork = FALSE
)
DATA_DIR <- file.path(TUTORIAL_DIR, "data")

# ---- Dataset loaders -------------------------------------------------------
load_deryugina <- function() {
    e <- new.env(); load(file.path(DATA_DIR, "deryugina_2017.RData"), envir = e); e$df
}
load_east <- function() {
    e <- new.env(); load(file.path(DATA_DIR, "east_et_al_2023.RData"), envir = e); e$df
}
load_gentzkow <- function() {
    e <- new.env(); load(file.path(DATA_DIR, "gentzkowetal_didtextbook.RData"), envir = e); e$df
}
load_hollingsworth <- function() {
    e <- new.env(); load(file.path(DATA_DIR, "hollingsworth_et_al_2022.RData"), envir = e); e$df
}
load_favara_imbs <- function() {
    # favara_imbs is bundled with the package.
    suppressPackageStartupMessages(library(DIDmultiplegtDYN))
    data("favara_imbs", package = "DIDmultiplegtDYN")
    as.data.frame(get("favara_imbs"))
}

# ---- Pre-processing helpers ------------------------------------------------
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

# Trimmed east-controls set: the full 16-control vector + effects=5/placebo=5
# combo from the original tutorial does not fit in 8GB during the baseline
# did_multiplegt_main() call (independent of bootstrap method). We trim to 5
# controls + effects=3/placebo=3 -- still exercises continuous + weight +
# multi-control regressions, but fits the test container.
EAST_CONTROLS <- c("stmarried", "sthsdrop", "pop25_44", "urate", "incpc")

# ---- Spec list -------------------------------------------------------------
# Each is a list with: name, loader, preprocess (optional), args (NAMED list
# of did_multiplegt_dyn() arguments EXCLUDING df / bootstrap / graph_off).
SPECS <- list(

    # -------- Chapter 01: deryugina_2017 ----------------------------------
    list(name = "ch01_T2_deryugina_basic",
         loader = "load_deryugina",
         args = list(
             outcome = "log_curr_trans_ind_gov_pc",
             group = "county_fips", time = "year", treatment = "hurricane",
             effects = 3, placebo = 3, cluster = "county_fips"
         )),

    list(name = "ch01_T3_deryugina_controls",
         loader = "load_deryugina",
         preprocess = "prep_deryugina_controls",
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
         loader = "load_east",
         args = list(
             outcome = "lbw_detrend81", group = "plborn",
             time = "dob_y_p", treatment = "newsimeli",
             effects = 3, placebo = 3, continuous = 1,
             controls = EAST_CONTROLS, weight = "births",
             effects_equal = "all"
         )),

    list(name = "ch02_T6_east_normalized",
         loader = "load_east",
         args = list(
             outcome = "lbw_detrend81", group = "plborn",
             time = "dob_y_p", treatment = "newsimeli",
             effects = 3, placebo = 3, continuous = 1,
             controls = EAST_CONTROLS, weight = "births",
             normalized = TRUE, effects_equal = "all"
         )),

    # -------- Chapter 03: gentzkow ----------------------------------------
    list(name = "ch03_T8_gentzkow_effects_equal",
         loader = "load_gentzkow",
         args = list(
             outcome = "prestout", group = "cnty90",
             time = "year", treatment = "numdailies",
             effects = 4, placebo = 4, effects_equal = "all"
         )),

    list(name = "ch03_T10_gentzkow_normalized",
         loader = "load_gentzkow",
         args = list(
             outcome = "prestout", group = "cnty90",
             time = "year", treatment = "numdailies",
             effects = 4, placebo = 4,
             normalized = TRUE, effects_equal = "all"
         )),

    list(name = "ch03_T11_gentzkow_same_switchers",
         loader = "load_gentzkow",
         preprocess = "prep_gentzkow_t11",
         args = list(
             outcome = "prestout", group = "cnty90",
             time = "year", treatment = "numdailies",
             effects = 2, effects_equal = "all", same_switchers = TRUE
         )),

    # -------- Chapter 04: hollingsworth -----------------------------------
    list(name = "ch04_T13_marijuana_medical",
         loader = "load_hollingsworth",
         preprocess = "prep_hollingsworth_t13",
         args = list(
             outcome = "ln_mj_use_365", group = "state",
             time = "year", treatment = "mm",
             effects = 3, placebo = 3, cluster = "state"
         )),

    list(name = "ch04_T14_marijuana_recreational_trends_nonparam",
         loader = "load_hollingsworth",
         preprocess = "prep_hollingsworth_t14",
         args = list(
             outcome = "ln_mj_use_365", group = "state",
             time = "year", treatment = "rm",
             effects = 3, placebo = 3,
             trends_nonparam = "fg1", cluster = "state"
         )),

    # -------- Option sweep (deryugina dataset, role substitutions) --------
    list(name = "stata_baseline",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane", effects = 5)),

    list(name = "stata_placebos",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2)),

    list(name = "stata_normalized",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, normalized = TRUE)),

    list(name = "stata_controls",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, controls = "log_pop1969")),

    list(name = "stata_trends_nonparam",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, trends_nonparam = "coastal")),

    list(name = "stata_trends_lin",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, trends_lin = TRUE)),

    list(name = "stata_continuous",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, continuous = 1)),

    list(name = "stata_weight",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, weight = "emp_rate1969")),

    list(name = "stata_cluster",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, cluster = "county_fips")),

    list(name = "stata_same_switchers",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, same_switchers = TRUE)),

    list(name = "stata_same_switchers_pl",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2,
                     same_switchers = TRUE, same_switchers_pl = TRUE)),

    list(name = "stata_switchers_in",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, switchers = "in")),

    list(name = "stata_switchers_out",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, switchers = "out")),

    list(name = "stata_only_never_switchers",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, only_never_switchers = TRUE)),

    list(name = "stata_ci_level_90",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, ci_level = 90)),

    list(name = "stata_ci_level_99",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, ci_level = 99)),

    list(name = "stata_less_conservative_se",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, less_conservative_se = TRUE)),

    list(name = "stata_dont_drop_larger_lower",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, dont_drop_larger_lower = TRUE)),

    list(name = "stata_effects_equal",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 5, placebo = 2, effects_equal = "all")),

    # =====================================================================
    #  Extra coverage of did_multiplegt_dyn() options that flow through the
    #  bootstrap path but were not yet exercised. Each verifies that NEW vs
    #  LEGACY produce the same result under the same draws.
    # =====================================================================

    # favara_imbs is the dataset bundled with the package. The user's bug
    # report used it; we exercise it directly.
    list(name = "favara_imbs_basic",
         loader = "load_favara_imbs",
         args = list(outcome = "Dl_hpi", group = "county",
                     time = "year", treatment = "inter_bra",
                     effects = 3, placebo = 2, cluster = "state_n")),

    list(name = "favara_imbs_weight",
         loader = "load_favara_imbs",
         args = list(outcome = "Dl_hpi", group = "county",
                     time = "year", treatment = "inter_bra",
                     effects = 3, placebo = 2, cluster = "state_n",
                     weight = "w1")),

    # Heterogeneity prediction; requires no controls and not normalized.
    list(name = "deryugina_predict_het",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 3, placebo = 2,
                     predict_het = list("coastal", -1))),

    # by-stratified bootstrap: the outer loop runs the bootstrap helper once
    # per stratum. Use coastal (binary) so we get two strata.
    list(name = "deryugina_by_coastal",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 3, placebo = 2, by = "coastal")),

    # normalized_weights = TRUE reports the implicit weights, computed after
    # bootstrap. Verifies normalized + bootstrap still match under both paths.
    list(name = "gentzkow_normalized_weights",
         loader = "load_gentzkow",
         args = list(outcome = "prestout", group = "cnty90",
                     time = "year", treatment = "numdailies",
                     effects = 3, placebo = 2,
                     normalized = TRUE, normalized_weights = TRUE)),

    # drop_if_d_miss_before_first_switch -- exercised on deryugina (no actual
    # missingness, but makes sure the option is passed through cleanly).
    list(name = "deryugina_drop_if_d_miss",
         loader = "load_deryugina",
         args = list(outcome = "log_curr_trans_ind_gov_pc", group = "county_fips",
                     time = "year", treatment = "hurricane",
                     effects = 3, placebo = 2,
                     drop_if_d_miss_before_first_switch = TRUE))
)

# Convenience: name -> spec
SPECS_BY_NAME <- setNames(SPECS, vapply(SPECS, function(s) s$name, character(1)))
