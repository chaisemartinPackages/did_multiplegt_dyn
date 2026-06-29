#' Compute the average number of time periods over which a treatment dose's
#' effect is accumulated, mirroring the Stata `avg_time_periods` option.
#'
#' Thin R wrapper around the Rcpp ports of `_compute_Tg` and `_avg_cumul`. The
#' polars dataframe returned by `did_multiplegt_main` is materialised into
#' R vectors via `pl_cols_to_r` (which works for both polars dataframes and
#' regular data.frames), then the heavy lifting is done in C++ via a
#' precomputed `(group, time) -> row_index` map.
#'
#' @param df_in  Polars dataframe (or data.frame) returned by
#'   `did_multiplegt_main`. Must contain `group_XX`, `time_XX`, `treatment_XX`,
#'   `outcome_XX`, `F_g_XX`, `T_g_XX`, `d_sq_int_XX`, `S_g_XX`, `N_gt_XX`,
#'   and either `d_sq_XX` or `d_sq_XX_orig` (when `continuous` is specified).
#' @param l_XX Number of event-study effects.
#' @param T_max_XX Maximum time period.
#' @param same_switchers Logical, from main call.
#' @param switchers Char ("", "in" or "out").
#' @param continuous Continuous polynomial order, or NULL.
#' @return A list with elements:
#'   * `avg_cumul`: weighted average horizon
#'   * `num`, `den`: the numerator / denominator behind `avg_cumul`
#'   * `nswitch`: numeric vector of length `l_XX` giving the per-effect
#'     switcher count.
#' @keywords internal
#' @noRd
did_multiplegt_avg_cumul <- function(df_in, l_XX, T_max_XX, same_switchers,
                                     switchers, continuous) {

  d1_col <- if (!is.null(continuous) && continuous > 0) "d_sq_XX_orig" else "d_sq_XX"

  keep <- c("group_XX", "time_XX", "treatment_XX", "outcome_XX",
            d1_col, "F_g_XX", "T_g_XX", "d_sq_int_XX", "S_g_XX", "N_gt_XX")

  # Materialise the requested columns as a list of R vectors regardless of
  # whether df_in is a polars dataframe, a data.table, or a plain data.frame.
  # Test for a polars-style `get_column` method (some polars class names have
  # shifted across versions, so probe the method directly).
  get_col <- tryCatch(df_in$get_column, error = function(e) NULL)
  cols <- if (is.function(get_col)) {
    lapply(stats::setNames(keep, keep), function(cn) as.vector(get_col(cn)))
  } else {
    stats::setNames(lapply(keep, function(cn) df_in[[cn]]), keep)
  }

  # Guard: every requested column must be present and non-empty.
  missing_cols <- keep[vapply(cols, function(v) is.null(v) || length(v) == 0L, logical(1))]
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "did_multiplegt_avg_cumul: required column(s) missing from df_in: %s (class: %s)",
      paste(missing_cols, collapse = ", "),
      paste(class(df_in), collapse = "/")
    ))
  }

  staged <- data.table::as.data.table(cols)
  data.table::setnames(staged, d1_col, "D1")
  data.table::setorderv(staged, c("group_XX", "time_XX"))

  G_int   <- as.integer(staged$group_XX)
  T_int   <- as.integer(staged$time_XX)
  D_num   <- as.numeric(staged$treatment_XX)
  Y_num   <- as.numeric(staged$outcome_XX)
  D1_num  <- as.numeric(staged$D1)
  F_num   <- as.numeric(staged$F_g_XX)
  TGC_num <- as.numeric(staged$T_g_XX)
  CLS_int <- as.integer(staged$d_sq_int_XX)
  SG_num  <- as.numeric(staged$S_g_XX)
  W_num   <- as.numeric(staged$N_gt_XX)

  EV_int  <- as.integer(!is.na(F_num) & F_num < (T_max_XX + 1L))
  NGT_int <- as.integer(!is.na(Y_num))

  Tg_ph <- compute_Tg_cpp(G_int, T_int, Y_num, F_num, TGC_num,
                          EV_int, CLS_int, NGT_int)
  Mg_ph <- pmin(as.numeric(l_XX), Tg_ph - F_num + 1)

  ssw_flag <- if (isTRUE(same_switchers)) 1L else 0L
  sw_dir   <- if (is.null(switchers)) "" else switchers

  out <- avg_cumul_cpp(G_int, T_int, D_num, Y_num, D1_num, F_num, Tg_ph,
                       EV_int, CLS_int, Mg_ph, NGT_int, SG_num, W_num,
                       as.integer(l_XX), ssw_flag, sw_dir)
  out
}
