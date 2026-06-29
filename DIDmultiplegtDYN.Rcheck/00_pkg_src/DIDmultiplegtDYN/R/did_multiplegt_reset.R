#' Reset option pre-processing for did_multiplegt_dyn
#'
#' Implements the Stata `reset(integer)` option: balances the panel to groups
#' with the maximum number of non-missing treatment observations, then splits
#' each group into sub-groups whenever the treatment has remained constant for
#' `reset` consecutive periods after a change. The original group identifier
#' is preserved as a clustering column when no `cluster` variable was supplied.
#'
#' @param df A data.frame containing the user's data.
#' @param group Name of the group column.
#' @param time Name of the time column.
#' @param treatment Name of the treatment column.
#' @param cluster Existing cluster column name (or NULL).
#' @param reset Positive integer: number of post-change periods of stability after which the group is "reset".
#' @return list(df, cluster) where df has the group column replaced and
#'         optionally an `old_group_XX` column added; cluster is the (possibly
#'         updated) cluster column name.
#' @keywords internal
#' @noRd
did_multiplegt_reset <- function(df, group, time, treatment, cluster, reset) {
  dt <- data.table::as.data.table(df)
  treatment_XX <- changed_XX <- spell_XX <- periods_since_change_XX <- NULL
  reset_dummy_XX <- reset_count_XX <- new_group_XX <- count_d_non_miss_XX <- NULL

  # Balance panel: keep only groups with the maximum count of non-missing treatment
  dt[, count_d_non_miss_XX := sum(!is.na(get(treatment))), by = group]
  max_count <- max(dt$count_d_non_miss_XX, na.rm = TRUE)
  dt <- dt[count_d_non_miss_XX >= max_count]
  dt[, count_d_non_miss_XX := NULL]

  data.table::setorderv(dt, c(group, time))

  # Indicator: treatment changed vs. previous obs within group (0 for first obs)
  dt[, changed_XX := as.integer(get(treatment) != data.table::shift(get(treatment))), by = group]
  dt[is.na(changed_XX), changed_XX := 0L]

  # Spell counter within group (cumulative number of changes so far)
  dt[, spell_XX := cumsum(changed_XX), by = group]

  # Periods since the latest change within (group, spell)
  dt[, periods_since_change_XX := seq_len(.N) - 1L, by = c(group, "spell_XX")]
  dt[spell_XX == 0L, periods_since_change_XX := 0L]

  # Hit the reset threshold? Cumulate within group to label new sub-groups.
  dt[, reset_dummy_XX := as.integer(periods_since_change_XX == reset)]
  dt[, reset_count_XX := cumsum(reset_dummy_XX), by = group]

  dt[, new_group_XX := .GRP, by = c(group, "reset_count_XX")]

  # When no cluster was passed, fall back to the original group ID
  if (is.null(cluster)) {
    cluster_out <- "old_group_XX"
    dt[, old_group_XX := get(group)]
  } else {
    cluster_out <- cluster
  }

  dt[[group]] <- dt$new_group_XX

  dt[, c("changed_XX", "spell_XX", "periods_since_change_XX",
         "reset_dummy_XX", "reset_count_XX", "new_group_XX") := NULL]

  list(df = as.data.frame(dt), cluster = cluster_out)
}
