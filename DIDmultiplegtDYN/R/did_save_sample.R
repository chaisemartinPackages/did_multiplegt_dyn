#' Internal function of did_multiplegt_dyn
#' @param data data
#' @param Gn Gn
#' @param Tn Tn
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom data.table setnames
#' @noRd
did_save_sample <- function(
    data, 
    Gn, 
    Tn
    ) {
  df <- data$df
  suppressWarnings({
  df_save <- subset(df, !is.na(df$G) & !is.na(df$T))
  df_save <- df_save %>% dplyr::select(.data$G, .data$T, .data$S_g_XX, .data$switchers_tag_XX)
  df_save <- data.table::setnames(df_save, old = c("G", "T", "S_g_XX", "switchers_tag_XX"), new = c(Gn, Tn, "did_sample", "did_effect"))
  df_save$did_sample <- ifelse(df_save$did_sample == 0, -1, df_save$did_sample)
  df_save$did_sample <- ifelse(is.na(df_save$did_sample), 0, df_save$did_sample)
  df_save$did_sample <- factor(df_save$did_sample, levels = c(0,1,-1), labels = c("Control", "Switcher-in", "Switchers-out"))
  })
  return(df_save)
}

#' Adjustment of save_sample output in case of by option
#' @param obj A did_multiplegt_dyn object
#' @noRd
adj_save_sample <- function(obj) {
  saved_sample <- obj$by_level_1$save_sample
  obj$by_level_1$save_sample <- NULL
  if (length(obj$by_levels) > 1) {
    for (j in 2:length(obj$by_levels)) {
      saved_sample <- rbind(saved_sample, obj[[paste0("by_level_", j)]]$save_sample)
      obj[[paste0("by_level_", j)]]$save_sample <- NULL
    }
  }
  saved_sample <- saved_sample[order(saved_sample[[obj$args$G]], saved_sample[[obj$args$T]]), , drop = FALSE]
  rownames(saved_sample) <- NULL
  obj <- append(obj, list(saved_sample))
  return(obj)
}


