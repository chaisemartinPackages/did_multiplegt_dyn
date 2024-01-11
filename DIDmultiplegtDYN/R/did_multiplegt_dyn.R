#' Core function for did_multiplegt_dyn
#' @md 
#' @param df Dataset
#' @param Y Outcome
#' @param G Group
#' @param T Time
#' @param D Treatment
#' @param effects N. of effects to be computed
#' @param placebo N. of placebos
#' @param ci_level Confidence level
#' @param switchers Switcher option
#' @param trends_nonparam Trends_nonparam
#' @param weight Weights
#' @param controls controls
#' @param dont_drop_larger_lower dont_drop_larger_lower
#' @param save_sample save_sample
#' @param drop_if_d_miss_before_first_switch drop_if_d_miss_before_first_switch
#' @param cluster cluster
#' @param same_switchers same_switchers
#' @param same_switchers_pl same_switchers_pl
#' @param effects_equal effects_equal
#' @param save_results save_results
#' @param normalized normalized
#' @param design design
#' @param date_first_switch date_first_switch
#' @param normalized_weights normalized_weights
#' @param graph_off graph_off
#' @param by by
#' @export 
did_multiplegt_dyn <- function(
    df, 
    Y, 
    G, 
    T, 
    D, 
    effects = 1, 
    placebo = 0, 
    ci_level = 95, 
    switchers = "", 
    trends_nonparam = NULL, 
    weight = NULL, 
    controls = NULL, 
    dont_drop_larger_lower = FALSE, 
    save_sample = FALSE, 
    drop_if_d_miss_before_first_switch = FALSE, 
    cluster = NULL, 
    same_switchers = FALSE, 
    same_switchers_pl = FALSE, 
    effects_equal = FALSE, 
    save_results = NULL, 
    normalized = FALSE, 
    design = NULL, 
    date_first_switch = NULL, 
    normalized_weights = NULL, 
    graph_off = FALSE,
    by = NULL
    ) { 
  
  args <- list()
  for (v in names(formals(did_multiplegt_dyn))) {
    if (v != "df") {
      args[[v]] <- get(v)
    }
  }
  did_multiplegt_dyn <- list(args)
  f_names <- c("args")

  by_levels <- c("_no_by")
  if (!is.null(by)) {
    if (!did_multiplegt_dyn_by_check(df, G, by)) {
      cat(sprintf("The variable %s is time-variant. The command will ignore the by option", by));cat("\n");
    } else {
      by_levels <- levels(factor(df[[by]]))
      did_multiplegt_dyn <- append(did_multiplegt_dyn, list(by_levels))
      f_names <- c(f_names, "by_levels")
    }
  }

  append_design <- FALSE
  append_dfs <- FALSE
  for (b in 1:length(by_levels)) {

    if (by_levels[b] != "_no_by") {
        df_main <- subset(df, df[[by]] == by_levels[b])
    } else {
      df_main <- df
    }

    df_est <- did_multiplegt_main(df_main, Y, G, T, D, 
    effects, placebo, ci_level, switchers, trends_nonparam, 
    weight, controls, dont_drop_larger_lower, 
    drop_if_d_miss_before_first_switch, cluster, 
    same_switchers, same_switchers_pl, effects_equal, 
    save_results, normalized)

    temp_obj <- list(df_est$did_multiplegt_dyn)
    names(temp_obj)[length(temp_obj)] <- "results"

    if (!is.null(design)) {
      temp_obj <- append(temp_obj, list(did_multiplegt_dyn_design(df_est, design, weight, by, by_levels[b], append_design)))
      names(temp_obj)[length(temp_obj)] <- "design"
      append_design <- TRUE
    }

    if (!is.null(date_first_switch)) {
      temp_obj <- append(temp_obj, list(did_multiplegt_dyn_dfs(df_est, date_first_switch, by, by_levels[b], append_dfs)))
      names(temp_obj)[length(temp_obj)] <- "date_first_switch"
      append_dfs <- TRUE
    }

    if (!is.null(normalized_weights)) {
      temp_obj <- append(temp_obj, 
          list(did_multiplegt_dyn_normweights(df_est, normalized, normalized_weights, same_switchers)))
      names(temp_obj)[length(temp_obj)] <- "normalized_weights"
    }

    temp_obj <- append(temp_obj, list(did_multiplegt_dyn_graph(df_est)))
    names(temp_obj)[length(temp_obj)] <- "plot"

    if (isTRUE(save_sample)) {
      df_save_XX <- did_save_sample(df_est, G, T)
      df_m <- merge(df, df_save_XX, by = c(G, T))
      temp_obj <- append(temp_obj, list(df_m))
      names(temp_obj)[length(temp_obj)] <- "save_sample"
    }

    if (by_levels[b] == "_no_by") {
      did_multiplegt_dyn <- c(did_multiplegt_dyn, temp_obj)
      f_names <- c(f_names, names(temp_obj))
    } else {
      temp_obj <- append(temp_obj, list(by_levels[b]))
      names(temp_obj)[length(temp_obj)] <- "level"
      did_multiplegt_dyn <- append(did_multiplegt_dyn, list(temp_obj))
      f_names <- c(f_names, paste0("by_level_",b))
    }
  }

  names(did_multiplegt_dyn) <- f_names

  if (!is.null(by)) {
    if (isTRUE(save_sample)) {
      did_multiplegt_dyn <- adj_save_sample(did_multiplegt_dyn)
      names(did_multiplegt_dyn)[length(did_multiplegt_dyn)] <- "save_sample"
    }
    did_multiplegt_dyn <- append(did_multiplegt_dyn, list(combine_plot(did_multiplegt_dyn)))
    names(did_multiplegt_dyn)[length(did_multiplegt_dyn)] <- "plot"
  }

  if (isFALSE(graph_off)) {
    print(did_multiplegt_dyn$plot)
  }
  
  class(did_multiplegt_dyn) <- c(class(did_multiplegt_dyn), "did_multiplegt_dyn")
  return(did_multiplegt_dyn)
}
