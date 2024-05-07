#' Internal function of did_multiplegt_dyn by treatment path
#' @param df df
#' @param by_path by_path
#' @param outcome outcome
#' @param group group
#' @param time time
#' @param treatment treatment
#' @param effects effects
#' @param placebo placebo
#' @param ci_level ci_level
#' @param switchers switchers
#' @param trends_nonparam trends_nonparam
#' @param weight weight
#' @param controls controls
#' @param dont_drop_larger_lower dont_drop_larger_lower
#' @param drop_if_d_miss_before_first_switch drop_if_d_miss_before_first_switch
#' @param cluster cluster
#' @param same_switchers same_switchers
#' @param same_switchers_pl same_switchers_pl
#' @param only_never_switchers only_never_switchers
#' @param effects_equal effects_equal
#' @param save_results save_results
#' @param normalized normalized
#' @param predict_het predict_het
#' @param trends_lin trends_lin
#' @param less_conservative_se less_conservative_se
#' @param continuous continuous
#' @returns A dataframe with the by_path classifier
#' @noRd

did_multiplegt_by_path <- function(
    df,
    outcome, 
    group, 
    time, 
    treatment, 
    effects, 
    placebo, 
    ci_level, 
    switchers, 
    trends_nonparam, 
    weight, 
    controls, 
    dont_drop_larger_lower, 
    drop_if_d_miss_before_first_switch, 
    cluster, 
    same_switchers, 
    same_switchers_pl,
    effects_equal, 
    only_never_switchers,
    save_results, 
    normalized,
    predict_het,
    trends_lin,
    less_conservative_se,
    continuous,
    by_path
    ) {
    
    data <- did_multiplegt_main(df = df, outcome = outcome, group =  group, time =  time, treatment = treatment, effects = effects, placebo = placebo, ci_level = ci_level,switchers = switchers, trends_nonparam = trends_nonparam, weight = weight, controls = controls, dont_drop_larger_lower = dont_drop_larger_lower, drop_if_d_miss_before_first_switch = drop_if_d_miss_before_first_switch, cluster = cluster, same_switchers = same_switchers, same_switchers_pl = same_switchers_pl, only_never_switchers = only_never_switchers, effects_equal = effects_equal, save_results = save_results, normalized = normalized, predict_het = predict_het, trends_lin = trends_lin, less_conservative_se = less_conservative_se, continuous = continuous, data_only = TRUE)
    class(data$df) <- "data.frame"

    design_base <- did_multiplegt_dyn_design(data = data, design_opt = list(1, "console"), weight = weight, by = NULL, by_index = "_no_by", append = FALSE)

    l_XX <- data$l_XX
    T_max_XX <- data$T_max_XX
    path_index <- data$df[c("group", "time", "time_XX", "treatment_XX", "F_g_XX")]
    data <- NULL
    if (by_path == -1) {
        by_path <- nrow(design_base$design_mat)
    }
    design_set <- design_base$design_mat[1:min(by_path, nrow(design_base$design_mat)), ]
    if (by_path > nrow(design_base$design_mat)) {
        message(sprintf("You requested %.0f treatment paths, but there are only %.0f paths in your data. The program will continue with the latter number of treatment paths."))
    }
    path <- design_set[,3]
    for (j in 1:l_XX) {
        path <- paste0(path,",",design_set[,3+j])
    }
    
    for (i in 0:l_XX) {
        path_index[[paste0("D_Fg",i)]] <- ifelse(path_index$time_XX == path_index$F_g_XX - 1 + i, path_index$treatment_XX, NA)
    }
    path_index$treatment_XX <- NULL
    for (j in 0:l_XX) {
        path_index <- path_index %>% group_by(.data$group) %>% 
            mutate(!!paste0("D_fg",j) := mean(.data[[paste0("D_Fg",j)]],na.rm = TRUE)) %>% ungroup()
        path_index[[paste0("D_Fg",j)]] <- NULL
    }

    path_index$path <- path_index$D_fg0
    path_index$D_fg0 <- NULL
    for (j in 1:l_XX) {
        path_index$path <- ifelse(!is.na(path_index[[paste0("D_fg",j)]]), paste0(path_index$path,",",path_index[[paste0("D_fg",j)]]), path_index$path)
        path_index[[paste0("D_fg",j)]] <- NULL
    }
    path_index$yet_to_switch_XX <- as.numeric(path_index$time_XX < path_index$F_g_XX)
    path_index$time_XX <- path_index$F_g_XX <- NULL
    path_index$baseline_XX <- substr(path_index$path,1,1)
    names(path_index) <- c(group, time, "path_XX", "yet_to_switch_XX", "baseline_XX")

    df <- merge(df, path_index, by = c(group, time))
    df <- df[order(df[[group]], df[[time]]), ]
    data <- list(df = df, path = path)
    return(data)
}

#' Auxiliary function to check if "," is in var
#' @param str str
#' @returns bool
#' @noRd
count_comma <- function(str) {
    tot = 0
    for (k in 1:nchar(str)) {
        if (substr(str,k,k) == ",") {
            tot = tot + 1
        }
    }
    return(tot)
}