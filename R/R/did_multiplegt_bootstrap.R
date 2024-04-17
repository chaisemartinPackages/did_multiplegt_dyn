#' Internal function of did_multiplegt_dyn - bootstrap se
#' @param df df
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
#' @param effects_equal effects_equal
#' @param save_results save_results
#' @param normalized normalized
#' @param predict_het predict_het
#' @param trends_lin trends_lin
#' @param less_conservative_se less_conservative_se
#' @param continuous continuous
#' @param bootstrap bootstrap
#' @param base base
#' @returns A list of the final results updated with the bootstrap standard errors
#' @noRd
did_multiplegt_bootstrap <- function(
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
  save_results, 
  normalized,
  predict_het,
  trends_lin,
  less_conservative_se,
  continuous,
  bootstrap, 
  base
){

    bresults_effects <- NULL
    bresults_ATE <- NULL
    bresults_placebo <- NULL

    n_effects <- nrow(base$Effects)
    bresults_effects <- matrix(NA, nrow = bootstrap, ncol = n_effects)
    if (isFALSE(trends_lin)) {
        bresults_ATE <- matrix(NA, nrow = bootstrap, ncol = 1)
    }
    if (placebo > 0) {
        n_placebo <- nrow(base$Placebos)
        bresults_placebo <- matrix(NA, nrow = bootstrap, ncol = n_placebo)
    }

    bs_group <- ifelse(!is.null(cluster), cluster, group)

    bdf <- data.frame(cbind(df[[bs_group]], 1:nrow(df)))
    colnames(bdf) <- c(bs_group, "id")
    xtset <- list()
    for (l in levels(factor(bdf[[bs_group]]))) {
        xtset[[l]] <- subset(bdf, bdf[[bs_group]] == l)[["id"]]
    }

    for (j in 1:bootstrap) {
        df_boot <- df[list_to_vec(xtset[sample(1:length(xtset), size = length(xtset), replace = TRUE)]), ]
        df_boot <- df_boot[order(df_boot[[group]], df_boot[[time]]), ]
        rownames(df_boot) <- 1:nrow(df_boot)
        suppressMessages({
        df_est <- did_multiplegt_main(df = df_boot, outcome = outcome, group =  group, time =  time, treatment = treatment, effects = effects, placebo = placebo, ci_level = ci_level,switchers = switchers, trends_nonparam = trends_nonparam, weight = weight, controls = controls, dont_drop_larger_lower = dont_drop_larger_lower, drop_if_d_miss_before_first_switch = drop_if_d_miss_before_first_switch, cluster = cluster, same_switchers = same_switchers, same_switchers_pl = same_switchers_pl, effects_equal = effects_equal, save_results = save_results, normalized = normalized, predict_het = predict_het, trends_lin = trends_lin, less_conservative_se = less_conservative_se, continuous = continuous)})

        res <- df_est$did_multiplegt_dyn
        for (i in 1:ncol(bresults_effects)) {
            if (i <= nrow(res$Effects)) {
                bresults_effects[j,i] <- res$Effects[i,1]
            }
        }

        if (!is.null(bresults_ATE)) {
            if (!is.null(res$ATE[1])) {
                bresults_ATE[j,1] <- res$ATE[1]
            }
        }

        if (!is.null(bresults_placebo)) {
            for (i in 1:ncol(bresults_placebo)) {
                if (i <= nrow(res$Placebos)) {
                    bresults_placebo[j,i] <- res$Placebos[i,1]
                }
            }
        }
        res <- df_est <- NULL
        progressBar(j, bootstrap)
    }

    ci_level <- ci_level / 100
    z_level <- qnorm(ci_level + (1 - ci_level)/2)

    for (i in 1:ncol(bresults_effects)) {
        if (!is.null(base$Effects[i,1])) {
            base$Effects[i,2] <- sd(bresults_effects[,i], na.rm = TRUE)
            base$Effects[i,3] <- base$Effects[i,1] - z_level * base$Effects[i,2]
            base$Effects[i,4] <- base$Effects[i,1] + z_level * base$Effects[i,2]
        }        
    }
    if (nrow(base$Effects) == 1) {
        class(base$Effects) <- "numeric"
    }

    if (!is.null(bresults_ATE)) {
        if (!is.null(base$ATE[1])) {
            base$ATE[2] <- sd(bresults_ATE, na.rm = TRUE)
            base$ATE[3] <- base$ATE[1] - z_level * base$ATE[2]
            base$ATE[4] <- base$ATE[1] + z_level * base$ATE[2]
        }
    }

    if (!is.null(bresults_placebo)) {
        for (i in 1:ncol(bresults_placebo)) {
            if (!is.null(base$Placebos[i,1])) {
                base$Placebos[i,2] <- sd(bresults_placebo[,i], na.rm = TRUE)
                base$Placebos[i,3] <- base$Placebos[i,1] - z_level * base$Placebos[i,2]
                base$Placebos[i,4] <- base$Placebos[i,1] + z_level * base$Placebos[i,2]
            }
        }
        if (nrow(base$Placebos) == 1) {
            class(base$Placebos) <- "numeric"
        }
    }
    return(base)    
}

#' Internal function to convert lists to vectors
#' @param lis A list
#' @returns A vector
#' @noRd
list_to_vec <- function(lis) {
    vec <- c()
    for (j in 1:length(lis)) {
        for (i in 1:length(lis[[j]])) {
            vec <- c(vec, lis[[j]][[i]])
        }
    }
    return(vec)
}
