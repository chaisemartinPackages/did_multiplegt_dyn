#' Core function of DIDmultiplegtDYNallPL
#' @md
#' @description Auxiliary function of DIDmultiplegtDYN to retrieve more placebos. 
#' @param df (dataframe) the estimation dataset. It has to be a <ins>balanced panel<ins> with a least <ins>one never-switcher<ins> per baseline treatment.
#' @param outcome (char) is the outcome variable. 
#' @param group (char) is the group variable. 
#' @param time (char) is the time period variable.
#' @param treatment (char) is the treatment variable.
#' @param effects (int) gives the number of event-study effects to be estimated.
#' @param placebo (int) gives the number of placebo estimators to be computed. 
#' @param switchers (char in ("", "in", "out")) one may be interested in estimating separately the treatment effect of switchers-in, whose average treatment after they switch is larger than their baseline treatment, and of switchers-out, whose average treatment after they switch is lower than their baseline treatment. In that case, one should run the command first with the \code{switchers = "in"} option, and then with the \code{switchers = "out"} option.
#' @param only_never_switchers (logical) if this option is specified, the command estimates the event-study effects using only never-switchers as control units.
#' @section Overview:
#' The following [vignette](https://github.com/chaisemartinPackages/did_multiplegt_dyn/blob/main/vignettes/vignette_3.md) provides an overview of the motivation and the usage of this function.
#' @import DIDmultiplegtDYN
#' @import dplyr
#' @examples
#' GG <- 10; TT <- 4;
#' set.seed(0)
#' df <- as.data.frame(matrix(nrow =GG * TT, ncol = 0))
#' df$G <- floor((1:nrow(df)-1)/TT) + 1
#' df$T <- ((1:nrow(df)-1) %% TT) + 1
#' df <- df[order(df$G, df$T), ]
#' df$D <- ifelse(df$T == TT & runif(n = nrow(df)) > 0.5, 1, 0)
#' df$Y <- runif(n = nrow(df)) * (1 + df$D)
#' did <- did_multiplegt_dyn_all_pl(df, "Y", "G", "T", "D", effects = 1, placebo = 2)
#' @returns The same output as did_multiplegt_dyn.
#' @export
did_multiplegt_dyn_all_pl <- function(
    df,
    outcome,
    group,
    time,
    treatment,
    effects = 1,
    placebo = 0,
    switchers = "",
    only_never_switchers = FALSE
) {

    did1 <- did_multiplegt_dyn(
        df = df, 
        outcome = outcome, 
        group = group,
        time = time,
        treatment = treatment,
        effects = effects,
        switchers = switchers,
        only_never_switchers = only_never_switchers,
        graph_off = T
        )

    data <- add_periods(
        df = df, 
        outcome = outcome, 
        group = group, 
        time = time, 
        treatment = treatment, 
        add = did1$results$max_pl_gap
    )

    did2 <- did_multiplegt_dyn(
        df = data, 
        outcome = outcome, 
        group = group,
        time = time,
        treatment = treatment,
        effects =  max(effects, placebo),
        placebo = placebo,
        switchers = switchers,
        only_never_switchers = only_never_switchers,
        graph_off = T
    )

    did2$args$effects <- effects
    did2$results$N_Effects <- nrow(did1$results$Effects)
    did2$results$Effects <- did1$results$Effects
    did2$results$ATE <- did1$results$ATE
    did2$results$delta_D_avg_total <- did1$results$delta_D_avg_total
    did2$max_pl <- did2$max_pl_gap <- NULL
    did2$coef <- NULL
    did2$plot <- NULL
    
    class(did2) <- "did_multiplegt_dyn_all_pl"
    return(did2)
}