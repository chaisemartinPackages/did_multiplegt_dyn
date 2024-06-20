#' Inner program to fill panel structure and retrieve placebos
#' @param df df
#' @param outcome outcome
#' @param group group
#' @param time time
#' @param treatment treatment
#' @param add add
#' @import dplyr
#' @noRd
add_periods <- function(
    df, 
    outcome, 
    group, 
    time, 
    treatment, 
    add
) {

    df <- df[order(df[[group]], df[[time]]), ]
    df <- df %>% group_by(.data[[group]]) %>% mutate(G_AP = cur_group_id())
    df <- df %>% group_by(.data[[time]]) %>% mutate(T_AP = cur_group_id()) 
    max_T <- max(df$T_AP, na.rm = T)
    max_G <- max(df$G_AP, na.rm = T)

    df <- df [order(df$T_AP, df$G_AP), ]
    for (j in 1:max_G) {
        N_temp <- nrow(df)
        to_add <- as.data.frame(matrix(NA,nrow=add,ncol=ncol(df)))
        colnames(to_add) <- colnames(df)
        df <- rbind(df,to_add)
        df[[outcome]][is.na(df[[outcome]])] <- 0
        df[[group]][is.na(df[[group]])] <- j
        df[[time]][is.na(df[[time]])] <- max_T + (1:add)
        df[[treatment]][is.na(df[[treatment]])] <- df[[treatment]][j]
    }
    df$T_AP <- df$G_AP <- NULL
    df <- df[order(df[[group]], df[[time]]), ]
    return(df)
}