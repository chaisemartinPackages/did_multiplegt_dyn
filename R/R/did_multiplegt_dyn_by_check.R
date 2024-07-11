#' Checks whether the variable specified in the by option is time-variant
#' The program allows only time-variant variables in the by option
#' @param df df
#' @param group group
#' @param by by
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @returns A logical value.
#' @noRd 
did_multiplegt_dyn_by_check <- function(
    df,
    group,
    by
) {
    df <- df %>% group_by(.data[[group]]) %>%
            mutate(sd_by = sd(.data[[by]], na.rm = TRUE))
    return(mean(df$sd_by, na.rm = TRUE) == 0)
}