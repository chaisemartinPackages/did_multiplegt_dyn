#' Internal function of did_multiplegt_dyn
#' @param df df
#' @param G G
#' @param by by
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @noRd 
did_multiplegt_dyn_by_check <- function(
    df,
    G,
    by
) {
    df <- df %>% group_by(.data[[G]]) %>%
            mutate(sd_by = sd(.data[[by]], na.rm = TRUE))
    return(mean(df$sd_by) == 0)
}