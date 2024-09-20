#' Checks whether the variable specified in the by option is time-variant
#' The program allows only time-variant variables in the by option
#' @param df df
#' @param group group
#' @param by by
#' @import data.table
#' @returns A logical value.
#' @noRd 
did_multiplegt_dyn_by_check <- function(
    df,
    group,
    by
) {
    sd_by <- NULL
    df <- data.table(df)
    df[, sd_by := sd(get(by), na.rm = TRUE), by = c(group)]
    return(mean(df$sd_by, na.rm = TRUE) == 0)
}