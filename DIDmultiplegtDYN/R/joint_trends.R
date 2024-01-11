#' Internal function of did_multiplegt_dyn
#' @param df df
#' @param var var
#' @param trends_nonparam trends_nonparam
#' @import dplyr
#' @importFrom magrittr %>%
#' @noRd
joint_trends <- function(
    df, 
    var, 
    trends_nonparam
    ) {
    df <- df %>% dplyr::select(-any_of(c("joint_trends_XX")))
    joint_vars <- c(var, trends_nonparam)
    df$joint_trends_XX <- df[joint_vars]
    df
}