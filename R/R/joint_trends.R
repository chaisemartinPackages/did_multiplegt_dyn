#' Internal function of did_multiplegt_dyn that adds to general purpose conditioning (dplyr) vector-wise conditioning (a necessary feature for treends_nonparam)
#' When we use group_by(...) in piped statements, the argument of the option is a set of vars in .data. The trends_nonparam option is specified as an atomic char or a vector of chars. The need for this function arises in this second case. As it is not uncommon that we need to group by trends_nonparam and other variables, this function handles this aggregation, yielding the same df with a new variable that groups any set of variables with the trends_nonparam argument.
#' @param df df
#' @param var var
#' @param trends_nonparam trends_nonparam
#' @import dplyr
#' @importFrom magrittr %>%
#' @returns A dataframe.
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