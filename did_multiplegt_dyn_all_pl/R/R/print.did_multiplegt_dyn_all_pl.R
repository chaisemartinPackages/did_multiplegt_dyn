#' Custom print method for did_multiplegt_dyn_all_pl objects
#' @name print.did_multiplegt_dyn_all_pl
#' @description A customized printed display for did_multiplegt_dyn_all_pl output
#' @param x A did_multiplegt_dyn_all_pl object
#' @param ... Undocumented
#' @returns No return, custom print method for did_multiplegt_dyn_all_pl objects. Estimation tables are fetched from the object and displayed in the same style as the Stata did_multiplegt_dyn_all_pl command.
#' @export
print.did_multiplegt_dyn_all_pl <- function(x, ...) {
    class(x) <- "did_multiplegt_dyn"
    print(x)
    cat("\nMod: Two-stage estimation of placebos.\n")
}