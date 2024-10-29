#' Option that shows the different treatment paths that switcher groups follow
#' @param data data
#' @param design_opt design_opt
#' @param weight weight
#' @param by by
#' @param by_index by_index
#' @param append append
#' @import data.table
#' @importFrom xlsx write.xlsx
#' @returns A list with the design option output.
#' @noRd
did_multiplegt_dyn_design <- function(
    data, 
    design_opt, 
    weight,
    by,
    by_index,
    append
    ) {

    # Inherited Globals #
    df <- data$df
    l_XX <- data$l_XX
    T_max_XX <- data$T_max_XX

    time_l_XX <- NULL
    group_XX <- NULL
    weight_XX <- NULL
    treatment_XX <- NULL
    N_XX <- NULL
    N_w_XX <- NULL
    treat_GRP <- NULL
    id_XX <- NULL
    in_table_XX <- NULL
    g_weight_XX <- NULL

	## Error message if the arguments in the option were specified wrong
  suppressWarnings({

	## Fetch the arguments 
  des_p <- as.numeric(design_opt[1])
  des_path <- design_opt[2]
  des_n <- l_XX
  des_per <- des_p * 100

	## keep periods up to ℓ periods after the first switch
  df$F_g_plus_n_XX <- df$F_g_XX + des_n - 1
  df$sel_XX <- df$time_XX >= df$F_g_XX - 1 & df$time_XX <= df$F_g_plus_n_XX
  df <- subset(df, df$time_XX >= df$F_g_XX - 1 & df$time_XX <= df$F_g_plus_n_XX)
  df <- df[order(df$group_XX, df$time_XX), ]
  df[, time_l_XX := seq_len(.N), by = group_XX]
  df <-  subset(df, select = c("group_XX", "time_l_XX", "weight_XX", "treatment_XX", "F_g_XX"))

	## Aggregate weights by group 
  if (!is.null(weight)) {
    df[, g_weight_XX := sum(weight_XX, na.rm = TRUE), by = group_XX]
  } else {
    df$g_weight_XX <- 1
  }
  df$weight_XX <- NULL

  max_time <- max(df$time_l_XX, na.rm = TRUE)
  treat_list <- c()
  treat_str <- ""
  for (i in 1:max_time) {
    df[, paste0("treatment_XX",i) := mean(treatment_XX[time_l_XX == i]), by = group_XX]
    treat_list <- c(treat_list, paste0("treatment_XX",i))
    treat_str <- paste0(treat_str,"treatment_XX",i,",")
  }
  treat_str <- substr(treat_str, 1, nchar(treat_str) - 1)
  df$time_l_XX <- df$treatment_XX <- NULL
  df <- unique(df)

	## Drop missing treatments 
  for (var in treat_list) {
    df <- subset(df, !is.na(df[[var]]))
  }

	## Creating varibale to store number of groups per treatment path and collapsing
  df$N_XX <- 1
  df$N_w_XX <- (df$g_weight_XX * df$N_XX) / sum(df$g_weight_XX, na.rm = TRUE)
  df$group_XX <- df$g_weight_XX <- NULL
  df[, N_XX := sum(N_XX, na.rm = TRUE), by = treat_list]
  df[, N_w_XX := sum(N_w_XX, na.rm = TRUE), by = treat_list]
  df$F_g_XX <- NULL
  df <- unique(df)
  tot_switch <- sum(df$N_XX, na.rm = TRUE)

	## Keep the observations amounting to p% of the detected treatment paths 
  df$neg_N_XX <- - df$N_XX
  df[, treat_GRP := .GRP, by = c(treat_list)]
  df <- df[order(df$neg_N_XX, df$treat_GRP), ]
  df$neg_N_XX <- df$treat_GRP <- NULL
  df$cum_sum_XX <- cumsum(df$N_w_XX)
  df$in_table_XX <- as.numeric(df$cum_sum_XX <= des_p)
  df <- df[order(df$in_table_XX, df$cum_sum_XX), ]
  df[, id_XX := seq_len(.N), by = in_table_XX]

	## Keep all observations up to the first exceeding the p%	
  df <- subset(df, df$in_table_XX == 1 | (df$in_table_XX == 0 & df$id_XX == 1))

	## Store the final % of groups included by the design option
  if (des_p < 1) {
    last_p <- 100 * min(df$cum_sum_XX[df$in_table_XX == 0])
  } else {
    last_p <- 100
  }
  df$neg_N_XX <- - df$N_XX
  df[, treat_GRP := .GRP, by = c(treat_list)]
  df <- df[order(df$neg_N_XX, df$treat_GRP), ]
  df <- subset(df, select = c("N_XX", "N_w_XX", treat_list))
  df$N_w_XX <- df$N_w_XX * 100

	## Prepare matrix for the output table
  coln <- c("N", "Share")
  rown <- c()
  desmat <- matrix(NA, nrow = dim(df)[1], ncol = 2 + 1 + l_XX)

	## Generate the column/row names and fill treatment path
  df <- data.frame(df)
  for (j in 1:(2 + 1 + l_XX)) {
    for (i in 1:dim(df)[1]) {
      if (j == 1) {
        rown <- c(rown, paste0("TreatPath",i))
      }
      desmat[i,j] <- as.numeric(df[i,j])
    }
    if (j > 2) {
      coln <- c(coln, paste0("\U2113","=",j - 2 - 1))
    }
  }
  df <- data.table(df)
  colnames(desmat) <- coln
  rownames(desmat) <- rown 
  
  desmat[, 2] <- noquote(sprintf("%s", format(round(desmat[,2], 2), big.mark=",", scientific=FALSE, trim=TRUE)))
  des_const <- c(l_XX, des_per, tot_switch, last_p)
  names(des_const) <- c("effects", "coverage_opt", "switchers", "detected_coverage")

  ## Save output as xlsx
  if (des_path != "console")  {
      by_add <- ""
      if (by_index != "_no_by") {
        by_add <- paste0(", ",abbreviate(by,5), "=", by_index)
      }
      write.xlsx(desmat, des_path, row.names = TRUE, col.names = TRUE, 
          sheetName = paste0("Design",by_add), append = append)
  }

  design <- list(
    design_path = des_path,
    design_mat = noquote(desmat),
    design_const = des_const
  )
  return(design)
  })
}