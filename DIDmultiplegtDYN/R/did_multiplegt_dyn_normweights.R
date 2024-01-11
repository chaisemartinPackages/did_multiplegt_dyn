#' Internal function of did_multiplegt_dyn
#' @param data data
#' @param normalized normalized
#' @param normopt normopt
#' @param same_switchers same_switchers
#' @import dplyr
#' @importFrom magrittr %>%
#' @noRd
did_multiplegt_dyn_normweights <- function(
  data, 
  normalized, 
  normopt, 
  same_switchers
  ) {
    # Inherited Globals #
    df <- data$df
    l_XX <- data$l_XX
    for (v in names(data$delta)) {
      assign(v, data$delta[[v]])
    }

  suppressWarnings({
  if (normalized == FALSE) {
    stop("normalized option required to compute normalized_weights")
  }

  lag_d <- normopt
  if (!(lag_d %in% c("by_k", "by_calendar"))) {
    stop("First argument of normalized_weights incorrectly specified. Normalized_weights requires by_k or by_calendar as arguments.")
  }

  weight_mat <- matrix(NA, nrow = l_XX, ncol = l_XX) 
  coln <- c()
  rown <- c()
  for (i in 1:l_XX) {
    coln <- c(coln, paste0("\U2113","=",i))
    for (k in 0:(i-1)) {
      if (lag_d == "by_k") {
        row <- k + 1
      } else {
        row <- i - k # Visualization by F_g - 1 + (l - k)
      }

      df[paste0("delta_",i,"_",k)] <- ifelse(df$time_XX == df$F_g_XX - 1 + i - k & df$F_g_XX - 1 + i <= df$T_g_XX, abs(df$treatment_XX - df$d_sq_XX), NA)

      if (same_switchers == TRUE) {
        df[[paste0("delta_",i,"_",k)]] <- ifelse(df$F_g_XX - 1 + l_XX > df$T_g_XX, 0, df[[paste0("delta_",i,"_",k)]])        
      }

      weight_mat[row, i] <- (sum(df[[paste0("delta_",i,"_",k)]], na.rm = TRUE) / get(paste0("delta_D_",i,"_global_XX"))) / data$mat_res_XX[i,6]
      #get(paste0("N_switchers_effect_",i,"_XX"))
    }
  }
  if (lag_d == "by_k") {
    for (j in 1:l_XX) {
      rown <- c(rown, paste0("k=",j-1))
    }     
  } else {
    rown <- c(rown, "D_Fg")
    for (j in 2:l_XX) {
      rown <- c(rown, paste0("D_Fg+",j-1))
    }
  }

  mat_total <- weight_mat
  mat_total[is.na(mat_total)] <- 0
  total <- matrix(1,nrow=1,ncol=l_XX) %*% mat_total
  weight_mat <- rbind(weight_mat, total)
  rownames(weight_mat) <- c(rown, "Total")
  colnames(weight_mat) <- coln
  by_opt_lag <- gsub("by_", "", lag_d)
  weight_mat[ , ] <- sprintf("%s", format(round(weight_mat[ , ], 3), big.mark=",", scientific=FALSE, trim=TRUE))

  return(list(by_opt_lag = by_opt_lag, norm_weight_mat = noquote(weight_mat)))
  })

}