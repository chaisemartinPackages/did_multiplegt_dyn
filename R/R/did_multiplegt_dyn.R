packages <- c("plm", "data.table", "dplyr", "logr", "estimatr", "huxtable", "magrittr", "fixest", "stringr", "matlib", "xlsx")

ds_libraries <- function(packages){
  for(package in packages){
    if(!require(package, character.only = TRUE)){
      install.packages(package, dependencies = TRUE,repos = "http://cran.us.r-project.org")
    }
    #Load package
    library(package, character.only = TRUE)
  }
}
ds_libraries(packages)

joint_trends <- function(df, var, trends_nonparam) {
    df <- df %>% dplyr::select(-any_of(c("joint_trends_XX")))
    joint_vars <- c(var, trends_nonparam)
    df$joint_trends_XX <- df[joint_vars]
    df
}

did_multiplegt_main <- function(df, Y, G, T, D, effects, placebo, ci_level, switchers, trends_nonparam, weight, controls, dont_drop_larger_lower, drop_if_d_miss_before_first_switch, cluster, same_switchers, same_switchers_pl,effects_equal, save_results, normalized) {

  suppressWarnings({
  # Renaming the variables in the dataset 
  original_names = c(c(Y, G, T, D), trends_nonparam, weight, controls, cluster)
  df <- data.frame(df)
  df <- df %>% select_at(original_names)
  df <- data.table::setnames(df, old = c(Y, G, T, D), new = c("Y", "G", "T", "D"))
  
  if (!is.null(trends_nonparam)) {
    df$trends_nonparam_XX <- df[trends_nonparam]
  }

  if (same_switchers == FALSE & same_switchers_pl == TRUE) {
    stop("The same_switchers_pl option only works if same_switchers is specified as well!")
  }

  # Patching the cluster variable
  if (!is.null(cluster)) {
    if (paste0(cluster) == paste0(G)) {
      cluster <- NULL
    }
    df$cluster_XX <- df[[cluster]]
  }

  # Dealing with missing data
  df <- df %>% filter(!is.na(.data$G) & !is.na(.data$T)) 
  if (!is.null(controls)) {
    for (var in controls) {
      df <- subset(df, !is.na(df[var]))
    }
  }

  if (!is.null(cluster)) {
    df <- subset(df, !is.na(df$cluster_XX))
  }
  df <- df %>% group_by(.data$G) %>% mutate(mean_D = mean(.data$D, na.rm = TRUE))
  df <- df %>% group_by(.data$G) %>% mutate(mean_Y = mean(.data$Y, na.rm = TRUE))
  df <- df %>% filter(!is.na(.data$mean_Y) & !is.na(.data$mean_D))  %>% 
      dplyr::select(-.data$mean_Y, -.data$mean_D) 

  # Creating the weight variable 
  if (is.null(weight)) {
    df$weight_XX <- 1
  } 
  else {
    df$weight_XX <- df[[weight]]
  }
  df$weight_XX <- ifelse(is.na(df$weight_XX), 0, df$weight_XX)

  # Checking if the data has to be collapsed
  df$counter_temp <- 1
  df <- df %>% group_by(.data$G, .data$T) %>% 
      mutate(counter = sum(.data$counter_temp))
  aggregated_data <- max(df$counter) == 1
  df <- df %>% dplyr::select(-.data$counter, -.data$counter_temp) 

  if (aggregated_data != 1) {
    df$weight_XX <- ifelse(is.na(df$D), 0, df$weight_XX)
    df <- df %>%  group_by(.data$G, .data$T) %>%
    summarise_at(vars(D, Y, trends_nonparam, weight, controls, cluster_XX, weight_XX), funs(weighted.mean(., w=weight_XX)))
    df$trends_nonparam_XX <- df[trends_nonparam]
  }

  # Generate factorized versions of Y, G, T and D
  df$outcome_XX <- df$Y
  df <- df %>% group_by(.data$G) %>% mutate(group_XX = cur_group_id()) %>% ungroup()
  df <- df %>% group_by(.data$T) %>% mutate(time_XX = cur_group_id()) %>% ungroup()
  df$treatment_XX <- df$D

  # Declaring that the dataset is a panel
  df <- pdata.frame(df, index = c("group_XX", "time_XX")) 
  df$time_XX <- as.numeric(as.character(df$time_XX))
  df$group_XX <- as.numeric(as.character(df$group_XX))

  # Preparation for the next blocks
  df$time_d_nonmiss_XX <- ifelse(!is.na(df$treatment_XX), df$time_XX, NA)
  df <- df %>% group_by(.data$group_XX) %>% 
      mutate(min_time_d_nonmiss_XX = min(.data$time_d_nonmiss_XX, na.rm = TRUE))
  df <- df %>% group_by(.data$group_XX) %>% 
      mutate(max_time_d_nonmiss_XX = max(.data$time_d_nonmiss_XX, na.rm = TRUE))

  df$time_y_nonmiss_XX <- ifelse(!is.na(df$outcome_XX), df$time_XX, NA)
  df <- df %>% group_by(.data$group_XX) %>%
     mutate(min_time_y_nonmiss_XX = min(.data$time_y_nonmiss_XX, na.rm = TRUE))

  df$time_d_miss_XX <- ifelse(is.na(df$treatment_XX) & df$time_XX >= df$min_time_y_nonmiss_XX,df$time_XX, NA)
  df <- df %>% group_by(.data$group_XX) %>%
       mutate(min_time_d_miss_aft_ynm_XX = min(.data$time_d_miss_XX, na.rm = TRUE))
  df <- df %>% dplyr::select(-.data$time_d_nonmiss_XX, -.data$time_y_nonmiss_XX, 
      -.data$time_d_miss_XX)

  # Status quo treatment D_{g,1}
  df$d_sq_XX <- ifelse(df$time_XX == df$min_time_d_nonmiss_XX,df$treatment_XX,NA)
  df <- df  %>% group_by(.data$group_XX)  %>% 
      mutate(d_sq_XX = mean(d_sq_XX, na.rm = TRUE)) 
  df$diff_from_sq_XX <- df$treatment_XX - df$d_sq_XX 

  df <- df[order(df$group_XX, df$time_XX), ]
  T_XX <- max(df$time_XX, na.rm = TRUE)

  if (!(dont_drop_larger_lower == TRUE)) {
    df <- df[order(df$group_XX, df$time_XX), ]

    df$ever_strict_increase_XX <- as.numeric(df$diff_from_sq_XX > 0 & !is.na(df$treatment_XX))
    df$ever_strict_decrease_XX <- as.numeric(df$diff_from_sq_XX < 0 & !is.na(df$treatment_XX))

    for (t in 1:(T_XX)) {
      for (var in c("group_XX", "ever_strict_increase_XX", "ever_strict_decrease_XX")) {
          df[paste0(var,"_lag")] <- as.numeric(lag(df[[var]], 1))
      }
      for (var in c("ever_strict_increase_XX", "ever_strict_decrease_XX")) {
        df[[var]] <- ifelse(df[[paste0(var,"_lag")]] == 1 & !is.na(df[[paste0(var,"_lag")]]) & df$group_XX == df$group_XX_lag, 1, df[[var]]) 
      }
      df <- df %>% select(-c("group_XX_lag","ever_strict_increase_XX_lag", "ever_strict_decrease_XX_lag" ))
    }
      df <- subset(df, !(df$ever_strict_increase_XX == 1 & df$ever_strict_decrease_XX == 1))
  }

  G_XX <- max(df$group_XX, na.rm = TRUE)


  # Ever changed treatment 
  df$ever_change_d_XX <- abs(df$diff_from_sq_XX) > 0 & !is.na(df$treatment_XX) 
  # Differently from Stata, R does not allow progressive replacement with [_n]
  for (i in 2:T_XX) {
    df$ever_change_d_XX[shift(df$ever_change_d_XX) == 1 & df$group_XX == shift(df$group_XX) & df$time_XX == i] <- 1
  }

  # Date of the first treatment change
  df <- df[order(df$group_XX, df$time_XX), ]
  df$temp_F_g_XX <- ifelse(df$ever_change_d_XX == 1 & shift(df$ever_change_d_XX) == 0,
     df$time_XX, 0)
  df <- df  %>% group_by(.data$group_XX)  %>% mutate(F_g_XX = max(.data$temp_F_g_XX, na.rm = TRUE))  %>% 
      dplyr::select(-.data$temp_F_g_XX)

  # Create a new value with integer levels of d_sq_XX
  df <- df %>% group_by(.data$d_sq_XX) %>% mutate(d_sq_int_XX = cur_group_id()) %>% ungroup()
  df$d_sq_int_XX <- as.numeric(as.character(df$d_sq_int_XX))

  # Dropping values of baseline treatment such that there is no variance in F_g within
  df <- joint_trends(df, "d_sq_XX", trends_nonparam)
  df <- df %>% group_by(.data$joint_trends_XX) %>% mutate(var_F_g_XX = sd(F_g_XX)) %>% ungroup()
  df <- subset(df, df$var_F_g_XX > 0)
  df <- df %>% dplyr::select(-.data$var_F_g_XX)

  if (nrow(df) == 0) {
      stop("No treatment effect can be estimated. This is because Assumption 1 in de Chaisemartin & D'Haultfoeuille (2023) is not satisfied in the data used for estimation, given the options requested. If this is caused by your baseline treatement being continuous you can try using the option continuous() which allows for a continous period-one treatement.")
  }
  
  # For each value of d_sq_XX, we drop time periods such that we do not have any control with the same baseline treatment afterwards
  df$never_change_d_XX <- 1 - df$ever_change_d_XX 
  df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
  df <- df %>% group_by(.data$joint_trends_XX) %>% mutate(controls_time_XX = max(never_change_d_XX))
  df <- subset(df, df$controls_time_XX > 0)

  # Computing t_min, T_max and adjusting F_g
  t_min_XX <- min(df$time_XX)
  T_max_XX <- max(df$time_XX)
  df$F_g_XX[df$F_g_XX == 0] <- T_max_XX + 1


  ### MISSING TREATMENT CONVENTIONS ### 
  # Dealing with missing treatments: most conservative option
  if (drop_if_d_miss_before_first_switch == TRUE) {
    df$outcome_XX <- ifelse(
      df$min_time_d_miss_aft_ynm_XX < df$F_g_XX & df$time_XX >= df$min_time_d_miss_aft_ynm_XX, NA, df$outcome_XX)
  }

  # Dealing with missing treatments: most liberal option
  # Let FD_g and LD_g respectively denote the first and last period where a group's treatment is non missing. Let FY_g denote the first period where a group's outcome is non missing.

  # For groups that experience at least one treatment change, let LDBF_g denote the last date before F_g where g's treatment is non missing. We have FD_g<=LDBF_g<F_g<=LD_g, and we will deal with missing treatments depending on when they occur with respect to those four dates. 

  df$last_obs_D_bef_switch_t_XX <- ifelse(df$time_XX < df$F_g_XX & !is.na(df$treatment_XX), df$time_XX, NA)
  df <- df %>% group_by(.data$group_XX) %>% 
  mutate(last_obs_D_bef_switch_XX = max(.data$last_obs_D_bef_switch_t_XX, na.rm = TRUE))

  # For groups that do not experience a treatment change, we just have FD_g<=LD_g, and we will deal with missing treatments depending on when they occur with respect to those two dates.

  # For t<FD_g, by default we are going to consider that g joins the panel at FD_g: any non-missing outcome before FD_g replaced as missing, but all non-missing outcomes after FD_g are kept. For groups such that FY_g<FD_g, this is a "liberal" convention: those groups exist before FD_g, so one could argue that their status quo treatment is missing and they should be dropped from the analysis. We give the user the option to do that, with drop_if_d_miss_before_first_switch option

  df$outcome_XX[df$time_XX < df$min_time_d_nonmiss_XX] <- NA

  # For groups that experience a treatment change, if D_gt missing at FD_g<t<LDBF_g, we replace their missing treatment by their status-quo treatment. Again, this is a liberal convention, so we give the user the option to not use those observations, with drop_if_d_miss_before_first_switch option.

  df$treatment_XX <- ifelse(df$F_g_XX < T_max_XX + 1 & is.na(df$treatment_XX) & df$time_XX < df$last_obs_D_bef_switch_XX & df$time_XX > df$min_time_d_nonmiss_XX, df$d_sq_XX, df$treatment_XX)

  # For groups that experience a treatment change, if D_gt missing at LDBF_g<t<F_g (equivalent to LDBF_g<F_g-1), we cannot know the exact date when their treatment has changed, even in a binary and staggered design. Therefore, we set their outcomes at missing starting at LDBF_g+1. We also redefine their F_g as T+1 because they are effectively control groups. We also define the trunc_control_XX as LDBF_g+1 for them, because they can only be used as controls till that date.

  df$outcome_XX <- ifelse(df$F_g_XX < T_max_XX + 1 & df$time_XX > df$last_obs_D_bef_switch_XX & df$last_obs_D_bef_switch_XX < df$F_g_XX - 1, NA, df$outcome_XX)
  df$trunc_control_XX <- ifelse(df$F_g_XX < T_max_XX + 1 & df$last_obs_D_bef_switch_XX < df$F_g_XX - 1, df$last_obs_D_bef_switch_XX + 1, NA)
  df$F_g_XX[df$F_g_XX < T_max_XX + 1 & df$last_obs_D_bef_switch_XX < df$F_g_XX - 1] <- T_max_XX + 1

  # For groups that experience a treatment change, if D_gt missing at F_g<t, we replace their missing treatment by D(g,F_g). This is again a liberal convention, but it is innocuous for the reduced-form parameters DID_l, so we do not give the user the option to overrule it (Note that overruling it could make the same_switchers option fail)

  df$d_F_g_temp_XX <- ifelse(df$time_XX == df$F_g_XX, df$treatment_XX, NA)
  df <- df %>% group_by(.data$group_XX) %>% mutate(d_F_g_XX = mean(.data$d_F_g_temp_XX, na.rm = TRUE))
  df$treatment_XX <- ifelse(df$F_g_XX < T_max_XX + 1 & is.na(df$treatment_XX) & df$time_XX > df$F_g_XX & df$last_obs_D_bef_switch_XX == df$F_g_XX - 1, df$d_F_g_XX, df$treatment_XX)

  # *For groups that do not experience a treatment change, if D_gt missing at FD_g<t<LD_g, we replace their missing treatment by D_g1. This is again a liberal convention, so we give the user the option to not use those observations, wi1th drop_if_d_miss_before_first_switch option.

  df$treatment_XX <- ifelse(df$F_g_XX == T_max_XX + 1 & is.na(df$treatment_XX) & df$time_XX > df$min_time_d_nonmiss_XX & df$time_XX < df$max_time_d_nonmiss_XX, df$d_sq_XX, df$treatment_XX)

  # For groups that do not experience a treatment change, we replace all their outcomes by missing at t>LD_g. Even in a binary and staggered design, we cannot infer their treatment at t>LD_g.

  df$outcome_XX[df$F_g_XX == T_max_XX + 1 &df$time_XX > df$max_time_d_nonmiss_XX] <- NA
  df$trunc_control_XX <- ifelse(df$F_g_XX == T_max_XX + 1, df$max_time_d_nonmiss_XX + 1, df$trunc_control_XX)

  ### END OF MISSING TREATMENT CONVENTIONS BLOCK ###

  # Balancing the panel
  df <- df %>% select(-any_of(c("joint_trends_XX", "trends_nonparam_XX")))

  df <- pdata.frame(df, index = c("group_XX", "time_XX")) 
  df <- make.pbalanced(df, balance.type = "fill")
  df$time_XX <- as.numeric(as.character(df$time_XX))
  df$group_XX <- as.numeric(as.character(df$group_XX))
  df <- df %>% group_by(.data$group_XX) %>% mutate(d_sq_XX = mean(.data$d_sq_XX, na.rm = TRUE))

  # Defining N_gt, th weights of each (g,t) cell
  df$N_gt_XX <- 1
  df$N_gt_XX <- ifelse(is.na(df$outcome_XX) | is.na(df$treatment_XX), 0, df$weight_XX * df$N_gt_XX)

  # Defining T_g, the last period whn there is still a not-yet-switcher group with the same treatment as group g
  df <- df %>% mutate(F_g_trunc_XX = min(.data$F_g_XX, .data$trunc_control_XX))
  df$F_g_trunc_XX <- ifelse(is.na(df$trunc_control_XX), df$F_g_XX, df$F_g_trunc_XX)

  df <- joint_trends(df, "d_sq_XX", trends_nonparam)
  df <- df %>% group_by(.data$joint_trends_XX) %>% 
        mutate(T_g_XX = max(.data$F_g_trunc_XX, na.rm = TRUE)) %>% ungroup()
  df$T_g_XX <- df$T_g_XX - 1

  # Defining S_g, an indicator variable for groups whose average post switch treatment value is larger than their initial value of treatment. They will be considered switchers in. If S_g==0, that means the group is a switcher out. For never-switchers, S_g will be undefined.

  df$treatment_XX_v1 <- ifelse(df$time_XX >= df$F_g_XX & df$time_XX <= df$T_g_XX, df$treatment_XX, NA)
  df <- df %>% group_by(.data$group_XX) %>% 
      mutate(avg_post_switch_treat_XX_temp = sum(.data$treatment_XX_v1, na.rm = TRUE))
  df <- df %>% dplyr::select(-.data$treatment_XX_v1)

  df$count_time_post_switch_XX_temp <- (df$time_XX >= df$F_g_XX & df$time_XX <= df$T_g_XX & !is.na( df$treatment_XX))
 
  df <- df %>% group_by(.data$group_XX) %>%
      mutate(count_time_post_switch_XX = sum(.data$count_time_post_switch_XX_temp, na.rm = TRUE)) 
  
  df$avg_post_switch_treat_XX_temp <- df$avg_post_switch_treat_XX_temp / df$count_time_post_switch_XX

  df <- df %>% group_by(.data$group_XX) %>% 
      mutate(avg_post_switch_treat_XX = mean(.data$avg_post_switch_treat_XX_temp, na.rm = TRUE)) %>%
      dplyr::select(-.data$avg_post_switch_treat_XX_temp)

  # When a group is a switching group, but its average post-treatment treatment value is exactly equal to its baseline treatment, we cannnot classify it as a swicher in or a switcher out, but it is not a control either. As such, we drop it from the estimation. Those groups are referred to as no-first-stage-switchers.

  df <- subset(df, !(df$avg_post_switch_treat_XX == df$d_sq_XX & df$F_g_XX != df$T_g_XX + 1))
  df$S_g_XX <- as.numeric(df$avg_post_switch_treat_XX > df$d_sq_XX)
  df$S_g_XX <- ifelse(df$F_g_XX != T_max_XX + 1, df$S_g_XX, NA)

  # Creating the variable L_g_XX = T_g_XX - F_g_XX
  df$L_g_XX <- df$T_g_XX - df$F_g_XX + 1

  # Creating the equivalent variable for placebos
  if (placebo > 0) {
    df <- df %>% rowwise() %>% 
        mutate(L_g_placebo_XX = min(.data$L_g_XX[.data$F_g_XX >= 3], 
          .data$F_g_XX[.data$F_g_XX >= 3] - 2)) 
    df$L_g_placebo_XX <- ifelse(df$L_g_placebo_XX == Inf, NA, df$L_g_placebo_XX)
  }

  # Flagging first observation of each group_XX
  df <- df[order(df$group_XX, df$time_XX), ]
  df <- df %>% group_by(.data$group_XX) %>% 
      mutate(first_obs_by_gp_XX = row_number() == 1)
  df$first_obs_by_gp_XX <- as.numeric(df$first_obs_by_gp_XX)

  # Flagging first obs in cluster and checking if the cluster variable is weakly coarser than the group one.
  if (!is.null(cluster)) {
    df <- df %>% group_by(.data$cluster_XX) %>%
        mutate(first_obs_by_clust_XX = row_number() == 1)
    df$first_obs_by_clust_XX <- as.numeric(df$first_obs_by_clust_XX)

    df <- df %>% group_by(.data$group_XX) %>%
        mutate(cluster_var_g_XX = sd(.data$cluster_XX))
    if (max(df$cluster_var_g_XX) > 0) {
      stop("The group variable should be nested within the clustering variable.")
    }
  }

  df <- pdata.frame(df, index = c("group_XX", "time_XX")) 
  df$time_XX <- as.numeric(as.character(df$time_XX))
  df$group_XX <- as.numeric(as.character(df$group_XX))
  df$diff_y_XX <- diff(df$outcome_XX)
  df$diff_d_XX <- diff(df$treatment_XX)

  if (!is.null(controls)) {
    count_controls <- 0
    df$fd_X_all_non_missing_XX <- 1
    for (var in controls) {
      count_controls <- count_controls + 1
      # Computing the first difference of control variables #
      df[paste0("diff_X", count_controls, "_XX")] <- diff(df[[var]])
      df$fd_X_all_non_missing_XX <- ifelse(
         is.na(df[[paste0("diff_X", count_controls, "_XX")]]), 0, df$fd_X_all_non_missing_XX)
    }

    count_controls <- 0
    mycontrols_XX <- c()
    prod_controls_y <- ""
    for (var in controls) {
      count_controls <- count_controls + 1
      df <- df %>% select(-any_of(c("sum_weights_contol_XX","avg_diff_temp_XX", "diff_y_wXX")))

      # First step of the residualization : computing \Delta X_{.,t}: average of controls' first difference at time t, among groups whose treatment has not changed yet.

      df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
      df <- df %>% group_by(.data$joint_trends_XX) %>% 
          mutate(sum_weights_control_XX = sum(.data$N_gt_XX[
            .data$ever_change_d_XX == 0 & !is.na(.data$diff_y_XX) &
            .data$fd_X_all_non_missing_XX == 1
          ])) %>% ungroup()
      df$sum_weights_control_XX <- ifelse(df$ever_change_d_XX == 0 & 
      !is.na(df$diff_y_XX) & df$fd_X_all_non_missing_XX == 1, df$sum_weights_control_XX, NA)

      df$avg_diff_temp_XX <- df$N_gt_XX * df[[paste0("diff_X",count_controls,"_XX")]]
      df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
      df <- df %>% group_by(.data$joint_trends_XX) %>% 
          mutate(!!paste0("avg_diff_X",count_controls,"_XX") := sum(.data$avg_diff_temp_XX[
            .data$ever_change_d_XX == 0 & !is.na(.data$diff_y_XX) &
            .data$fd_X_all_non_missing_XX == 1
          ])) %>% ungroup()
      df[[paste0("avg_diff_X", count_controls,"_XX")]] <- ifelse(df$ever_change_d_XX == 0 & 
      !is.na(df$diff_y_XX) & df$fd_X_all_non_missing_XX == 1, df[[paste0("avg_diff_X", count_controls,"_XX")]], NA)
      df[[paste0("avg_diff_X", count_controls,"_XX")]] <-df[[paste0("avg_diff_X", count_controls,"_XX")]]  / df$sum_weights_control_XX

      # Computing \Delta\Dot{X}_{g,t}, the difference between the first differences of covariates and the average of their first-difference, which gives us the residuals of a regression of covariates on time fixed effects. 
      # Multiply by sqrt(N_gt_XX) to replicate weighted regression

      df[paste0("resid_X", count_controls, "_time_FE_XX")] <- sqrt(df$N_gt_XX) * 
          (df[[paste0("diff_X", count_controls,"_XX")]] - 
            df[[paste0("avg_diff_X",count_controls, "_XX")]])
      df[[paste0("resid_X", count_controls, "_time_FE_XX")]] <- ifelse(
        is.na(df[paste0("resid_X", count_controls, "_time_FE_XX")]), 0, 
        df[[paste0("resid_X", count_controls, "_time_FE_XX")]])

      # Storing the obtained residuals for the computation of theta_d
      mycontrols_XX <- c(mycontrols_XX, paste0("resid_X", count_controls, "_time_FE_XX"))

      # Generating the product between \Delta\Dot{X}_{g,t} and \Delta Y_{g,t}
      # Multiply by sqrt(N_gt_XX) to replicate weighted regression
      df$diff_y_wXX <- sqrt(df$N_gt_XX) * df$diff_y_XX
      df[[paste0("prod_X", count_controls, "_diff_y_temp_XX")]] <- ifelse(
        df$time_XX >= 2 & df$time_XX < df$F_g_XX, 
          df[[paste0("resid_X", count_controls,"_time_FE_XX")]] * df$diff_y_XX, NA)
      df[[paste0("prod_X", count_controls, "_diff_y_temp_XX")]] <- ifelse(
        is.na(df[[paste0("prod_X", count_controls, "_diff_y_temp_XX")]]), 0,
          df[[paste0("prod_X", count_controls, "_diff_y_temp_XX")]])

      # Computing the sum for each group to obtain the term \sum_{t=2}^{F_g-1}*N_{g,t}*\Delta \Dot{X}_{g,t}* \Delta Y_{g,t}
      df <- df %>% group_by(.data$group_XX) %>% 
          mutate(!!paste0("prod_X",count_controls,"_diff_y_XX") := 
            sum(.data[[paste0("prod_X", count_controls,"_diff_y_temp_XX")]]))
    }

    store_singular_XX <- ""
    store_noresidualization_XX <- c()
    levels_d_sq_XX <- levels(as.factor(df$d_sq_int_XX))

    for (l in levels_d_sq_XX) {
        data_XX <- df
        assign(paste0("store_singular_", l, "_XX"), 0)
        assign(paste0("useful_res_", l, "_XX"), 
            length(levels(as.factor(data_XX$F_g_XX[data_XX$d_sq_int_XX == l]))))
        if (get(paste0("useful_res_",l,"_XX")) > 1) {

          # Isolate the observations used for the computation of theta_d
          data_XX <- subset(data_XX, data_XX$ever_change_d_XX == 0 & !is.na(data_XX$diff_y_XX) & data_XX$fd_X_all_non_missing_XX == 1 & data_XX$d_sq_int_XX == l)

          # R version of the accum function
          #-- The final matrix should be order k + 1 with k n. of controls
          Y_vec <- as.matrix(data_XX$diff_y_wXX)
          X_vec <- as.matrix(data_XX[mycontrols_XX])
          W_vec <- as.matrix(data_XX$weight_XX)
          YX_vec <- cbind(Y_vec, X_vec, matrix(1, nrow = length(Y_vec), ncol = 1))
#          YX_vec <- cbind(Y_vec, X_vec, W_vec)
          overall_XX <- matrix(NA, nrow = count_controls + 2, ncol = count_controls + 2)
          for (i in 1:(2 + count_controls)) {
            for (j in 1:(2 + count_controls)) {
              overall_XX[i,j] <- t(YX_vec[,i]) %*% YX_vec[,j]
            }
          }
          # Check if the corresponding matrix exists
          e_vec <- matrix(1, nrow = 1, ncol = count_controls + 2)
          if (is.na(e_vec %*% overall_XX %*% t(e_vec))) {
              assign(paste0("store_singular_", l, "_XX"), 1)
              store_noresidualization_XX <- c(store_noresidualization_XX, l)
              assign(paste0("useful_res_", l, "_XX"), 1)
          }
          else {
            # Isolate the part of the matrix we need
            didmgt_XX <- overall_XX[2:(count_controls + 1), 2:(count_controls + 1)]
            didmgt_XY <- overall_XX[2:(count_controls + 1), 1]

            # Computing the vectors of coefficients \theta_d for each l
            assign(paste0("coefs_sq_", l, "_XX"), Ginv(as.matrix(didmgt_XX), tol = 10^(-16)) %*% didmgt_XY)

            # Computing the matrix Denom^{-1}
            # Check if the matrix is invertible
            det_XX <- det(as.matrix(didmgt_XX))
            if (abs(det_XX) <= 10^(-16)) {
              assign(paste0("store_singular_", l, "_XX"), 1)
            } 

            assign(paste0("inv_Denom_",l,"_XX"), Ginv(as.matrix(didmgt_XX), tol = 10^(-16)) * G_XX)
            }
            
          }
        assign(paste0("coefs_sq_", l, "_XX"), get(paste0("coefs_sq_", l, "_XX")), envir = globalenv())
        assign(paste0("useful_res_", l, "_XX"), get(paste0("useful_res_", l, "_XX")), envir = globalenv())
        }
    levels_d_sq_bis_XX <- levels(as.factor(df$d_sq_XX))
    index_sing_XX <- 0 
    for (l in levels_d_sq_bis_XX) {
      index_sing_XX <- index_sing_XX + 1
      if(get(paste0("store_singular_", index_sing_XX,"_XX")) == 1) {
        store_singular_XX <- paste(store_singular_XX, l)
      }
    }

    if (store_singular_XX != "") {
      cat(sprintf("Some control variables are not taken into account for groups with baseline treatment equal to: %s. This may occur in the following situations:", store_singular_XX))
      cat("\n")
      cat("1. For groups with those values of the baseline treatment, the regression of the outcome first difference on the controls' first differences and time fixed effects has fewer observations than variables. Note that for each value of the baseline treatment, those regressions are estimated among (g,t)s such that g has not changed treatment yet at t.")
      cat("\n")
      cat("2. For groups with those values of the baseline treatment, two or more of your control variables are perfectly collinear in the sample where the regression is run, for instance because those control variables do not vary over time.")
      cat("\n")
    }

    for (l in store_noresidualization_XX) {
      df <- subset(df, df$d_sq_int_XX != l)
    }

  }

  # Initialize L_u_XX/L_a_XX
  L_u_XX <- NA
  L_a_XX <- NA

  if (switchers == "" | switchers == "in") {
    L_u_XX <- max(df$L_g_XX[df$S_g_XX == 1], na.rm = TRUE)  
    if (placebo != 0) {
      L_placebo_u_XX <- max(df$L_g_placebo_XX[df$S_g_XX == 1], na.rm = TRUE)  
    }
  }

  if (switchers == "" | switchers == "out") {
    L_a_XX <- max(df$L_g_XX[df$S_g_XX == 0], na.rm = TRUE)  
    if (placebo != 0) {
      L_placebo_a_XX <- max(df$L_g_placebo_XX[df$S_g_XX == 0], na.rm = TRUE)  
    }
  }

  if (
    (switchers == "in" & (is.na(L_u_XX) | L_u_XX == 0)) | 
    (switchers == "out" & (is.na(L_a_XX) | L_a_XX == 0)) | 
    (switchers == "" &  ((is.na(L_u_XX) | L_u_XX == 0) & (is.na(L_a_XX) | L_a_XX == 0)))
    ) {
    stop("No treatment effect can be estimated. This is because Assumption 1 in de Chaisemartin & D'Haultfoeuille (2023) is not satisfied in the data used forestimation, given the options requested. If this is caused by your baseline treatement being continuous you can try using the option continuous() which allows for a continous period-one treatement.")
  }

  if (switchers == "" ) {
    l_XX <- max(L_a_XX, L_u_XX, na.rm = TRUE)
    l_XX <- min(l_XX, effects)
    if (placebo != 0) {
      l_placebo_XX <- max(L_placebo_a_XX, L_placebo_u_XX)
      l_placebo_XX <- min(l_placebo_XX, placebo)
      # The number of placebos cannot be greater than the number of effects computed:
      l_placebo_XX <- min(l_placebo_XX, effects)
    } else {
      l_placebo_XX <- 0
    }
  }

  if (switchers == "in") {
    l_XX <- min(effects, L_u_XX)
    if (placebo != 0) {
      l_placebo_XX <- min(placebo, L_placebo_u_XX)
      # The number of placebos cannot be greater than the number of effects computed:
      l_placebo_XX <- min(l_placebo_XX, effects)
    }
    else {
      l_placebo_XX <- 0
    }
  }

  if (switchers == "out") {
    l_XX <- min(effects, L_a_XX)
    if (placebo != 0) {
      l_placebo_XX <- min(placebo, L_placebo_a_XX)
      # The number of placebos cannot be greater than the number of effects computed:
      l_placebo_XX <- min(l_placebo_XX, effects)
    }
    else {
      l_placebo_XX <- 0
    }
  }

  if (l_XX < effects) {
    print(warning(sprintf("The number of effects requested is too large. The number of effects which can be estimated is at most %.0f. The command will therefore try to estimante %.0f effect(s)", l_XX, l_XX)))
  }

  if (placebo != 0) {
    if (l_placebo_XX < placebo & effects >= placebo) {
      print(warning(sprintf("The number of placebos which can be estimated is at most %.0f.The command will therefore try to estimate %.0f placebo(s).", l_placebo_XX, l_placebo_XX)))
    }
    if (effects < placebo) {
      print(warning(sprintf("The number of placebo requested cannot be larger than the number of effects requested. The command cannot compute more than %.0f placebo(s).", l_placebo_XX)))
    }
  }

  # Generating default values for the variables which will be aggregated.
  df[paste0("U_Gg", 1:l_XX, "_plus_XX")] <- lapply(1:l_XX, function(i) 0)
  df[paste0("U_Gg", 1:l_XX, "_minus_XX")] <- lapply(1:l_XX, function(i) 0)
  df[paste0("count", 1:l_XX, "_plus_XX")] <- lapply(1:l_XX, function(i) 0)
  df[paste0("count", 1:l_XX, "_minus_XX")] <- lapply(1:l_XX, function(i) 0)
  df[paste0("U_Gg_var_", 1:l_XX, "_in_XX")] <- lapply(1:l_XX, function(i) 0)
  df[paste0("U_Gg_var_", 1:l_XX, "_out_XX")] <- lapply(1:l_XX, function(i) 0)
  assign("sum_for_var_in_XX", 0, envir = globalenv())
  assign("sum_for_var_out_XX", 0, envir = globalenv())
  if (placebo != 0) {
    df[paste0("U_Gg_pl_", 1:l_XX, "_plus_XX")] <- lapply(1:l_XX, function(i) 0)
    df[paste0("U_Gg_pl_", 1:l_XX, "_minus_XX")] <- lapply(1:l_XX, function(i) 0)
    df[paste0("count", 1:l_XX, "_pl_plus_XX")] <- lapply(1:l_XX, function(i) 0)
    df[paste0("count", 1:l_XX, "_pl_minus_XX")] <- lapply(1:l_XX, function(i) 0)
    df[paste0("U_Gg_var_pl_", 1:l_XX, "_in_XX")] <- lapply(1:l_XX, function(i) 0)
    df[paste0("U_Gg_var_pl_", 1:l_XX, "_out_XX")] <- lapply(1:l_XX, function(i) 0)
    assign("sum_for_var_placebo_in_XX", 0, envir = globalenv())
    assign("sum_for_var_placebo_out_XX", 0, envir = globalenv())
  }

  for (i in 1:l_XX) {
    assign(paste0("N1_", i, "_XX"), 0, envir = globalenv()) 
    assign(paste0("N1_", i, "_XX_new"), 0, envir = globalenv()) 
    assign(paste0("N0_", i, "_XX"), 0, envir = globalenv()) 
    assign(paste0("N0_", i, "_XX_new"), 0, envir = globalenv()) 
    if (normalized == TRUE) {
      assign(paste0("delta_D_",i,"_in_XX"), 0, envir = globalenv())
      assign(paste0("delta_D_",i,"_out_XX"), 0, envir = globalenv())      
    }
    if (placebo != 0) {
      assign(paste0("N1_placebo_", i, "_XX"), 0, envir = globalenv()) 
      assign(paste0("N1_placebo_", i, "_XX_new"), 0, envir = globalenv()) 
      assign(paste0("N0_placebo_", i, "_XX"), 0, envir = globalenv()) 
      assign(paste0("N0_placebo_", i, "_XX_new"), 0, envir = globalenv()) 
      if (normalized == TRUE) {
        assign(paste0("delta_D_pl_",i,"_in_XX"), 0, envir = globalenv())
        assign(paste0("delta_D_pl_",i,"_out_XX"), 0, envir = globalenv())      
      }
    }
  }

  df$U_Gg_plus_XX <- 0
  df$U_Gg_minus_XX <- 0
  assign("U_Gg_den_plus_XX", 0, envir = globalenv())
  assign("U_Gg_den_minus_XX", 0, envir = globalenv())
  assign("sum_N1_l_XX", 0, envir = globalenv())
  assign("sum_N0_l_XX", 0, envir = globalenv())
  df$U_Gg_var_plus_XX <- 0
  df$U_Gg_var_minus_XX <- 0  

  # Saving useful scalars to the Global Environment
  assign("L_u_XX",L_u_XX, envir = globalenv())
  assign("L_a_XX",L_a_XX, envir = globalenv())
  assign("L_placebo_u_XX",L_u_XX, envir = globalenv())
  assign("L_placebo_a_XX",L_a_XX, envir = globalenv())
  assign("l_XX",l_XX, envir = globalenv())
  assign("t_min_XX",t_min_XX, envir = globalenv())
  assign("T_max_XX",T_max_XX, envir = globalenv())
  assign("G_XX",G_XX, envir = globalenv())

  if (switchers == "" | switchers == "in") {
    if (!is.na(L_u_XX) & L_u_XX != 0) {
      df <- did_multiplegt_dyn_core(df, outcome_XX, group_XX, time_XX, treatment_XX, effects = l_XX, placebo = l_placebo_XX, switchers_core = "in", trends_nonparam = trends_nonparam, controls = controls, same_switchers, same_switchers_pl, normalized)
      for (i in 1:l_XX) {
        if (get(paste0("N1_",i,"_XX")) != 0) {
          df[[paste0("U_Gg",i,"_plus_XX")]] <- df[[paste0("U_Gg",i,"_XX")]]
          df[[paste0("count",i,"_plus_XX")]] <- df[[paste0("count",i,"_core_XX")]]
          df[[paste0("U_Gg_var_",i,"_in_XX")]] <- df[[paste0("U_Gg",i,"_var_XX")]]
          assign(paste0("N1_",i,"_XX_new"), get(paste0("N1_",i,"_XX")))

          if (normalized == TRUE) {
            assign(paste0("delta_D_",i,"_in_XX"), get(paste0("delta_norm_",i,"_XX")))
          }
        }
      }

      if (l_placebo_XX != 0) {
        for (i in 1:l_placebo_XX) {
          if (get(paste0("N1_placebo_",i,"_XX")) != 0) {
            df[[paste0("U_Gg_pl_",i,"_plus_XX")]] <- df[[paste0("U_Gg_placebo_",i,"_XX")]]
            df[[paste0("count",i,"_pl_plus_XX")]] <- df[[paste0("count",i,"_pl_core_XX")]]
            df[[paste0("U_Gg_var_pl_",i,"_in_XX")]] <- df[[paste0("U_Gg_pl_",i,"_var_XX")]]
            assign(paste0("N1_placebo_",i,"_XX_new"), get(paste0("N1_placebo_",i,"_XX")))

            if (normalized == TRUE) {
              assign(paste0("delta_D_pl_",i,"_in_XX"), 
                  get(paste0("delta_norm_pl_",i,"_XX")))
            }
          }

        }
      }

      if (sum_N1_l_XX != 0) {
        df$U_Gg_plus_XX <- df$U_Gg_XX
        df$U_Gg_den_plus_XX <- df$U_Gg_den_XX
        df$U_Gg_var_plus_XX <- df$U_Gg_var_XX
      }
    }
  }

  if (switchers == "" | switchers == "out") {
    if (!is.na(L_a_XX) & L_a_XX != 0) {
      df <- did_multiplegt_dyn_core(df, outcome_XX, group_XX, time_XX, treatment_XX, effects = l_XX, placebo = l_placebo_XX, switchers_core = "out", trends_nonparam = trends_nonparam, controls = controls, same_switchers, same_switchers_pl, normalized)
      for (i in 1:l_XX) {
        if (get(paste0("N0_",i,"_XX")) != 0) {
          df[[paste0("U_Gg",i,"_minus_XX")]] <- - df[[paste0("U_Gg",i,"_XX")]]
          df[[paste0("count",i,"_minus_XX")]] <- df[[paste0("count",i,"_core_XX")]]
          df[[paste0("U_Gg_var_",i,"_out_XX")]] <- df[[paste0("U_Gg",i,"_var_XX")]]
          assign(paste0("N0_",i,"_XX_new"), get(paste0("N0_",i,"_XX")))

          if (normalized == TRUE) {
            assign(paste0("delta_D_",i,"_out_XX"), get(paste0("delta_norm_",i,"_XX")))
          }
        }
      }

      if (l_placebo_XX != 0) {
        for (i in 1:l_placebo_XX) {
          if (get(paste0("N0_placebo_",i,"_XX")) != 0) {
            df[[paste0("U_Gg_pl_",i,"_minus_XX")]] <- - df[[paste0("U_Gg_placebo_",i,"_XX")]]
            df[[paste0("count",i,"_pl_minus_XX")]] <- df[[paste0("count",i,"_pl_core_XX")]]
            df[[paste0("U_Gg_var_pl_",i,"_out_XX")]] <- df[[paste0("U_Gg_pl_",i,"_var_XX")]]
            assign(paste0("N0_placebo_",i,"_XX_new"), get(paste0("N0_placebo_",i,"_XX")))

            if (normalized == TRUE) {
              assign(paste0("delta_D_pl_",i,"_out_XX"), get(paste0("delta_norm_pl_",i,"_XX")))
            }
        }
        }
      }

      if (sum_N0_l_XX != 0) {
        df$U_Gg_minus_XX <- - df$U_Gg_XX
        df$U_Gg_den_minus_XX <- df$U_Gg_den_XX
        df$U_Gg_var_minus_XX <- df$U_Gg_var_XX
      }
    }
  }
  rownames <- c()
  # Aggregating the results for switchers in and out. #####

  mat_res_XX <- matrix(NA, nrow = l_XX + l_placebo_XX + 1, ncol = 7)

  #------ DID_l ------------------------------------------#

  for (i in 1:l_XX) {
    df[paste0("U_Gg",i,"_global_XX")] <- get(paste0("N1_",i,"_XX_new")) / (get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new"))) * df[[paste0("U_Gg",i,"_plus_XX")]] + get(paste0("N0_",i,"_XX_new")) / (get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new"))) * df[[paste0("U_Gg",i,"_minus_XX")]]
    df[[paste0("U_Gg",i,"_global_XX")]][df$first_obs_by_gp_XX == 0] <- NA

    df <- df %>% rowwise() %>% 
        mutate(!!paste0("count",i,"_global_XX") :=  max(.data[[paste0("count",i,"_plus_XX")]], .data[[paste0("count",i,"_minus_XX")]], na.rm = TRUE)) 
    df[[paste0("count",i,"_global_XX")]][df[[paste0("count",i,"_global_XX")]] == -Inf] <- NA

    if (normalized == TRUE) {
      assign(paste0("delta_D_",i,"_global_XX"), 
      ( get(paste0("N1_",i,"_XX_new")) /( get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new")) ) ) * get(paste0("delta_D_",i,"_in_XX")) +  (get(paste0("N0_",i,"_XX_new"))/(get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new")))) * get(paste0("delta_D_",i,"_out_XX")), envir = globalenv())
    }

    # Number of switchers
    assign(paste0("N_switchers_effect_",i,"_XX"), get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new")), envir = globalenv())
    mat_res_XX[i,6] <- get(paste0("N_switchers_effect_",i,"_XX"))
    mat_res_XX[i,7] <- i
    assign(paste0("N_switchers_effect_",i), get(paste0("N_switchers_effect_",i,"_XX")), envir = globalenv())
    # Number of observations used in the estimation 
    df[paste0("N_effect_",i,"_XX")] <- sum(df[[paste0("count",i,"_global_XX")]], na.rm = TRUE)
    assign(paste0("N_effect_",i,"_XX"), mean(df[[paste0("N_effect_",i,"_XX")]]))
    assign(paste0("N_effect_",i), get(paste0("N_effect_",i,"_XX")), envir = globalenv())
    mat_res_XX[i,5] <- get(paste0("N_effect_",i,"_XX"))

    if (get(paste0("N_switchers_effect_",i,"_XX")) == 0 | get(paste0("N_effect_",i,"_XX")) == 0) {
      warning(paste0("Effect_",i,"cannot be estimated. There is no switcher or no control for this effect."))
    }

    # DID_l
    df[paste0("DID_",i,"_XX")] <- sum(df[[paste0("U_Gg",i,"_global_XX")]], na.rm = TRUE)
    df[[paste0("DID_",i,"_XX")]] <- df[[paste0("DID_",i,"_XX")]] / G_XX

    if (normalized == TRUE) {
        df[[paste0("DID_",i,"_XX")]] <- df[[paste0("DID_",i,"_XX")]] / get(paste0("delta_D_",i,"_global_XX"))
    }
    assign(paste0("DID_",i,"_XX"), mean(df[[paste0("DID_",i,"_XX")]]))

    if ((switchers == "" & get(paste0("N1_",i,"_XX_new")) == 0 & get(paste0("N0_",i,"_XX_new")) == 0) | (switchers == "out" & get(paste0("N0_",i,"_XX_new")) == 0 ) |
    (switchers == "in" & get(paste0("N1_",i,"_XX_new")) == 0 )) {
        assign(paste0("DID_",i,"_XX"), NA)
    }

    assign(paste0("Effect_",i), get(paste0("DID_",i,"_XX")), envir = globalenv())
    mat_res_XX[i,1] <- get(paste0("DID_",i,"_XX")) 
    rownames <- append(rownames, paste0("Effect_",i, strrep(" ",(12 - nchar(paste0("Effect_",i))))))
  }

  #------ ATE --------------------------------------------#
  U_Gg_den_plus_XX <- mean(df$U_Gg_den_plus_XX, na.rm = TRUE)
  U_Gg_den_minus_XX <- mean(df$U_Gg_den_minus_XX, na.rm = TRUE)

  if (switchers == "") {
    w_plus_XX <- U_Gg_den_plus_XX * sum_N1_l_XX / (U_Gg_den_plus_XX * sum_N1_l_XX + U_Gg_den_minus_XX * sum_N0_l_XX)
  }
  if (switchers == "out") {
    w_plus_XX <- 0
  }
  if (switchers == "in") {
    w_plus_XX <- 1
  }

  df$U_Gg_global_XX <- w_plus_XX * df$U_Gg_plus_XX + (1 - w_plus_XX) * df$U_Gg_minus_XX
  df$U_Gg_global_XX[df$first_obs_by_gp == 0] <- NA

  df$delta_XX <- sum(df$U_Gg_global_XX, na.rm = TRUE) / G_XX
  delta_XX <- mean(df$delta_XX, na.rm = TRUE)
  assign("Av_tot_effect", delta_XX, envir = globalenv())
  mat_res_XX[l_XX+1,1] <- delta_XX
  N_switchers_effect_XX <- 0
  for (i in 1:l_XX) {
    N_switchers_effect_XX <- N_switchers_effect_XX + get(paste0("N_switchers_effect_",i,"_XX"))
  }
  mat_res_XX[l_XX+1,6] <- N_switchers_effect_XX
  mat_res_XX[l_XX+1,7] <- 0
  assign("N_switchers_effect_average", N_switchers_effect_XX)
  df$count_global_XX <- 0
  for (i in 1:l_XX) {
    df <- df %>% rowwise() %>% 
    mutate(count_global_XX = max(.data$count_global_XX, .data[[paste0("count",i,"_global_XX")]], na.rm = TRUE))
  }
  N_effect_XX <- sum(df$count_global_XX, na.rm = TRUE)
  mat_res_XX[l_XX+1,5] <- N_effect_XX
  assign("N_avg_total_effect", N_effect_XX, envir = globalenv())
  rownames <- append(rownames, paste0("Av_tot_eff", strrep(" ",(12 - nchar("Av_tot_eff")))))

  #-- Placebos ---------------------------------------------#

  if (l_placebo_XX != 0) {
    for (i in 1:l_placebo_XX) {
      df[paste0("U_Gg_pl_",i,"_global_XX")] <- get(paste0("N1_placebo_",i,"_XX_new")) / (get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new"))) * df[[paste0("U_Gg_pl_",i,"_plus_XX")]] + get(paste0("N0_placebo_",i,"_XX_new")) / (get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new"))) * df[[paste0("U_Gg_pl_",i,"_minus_XX")]]
      df[[paste0("U_Gg_pl_",i,"_global_XX")]][df$first_obs_by_gp_XX == 0] <- NA

      df <- df %>% rowwise() %>% 
          mutate(!!paste0("count",i,"_pl_global_XX") :=  max(.data[[paste0("count",i,"_pl_plus_XX")]], .data[[paste0("count",i,"_pl_minus_XX")]], na.rm = TRUE)) 
      df[[paste0("count",i,"_pl_global_XX")]][df[[paste0("count",i,"_pl_global_XX")]] == -Inf] <- NA

      if (normalized == TRUE) {
        assign(paste0("delta_D_pl_",i,"_global_XX"), 
        ( get(paste0("N1_placebo_",i,"_XX_new")) /( get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new")) ) ) * get(paste0("delta_D_pl_",i,"_in_XX")) +  (get(paste0("N0_placebo_",i,"_XX_new"))/(get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new")))) * get(paste0("delta_D_pl_",i,"_out_XX")))
      }

      df[paste0("DID_placebo_",i,"_XX")] <- sum(df[[paste0("U_Gg_pl_",i,"_global_XX")]], na.rm = TRUE)
      df[[paste0("DID_placebo_",i,"_XX")]] <- df[[paste0("DID_placebo_",i,"_XX")]] / G_XX

      if (normalized == TRUE) {
          df[[paste0("DID_placebo_",i,"_XX")]] <- df[[paste0("DID_placebo_",i,"_XX")]] / get(paste0("delta_D_pl_",i,"_global_XX"))
      }
      assign(paste0("DID_placebo_",i,"_XX"), mean(df[[paste0("DID_placebo_",i,"_XX")]]))

      if ((switchers == "" & get(paste0("N1_placebo_",i,"_XX_new")) == 0 & get(paste0("N0_placebo_",i,"_XX_new")) == 0) | (switchers == "out" & get(paste0("N0_placebo_",i,"_XX_new")) == 0 ) | (switchers == "in" & get(paste0("N1_placebo_",i,"_XX_new")) == 0 )) {
          assign(paste0("DID_placebo_",i,"_XX"), NA)
      }

      assign(paste0("Placebo_",i), get(paste0("DID_placebo_",i,"_XX")), envir = globalenv())
      mat_res_XX[l_XX+1+i,1] <- get(paste0("DID_placebo_",i,"_XX")) 
      rownames <- append(rownames, paste0("Placebo_",i, strrep(" ",(12 - nchar(paste0("Placebo_",i))))))

      # Number of switchers
      assign(paste0("N_switchers_placebo_",i,"_XX"), get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new")))
      mat_res_XX[l_XX+1+i,6] <- get(paste0("N_switchers_placebo_",i,"_XX"))
      mat_res_XX[l_XX+1+i,7] <- -i
      assign(paste0("N_switchers_placebo_",i), get(paste0("N_switchers_placebo_",i,"_XX")), envir = globalenv())
      # Number of observations used in the estimation 
      df[paste0("N_placebo_",i,"_XX")] <- sum(df[[paste0("count",i,"_pl_global_XX")]], na.rm = TRUE)
      assign(paste0("N_placebo_",i,"_XX"), mean(df[[paste0("N_placebo_",i,"_XX")]]))
      assign(paste0("N_placebo_",i), get(paste0("N_placebo_",i,"_XX")), envir = globalenv())
      mat_res_XX[l_XX + 1 + i,5] <- get(paste0("N_placebo_",i,"_XX"))

      if (get(paste0("N_switchers_placebo_",i,"_XX")) == 0 | get(paste0("N_placebo_",i,"_XX")) == 0) {
        warning(paste0("Placebo_",i,"cannot be estimated. There is no switcher or no control for this placebo."))
      }

    }
  }

  # Estimating the asymptotic variances #####################################

  ci_level <- ci_level / 100
  z_level <- qnorm(ci_level + (1 - ci_level)/2)

  #-- Estimating \hat{\sigma}^2_l ------------------------------------------#
  for (i in 1:l_XX) {
    if ((switchers == "" & get(paste0("N1_",i,"_XX_new")) != 0 & get(paste0("N0_",i,"_XX_new")) != 0) | (switchers == "out" & get(paste0("N0_",i,"_XX_new")) != 0 ) | (switchers == "in" & get(paste0("N1_",i,"_XX_new")) != 0 )) {
        
        df[paste0("U_Gg_var_glob_",i,"_XX")] <- df[[paste0("U_Gg_var_",i,"_in_XX")]] * (get(paste0("N1_",i,"_XX_new")) / (get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new")))) + df[[paste0("U_Gg_var_",i,"_out_XX")]] * (get(paste0("N0_",i,"_XX_new")) / (get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new"))))

        if (is.null(cluster)) {
        df[paste0("U_Gg_var_glob_eff",i,"_sqrd_XX")] <- df[[paste0("U_Gg_var_glob_",i,"_XX")]]^2 * df$first_obs_by_gp_XX
        assign(paste0("sum_for_var_",i,"_XX"), sum(df[[paste0("U_Gg_var_glob_eff",i,"_sqrd_XX")]], na.rm = TRUE) / G_XX^2) 
        } else {
          df[[paste0("U_Gg_var_glob_",i,"_XX")]] <- df[[paste0("U_Gg_var_glob_",i,"_XX")]] * df$first_obs_by_gp_XX
          df <- df %>% group_by(.data$cluster_XX) %>%
              mutate(!!paste0("clust_U_Gg_var_glob_",i,"_XX") 
                  := sum(.data[[paste0("U_Gg_var_glob_",i,"_XX")]], na.rm = TRUE))
          df[paste0("clust_U_Gg_var_glob_",i,"_2_XX")] <-
              df[[paste0("clust_U_Gg_var_glob_", i, "_XX")]]^2 * df$first_obs_by_clust_XX
          assign(paste0("sum_for_var_",i,"_XX"), 
              sum(df[[paste0("clust_U_Gg_var_glob_", i,"_2_XX")]], na.rm = TRUE) / G_XX^2)
          df[[paste0("U_Gg_var_glob_",i,"_XX")]] <- 
              df[[paste0("clust_U_Gg_var_glob_",i,"_XX")]]
        }

        assign(paste0("se_",i,"_XX"), sqrt(get(paste0("sum_for_var_",i,"_XX"))))
        if (normalized == TRUE) {
          assign(paste0("se_",i,"_XX"), get(paste0("se_",i,"_XX")) / get(paste0("delta_D_",i,"_global_XX")))
        }
        mat_res_XX[i,2] <- get(paste0("se_",i,"_XX"))
        assign(paste0("se_effect_",i),get(paste0("se_",i,"_XX")), envir = globalenv())

        # CI level
        assign(paste0("LB_CI_",i,"_XX"), get(paste0("DID_",i,"_XX")) - z_level * get(paste0("se_",i,"_XX")))
        mat_res_XX[i,3] <- get(paste0("LB_CI_",i,"_XX"))
        assign(paste0("UB_CI_",i,"_XX"), get(paste0("DID_",i,"_XX")) + z_level * get(paste0("se_",i,"_XX")))
        mat_res_XX[i,4] <- get(paste0("UB_CI_",i,"_XX"))
    }
  }

  #-- Estimating \hat{\sigma}^2_pl ------------------------------------------#
  if (l_placebo_XX != 0) {
    for (i in 1:l_placebo_XX) {
      if ((switchers == "" & get(paste0("N1_placebo_",i,"_XX_new")) != 0 & get(paste0("N0_placebo_",i,"_XX_new")) != 0) | (switchers == "out" & get(paste0("N0_placebo_",i,"_XX_new")) != 0 ) | (switchers == "in" & get(paste0("N1_placebo_",i,"_XX_new")) != 0 )) {
          
          df[paste0("U_Gg_var_glob_pl_",i,"_XX")] <- df[[paste0("U_Gg_var_pl_",i,"_in_XX")]] * (get(paste0("N1_placebo_",i,"_XX_new")) / (get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new")))) + df[[paste0("U_Gg_var_pl_",i,"_out_XX")]] * (get(paste0("N0_placebo_",i,"_XX_new")) / (get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new"))))

          if (is.null(cluster)) {
          df[paste0("U_Gg_var_glob_pl_",i,"_2_XX")] <- df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]]^2 * df$first_obs_by_gp_XX
          assign(paste0("sum_for_var_placebo_",i,"_XX"), sum(df[[paste0("U_Gg_var_glob_pl_",i,"_2_XX")]], na.rm = TRUE) / G_XX^2) 
          } else {
            df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]] <- df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]] * df$first_obs_by_gp_XX
            df <- df %>% group_by(.data$cluster_XX) %>%
                mutate(!!paste0("clust_U_Gg_var_glob_pl_",i,"_XX") 
                    := sum(.data[[paste0("U_Gg_var_glob_pl_",i,"_XX")]], na.rm = TRUE))
            df[paste0("clust_U_Gg_var_glob_pl_",i,"_2_XX")] <-
                df[[paste0("clust_U_Gg_var_glob_pl_", i, "_XX")]]^2 * df$first_obs_by_clust_XX
            assign(paste0("sum_for_var_placebo_",i,"_XX"), 
                sum(df[[paste0("clust_U_Gg_var_glob_pl_", i,"_2_XX")]], na.rm = TRUE) / G_XX^2)
            df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]] <- 
                df[[paste0("clust_U_Gg_var_glob_pl_",i,"_XX")]]
          }

          assign(paste0("se_placebo_",i,"_XX"), sqrt(get(paste0("sum_for_var_placebo_",i,"_XX"))))
          if (normalized == TRUE) {
            assign(paste0("se_placebo_",i,"_XX"), get(paste0("se_placebo_",i,"_XX")) / get(paste0("delta_D_pl_",i,"_global_XX")))
          }
          mat_res_XX[l_XX + 1 + i,2] <- get(paste0("se_placebo_",i,"_XX"))
          assign(paste0("se_placebo_",i),get(paste0("se_placebo_",i,"_XX")), envir = globalenv())

          # CI level
          assign(paste0("LB_CI_placebo_",i,"_XX"), get(paste0("DID_placebo_",i,"_XX")) - z_level * get(paste0("se_placebo_",i,"_XX")))
          mat_res_XX[l_XX + 1 + i,3] <- get(paste0("LB_CI_placebo_",i,"_XX"))
          assign(paste0("UB_CI_placebo_",i,"_XX"), get(paste0("DID_placebo_",i,"_XX")) + z_level * get(paste0("se_placebo_",i,"_XX")))
          mat_res_XX[l_XX + 1 + i,4] <- get(paste0("UB_CI_placebo_",i,"_XX"))
      }
    }
  }

  #-- Estimating \hat{\sigma}^2 __------------------------------------------#
if ((switchers=="" & (sum_N1_l_XX!=0|sum_N0_l_XX!=0))|(switchers =="out" & sum_N0_l_XX!=0)|(switchers=="in" & sum_N1_l_XX !=0)) {

  df$U_Gg_var_global_XX <- w_plus_XX * df$U_Gg_var_plus_XX + (1 - w_plus_XX) * df$U_Gg_var_minus_XX

  if (is.null(cluster)) {
  df$U_Gg_var_global_2_XX <- df$U_Gg_var_global_XX^2 * df$first_obs_by_gp_XX
  assign("sum_for_var_XX", sum(df$U_Gg_var_global_2_XX, na.rm = TRUE) / G_XX^2) 
  } else {
    df$U_Gg_var_global_XX <- df$U_Gg_var_global_XX * df$first_obs_by_gp_XX
    df <- df %>% group_by(.data$cluster_XX) %>%
        mutate(clust_U_Gg_var_global_XX = sum(.data$U_Gg_var_global_XX, na.rm = TRUE))
    df$clust_U_Gg_var_global_XX <- df$clust_U_Gg_var_global_XX^2 * df$first_obs_by_clust_XX
    assign("sum_for_var_XX", sum(df$clust_U_Gg_var_global_XX)/G_XX^2)
  }

  assign("se_XX", sqrt(sum_for_var_XX))
  mat_res_XX[l_XX + 1,2] <- se_XX
  assign("se_avg_total_effect",se_XX, envir = globalenv())

  # CI level
  LB_CI_XX <- delta_XX - z_level * se_XX
  mat_res_XX[l_XX + 1,3] <- LB_CI_XX
  UB_CI_XX <- delta_XX + z_level * se_XX
  mat_res_XX[l_XX + 1,4] <- UB_CI_XX
}

#-- F tests ---------------------------------------------------#
# If the option cluster is specified, we have previously replaced U_Gg_var_glob_pl_`i'_XX by clust_U_Gg_var_glob_pl_`i'_XX, and U_Gg_var_glob_`i'_XX by clust_U_Gg_var_glob_`i'_XX. 
# Now, we must also replace first_obs_by_gp_XX by first_obs_by_clust_XX
if (!is.null(cluster)) {
  df$first_obs_by_gp_XX <- df$first_obs_by_clust_XX
}

# Performing a test to see whether all placebos are jointly equal to 0
all_Ns_pl_not_zero <- NA
if (l_placebo_XX != 0 & l_placebo_XX > 1) {
  all_Ns_pl_not_zero <- 0 

  for (i in 1:l_placebo_XX) {
    if ( (switchers == "" & (get(paste0("N1_placebo_",i,"_XX_new"))!= 0 | get(paste0("N0_placebo_",i,"_XX_new"))!= 0 )) | (switchers == "out" & get(paste0("N0_placebo_",i,"_XX_new")) != 0) | (switchers == "in" & get(paste0("N1_placebo_",i,"_XX_new")) != 0) ) {
      all_Ns_pl_not_zero = all_Ns_pl_not_zero + 1
    }
  }

  if (all_Ns_pl_not_zero == l_placebo_XX) {
    didmgt_Placebo <- matrix(0, nrow = l_placebo_XX, ncol = 1)
    didmgt_Var_Placebo <- matrix(0, nrow = l_placebo_XX, ncol = l_placebo_XX)

    for (i in 1:l_placebo_XX) {
      didmgt_Placebo[i,1] <- get(paste0("DID_placebo_",i,"_XX"))
      didmgt_Var_Placebo[i,i] <- get(paste0("se_placebo_",i,"_XX"))^2

      if (i < l_placebo_XX) {
        for (j in (i+1):l_placebo_XX) {
          if (normalized == FALSE) {
            df[[paste0("U_Gg_var_pl_",i,"_",j,"_XX")]] <- df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]] +  df[[paste0("U_Gg_var_glob_pl_",j,"_XX")]] 
          } else {
            df[[paste0("U_Gg_var_pl_",i,"_",j,"_XX")]] <- df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]] / get(paste0("delta_D_pl_",i,"_global_XX")) +  df[[paste0("U_Gg_var_glob_pl_",j,"_XX")]] / get(paste0("delta_D_pl_",j,"_global_XX"))
          }

          df[[paste0("U_Gg_var_pl_",i,"_",j,"_2_XX")]] <- df[[paste0("U_Gg_var_pl_",i,"_",j,"_XX")]]^2 * df$first_obs_by_gp_XX
          assign(paste0("var_sum_pl_",i,"_",j,"_XX"), sum( df[[paste0("U_Gg_var_pl_",i,"_",j,"_2_XX")]], na.rm = TRUE) / G_XX^2)
          assign(paste0("cov_pl_",i,"_",j,"_XX"), (get(paste0("var_sum_pl_",i,"_",j,"_XX")) - get(paste0("se_placebo_",i,"_XX"))^2 - get(paste0("se_placebo_",j,"_XX"))^2)/2) 

          didmgt_Var_Placebo[i,j] <- get(paste0("cov_pl_",i,"_",j,"_XX"))
          didmgt_Var_Placebo[j,i] <- get(paste0("cov_pl_",i,"_",j,"_XX"))
        }
      }
    }

    didmgt_Var_Placebo_inv <- Ginv(didmgt_Var_Placebo)
    didmgt_chi2placebo <- t(didmgt_Placebo) %*% didmgt_Var_Placebo_inv  %*% didmgt_Placebo
    p_jointplacebo <- 1 - pchisq(didmgt_chi2placebo[1,1], df = l_placebo_XX)
    assign("p_jointplacebo", p_jointplacebo, envir = globalenv())
  } else {
    warning("Some placebos could not be estimated. Therefore, the test of joint nullity of the placebos could not be computed.")
  }
}

# Performing a test that all the DID_l are equal

if (effects_equal == TRUE & l_XX > 1) {
  all_Ns_not_zero <- 0
  for (i in 1:l_XX) {
    if (((switchers == "" & (get(paste0("N1_",i,"_XX_new")) != 0 | get(paste0("N0_",i,"_XX_new")) != 0)) | 
    (switchers == "out" & get(paste0("N0_",i,"_XX_new")) != 0) | 
    (switchers == "in" & get(paste0("N1_",i,"_XX_new")) != 0))) {
      all_Ns_not_zero <- all_Ns_not_zero + 1
    }
  }
  if (all_Ns_not_zero == l_XX) {
    didmgt_Effects <- mat_res_XX[1:l_XX, 1]
    didmgt_Var_Effects <- matrix(0, nrow = l_XX, ncol = l_XX)
    didmgt_identity <- matrix(0, nrow = l_XX - 1, ncol = l_XX)

  # Computing and storing the covariances.
  for (i in 1:l_XX) {
    if (((switchers == "" & (get(paste0("N1_",i,"_XX_new")) != 0 | get(paste0("N0_",i,"_XX_new")) != 0)) | 
    (switchers == "out" & get(paste0("N0_",i,"_XX_new")) != 0) | 
    (switchers == "in" & get(paste0("N1_",i,"_XX_new")) != 0))) {

      didmgt_Var_Effects[i,i] <- get(paste0("se_", i, "_XX")) ^ 2
      if (i < l_XX) {
        didmgt_identity[i,i] <- 1
      }

      if (i < l_XX) {
        for (j in (i + 1):l_XX) {
          if (normalized == FALSE) {
            df[[paste0("U_Gg_var_", i, "_", j,"_XX")]] <- df[[paste0("U_Gg_var_glob_", i, "_XX")]] +  df[[paste0("U_Gg_var_glob_", j, "_XX")]] 
          } else {
            df[[paste0("U_Gg_var_", i, "_", j,"_XX")]] <- 
            (df[[paste0("U_Gg_var_glob_", i, "_XX")]] / get(paste0("delta_D_",i,"_global_XX"))) +  
            (df[[paste0("U_Gg_var_glob_", j, "_XX")]] / get(paste0("delta_D_",j,"_global_XX")))
          }

          df[[paste0("U_Gg_var_", i, "_", j, "_2_XX")]] <- df[[paste0("U_Gg_var_", i, "_", j, "_XX")]]^2 * df$first_obs_by_gp_XX
          assign(paste0("var_sum_",i,"_",j,"_XX"), 
              sum(df[[paste0("U_Gg_var_", i, "_", j, "_2_XX")]], na.rm = TRUE)/
                  G_XX^2)
          assign(paste0("cov_",i,"_",j,"_XX"), 
          (get(paste0("var_sum_",i,"_",j,"_XX")) - get(paste0("se_",i,"_XX"))^2 - get(paste0("se_",j,"_XX"))^2)/2)

          didmgt_Var_Effects[i,j] <- get(paste0("cov_",i,"_",j,"_XX"))
          didmgt_Var_Effects[j,i] <- get(paste0("cov_",i,"_",j,"_XX"))
        }
      }

    }
  }

  # Computing the vector of recentered effects and its variance matrix
  didmgt_D <- didmgt_identity - matrix(1/l_XX, nrow = l_XX - 1, ncol = l_XX)
  didmgt_test_effects <- didmgt_D %*% didmgt_Effects
  didmgt_test_var <- didmgt_D %*% didmgt_Var_Effects %*% t(didmgt_D)
  didmgt_test_var <- (didmgt_test_var + t(didmgt_test_var)) / 2

  # Performing the test
  didmgt_chi2_equal_ef <- t(didmgt_test_effects) %*% Ginv(didmgt_test_var) %*% didmgt_test_effects
  p_equality_effects <- 
      1 - pchisq(didmgt_chi2_equal_ef[1,1], df = l_XX - 1)
  assign("p_equality_effects", p_equality_effects, envir = globalenv())
} else {
  warning("Some effects could not be estimated. Therefore, the test of equality of effects could not be computed.")
}
}
  
# Returning the results of the estimation ##################################

dis_mat_XX <- matrix(data = 0, nrow = nrow(mat_res_XX) , ncol = ncol(mat_res_XX) - 1)
dis_mat_XX[,1:4] <- sprintf("%s", format(round(mat_res_XX[,1:4], 5), big.mark=",", scientific=FALSE, trim=TRUE))
dis_mat_XX[,5:6] <- sprintf("%s", format(round(mat_res_XX[,5:6], 0), big.mark=",", scientific=FALSE, trim=TRUE))
rownames(dis_mat_XX) <- rownames
colnames(dis_mat_XX) <- c("Estimate", "SE", "LB CI", "UB CI", "N", "Switchers")


cat("\n")
cat(noquote(strrep("-", 80)));cat("\n");
cat(strrep(" ", 13));cat("Estimation of treatment effects: Event-study effects");cat("\n");
cat(noquote(strrep("-", 80)));cat("\n");
print(noquote(dis_mat_XX[1:l_XX, , drop = FALSE]))
cat("\n")
if (effects_equal == TRUE) {
    cat(sprintf("Test of equality of the effects : p-value = %.4f", p_equality_effects))
    cat("\n")
}

cat("\n")
cat(noquote(strrep("-", 80)));cat("\n");
cat(strrep(" ", 4));cat("Estimation of treatment effects: Average total effect per treatment unit");cat("\n");
cat(noquote(strrep("-", 80)));cat("\n");
print(noquote(dis_mat_XX[l_XX+1, , drop = FALSE]))
cat("\n")

if (l_placebo_XX != 0) {
  cat("\n")
  cat(noquote(strrep("-", 80)));cat("\n");
  cat(strrep(" ", 4));cat(" Testing the parallel trends and no anticipation assumptions");cat("\n");
  cat(noquote(strrep("-", 80)));cat("\n");
  print(noquote(dis_mat_XX[(l_XX+2):(dim(dis_mat_XX)[1]), , drop = FALSE]))
  cat("\n")
  cat(sprintf("Test of joint nullity of the placebos : p-value = %.4f", p_jointplacebo))
  cat("\n")
}

# Saving the results if requested 
if (!is.null(save_results)) {
  write.csv(dis_mat_XX, save_results, row.names = TRUE, col.names = TRUE)
}
df
})
}


#########################################################################################

did_multiplegt_dyn_core <- function(df, Y, G, T, D, effects, placebo, switchers_core = NULL, trends_nonparam, controls, same_switchers, same_switchers_pl, normalized) {

  if (switchers_core == "in") {
    l_u_a_XX <- min(L_u_XX, effects, na.rm = TRUE)
    if (placebo != 0) {
      if (!is.na(L_placebo_u_XX) & L_placebo_u_XX != 0) {
        l_placebo_u_a_XX <- min(placebo, L_placebo_u_XX)
      }
    }
    increase_XX <- 1
  }

  if (switchers_core == "out") {
    l_u_a_XX <- min(L_a_XX, effects, na.rm = TRUE)
    if (placebo != 0) {
      if (!is.na(L_placebo_u_XX) & L_placebo_u_XX != 0) {
        l_placebo_u_a_XX <- min(placebo, L_placebo_u_XX)
      }
    }
    increase_XX <- 0
  }

  levels_d_sq_XX <- levels(as.factor(df$d_sq_int_XX))

  # Estimating the DID_{+,l} or DID_{-,l}
  for (i in 1:l_u_a_XX) {

     df <- df %>% dplyr::select(-dplyr::any_of(c(
      paste0("distance_to_switch_",i,"_XX"), paste0("never_change_d_",i,"_XX"),
      paste0("N",increase_XX,"_t_",i,"_XX"), paste0("N",increase_XX,"_t_",i,"_g_XX"),
      paste0("N_gt_control_",i,"_XX"), paste0("diff_y_",i,"_XX"),
      paste0("diff_y_",i,"_XX_temp"), paste0("dummy_U_Gg",i,"_XX"),
      paste0("U_Gg",i,"_temp_XX"), paste0("U_Gg",i,"_XX"),
      paste0("count",i,"_core_XX"), paste0("mean_diff_y_",i,"nd_sq_t_XX"),
      paste0("mean_diff_y_",i,"d_sq_t_XX"), paste0("U_Gg",i,"_temp_var_XX"),
      paste0("U_Gg",i,"var_XX"),paste0("U_Gg",i,"var_2_XX"),
      paste0("count_var_",i,"_ntreat_XX_temp"), paste0("count_var_",i,"ntreat_XX"),
      paste0("count_var_",i,"_treat_XX_temp"), paste0("count_var_",i,"treat_XX"),
      paste0("avg_diff_y_",i,"tnp_XX"), paste0("count_diff_y_",i,"_nd_sq_t_XX"),
      paste0("count_diff_y_",i,"_d_sq_t_XX"), paste0("never_change_d_",i,"_wXX"),
      paste0("distance_to_switch_",i,"_wXX")
     ))) 

    df <- df %>% 
      group_by(.data$group_XX) %>% 
      mutate(!!paste0("diff_y_", i, "_XX") := .data$outcome_XX - lag(.data$outcome_XX, i)) %>% 
      ungroup()

    if (!is.null(controls)) {
      count_controls <- 0
      for (var in controls) {
        count_controls <- count_controls + 1
        df[paste0("diff_X", count_controls, "_", i, "_XX")]  <- df[[var]] - lag(df[[var]], i)
        for (l in levels_d_sq_XX) {
          if (get(paste0("useful_res_", l, "_XX")) > 1) {
            df[[paste0("diff_y_", i, "_XX")]] <- ifelse(df$d_sq_int_XX == l,              
              df[[paste0("diff_y_", i, "_XX")]] - get(paste0("coefs_sq_", l, "_XX"))[count_controls, 1] * df[[paste0("diff_X", count_controls, "_", i, "_XX")]]
              , df[[paste0("diff_y_", i, "_XX")]])              
          }
        }
      }
    }

    df[paste0("never_change_d_", i, "_XX")] <- as.numeric(df$F_g_XX > df$time_XX)
    df[[paste0("never_change_d_", i, "_XX")]] <- ifelse(is.na(df[[paste0("diff_y_", i, "_XX")]]), NA,  df[[paste0("never_change_d_", i, "_XX")]]) 
    df[paste0("never_change_d_", i, "_wXX")] <- df[[paste0("never_change_d_", i, "_XX")]] * df$N_gt_XX

    df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
    df <- df %>% group_by(.data$joint_trends_XX) %>%
    mutate(!!paste0("N_gt_control_", i, "_XX") := sum(.data[[paste0("never_change_d_", i, "_wXX")]], na.rm = TRUE)) %>% ungroup()

    if (same_switchers == TRUE) {
      df <- df[order(df$group_XX, df$time_XX), ]

      if (same_switchers_pl == TRUE) {
        df$relevant_y_missing_XX <- (is.na(df$outcome_XX) & df$time_XX >= df$F_g_XX - 1 - placebo & df$time_XX <= df$F_g_XX - 1 + effects)
        if (!is.null(controls)) {
          df$relevant_y_missing_XX <- ifelse(df$fd_X_all_non_missing_XX == 0 & df$time_XX >= df$F_g_XX - placebo & df$time_XX <= df$F_g_XX - 1 + effects, 1, df$relevant_y_missing_XX)
        }
        df <- df %>% group_by(.data$group_XX) %>%
            mutate(cum_fillin_XX = sum(.data$relevant_y_missing_XX, na.rm = TRUE))
        df$dum_fillin_temp_XX <- 
            df$cum_fillin_XX == 0 & df$time_XX == df$F_g_XX - 1 + effects
        df <- df %>% group_by(.data$group_XX) %>% 
            mutate(fillin_g_XX = sum(.data$dum_fillin_temp_XX, na.rm = TRUE))
        df[paste0("still_switcher_",i,"_XX")] <- 
            df$F_g_XX - 1 + effects <= df$T_g_XX & df$fillin_g_XX > 0
        df[paste0("distance_to_switch_", i, "_XX")] <- 
        ifelse(!is.na(df[[paste0("diff_y_", i, "_XX")]]),
        df[[paste0("still_switcher_", i, "_XX")]] == 1 &
        df$time_XX == df$F_g_XX- 1 + i & i <= df$L_g_XX & df$S_g_XX == increase_XX &
        df[[paste0("N_gt_control_", i, "_XX")]] > 0 & !is.na(df[[paste0("N_gt_control_", i, "_XX")]]), NA)
      } else {
        df$relevant_y_missing_XX <- (is.na(df$outcome_XX) & df$time_XX >= df$F_g_XX - 1 & df$time_XX <= df$F_g_XX - 1 + effects)
        if (!is.null(controls)) {
          df$relevant_y_missing_XX <- ifelse(df$fd_X_all_non_missing_XX == 0 & df$time_XX >= df$F_g_XX & df$time_XX <= df$F_g_XX - 1 + effects, 1, df$relevant_y_missing_XX)
        }
        df <- df %>% group_by(.data$group_XX) %>%
            mutate(cum_fillin_XX = sum(.data$relevant_y_missing_XX, na.rm = TRUE))
        df$dum_fillin_temp_XX <- 
            df$cum_fillin_XX == 0 & df$time_XX == df$F_g_XX - 1 + effects
        df <- df %>% group_by(.data$group_XX) %>% 
            mutate(fillin_g_XX = sum(.data$dum_fillin_temp_XX, na.rm = TRUE))
        df[paste0("still_switcher_",i,"_XX")] <- 
            df$F_g_XX - 1 + effects <= df$T_g_XX & df$fillin_g_XX > 0
        df[paste0("distance_to_switch_", i, "_XX")] <- 
        ifelse(!is.na(df[[paste0("diff_y_", i, "_XX")]]),
        df[[paste0("still_switcher_", i, "_XX")]] == 1 &
        df$time_XX == df$F_g_XX- 1 + i & i <= df$L_g_XX & df$S_g_XX == increase_XX &
        df[[paste0("N_gt_control_", i, "_XX")]] > 0 & !is.na(df[[paste0("N_gt_control_", i, "_XX")]]), NA)
      }
    }
    else {
      df[paste0("distance_to_switch_", i, "_XX")] <- 
      ifelse(!is.na(df[[paste0("diff_y_", i, "_XX")]]),
      df$time_XX == df$F_g_XX - 1 + i & i <= df$L_g_XX & df$S_g_XX == increase_XX & df[[paste0("N_gt_control_", i, "_XX")]] > 0 & !is.na(df[[paste0("N_gt_control_", i, "_XX")]]),
      NA)
    }
    df[paste0("distance_to_switch_", i, "_XX")] <- as.numeric(df[[paste0("distance_to_switch_", i, "_XX")]])

    df[paste0("distance_to_switch_", i, "_wXX")] <- 
    df[[paste0("distance_to_switch_", i, "_XX")]] * df$N_gt_XX
    df <- df %>% group_by(.data$time_XX) %>% mutate(!!paste0("N", increase_XX,"_t_", i, "_XX") := sum(.data[[paste0("distance_to_switch_", i, "_wXX")]], na.rm = TRUE))

    # Computing N^1_l and N^0_l
    assign(paste0("N",increase_XX,"_",i,"_XX"), 0, envir = globalenv())
    for (t in t_min_XX:T_max_XX) {
      assign(paste0("N",increase_XX,"_",i,"_XX"), 
        get(paste0("N",increase_XX,"_",i,"_XX")) + mean(df[[paste0("N", increase_XX,"_t_", i, "_XX")]][df$time_XX == t], na.rm = TRUE))
    }

    # Computing N^0_{t,l,g} and N^1_{t,l,g}
    df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
    df <- df %>% group_by(.data$joint_trends_XX) %>%
        mutate(!!paste0("N",increase_XX,"_t_",i,"_g_XX") := sum(.data[[paste0("distance_to_switch_",i,"_wXX")]], na.rm = TRUE)) %>% ungroup()

    # Computing the mean of differences of outcomes for non-treated and treated separately
    df[paste0("diff_y_", i,"_N_gt_XX")] <- df[[paste0("diff_y_", i,"_XX")]] * df$N_gt_XX
    df[paste0("dof_y_", i,"_N_gt_XX")]  <- as.numeric(df$N_gt_XX != 0 & !is.na(df[[paste0("diff_y_", i,"_XX")]]))

    df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
    df <- df %>% group_by(.data$joint_trends_XX) %>%
        mutate(!!paste0("mean_diff_y_",i,"_nd_sq_t_XX") := 
        mean(.data[[paste0("diff_y_",i,"_N_gt_XX")]][
          .data[[paste0("never_change_d_",i,"_XX")]] == 1 &
          .data[[paste0("N", increase_XX,"_t_", i, "_XX")]] > 0 &
          !is.na(.data[[paste0("N", increase_XX,"_t_", i, "_XX")]])
        ], na.rm = TRUE)) %>% ungroup()
    df[[paste0("mean_diff_y_",i,"_nd_sq_t_XX")]][!(
          df[[paste0("never_change_d_",i,"_XX")]] == 1 &
          df[[paste0("N", increase_XX,"_t_", i, "_XX")]] > 0 &
          !is.na(df[[paste0("N", increase_XX,"_t_", i, "_XX")]]))] <- NA

    df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
    df <- df %>% group_by(.data$joint_trends_XX) %>%
        mutate(!!paste0("count_diff_y_",i,"_nd_sq_t_XX") := 
        sum(.data[[paste0("dof_y_",i,"_N_gt_XX")]][
          .data[[paste0("never_change_d_",i,"_XX")]] == 1 &
          .data[[paste0("N", increase_XX,"_t_", i, "_XX")]] > 0 &
          !is.na(.data[[paste0("N", increase_XX,"_t_", i, "_XX")]])
        ], na.rm = TRUE)) %>% ungroup()
    df[[paste0("count_diff_y_",i,"_nd_sq_t_XX")]][!(
          df[[paste0("never_change_d_",i,"_XX")]] == 1 &
          df[[paste0("N", increase_XX,"_t_", i, "_XX")]] > 0 &
          !is.na(df[[paste0("N", increase_XX,"_t_", i, "_XX")]]))] <- NA

    df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
    df <- df %>% group_by(.data$joint_trends_XX) %>%
        mutate(!!paste0("mean_diff_y_",i,"_d_sq_t_XX") := 
        mean(.data[[paste0("diff_y_",i,"_N_gt_XX")]][
          .data[[paste0("distance_to_switch_",i,"_XX")]] == 1 
        ], na.rm = TRUE)) %>% ungroup()
    df[[paste0("mean_diff_y_",i,"_d_sq_t_XX")]][!(
          df[[paste0("distance_to_switch_",i,"_XX")]] == 1) |
          is.na(df[[paste0("distance_to_switch_",i,"_XX")]])] <- NA

    df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
    df <- df %>% group_by(.data$joint_trends_XX) %>%
        mutate(!!paste0("count_diff_y_",i,"_d_sq_t_XX") := 
        sum(.data[[paste0("dof_y_",i,"_N_gt_XX")]][
          .data[[paste0("distance_to_switch_",i,"_XX")]] == 1 
        ], na.rm = TRUE)) %>% ungroup()
    df[[paste0("count_diff_y_",i,"_d_sq_t_XX")]][!(
          df[[paste0("distance_to_switch_",i,"_XX")]] == 1) |
           is.na(df[[paste0("distance_to_switch_",i,"_XX")]])] <- NA

    if (get(paste0("N",increase_XX,"_",i,"_XX")) != 0) {
      df[paste0("dummy_U_Gg",i,"_XX")] <- as.numeric(i <= df$T_g_XX - 1)
    
      df[paste0("U_Gg",i,"_temp_XX")] <- df[[paste0("dummy_U_Gg",i,"_XX")]] * G_XX / get(paste0("N",increase_XX,"_",i,"_XX")) * as.numeric(df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df$N_gt_XX * (df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]])
      df[paste0("U_Gg",i,"_temp_XX")] <- df[[paste0("U_Gg",i,"_temp_XX")]] *  df[[paste0("diff_y_",i,"_XX")]]
      df <- df %>% group_by(.data$group_XX) %>% 
          mutate(!!paste0("U_Gg",i,"_XX") := sum(.data[[paste0("U_Gg",i,"_temp_XX")]], na.rm = TRUE))
      df[[paste0("U_Gg",i,"_XX")]] <- df[[paste0("U_Gg",i,"_XX")]] * df$first_obs_by_gp_XX

      df[paste0("count",i,"_core_XX")] <- ifelse(!is.na(df[[paste0("U_Gg",i,"_temp_XX")]]) & df[[paste0("U_Gg",i,"_temp_XX")]] != 0 | df[[paste0("U_Gg",i,"_temp_XX")]] == 0 & df[[paste0("diff_y_",i,"_XX")]]==0 & (df[[paste0("distance_to_switch_",i,"_XX")]] != 0 | df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]] != 0 & df[[paste0("never_change_d_",i,"_XX")]] != 0), df$N_gt_XX, 0)
      df[paste0("count",i,"_core_XX")] <- as.numeric(df[[paste0("count",i,"_core_XX")]])

      # Computing the "alternative" U_{G,g,l}

      df[paste0("U_Gg",i,"_temp_var_XX")] <- 0

      # For controls
      df[paste0("U_Gg",i,"_temp_var_XX")] <- ifelse(
        df[[paste0("never_change_d_",i,"_XX")]] == 1 & df[[paste0("count_diff_y_",i,"_nd_sq_t_XX")]] > 1 & !is.na(df[[paste0("count_diff_y_",i,"_nd_sq_t_XX")]]), 
      df[[paste0("dummy_U_Gg",i,"_XX")]] * G_XX / get(paste0("N",increase_XX,"_",i,"_XX")) *(df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * (df[[paste0("diff_y_",i,"_N_gt_XX")]] - df[[paste0("mean_diff_y_",i,"_nd_sq_t_XX")]]) * sqrt(df[[paste0("count_diff_y_",i,"_nd_sq_t_XX")]] /(df[[paste0("count_diff_y_",i,"_nd_sq_t_XX")]]-1)), 
      df[[paste0("U_Gg",i,"_temp_var_XX")]])

      df[paste0("U_Gg",i,"_temp_var_XX")] <- ifelse(
        df[[paste0("never_change_d_",i,"_XX")]] == 1 & df[[paste0("count_diff_y_",i,"_nd_sq_t_XX")]] == 1, 
      df[[paste0("dummy_U_Gg",i,"_XX")]] * G_XX / get(paste0("N",increase_XX,"_",i,"_XX")) *(df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df[[paste0("diff_y_",i,"_N_gt_XX")]], 
      df[[paste0("U_Gg",i,"_temp_var_XX")]])

      df[paste0("U_Gg",i,"_temp_var_XX")] <- ifelse(
        df[[paste0("distance_to_switch_",i,"_XX")]] == 1 & df[[paste0("count_diff_y_",i,"_d_sq_t_XX")]] > 1 & !is.na(df[[paste0("count_diff_y_",i,"_d_sq_t_XX")]]), 
      df[[paste0("dummy_U_Gg",i,"_XX")]] * G_XX / get(paste0("N",increase_XX,"_",i,"_XX")) *(df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * (df[[paste0("diff_y_",i,"_N_gt_XX")]] - df[[paste0("mean_diff_y_",i,"_d_sq_t_XX")]]) * sqrt(df[[paste0("count_diff_y_",i,"_d_sq_t_XX")]] /(df[[paste0("count_diff_y_",i,"_d_sq_t_XX")]]-1)), 
      df[[paste0("U_Gg",i,"_temp_var_XX")]])

      df[paste0("U_Gg",i,"_temp_var_XX")] <- ifelse(
        df[[paste0("distance_to_switch_",i,"_XX")]] == 1 & df[[paste0("count_diff_y_",i,"_d_sq_t_XX")]] == 1, 
      df[[paste0("dummy_U_Gg",i,"_XX")]] * G_XX / get(paste0("N",increase_XX,"_",i,"_XX")) *(df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df[[paste0("diff_y_",i,"_N_gt_XX")]], 
      df[[paste0("U_Gg",i,"_temp_var_XX")]])

      df[paste0("U_Gg",i,"_temp_var_XX")] <- as.numeric(df[[paste0("U_Gg",i,"_temp_var_XX")]])

      # Summing the U_{G,g,l} over time periods for each group
      df <- df %>% group_by(.data$group_XX) %>% 
          mutate(!!paste0("U_Gg",i,"_var_XX"):= sum(.data[[paste0("U_Gg", i,"_temp_var_XX")]], na.rm = TRUE))
    }

    if (normalized == TRUE) {
      #if (continuous == TRUE){
        df$sum_temp_XX <- ifelse(df$time_XX >= df$F_g_XX & df$time_XX <= df$F_g_XX - 1 + i & df$S_g_XX == increase_XX, df$treatment_XX - df$d_sq_XX, NA)
      #} else {}
      df <- df %>% group_by(.data$group_XX) %>% 
          mutate(!!paste0("sum_treat_until_",i,"_XX") := sum(.data$sum_temp_XX, na.rm = TRUE)) %>%
          dplyr::select(-.data$sum_temp_XX)
      df[[paste0("delta_D_",i,"_cum_temp_XX")]] <- ifelse(
        df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
        (df$N_gt_XX/get(paste0("N",increase_XX,"_",i,"_XX"))) * (
          df$S_g_XX * df[[paste0("sum_treat_until_",i,"_XX")]] +
          (1 - df$S_g_XX) * (-df[[paste0("sum_treat_until_",i,"_XX")]]) 
        ), NA)
      assign(paste0("delta_norm_",i,"_XX"), sum(df[[paste0("delta_D_",i,"_cum_temp_XX")]], na.rm = TRUE), envir = globalenv())    
    }
  }

  ## Computation of the Placebos 
  if (placebo != 0) {
    if (l_placebo_u_a_XX >= 1) {
      for (i in 1:l_placebo_u_a_XX) {

        df <- df %>% dplyr::select(-dplyr::any_of(c(
          paste0("diff_y_pl_",i,"_XX"),
          paste0("U_Gg_pl_",i,"_temp_XX"), paste0("U_Gg_placebo_",i,"_XX"),
          paste0("U_Gg_pl_",i,"_temp_var_XX"), paste0("U_Gg_pl_",i,"_var_XX"),
          paste0("mean_diff_y_pl_",i,"_nd_sq_t_XX"), 
          paste0("mean_diff_y_pl_",i,"_d_sq_t_XX"), 
          paste0("count_diff_y_pl_",i,"_nd_sq_t_XX"), 
          paste0("count_diff_y_pl_",i,"_d_sq_t_XX"), 
          paste0("dist_to_switch_pl_",i,"_XX"), paste0("never_change_d_pl_",i,"_XX"),
          paste0("N", increase_XX,"_t_placebo_",i,"_XX"),
          paste0("N", increase_XX,"_t_placebo_",i,"_g_XX"),
          paste0("N_gt_control_placebo_",i,"_XX"),
          paste0("dummy_U_Gg_pl_",i,"_XX"),
          paste0("never_change_d_pl_",i,"_wXX"),
          paste0("dist_to_switch_pl_",i,"_wXX")
        ))) 

    # The main trick to computing the placebo point estimates is:
    # 1. to place the corresponding outcome (y_{F_g-1} - y_{F_g - l - 1})) values in the same row of that (y_{F_g + l -1} - y_{F_g - 1}) of the symmetric DID_l. 
    # 2. The other variables, such as N_gt, N0_l or N1_l, remain unchanged, except that we have to check if diff_y_placebo ( = y_{F_g - 2l -2}- y_{F_g - l -1}) exists. 
    # 3. If y_{F_g - l -1} does not exist for a specific group, that group is excluded from the calculation, hence, for example, one always has #Switchers for DID_l>= #Switchers for DID_pl.

        df <- df %>% 
          group_by(.data$group_XX) %>% 
          mutate(!!paste0("diff_y_pl_", i, "_XX") := lag(.data$outcome_XX, 2*i) - lag(.data$outcome_XX, i)) %>% 
          ungroup()

        if (!is.null(controls)) {
          count_controls <- 0
          for (var in controls) {
            count_controls <- count_controls + 1
            df[paste0("diff_X", count_controls, "_placebo_", i, "_XX")]  <- lag(df[[var]], 2*i) - lag(df[[var]], i)
            for (l in levels_d_sq_XX) {
              if (get(paste0("useful_res_", l, "_XX")) > 1) {
                df[[paste0("diff_y_pl_", i, "_XX")]] <- ifelse(df$d_sq_int_XX == l,              
                  df[[paste0("diff_y_pl_", i, "_XX")]] - get(paste0("coefs_sq_", l, "_XX"))[count_controls, 1] * df[[paste0("diff_X", count_controls, "_placebo_", i, "_XX")]]
                  , df[[paste0("diff_y_pl_", i, "_XX")]])              
              }
            }
          }
        }

        df[paste0("never_change_d_pl_", i, "_XX")] <- df[[paste0("never_change_d_", i, "_XX")]] * (!is.na(df[[paste0("diff_y_pl_", i,"_XX")]]))
        df[paste0("never_change_d_pl_", i, "_wXX")] <- df[[paste0("never_change_d_pl_", i, "_XX")]] * df$N_gt_XX

        df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
        df <- df %>% group_by(.data$joint_trends_XX) %>%
        mutate(!!paste0("N_gt_control_placebo_", i, "_XX") := sum(.data[[paste0("never_change_d_pl_", i, "_wXX")]], na.rm = TRUE)) %>% ungroup()

        df[paste0("dist_to_switch_pl_", i, "_XX")] <- df[[paste0("distance_to_switch_", i, "_XX")]] * (!is.na(df[[paste0("diff_y_pl_",i,"_XX")]]))
        df[paste0("dist_to_switch_pl_", i, "_wXX")] <- 
        df[[paste0("dist_to_switch_pl_", i, "_XX")]] * df$N_gt_XX

        df <- df %>% group_by(.data$time_XX) %>% mutate(!!paste0("N", increase_XX,"_t_placebo_", i, "_XX") := sum(.data[[paste0("dist_to_switch_pl_", i, "_wXX")]], na.rm = TRUE))

        # Computing N^1_l and N^0_l
        assign(paste0("N",increase_XX,"_placebo_",i,"_XX"), 0, envir = globalenv())
        for (t in t_min_XX:T_max_XX) {
          assign(paste0("N",increase_XX,"_placebo_",i,"_XX"), 
            get(paste0("N",increase_XX,"_placebo_",i,"_XX")) + mean(df[[paste0("N", increase_XX,"_t_placebo_", i, "_XX")]][df$time_XX == t], na.rm = TRUE))
        }

        # Computing N^0_{t,l,g} and N^1_{t,l,g}
        df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
        df <- df %>% group_by(.data$joint_trends_XX) %>%
            mutate(!!paste0("N",increase_XX,"_t_placebo_",i,"_g_XX") := sum(.data[[paste0("dist_to_switch_pl_",i,"_wXX")]], na.rm = TRUE)) %>% ungroup()

        # Computing the mean of differences of outcomes for non-treated and treated separately
        df[paste0("diff_y_pl_", i,"_N_gt_XX")] <- df[[paste0("diff_y_pl_", i,"_XX")]] * df$N_gt_XX
        df[paste0("dof_y_pl_", i,"_N_gt_XX")]  <- as.numeric(df$N_gt_XX != 0 & !is.na(df[[paste0("diff_y_pl_", i,"_XX")]]))

        df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
        df <- df %>% group_by(.data$joint_trends_XX) %>%
            mutate(!!paste0("mean_diff_y_pl_",i,"_nd_sq_t_XX") := 
            mean(.data[[paste0("diff_y_pl_",i,"_N_gt_XX")]][
              .data[[paste0("never_change_d_pl_",i,"_XX")]] == 1 &
              .data[[paste0("N", increase_XX,"_t_placebo_", i, "_XX")]] > 0 &
              !is.na(.data[[paste0("N", increase_XX,"_t_placebo_", i, "_XX")]])
            ], na.rm = TRUE)) %>% ungroup()
        df[[paste0("mean_diff_y_pl_",i,"_nd_sq_t_XX")]][!(
              df[[paste0("never_change_d_pl_",i,"_XX")]] == 1 &
              df[[paste0("N", increase_XX,"_t_placebo_", i, "_XX")]] > 0 &
              !is.na(df[[paste0("N", increase_XX,"_t_placebo_", i, "_XX")]]))] <- NA

        df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
        df <- df %>% group_by(.data$joint_trends_XX) %>%
            mutate(!!paste0("count_diff_y_pl_",i,"_nd_sq_t_XX") := 
            sum(.data[[paste0("dof_y_pl_",i,"_N_gt_XX")]][
              .data[[paste0("never_change_d_pl_",i,"_XX")]] == 1 &
              .data[[paste0("N", increase_XX,"_t_placebo_", i, "_XX")]] > 0 &
              !is.na(.data[[paste0("N", increase_XX,"_t_placebo_", i, "_XX")]])
            ], na.rm = TRUE)) %>% ungroup()
        df[[paste0("count_diff_y_pl_",i,"_nd_sq_t_XX")]][!(
              df[[paste0("never_change_d_pl_",i,"_XX")]] == 1 &
              df[[paste0("N", increase_XX,"_t_placebo_", i, "_XX")]] > 0 &
              !is.na(df[[paste0("N", increase_XX,"_t_placebo_", i, "_XX")]]))] <- NA

        df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
        df <- df %>% group_by(.data$joint_trends_XX) %>%
            mutate(!!paste0("mean_diff_y_pl_",i,"_d_sq_t_XX") := 
            mean(.data[[paste0("diff_y_pl_",i,"_N_gt_XX")]][
              .data[[paste0("dist_to_switch_pl_",i,"_XX")]] == 1 
            ], na.rm = TRUE)) %>% ungroup()
        df[[paste0("mean_diff_y_pl_",i,"_d_sq_t_XX")]][!(
              df[[paste0("distance_to_switch_pl_",i,"_XX")]] == 1) |
              is.na(df[[paste0("dist_to_switch_pl_",i,"_XX")]])] <- NA

        df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
        df <- df %>% group_by(.data$joint_trends_XX) %>%
            mutate(!!paste0("count_diff_y_pl_",i,"_d_sq_t_XX") := 
            sum(.data[[paste0("dof_y_pl_",i,"_N_gt_XX")]][
              .data[[paste0("dist_to_switch_pl_",i,"_XX")]] == 1 
            ], na.rm = TRUE)) %>% ungroup()
        df[[paste0("count_diff_y_pl_",i,"_d_sq_t_XX")]][!(
              df[[paste0("distance_to_switch_",i,"_XX")]] == 1) |
              is.na(df[[paste0("dist_to_switch_pl_",i,"_XX")]])] <- NA

        df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] <- i <= df$T_g_XX - 1

        if (get(paste0("N",increase_XX,"_placebo_",i,"_XX")) != 0) {
        
          df[paste0("U_Gg_pl_",i,"_temp_XX")] <- 
          df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] *
           (G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) * df$N_gt_XX * (df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) *  df[[paste0("diff_y_pl_",i,"_XX")]]

          df <- df %>% group_by(.data$group_XX) %>% 
              mutate(!!paste0("U_Gg_placebo_",i,"_XX") := sum(.data[[paste0("U_Gg_pl_",i,"_temp_XX")]], na.rm = TRUE))
          df[[paste0("U_Gg_placebo_",i,"_XX")]] <- df[[paste0("U_Gg_placebo_",i,"_XX")]] * df$first_obs_by_gp_XX

          df[paste0("count",i,"_pl_core_XX")] <- ifelse(!is.na(df[[paste0("U_Gg_pl_",i,"_temp_XX")]]) & df[[paste0("U_Gg_pl_",i,"_temp_XX")]] != 0 | df[[paste0("U_Gg_pl_",i,"_temp_XX")]] == 0 & df[[paste0("diff_y_pl_",i,"_XX")]]==0 & (df[[paste0("dist_to_switch_pl_",i,"_XX")]] != 0 | df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]] != 0 & df[[paste0("never_change_d_pl_",i,"_XX")]] != 0), df$N_gt_XX, 0)
          df[paste0("count",i,"_pl_core_XX")] <- as.numeric(df[[paste0("count",i,"_pl_core_XX")]])

          # Computing the "alternative" U_{G,g,l}

          df[paste0("U_Gg_pl_",i,"_temp_var_XX")] <- 0

          # For controls
          df[paste0("U_Gg_pl_",i,"_temp_var_XX")] <- ifelse(
            df[[paste0("never_change_d_pl_",i,"_XX")]] == 1 & df[[paste0("count_diff_y_pl_",i,"_nd_sq_t_XX")]] > 1 & !is.na(df[[paste0("count_diff_y_pl_",i,"_nd_sq_t_XX")]]), 
          df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] * G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX")) *(df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) * (df$time_XX >= i + 2 & df$time_XX <= df$T_g_XX) * (df[[paste0("diff_y_pl_",i,"_N_gt_XX")]] - df[[paste0("mean_diff_y_pl_",i,"_nd_sq_t_XX")]]) * sqrt(df[[paste0("count_diff_y_pl_",i,"_nd_sq_t_XX")]] /(df[[paste0("count_diff_y_pl_",i,"_nd_sq_t_XX")]]-1)), 
          df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]])

          df[paste0("U_Gg_pl_",i,"_temp_var_XX")] <- ifelse(
            df[[paste0("never_change_d_pl_",i,"_XX")]] == 1 & df[[paste0("count_diff_y_pl_",i,"_nd_sq_t_XX")]] == 1, 
          df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] * G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX")) *(df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) * (df$time_XX >= i + 2 & df$time_XX <= df$T_g_XX) * df[[paste0("diff_y_pl_",i,"_N_gt_XX")]], 
          df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]])

          df[paste0("U_Gg_pl_",i,"_temp_var_XX")] <- ifelse(
            df[[paste0("dist_to_switch_pl_",i,"_XX")]] == 1 & df[[paste0("count_diff_y_pl_",i,"_d_sq_t_XX")]] > 1 & !is.na(df[[paste0("count_diff_y_pl_",i,"_d_sq_t_XX")]]), 
          df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] * G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX")) *(df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) * (df$time_XX >= i + 2 & df$time_XX <= df$T_g_XX) * (df[[paste0("diff_y_pl_",i,"_N_gt_XX")]] - df[[paste0("mean_diff_y_pl_",i,"_d_sq_t_XX")]]) * sqrt(df[[paste0("count_diff_y_pl_",i,"_d_sq_t_XX")]] /(df[[paste0("count_diff_y_pl_",i,"_d_sq_t_XX")]]-1)), 
          df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]])

          df[paste0("U_Gg_pl_",i,"_temp_var_XX")] <- ifelse(
            df[[paste0("dist_to_switch_pl_",i,"_XX")]] == 1 & df[[paste0("count_diff_y_pl_",i,"_d_sq_t_XX")]] == 1, 
          df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] * G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX")) *(df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) * (df$time_XX >= i + 2 & df$time_XX <= df$T_g_XX) * df[[paste0("diff_y_pl_",i,"_N_gt_XX")]], 
          df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]])

          df[paste0("U_Gg_pl_",i,"_temp_var_XX")] <- 
              as.numeric(df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]])

          # Summing the U_{G,g,l} over time periods for each group
          df <- df %>% group_by(.data$group_XX) %>% 
              mutate(!!paste0("U_Gg_pl_",i,"_var_XX"):= sum(.data[[paste0("U_Gg_pl_", i,"_temp_var_XX")]], na.rm = TRUE))
        }

        if (normalized == TRUE) {
          df$sum_temp_pl_XX <- ifelse(df$time_XX >= df$F_g_XX & df$time_XX <= df$F_g_XX - 1 + i & !is.na(df[[paste0("diff_y_pl_",i,"_XX")]]) & df$S_g_XX == increase_XX, df$treatment_XX - df$d_sq_XX, NA)
          df <- df %>% group_by(.data$group_XX) %>% 
              mutate(!!paste0("sum_treat_until_",i,"_pl_XX") := sum(.data$sum_temp_pl_XX, na.rm = TRUE)) %>%
              dplyr::select(-.data$sum_temp_pl_XX)
          df[[paste0("delta_D_pl_",i,"_cum_temp_XX")]] <- ifelse(
            df[[paste0("dist_to_switch_pl_",i,"_XX")]] == 1,
            (df$N_gt_XX/get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) * (
              df$S_g_XX * df[[paste0("sum_treat_until_",i,"_pl_XX")]] +
              (1 - df$S_g_XX) * (-df[[paste0("sum_treat_until_",i,"_pl_XX")]]) 
            ), NA)
          assign(paste0("delta_norm_pl_",i,"_XX"), sum(df[[paste0("delta_D_pl_",i,"_cum_temp_XX")]], na.rm = TRUE), envir = globalenv())    
        }
      }
    }
  }

  # End of the Placebo computation 

  assign(paste0("sum_N",increase_XX,"_l_XX"), 0, envir = globalenv())
  for (i in 1:l_u_a_XX) {
    assign(paste0("sum_N",increase_XX,"_l_XX"), get(paste0("sum_N",increase_XX,"_l_XX"))+
    get(paste0("N",increase_XX,"_",i,"_XX")))
  }

  df <- df %>% dplyr::select(-any_of(c(
    "U_Gg_XX", "U_Gg_num_XX", "U_Gg_den_XX", "U_Gg_num_var_XX", "U_Gg_var_XX"
  )))

  df$U_Gg_num_XX <- 0
  df$U_Gg_den_XX <- 0
  df$U_Gg_num_var_XX <- 0
  for (i in 1:l_u_a_XX) {
    if (get(paste0("N",increase_XX,"_",i,"_XX")) != 0) {

      df <- df %>% dplyr::select(-any_of(c(
        paste0("delta_D_",i,"_temp_XX"), paste0("delta_D_",i,"_XX")
        )))

      assign(paste0("w_",i,"_XX"), get(paste0("N",increase_XX,"_",i,"_XX")) / get(paste0("sum_N",increase_XX,"_l_XX")))

      df[paste0("delta_D_",i,"_temp_XX")] <- df$N_gt_XX/get(paste0("N",increase_XX,"_",i,"_XX")) * ((df$treatment_XX - df$d_sq_XX) * df$S_g_XX + (1 - df$S_g_XX) * (df$d_sq_XX - df$treatment_XX))
      df[[paste0("delta_D_",i,"_temp_XX")]] <- ifelse(df[[paste0("distance_to_switch_",i,"_XX")]] == 1, 
          df[[paste0("delta_D_",i,"_temp_XX")]], NA)
      df[[paste0("delta_D_",i,"_temp_XX")]][is.na(
        df[[paste0("delta_D_",i,"_temp_XX")]])] <- 0

      df[paste0("delta_D_",i,"_XX")] <- sum(df[[paste0("delta_D_",i,"_temp_XX")]], na.rm = TRUE)
      df <- df %>% dplyr::select(-.data[[paste0("delta_D_",i,"_temp_XX")]])

      df$U_Gg_num_XX <- df$U_Gg_num_XX + get(paste0("w_",i,"_XX")) * df[[paste0("U_Gg",i,"_XX")]]
      df$U_Gg_num_var_XX <- df$U_Gg_num_var_XX + get(paste0("w_",i,"_XX")) * df[[paste0("U_Gg",i,"_var_XX")]]
      df$U_Gg_den_XX <- df$U_Gg_den_XX + get(paste0("w_",i,"_XX")) * df[[paste0("delta_D_",i,"_XX")]]
    }
  }

  # Computing the U^+_{G,g}s.
  df$U_Gg_XX <- df$U_Gg_num_XX/df$U_Gg_den_XX
  df$U_Gg_var_XX <- df$U_Gg_num_var_XX/df$U_Gg_den_XX

  assign(paste0("sum_N",increase_XX,"_l_XX"), get(paste0("sum_N",increase_XX,"_l_XX")), envir = globalenv())
  for (i in 1:l_u_a_XX) {
    assign(paste0("N",increase_XX,"_",i,"_XX") ,get(paste0("N",increase_XX,"_",i,"_XX")), envir = globalenv())
  }
  if (placebo != 0) {
    if (l_placebo_u_a_XX >= 1) {
      for (i in 1:l_placebo_u_a_XX) {
        assign(paste0("N",increase_XX,"_placebo_",i,"_XX") ,get(paste0("N",increase_XX,"_placebo_",i,"_XX")), envir = globalenv())
      }
    }
  }

  df
}

did_save_sample <- function(df, Gn, Tn) {
  df_save <- subset(df, !is.na(df$G) & !is.na(df$T))
  df_save <- df_save %>% dplyr::select(.data$G, .data$T, .data$S_g_XX)
  df_save <- data.table::setnames(df_save, old = c("G", "T", "S_g_XX"), new = c(Gn, Tn, "did_sample"))
  df_save$did_sample <- ifelse(df_save$did_sample == 0, -1, df_save$did_sample)
  df_save$did_sample <- ifelse(is.na(df_save$did_sample), 0, df_save$did_sample)
  df_save$did_sample <- factor(df_save$did_sample, levels = c(0,1,-1), labels = c("Never-switcher", "Switcher-in", "Switchers-out"))
  df_save
}

did_multiplegt_dyn_design <- function(df, design_opt, weight) {
  if (length(design_opt) != 2) {
    stop("Syntax error in design option.")
  }
  des_p <- as.numeric(design_opt[1])
  des_path <- design_opt[2]
  des_n <- l_XX
  des_per <- des_p * 100

  df$F_g_plus_n_XX <- df$F_g_XX + des_n - 1
  df <- subset(df, df$time_XX >= df$F_g_XX - 1 & df$time_XX <= F_g_plus_n_XX)
  df <- df[order(df$group_XX, df$time_XX), ]
  df <- df %>% group_by(.data$group_XX) %>% 
      mutate(time_l_XX = row_number())  %>% ungroup()
  df <- df %>% dplyr::select(.data$group_XX, .data$time_l_XX, .data$weight_XX, .data$treatment_XX, .data$F_g_XX)

  if (!is.null(weight)) {
    df <- df %>% group_by(.data$group_XX) %>%
        mutate(g_weight = sum(weight_XX), na.rm = TRUE) 
  } else {
    df$g_weight_XX <- 1
  }
  df <- df %>% dplyr::select(-.data$weight_XX)

  max_time <- max(df$time_l_XX, na.rm = TRUE)
  treat_list <- c()
  treat_str <- ""
  for (i in 1:max_time) {
    df <- df %>% group_by(.data$group_XX) %>%
        mutate(!!paste0("treatment_XX",i) := mean(.data$treatment_XX[.data$time_l_XX == i])) %>%
        ungroup()
    treat_list <- c(treat_list, paste0("treatment_XX",i))
    treat_str <- paste0(treat_str,"treatment_XX",i,",")
  }
  treat_str <- substr(treat_str, 1, nchar(treat_str) - 1)
  df <- df %>% dplyr::select(-.data$time_l_XX, -.data$treatment_XX)
  df <- unique(df)
  for (var in treat_list) {
    df <- subset(df, !is.na(df[[var]]))
  }

  df$N_XX <- 1
  df$N_w_XX <- (df$g_weight_XX * df$N_XX) / sum(df$g_weight_XX, na.rm = TRUE)
  df <- df %>% dplyr::select(-.data$group_XX, -.data$g_weight_XX)
  df$treatments <- df[treat_list]
  df <- df %>% group_by(.data$treatments) %>%
      mutate(N_XX = sum(.data$N_XX, na.rm = TRUE)) %>%
      mutate(N_w_XX = sum(.data$N_w_XX, na.rm = TRUE))
  df <- df %>% dplyr::select(-.data$F_g_XX)
  df <- unique(df)
  tot_switch <- sum(df$N_XX, na.rm = TRUE)

  df <- df %>% dplyr::arrange(-.data$N_XX, .data$treatments) 
  df$cum_sum_XX <- cumsum(df$N_w_XX)
  df$in_table_XX <- as.numeric(df$cum_sum_XX <= des_p)
  df <- df[order(df$in_table_XX, df$cum_sum_XX), ]
  df <- df %>% group_by(.data$in_table_XX) %>% mutate(id_XX = row_number())
  df <- subset(df, df$in_table_XX == 1 | (df$in_table_XX == 0 & df$id_XX == 1))

  if (des_p < 1) {
    last_p <- 100 * min(df$cum_sum_XX[df$in_table_XX == 0])
  } else {
    last_p <- 100
  }
  df <- df %>% dplyr::arrange(-.data$N_XX, .data$treatments) 
  df <- df[c("N_XX", "N_w_XX", treat_list)]
  df$N_w_XX <- df$N_w_XX * 100

  coln <- c("N", "Share")
  rown <- c()
  desmat <- matrix(NA, nrow = dim(df)[1], ncol = 2 + 1 + l_XX)
  for (j in 1:(2 + 1 + l_XX)) {
    for (i in 1:dim(df)[1]) {
      if (j == 1) {
        rown <- c(rown, paste0("TreatPath",i))
      }
      desmat[i,j] <- as.numeric(df[i,j])
    }
    if (j > 2) {
      coln <- c(coln, paste0("ℓ=",j - 2 - 1))
    }
  }
  colnames(desmat) <- coln
  rownames(desmat) <- rown 
  
  desmat[, 2] <- sprintf("%s", format(round(desmat[,2], 2), big.mark=",", scientific=FALSE, trim=TRUE))
  
  if (des_path == "console") {
    cat("\n")
    cat(noquote(strrep("-", 80)));cat("\n");
    cat(strrep(" ", 13));cat(sprintf("Detection of treatment paths - %.0f periods after first switch", l_XX));cat("\n");
    cat(noquote(strrep("-", 80)));cat("\n");
    print(noquote(desmat[ , , drop = FALSE]))
    cat(sprintf("Treatment paths detected in at least %.2f%% of the %.0f switching groups for which %.0f effects could be estimated.", des_per, tot_switch, l_XX))
    cat(sprintf(" Total %% = %.2f%%", last_p));cat("\n");
    cat("Design interpretation (first row):")
    n_groups <- desmat[1,1]
    d_start <- desmat[1,3]
    d_vec <- "("
    for (i in 1:l_XX) {
      d_vec <- paste0(d_vec, desmat[1, 3 + i],",")
    }
    d_vec <- paste0(substr(d_vec,1,nchar(d_vec)-1),")")
    cat(sprintf(" %s groups started with treatment %s and then experienced treatment vector %s", n_groups, d_start, d_vec))
    cat("\n")
  }
  else {
    cat(sprintf("Design exported to %s", des_path)); cat("\n")
    write.csv(desmat, des_path, row.names = TRUE, col.names = TRUE)
  }
}

did_multiplegt_dyn_dfs <- function(df, dfs) {
  if (length(dfs) != 2) {
    stop("Syntax error in date_first_switch option.")
  }
  dfs_opt <- dfs[1]
  dfs_path <- dfs[2]

  if (dfs_opt != "" & dfs_opt != "by_baseline_treat") {
    stop("Only option by_baseline_treat allowed.")
  }

  df <- subset(df, !(df$F_g_XX == T_max_XX + 1 | is.na(F_g_XX)))
  df <- subset(df, df$time_XX == df$F_g_XX)
  df <- df %>% dplyr::select(.data$G, .data$T, .data$F_g_XX, .data$d_sq_XX)

  if (dfs_opt == "") {
    df$tot_s <- 1
    df <- df %>% group_by(.data$T) %>% 
    mutate(tot_s = sum(.data$tot_s, na.rm = TRUE)) %>% 
    dplyr::select(-.data$G, -.data$F_g_XX, -.data$d_sq_XX)
    df <- unique(df)
    df <- df[order(df$T), ]
    df$share_XX <- (df$tot_s / sum(df$tot_s, na.rm = TRUE)) * 100
    df <- df[c("tot_s", "share_XX", "T")]
    dfsmat <- matrix(NA, ncol = 2, nrow = dim(df)[1])
    rown <- c()
    coln <- c("N", "Share")
    for (j in 1:2) {
      for (i in 1:dim(df)[1]) {
        if (j == 1) {
          rown <- c(rown, df$T[i])
        }
        dfsmat[i,j] <- as.numeric(df[i,j])
      }
    }
    colnames(dfsmat) <- coln
    rownames(dfsmat) <- rown 
    
    dfsmat[, 2] <- sprintf("%s", format(round(dfsmat[,2], 2), big.mark=",", scientific=FALSE, trim=TRUE))
    if (dfs_path == "console") {
      cat("\n")
      cat(noquote(strrep("-", 40)));cat("\n");
      cat(strrep(" ", 7));cat("Switching dates");cat("\n");
      cat(noquote(strrep("-", 40)));cat("\n");
      cat("By any status quo treatment");cat("\n");
      print(noquote(dfsmat[ , , drop = FALSE]))
      cat("\n")
    }
    else {
      cat(sprintf("Switching dates exported to %s", dfs_path)); cat("\n")
      write.xlsx(dfsmat, dfs_path, row.names = TRUE, col.names = TRUE, sheet = "Switching Dates")
    }
  }
  if (dfs_opt == "by_baseline_treat") {
    df$tot_s <- 1
    df <- df %>% group_by(.data$T, .data$d_sq_XX) %>% 
    mutate(tot_s = sum(.data$tot_s, na.rm = TRUE)) %>% 
    dplyr::select(-.data$G, -.data$F_g_XX)
    df <- unique(df)
    df <- df[order(df$d_sq_XX, df$T), ]

    if (dfs_path == "console") {
      cat("\n")
      cat(noquote(strrep("-", 40)));cat("\n");
      cat(strrep(" ", 7));cat("Switching dates");cat("\n");
      cat(noquote(strrep("-", 40)));cat("\n");
    }

    levels_d_sq_XX <- levels(factor(df$d_sq_XX))
    new_file <- 1
    for (l in levels_d_sq_XX) {
      df_by <- subset(df, df$d_sq_XX == l)
      df_by$share_XX <- (df_by$tot_s / sum(df_by$tot_s, na.rm = TRUE)) * 100
      df_by <- df_by[c("tot_s", "share_XX", "T")]
      dfsmat <- matrix(NA, ncol = 2, nrow = dim(df_by)[1])
      rown <- c()
      coln <- c("N", "Share")
      for (j in 1:2) {
        for (i in 1:dim(df_by)[1]) {
          if (j == 1) {
            rown <- c(rown, df_by$T[i])
          }
          dfsmat[i,j] <- as.numeric(df_by[i,j])
        }
      }
      colnames(dfsmat) <- coln
      rownames(dfsmat) <- rown 
      
      dfsmat[, 2] <- sprintf("%s", format(round(dfsmat[,2], 2), big.mark=",", scientific=FALSE, trim=TRUE))
      if (dfs_path == "console") {
        cat(sprintf("Status quo treatment = %s", l));cat("\n");
        print(noquote(dfsmat[ , , drop = FALSE]))
        cat("\n")
      }
      else {
        sheetn <- paste0("Switch. dates base treat. ", l)
        if (new_file == 1) {
          write.xlsx(dfsmat, dfs_path, row.names = TRUE, col.names = TRUE, sheet = sheetn, append = FALSE)
        } else {
          write.xlsx(dfsmat, dfs_path, row.names = TRUE, col.names = TRUE, sheet = sheetn, append = TRUE)
        }
        new_file <- 0
      }
    }
    if (dfs_path != "console") {
      cat(sprintf("Switching dates exported to %s", dfs_path)); cat("\n")
    }
  }
}

did_multiplegt_dyn_normweights <- function(df, normalized, normopt) {
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
    coln <- c(coln, paste0("ℓ=",i))
    for (k in 0:(i-1)) {
      if (lag_d == "by_k") {
        row <- k + 1
      } else {
        row <- i - k # Visualization by F_g - 1 + (l - k)
      }

      df[paste0("delta_",i,"_",k)] <- ifelse(df$time_XX == df$F_g_XX - 1 + i - k & df$F_g_XX - 1 + i <= df$T_g_XX, abs(df$treatment_XX - df$d_sq_XX), NA)
      weight_mat[row, i] <- (sum(df[[paste0("delta_",i,"_",k)]], na.rm = TRUE) / get(paste0("delta_D_",i,"_global_XX"))) / get(paste0("N_switchers_effect_",i,"_XX"))
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

  cat("\n")
  cat(noquote(strrep("-", 80)));cat("\n");
  cat(strrep(" ", 7));cat(sprintf("Weights on treatment lags - by %s", by_opt_lag));cat("\n");
  cat(noquote(strrep("-", 80)));cat("\n");
  print(noquote(weight_mat))
  cat("\n")


}

env_clean <- function() {
  rm(list = ls(envir = globalenv())[grep("_XX$", ls(envir = globalenv()))], envir = globalenv())
  rm(list = ls(envir = globalenv())[grep("_XX_new$", ls(envir = globalenv()))], envir = globalenv())
}

did_multiplegt_dyn <- function(df, Y, G, T, D, effects = 1, placebo = 0, ci_level = 95, switchers = "", trends_nonparam = NULL, weight = NULL, controls = NULL, dont_drop_larger_lower = FALSE, save_sample = FALSE, drop_if_d_miss_before_first_switch = FALSE, cluster = NULL, same_switchers = FALSE, same_switchers_pl = FALSE, effects_equal = FALSE, save_results = NULL, normalized = FALSE, design = NULL, date_first_switch = NULL, normalized_weights = NULL)
{ 
  df_main <- df
  df_est <- did_multiplegt_main(df_main, Y, G, T, D, effects, placebo,
  ci_level, switchers, trends_nonparam, 
  weight, controls, dont_drop_larger_lower, 
  drop_if_d_miss_before_first_switch, cluster, 
  same_switchers, same_switchers_pl, effects_equal, 
  save_results, normalized)

  if (!is.null(design)) {
    did_multiplegt_dyn_design(df_est, design, weight)
  }

  if (!is.null(date_first_switch)) {
    did_multiplegt_dyn_dfs(df_est, date_first_switch)
  }

  if (!is.null(normalized_weights)) {
    did_multiplegt_dyn_normweights(df_est, normalized, normalized_weights)
  }

  if (save_sample == TRUE) {
    df_save_XX <- did_save_sample(df_est, G, T)
    df <- merge(df, df_save_XX, by = c(G, T))
  }
  env_clean()
  
  df_est
  #df

  cat("\n")
  cat("The development of this package was funded by the European Union (ERC, REALLYCREDIBLE,GA N°101043899).")
  cat("\n")
}
