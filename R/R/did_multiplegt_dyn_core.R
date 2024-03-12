#' Internal function of did_multiplegt_dyn that computes U_Gg_plus_XX, U_Gg_minus_XX, U_Gg_var_plus_XX, and U_Gg_var_minus_XX. These are essential variables for the computation of the DID_\ell estimators and their variances.
#' @param df df
#' @param outcome outcome
#' @param group group
#' @param time time
#' @param treatment treatment
#' @param effects effects
#' @param placebo placebp
#' @param switchers_core switchers_core
#' @param trends_nonparam trends_nonparam
#' @param controls controls
#' @param same_switchers same_switchers
#' @param same_switchers_pl same_switchers_pl
#' @param normalized normalized
#' @param globals globals
#' @param const constants
#' @param trends_lin trends_lin
#' @param controls_globals controls_globals
#' @param less_conservative_se less_conservative_se
#' @param continuous continuous
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang := 
#' @importFrom rlang .data
#' @returns A dataframe for the values of U_Gg_plus_XX, U_Gg_minus_XX, U_Gg_var_plus_XX, and U_Gg_var_minus_XX.
#' @noRd
did_multiplegt_dyn_core <- function(
    df, 
    outcome, 
    group, 
    time, 
    treatment, 
    effects, 
    placebo, 
    switchers_core = NULL, 
    trends_nonparam, 
    controls, 
    same_switchers, 
    same_switchers_pl, 
    normalized,
    globals,
    const,
    trends_lin,
    controls_globals,
    less_conservative_se,
    continuous
    ) {
  
  # Inherited Globals #
  L_u_XX <- globals$L_u_XX
  L_placebo_u_XX <- globals$L_placebo_u_XX
  L_placebo_a_XX <- globals$L_placebo_a_XX
  L_a_XX <-  globals$L_a_XX
  t_min_XX <- globals$t_min_XX
  T_max_XX <- globals$T_max_XX
  G_XX <- globals$G_XX
  t_min_XX <- globals$t_min_XX
  T_max_XX <- globals$T_max_XX

  for (e in names(const)) {
    assign(e, const[[e]])
  }

  if (!is.null(controls)) {
    for (e in names(controls_globals)) {
      assign(e, controls_globals[[e]])
    }
  }

  suppressWarnings({
  ####### 1. Scalars initialization

  ## Initializing the number of effects and placebos to estimate, depending on whether we consider switchers in or switchers out. 
  ## Initializing also a scalar for whether the estimation is for switchers in or out.

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
      if (!is.na(L_placebo_a_XX) & L_placebo_a_XX != 0) {
        l_placebo_u_a_XX <- min(placebo, L_placebo_a_XX)
      }
    }
    increase_XX <- 0
  }

  ## Initializing values of baseline treatment
  levels_d_sq_XX <- levels(as.factor(df$d_sq_int_XX))
  df <- df %>% dplyr::select(-dplyr::any_of(c("num_g_paths_0_XX", "cohort_fullpath_0_XX")))

  ## 2. Data preparation steps to generate variables necessary for computation of event-study effects	
  ## Loop over the number of dynamic effects

  for (i in 1:l_u_a_XX) {

     df <- df %>% dplyr::select(-dplyr::any_of(c(
      paste0("distance_to_switch_",i,"_XX"), paste0("never_change_d_",i,"_XX"),
      paste0("N",increase_XX,"_t_",i,"_XX"), paste0("N",increase_XX,"_t_",i,"_g_XX"),
      paste0("N_gt_control_",i,"_XX"), paste0("diff_y_",i,"_XX"),
      paste0("diff_y_",i,"_XX_temp"), paste0("dummy_U_Gg",i,"_XX"),
      paste0("U_Gg",i,"_temp_XX"), paste0("U_Gg",i,"_XX"),
      paste0("count",i,"_core_XX"), paste0("mean_diff_y_",i,"_nd_sq_t_XX"),
      paste0("mean_diff_y_",i,"_d_sq_t_XX"), paste0("U_Gg",i,"_temp_var_XX"),
      paste0("U_Gg",i,"_var_XX"),paste0("U_Gg",i,"_var_2_XX"),
      paste0("count_var_",i,"_ntreat_XX_temp"), paste0("count_var_",i,"_ntreat_XX"),
      paste0("count_var_",i,"_treat_XX_temp"), paste0("count_var_",i,"_treat_XX"),
      paste0("avg_diff_y_",i,"_tnp_XX"), paste0("count_diff_y_",i,"_nd_sq_t_XX"),
      paste0("count_diff_y_",i,"_d_sq_t_XX"), paste0("never_change_d_",i,"_wXX"),
      paste0("distance_to_switch_",i,"_wXX"),

      paste0("dof_cohort_",i,"_ns_t_XX"), paste0("dof_cohort_",i,"_s_t_XX"),
      paste0("dof_cohort_",i,"_s0_t_XX"), paste0("dof_cohort_",i,"_s1_t_XX"),
      paste0("dof_cohort_",i,"_s2_t_XX"),
      paste0("count_cohort_",i,"_ns_t_XX"), paste0("count_cohort_",i,"_s_t_XX"),
      paste0("count_cohort_",i,"_s0_t_XX"), paste0("count_cohort_",i,"_s1_t_XX"),
      paste0("count_cohort_",i,"_s2_t_XX"),
      paste0("total_cohort_",i,"_ns_t_XX"), paste0("total_cohort_",i,"_s_t_XX"),
      paste0("total_cohort_",i,"_s0_t_XX"), paste0("total_cohort_",i,"_s1_t_XX"),
      paste0("total_cohort_",i,"_s2_t_XX"),
      paste0("mean_cohort_",i,"_ns_t_XX"), paste0("mean_cohort_",i,"_s_t_XX"),
      paste0("mean_cohort_",i,"_s0_t_XX"), paste0("mean_cohort_",i,"_s1_t_XX"),
      paste0("mean_cohort_",i,"_s2_t_XX")
     ))) 

    ## Creating long difference of outcome
    df <- df[order(df$group_XX, df$time_XX), ]
    df <- df %>% 
      group_by(.data$group_XX) %>% 
      mutate(!!paste0("diff_y_", i, "_XX") := .data$outcome_XX - lag(.data$outcome_XX, i)) %>% 
      ungroup()

    ## Creating treatment paths if less_conservative_se option specified
    if (isTRUE(less_conservative_se)) {

      ## Creating a time-invariant, group-level variable, containing g's treatment at F_g-1+\ell
      df$d_fg_XX_temp <- ifelse(df$time_XX == df$F_g_XX +i-1,
          df$treatment_XX, NA)
      df <- df %>% group_by(.data$group_XX) %>% 
          mutate(!!paste0("d_fg",i,"_XX") := mean(.data$d_fg_XX_temp, na.rm = TRUE))

      ## This variable might be missing, for groups whose treatment never changes, and for groups not observed \ell periods after treatment change. Then, we impute their treatment at F_g-1+\ell-1. Inconsequential, just to avoid missing values. We also need to initialize a variable d_fg0_XX, when \ell=1.      
      if (i == 1) {
        df$d_fg0_XX <- df$d_sq_XX
        df <- df %>% group_by(.data$d_fg0_XX, .data$F_g_XX) %>%
          mutate(path_0_XX = cur_group_id()) %>% ungroup()
      }

      df[[paste0("d_fg",i,"_XX")]] <- ifelse(is.na(df[[paste0("d_fg",i,"_XX")]]),
          df[[paste0("d_fg",i-1,"_XX")]], df[[paste0("d_fg",i,"_XX")]])
      df <- df %>% group_by(.data[[paste0("path_",i-1,"_XX")]],
          .data[[paste0("d_fg",i,"_XX")]]) %>%
          mutate(!!paste0("path_",i,"_XX") := cur_group_id()) %>% ungroup()
      
      df$d_fg_XX_temp <- NULL

      ## For each group, define a variable counting how many groups belong to the same cohort, with cohorts defined as d_fg0_XX F_g_XX, as well as the full path.
      if (i == 1) {
        df <- df %>% group_by(.data$path_0_XX) %>% 
            mutate(num_g_paths_0_XX = n_distinct(.data$group_XX))
        df$cohort_fullpath_0_XX <- as.numeric(df$num_g_paths_0_XX > 1)
      }
      ## For each group, generate a dummy for whether at least two groups in their cohort.
      df <- df %>% group_by(.data[[paste0("path_",i,"_XX")]]) %>% 
          mutate(!!paste0("num_g_paths_",i,"_XX") := n_distinct(.data$group_XX))
      df[[paste0("cohort_fullpath_",i,"_XX")]] <- as.numeric(df[[paste0("num_g_paths_",i,"_XX")]] > 1)
    }

    ## Identifying the control (g,t)s in the estimation of dynamic effect i 
    df[paste0("never_change_d_", i, "_XX")] <- as.numeric(df$F_g_XX > df$time_XX)
    df[[paste0("never_change_d_", i, "_XX")]] <- ifelse(is.na(df[[paste0("diff_y_", i, "_XX")]]), NA,  df[[paste0("never_change_d_", i, "_XX")]]) 

    ## Creating N^g_t:
    ## number of control groups for g at t
    df[paste0("never_change_d_", i, "_wXX")] <- df[[paste0("never_change_d_", i, "_XX")]] * df$N_gt_XX
    df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
    df <- df %>% group_by(.data$joint_trends_XX) %>%
    mutate(!!paste0("N_gt_control_", i, "_XX") := sum(.data[[paste0("never_change_d_", i, "_wXX")]], na.rm = TRUE)) %>% ungroup()

    ###### Creating binary variable indicating whether g is \ell periods away from switch
    ## If the same_switchers option is specified:
    if (same_switchers == TRUE) {
      df <- df[order(df$group_XX, df$time_XX), ]
      ## If the same_switchers_pl option is specified:
      if (same_switchers_pl == TRUE) {
        ## Generate a variable tagging the switchers that should be dropped
        ## Is the case if at least one of the placebos or effects we try to estimate is missing:
        df$relevant_y_missing_XX <- (is.na(df$outcome_XX) & df$time_XX >= df$F_g_XX - 1 - placebo & df$time_XX <= df$F_g_XX - 1 + effects)
        ## Or if some of the controls are missing:
        if (!is.null(controls)) {
          df$relevant_y_missing_XX <- ifelse(df$fd_X_all_non_missing_XX == 0 & df$time_XX >= df$F_g_XX - 1 - placebo & df$time_XX <= df$F_g_XX - 1 + effects, 1, df$relevant_y_missing_XX)
        }
        df <- df %>% group_by(.data$group_XX) %>%
            mutate(cum_fillin_XX = sum(.data$relevant_y_missing_XX, na.rm = TRUE))
        df$dum_fillin_temp_XX <- 
            df$cum_fillin_XX == 0 & df$time_XX == df$F_g_XX - 1 + effects
        df <- df %>% group_by(.data$group_XX) %>% 
            mutate(fillin_g_XX = sum(.data$dum_fillin_temp_XX, na.rm = TRUE))
        #v1.0.1
        if (!is.null(placebo)) {
          df$dum_fillin_temp_pl_XX <- 
              df$cum_fillin_XX == 0 & df$time_XX == df$F_g_XX - 1 - placebo
          df <- df %>% group_by(.data$group_XX) %>% 
              mutate(fillin_g_pl_XX = sum(.data$dum_fillin_temp_pl_XX, na.rm = TRUE))
        }

        ## tag switchers who have no missings from F_g_XX-1-placebo to F_g_XX-1+effects
        df[paste0("still_switcher_",i,"_XX")] <- 
            df$F_g_XX - 1 + effects <= df$T_g_XX & df$fillin_g_XX > 0
        df[paste0("distance_to_switch_", i, "_XX")] <- 
        ifelse(!is.na(df[[paste0("diff_y_", i, "_XX")]]),
        df[[paste0("still_switcher_", i, "_XX")]] == 1 &
        df$time_XX == df$F_g_XX- 1 + i & i <= df$L_g_XX & df$S_g_XX == increase_XX &
        df[[paste0("N_gt_control_", i, "_XX")]] > 0 & !is.na(df[[paste0("N_gt_control_", i, "_XX")]]), NA)
      } else {
        ## If the same_switchers_pl option is not specified:
        ## Generate a variable tagging the switchers that should be dropped
        ## Is the case if at least one of the effects we try to estimate is missing:
        df$relevant_y_missing_XX <- (is.na(df$outcome_XX) & df$time_XX >= df$F_g_XX - 1 & df$time_XX <= df$F_g_XX - 1 + effects)
        ## Or if some of the controls are missing:
        if (!is.null(controls)) {
          df$relevant_y_missing_XX <- ifelse(df$fd_X_all_non_missing_XX == 0 & df$time_XX >= df$F_g_XX & df$time_XX <= df$F_g_XX - 1 + effects, 1, df$relevant_y_missing_XX)
        }
        df <- df %>% group_by(.data$group_XX) %>%
            mutate(cum_fillin_XX = sum(.data$relevant_y_missing_XX, na.rm = TRUE))
        df$dum_fillin_temp_XX <- 
            df$cum_fillin_XX == 0 & df$time_XX == df$F_g_XX - 1 + effects
        df <- df %>% group_by(.data$group_XX) %>% 
            mutate(fillin_g_XX = sum(.data$dum_fillin_temp_XX, na.rm = TRUE))
        ## tag switchers who have no missings from F_g_XX-1 to F_g_XX-1+effects
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
      ## If the same_switchers option is not specified
      df[paste0("distance_to_switch_", i, "_XX")] <- 
      ifelse(!is.na(df[[paste0("diff_y_", i, "_XX")]]),
      df$time_XX == df$F_g_XX - 1 + i & i <= df$L_g_XX & df$S_g_XX == increase_XX & df[[paste0("N_gt_control_", i, "_XX")]] > 0 & !is.na(df[[paste0("N_gt_control_", i, "_XX")]]),
      NA)
    }
    df[paste0("distance_to_switch_", i, "_XX")] <- as.numeric(df[[paste0("distance_to_switch_", i, "_XX")]])

    #### Creating a variable counting the number of groups \ell periods away from switch at t
    df[paste0("distance_to_switch_", i, "_wXX")] <- 
    df[[paste0("distance_to_switch_", i, "_XX")]] * df$N_gt_XX
    df <- df %>% group_by(.data$time_XX) %>% mutate(!!paste0("N", increase_XX,"_t_", i, "_XX") := sum(.data[[paste0("distance_to_switch_", i, "_wXX")]], na.rm = TRUE))

    #### Computing N^1_\ell/N^0_\ell.
    assign(paste0("N",increase_XX,"_",i,"_XX"), 0)
    for (t in t_min_XX:T_max_XX) {
      assign(paste0("N",increase_XX,"_",i,"_XX"), 
        get(paste0("N",increase_XX,"_",i,"_XX")) + mean(df[[paste0("N", increase_XX,"_t_", i, "_XX")]][df$time_XX == t], na.rm = TRUE))
    }

    #### Creating N^1_{t,\ell,g}/N^0_{t,\ell,g}: Variable counting number of groups \ell periods away from switch at t, and with same D_{g,1} and trends_nonparam.
    df <- df %>% group_by(.data$joint_trends_XX) %>%
        mutate(!!paste0("N",increase_XX,"_t_",i,"_g_XX") := sum(.data[[paste0("distance_to_switch_",i,"_wXX")]], na.rm = TRUE)) %>% ungroup()

    #### Creating all the adjustment terms to compute estimators with controls, and their variances
    if (!is.null(controls)) {

      #### Initialize intermediate Variable needed later	
      df[[paste0("part2_switch",increase_XX,"_",i,"_XX")]] <- 0

      ## generation of the T_d variable = max_{g:D_g,1 = d} F_g - 1: 
      ## last period when treatment effects can still be estimated for groups with baseline treatment equal to d
      df <- df %>% group_by(.data$d_sq_int_XX) %>% 
          mutate(T_d_XX = max(.data$F_g_XX, na.rm = TRUE))
      df$T_d_XX <- df$T_d_XX - 1

      ## Computing the long differences of the control variables (X_g_t - X_g_t-l)
      count_controls <- 0
      for (var in controls) {
        count_controls <- count_controls + 1
        df[paste0("diff_X",count_controls,"_", i, "_XX")]  <- df[[var]] - lag(df[[var]], i)

        ## Computing N_g_t * (X_g_t - X_g_t-l)
        df[paste0("diff_X",count_controls,"_",i,"_N_XX")] <- df$N_gt_XX * 
            df[[paste0("diff_X",count_controls,"_",i,"_XX")]]

        ## index l corresponds to d in the paper
        for (l in levels_d_sq_XX) { 

          ## intermediate variable to count the number of groups within each not yet switched cohort          
          df <- df %>% dplyr::select(-dplyr::any_of("dummy_XX"))
          df$dummy_XX <- as.numeric(df$F_g_XX > df$time_XX & df$d_sq_int_XX == l)

          ## Computing coordinates of vectors m^+_{g,d,\ell} and m^-_{g,d,\ell}
          # small m
          ## Creating variable inside the summation across t in m^+_{g,d,\ell}/m^-_{g,d,\ell}
          df[[paste0("m",increase_XX,"_g_",l,"_",count_controls,"_",i,"_XX")]] <- 
            (i <= df$T_g_XX - 2 & df$d_sq_int_XX == l) * (G_XX / get(paste0("N",increase_XX,"_",i,"_XX"))) * ((df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * (df[[paste0("diff_X",count_controls,"_",i,"_N_XX")]]))

          ## Summing that variable across t, and leaving one non missing observation per g	
          df <- df %>% group_by(.data$group_XX) %>% 
              mutate(!!paste0("m",increase_XX,"_",l,"_", count_controls,"_",i,"_XX") :=
                  sum(.data[[paste0("m",increase_XX,"_g_",l,"_", count_controls,"_",i,"_XX")]], na.rm = TRUE))
          df[[paste0("m",increase_XX,"_",l,"_",count_controls,"_",i,"_XX")]] <- ifelse(
            df$first_obs_by_gp_XX == 1, df[[paste0("m",increase_XX,"_",l,"_", count_controls,"_",i,"_XX")]], NA)
          
          ## Computing coordinates of vectors M^+_{d,\ell} and M^-_{d,\ell}
          df[paste0("M",increase_XX,"_",l,"_", count_controls,"_",i,"_XX")] <- 
              sum(df[[paste0("m",increase_XX,"_",l,"_", count_controls,"_",i,"_XX")]], na.rm = TRUE) / G_XX          

          ## number of groups within each not yet switched cohort
          df <- df %>% group_by(.data$time_XX, .data$d_sq_int_XX) %>% 
              mutate(!!paste0("E_hat_denom_", count_controls,"_", l, "_XX") :=
                sum(.data$dummy_XX[.data$d_sq_int_XX == l], na.rm = TRUE))
          df[[paste0("E_hat_denom_", count_controls,"_", l, "_XX")]] <- ifelse(
            df$d_sq_int_XX == l, df[[paste0("E_hat_denom_", count_controls,"_", l, "_XX")]], NA)

          ## Add the indicator for at least two groups in the cohort to E_y_hat_gt_`l'_XX (demeaning is possible)
          df[[paste0("E_y_hat_gt_",l,"_XX")]] <- df[[paste0("E_y_hat_gt_int_",l,"_XX")]] *
              (df[[paste0("E_hat_denom_",count_controls,"_",l,"_XX")]] >= 2)

          ## Computing the summation from t=2 to F_g-1 that appears in the last term 
          ## of U^{+,var,X}_{g,\ell} and U^{-,var,X}_{g,\ell}, defined in the companion paper.          
          df <- df %>% dplyr::select(-dplyr::any_of(c(
            paste0("N_c_",l,"_temp_XX"), paste0("N_c_",l,"_XX"),
            paste0("in_sum_temp_",count_controls,"_",l,"_XX"))))

          df[paste0("N_c_",l,"_temp_XX")] <- df$d_sq_int_XX == l & df$time_XX >= 2 &
              df$time_XX <= df$T_d_XX & df$time_XX < df$F_g_XX 
          df[paste0("N_c_",l,"_XX")] <- sum(df[[paste0("N_c_",l,"_temp_XX")]], na.rm = TRUE)

          df[paste0("in_sum_temp_",count_controls,"_",l,"_XX")] <- 
            (df[[paste0("prod_X",count_controls,"_Ngt_XX")]] * (1 + 
              (df[[paste0("E_hat_denom_", count_controls,"_",l,"_XX")]] >= 2) *
              (sqrt(df[[paste0("E_hat_denom_", count_controls,"_",l,"_XX")]] / 
              (df[[paste0("E_hat_denom_", count_controls,"_",l,"_XX")]] - 1)) -1)) *
              (df$diff_y_XX - df[[paste0("E_y_hat_gt_",l,"_XX")]]) *
              (df$time_XX >= 2 & df$time_XX <= df$F_g_XX - 1)) / df[[paste0("N_c_",l,"_XX")]]

          df <- df %>% group_by(.data$group_XX) %>% 
              mutate(!!paste0("in_sum_",count_controls,"_",l,"_XX") := 
              sum(.data[[paste0("in_sum_temp_",count_controls,"_",l,"_XX")]], na.rm=TRUE))

          ## Residualize the outcome difference wrt control differences:
          ## Yg,t − Yg,t−l − (Xg,t − Xg,t−l)*\theta_{Dg,1}
          if (get(paste0("useful_res_", l, "_XX")) > 1) {
            df[[paste0("diff_y_", i, "_XX")]] <- ifelse(df$d_sq_int_XX == l,              
              df[[paste0("diff_y_", i, "_XX")]] - get(paste0("coefs_sq_", l, "_XX"))[count_controls, 1] * df[[paste0("diff_X", count_controls, "_", i, "_XX")]]
              , df[[paste0("diff_y_", i, "_XX")]])              
              ## N.B. : in the above line, we do not add "&diff_X`count_controls'_`i'_XX!=." because we want to exclude from the estimation any first/long-difference for which the covariates are missing.

            ## Initialize intermediate Variable needed later
            df[[paste0("in_brackets_",l,"_",count_controls,"_XX")]] <- 0
          }
        }
      }
    }

    #### Computing the variables used for the demeaning of outcome's long difference diff_y_`i',
    #### which we will use to create the U_g^{var} variables, which are used to compute 
    #### the estimators' variances. 

    ## Generate long-differences of outcomes time N_{g,t}, and dummy for (g,t) such that
    ## diff_y_`i' non missing and N_gt non missing.

    df[paste0("diff_y_", i,"_N_gt_XX")] <- df[[paste0("diff_y_", i,"_XX")]] * df$N_gt_XX
    df[paste0("dof_y_", i,"_N_gt_XX")]  <- as.numeric(df$N_gt_XX != 0 & !is.na(df[[paste0("diff_y_", i,"_XX")]]))

    ##  For never switchers, demeaning wrt to cohorts defined by D_{g,1}, `trends_nonparam' (\mathcal{D}_k in companion paper)

    df <- joint_trends(df, c("d_sq_XX", "time_XX"), trends_nonparam)

    # Denominator of the mean
    df <- df %>% group_by(.data$joint_trends_XX) %>% 
        mutate(!!paste0("count_cohort_",i,"_ns_t_XX") := 
        sum(.data[["N_gt_XX"]][
          !is.na(.data[[paste0("diff_y_",i,"_XX")]]) &
          .data[[paste0("never_change_d_",i,"_XX")]] == 1 &
          .data[[paste0("N",increase_XX,"_t_",i,"_XX")]] > 0 &
          !is.na(.data[[paste0("N",increase_XX,"_t_",i,"_XX")]])
        ],na.rm = TRUE)) %>% ungroup()
    df[[paste0("count_cohort_",i,"_ns_t_XX")]] <- ifelse(
          !is.na(df[[paste0("diff_y_",i,"_XX")]]) &
          df[[paste0("never_change_d_",i,"_XX")]] == 1 &
          df[[paste0("N",increase_XX,"_t_",i,"_XX")]] > 0 &
          !is.na(df[[paste0("N",increase_XX,"_t_",i,"_XX")]]),
          df[[paste0("count_cohort_",i,"_ns_t_XX")]], NA)

    # Numerator of the mean
    df <- df %>% group_by(.data$joint_trends_XX) %>% 
        mutate(!!paste0("total_cohort_",i,"_ns_t_XX") := 
        sum(.data[[paste0("diff_y_",i,"_N_gt_XX")]][
          .data[[paste0("never_change_d_",i,"_XX")]] == 1 &
          .data[[paste0("N",increase_XX,"_t_",i,"_XX")]] > 0 &
          !is.na(.data[[paste0("N",increase_XX,"_t_",i,"_XX")]])
        ],na.rm = TRUE)) %>% ungroup()
    df[[paste0("total_cohort_",i,"_ns_t_XX")]] <- ifelse(
          df[[paste0("never_change_d_",i,"_XX")]] == 1 &
          df[[paste0("N",increase_XX,"_t_",i,"_XX")]] > 0 &
          !is.na(df[[paste0("N",increase_XX,"_t_",i,"_XX")]]),
          df[[paste0("total_cohort_",i,"_ns_t_XX")]], NA)

    ## Mean 
    df[[paste0("mean_cohort_",i,"_ns_t_XX")]] <- df[[paste0("total_cohort_",i,"_ns_t_XX")]] /
       df[[paste0("count_cohort_",i,"_ns_t_XX")]]

    # Counting number of groups for DOF adjustment
    df <- df %>% group_by(.data$joint_trends_XX) %>% 
        mutate(!!paste0("dof_cohort_",i,"_ns_t_XX") := 
        sum(.data[[paste0("dof_y_",i,"_N_gt_XX")]][
          !is.na(.data[[paste0("diff_y_",i,"_XX")]]) &
          .data[[paste0("never_change_d_",i,"_XX")]] == 1 &
          .data[[paste0("N",increase_XX,"_t_",i,"_XX")]] > 0 &
          !is.na(.data[[paste0("N",increase_XX,"_t_",i,"_XX")]])
        ],na.rm = TRUE)) %>% ungroup()
    df[[paste0("dof_cohort_",i,"_ns_t_XX")]] <- ifelse(
          !is.na(df[[paste0("diff_y_",i,"_XX")]]) &
          df[[paste0("never_change_d_",i,"_XX")]] == 1 &
          df[[paste0("N",increase_XX,"_t_",i,"_XX")]] > 0 &
          !is.na(df[[paste0("N",increase_XX,"_t_",i,"_XX")]]),
          df[[paste0("dof_cohort_",i,"_ns_t_XX")]], NA)

    ## For switchers, if option less_conservative_se not specified, demeaning wrt to cohorts defined by D_{g,1}, F_g, D_{g,F_g}, `trends_nonparam' (\mathcal{C}_k in companion paper).
    if (isFALSE(less_conservative_se)) {

    df <- joint_trends(df, c("d_sq_XX", "F_g_XX", "d_fg_XX", paste0("distance_to_switch_",i,"_XX")), trends_nonparam)

    # Denominator of the mean
    df <- df %>% group_by(.data$joint_trends_XX) %>% 
        mutate(!!paste0("count_cohort_",i,"_s_t_XX") := 
        sum(.data[["N_gt_XX"]],na.rm = TRUE)) %>% ungroup()
    df[[paste0("count_cohort_",i,"_s_t_XX")]] <- ifelse(
          df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
          df[[paste0("count_cohort_",i,"_s_t_XX")]], NA)

    # Numerator of the mean
    df <- df %>% group_by(.data$joint_trends_XX) %>% 
        mutate(!!paste0("total_cohort_",i,"_s_t_XX") := 
        sum(.data[[paste0("diff_y_",i,"_N_gt_XX")]],na.rm = TRUE)) %>% ungroup()
    df[[paste0("total_cohort_",i,"_s_t_XX")]] <- ifelse(
          df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
          df[[paste0("total_cohort_",i,"_s_t_XX")]], NA)

    ## Mean
    df[paste0("mean_cohort_",i,"_s_t_XX")] <- df[[paste0("total_cohort_",i,"_s_t_XX")]] /
         df[[paste0("count_cohort_",i,"_s_t_XX")]]

    # Counting number of groups for DOF adjustment
    df <- df %>% group_by(.data$joint_trends_XX) %>% 
        mutate(!!paste0("dof_cohort_",i,"_s_t_XX") := 
        sum(.data[[paste0("dof_y_",i,"_N_gt_XX")]],na.rm = TRUE)) %>% ungroup()
    df[[paste0("dof_cohort_",i,"_s_t_XX")]] <- ifelse(
          df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
          df[[paste0("dof_cohort_",i,"_s_t_XX")]], NA)

    } else {
      ## For switchers, if option less_conservative_se specified, demeaning wrt to cohorts defined by D_{g,1} F_g, D_{g,F_g},..., D_{g,F_g+\ell}, if that cohort has at least two groups, if not: demeaning wrt to cohorts defined by D_{g,1} F_g, D_{g,F_g}, if that cohort has at least two groups, if not: demeaning wrt D_{g,1} F_g.

      # by D_{g,1}, F_g, `trends_nonparam':
      df <- joint_trends(df, "path_0_XX", trends_nonparam)

      ### Denominator of the mean
      df <- df %>% group_by(.data$joint_trends_XX) %>% 
          mutate(!!paste0("count_cohort_",i,"_s0_t_XX") := 
          sum(.data[["N_gt_XX"]][
            .data[[paste0("distance_to_switch_",i,"_XX")]] == 1
          ],na.rm = TRUE)) %>% ungroup()
      df[[paste0("count_cohort_",i,"_s0_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
            df[[paste0("count_cohort_",i,"_s0_t_XX")]], NA)

      ### Numerator of the mean
      df <- df %>% group_by(.data$joint_trends_XX) %>% 
          mutate(!!paste0("total_cohort_",i,"_s0_t_XX") := 
          sum(.data[[paste0("diff_y_",i,"_N_gt_XX")]][
            .data[[paste0("distance_to_switch_",i,"_XX")]] == 1
          ],na.rm = TRUE)) %>% ungroup()
      df[[paste0("total_cohort_",i,"_s0_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
            df[[paste0("total_cohort_",i,"_s0_t_XX")]], NA)

      ### DOF
      df <- df %>% group_by(.data$joint_trends_XX) %>% 
          mutate(!!paste0("dof_cohort_",i,"_s0_t_XX") := 
          sum(.data[[paste0("dof_y_",i,"_N_gt_XX")]][
            .data[[paste0("distance_to_switch_",i,"_XX")]] == 1
          ],na.rm = TRUE)) %>% ungroup()
      df[[paste0("dof_cohort_",i,"_s0_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
        df[[paste0("dof_cohort_",i,"_s0_t_XX")]], NA)

      # by D_{g,1}, F_g, D_{g,F_g}, `trends_nonparam':
      df <- joint_trends(df, "path_1_XX", trends_nonparam)

      ### Denominator of the mean
      df <- df %>% group_by(.data$joint_trends_XX) %>% 
          mutate(!!paste0("count_cohort_",i,"_s1_t_XX") := 
          sum(.data[["N_gt_XX"]][
            .data[[paste0("distance_to_switch_",i,"_XX")]] == 1
          ],na.rm = TRUE))
      df[[paste0("count_cohort_",i,"_s1_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
            df[[paste0("count_cohort_",i,"_s1_t_XX")]], NA)

      ### Numerator of the mean
      df <- df %>% group_by(.data$joint_trends_XX) %>% 
          mutate(!!paste0("total_cohort_",i,"_s1_t_XX") := 
          sum(.data[[paste0("diff_y_",i,"_N_gt_XX")]][
            .data[[paste0("distance_to_switch_",i,"_XX")]] == 1
          ],na.rm = TRUE)) %>% ungroup()
      df[[paste0("total_cohort_",i,"_s1_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
            df[[paste0("total_cohort_",i,"_s1_t_XX")]], NA)

      ### DOF
      df <- df %>% group_by(.data$joint_trends_XX) %>% 
          mutate(!!paste0("dof_cohort_",i,"_s1_t_XX") := 
          sum(.data[[paste0("dof_y_",i,"_N_gt_XX")]][
            .data[[paste0("distance_to_switch_",i,"_XX")]] == 1
          ],na.rm = TRUE)) %>% ungroup()
      df[[paste0("dof_cohort_",i,"_s1_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
        df[[paste0("dof_cohort_",i,"_s1_t_XX")]], NA)

      # by D_{g,1}, F_g, D_{g,F_g},..., D_{g,F_g+\ell}, `trends_nonparam':
      df <- joint_trends(df, paste0("path_",i,"_XX"), trends_nonparam)

      ### Denominator of the mean
      df <- df %>% group_by(.data$joint_trends_XX) %>% 
          mutate(!!paste0("count_cohort_",i,"_s2_t_XX") := 
          sum(.data[["N_gt_XX"]][
            .data[[paste0("distance_to_switch_",i,"_XX")]] == 1
          ],na.rm = TRUE)) %>% ungroup()
      df[[paste0("count_cohort_",i,"_s2_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
            df[[paste0("count_cohort_",i,"_s2_t_XX")]], NA)

      ### Numerator of the mean
      df <- df %>% group_by(.data$joint_trends_XX) %>% 
          mutate(!!paste0("total_cohort_",i,"_s2_t_XX") := 
          sum(.data[[paste0("diff_y_",i,"_N_gt_XX")]][
            .data[[paste0("distance_to_switch_",i,"_XX")]] == 1
          ],na.rm = TRUE)) %>% ungroup()
      df[[paste0("total_cohort_",i,"_s2_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
            df[[paste0("total_cohort_",i,"_s2_t_XX")]], NA)

      ### DOF
      df <- df %>% group_by(.data$joint_trends_XX) %>% 
          mutate(!!paste0("dof_cohort_",i,"_s2_t_XX") := 
          sum(.data[[paste0("dof_y_",i,"_N_gt_XX")]][
            .data[[paste0("distance_to_switch_",i,"_XX")]] == 1
          ],na.rm = TRUE)) %>% ungroup()
      df[[paste0("dof_cohort_",i,"_s2_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
        df[[paste0("dof_cohort_",i,"_s2_t_XX")]], NA)

      ## Mean
      df[paste0("mean_cohort_",i,"_s_t_XX")] <- ifelse( df[[paste0("cohort_fullpath_",i,"_XX")]] == 1, df[[paste0("total_cohort_",i,"_s2_t_XX")]] / df[[paste0("count_cohort_",i,"_s2_t_XX")]], NA)
      df[[paste0("mean_cohort_",i,"_s_t_XX")]] <- ifelse(df[[paste0("cohort_fullpath_",i,"_XX")]] == 0 & df$cohort_fullpath_1_XX == 1, df[[paste0("total_cohort_",i,"_s1_t_XX")]] /df[[paste0("count_cohort_",i,"_s1_t_XX")]], df[[paste0("mean_cohort_",i,"_s_t_XX")]])
      df[[paste0("mean_cohort_",i,"_s_t_XX")]] <- ifelse(df$cohort_fullpath_1_XX == 0, df[[paste0("total_cohort_",i,"_s0_t_XX")]] /df[[paste0("count_cohort_",i,"_s0_t_XX")]], df[[paste0("mean_cohort_",i,"_s_t_XX")]])

      ## Counting number of groups for DOF adjustment
      df[paste0("dof_cohort_",i,"_s_t_XX")] <- ifelse( df[[paste0("cohort_fullpath_",i,"_XX")]] == 1, df[[paste0("dof_cohort_",i,"_s2_t_XX")]], NA)
      df[[paste0("dof_cohort_",i,"_s_t_XX")]] <- ifelse(df[[paste0("cohort_fullpath_",i,"_XX")]] == 0 & df$cohort_fullpath_1_XX == 1, df[[paste0("dof_cohort_",i,"_s1_t_XX")]], df[[paste0("dof_cohort_",i,"_s_t_XX")]])
      df[[paste0("dof_cohort_",i,"_s_t_XX")]] <- ifelse(df$cohort_fullpath_1_XX == 0, df[[paste0("dof_cohort_",i,"_s0_t_XX")]], df[[paste0("dof_cohort_",i,"_s_t_XX")]])
    }

    #### From those parts, generate variables for the demeaning and the DOF adjustment 
    ## E_hat_(g,t), defined from parts depending on the cohort definition 
    df[[paste0("E_hat_g_t_",i,"_XX")]] <- ifelse(df$F_g_XX > df$time_XX, df[[paste0("mean_cohort_",i,"_ns_t_XX")]] * (df[[paste0("dof_cohort_",i,"_ns_t_XX")]] >= 2), NA)
    df[[paste0("E_hat_g_t_",i,"_XX")]] <- ifelse(df$time_XX == df$F_g_XX -1 + i,  df[[paste0("mean_cohort_",i,"_s_t_XX")]] * (df[[paste0("dof_cohort_",i,"_s_t_XX")]] >= 2), df[[paste0("E_hat_g_t_",i,"_XX")]])

    ## DOF_(g,t) for DOF adjustement, defined from parts depending on the cohort definition 
    df[[paste0("DOF_g_t_",i,"_XX")]] <- ifelse(df$F_g_XX-1+i==df$time_XX,1+(df[[paste0("dof_cohort_",i,"_s_t_XX")]] >= 2) * ((sqrt( df[[paste0("dof_cohort_",i,"_s_t_XX")]]/(df[[paste0("dof_cohort_",i,"_s_t_XX")]] -1)))-1), NA)
    df[[paste0("DOF_g_t_",i,"_XX")]] <- ifelse(df$F_g_XX > df$time_XX, 1 + (df[[paste0("dof_cohort_",i,"_ns_t_XX")]] >= 2) * ((sqrt( df[[paste0("dof_cohort_",i,"_ns_t_XX")]]/(df[[paste0("dof_cohort_",i,"_ns_t_XX")]] -1)))-1), df[[paste0("DOF_g_t_",i,"_XX")]])

    ###### 3. Computing U_Gg_\ell variables
    #### If the dynamic effect can be estimated (as there are switchers), we compute the U_Gg_\ell variables etc.

    if (get(paste0("N",increase_XX,"_",i,"_XX")) != 0) {

      ## Creating a dummy variable indicating whether l<=T_g_XX-1
      df[paste0("dummy_U_Gg",i,"_XX")] <- as.numeric(i <= df$T_g_XX - 1)
    
      ## Computing U^+_{G,g,l}
      df[paste0("U_Gg",i,"_temp_XX")] <- df[[paste0("dummy_U_Gg",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_",i,"_XX"))) * as.numeric(df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df$N_gt_XX * (df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]])
      df[paste0("U_Gg",i,"_temp_XX")] <- df[[paste0("U_Gg",i,"_temp_XX")]] *  df[[paste0("diff_y_",i,"_XX")]]

      df <- df %>% group_by(.data$group_XX) %>% 
          mutate(!!paste0("U_Gg",i,"_XX") := sum(.data[[paste0("U_Gg",i,"_temp_XX")]], na.rm = TRUE))
      df[[paste0("U_Gg",i,"_XX")]] <- df[[paste0("U_Gg",i,"_XX")]] * df$first_obs_by_gp_XX

      # Counting the number of groups for which we can estimate U_Gg`i'_temp_XX - to help compute the "N" displayed by the command 
      df[paste0("count",i,"_core_XX")] <- ifelse(
        (!is.na(df[[paste0("U_Gg",i,"_temp_XX")]]) & df[[paste0("U_Gg",i,"_temp_XX")]] != 0) | 
        (df[[paste0("U_Gg",i,"_temp_XX")]] == 0 & df[[paste0("diff_y_",i,"_XX")]]==0 & (df[[paste0("distance_to_switch_",i,"_XX")]] != 0 | df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]] != 0 & df[[paste0("never_change_d_",i,"_XX")]] != 0)), df$N_gt_XX, 0)
      df[paste0("count",i,"_core_XX")] <- as.numeric(df[[paste0("count",i,"_core_XX")]])

      ## Computing U^(+,var)_{G,g,l}
      df[paste0("U_Gg",i,"_temp_var_XX")] <- 0

      df[[paste0("U_Gg",i,"_temp_var_XX")]] <- df[[paste0("dummy_U_Gg",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_",i,"_XX"))) * (df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df$N_gt_XX * df[[paste0("DOF_g_t_",i,"_XX")]] * (df[[paste0("diff_y_",i,"_XX")]] - df[[paste0("E_hat_g_t_",i,"_XX")]])

      ## Adding the additional part of U^(+,var,X)_{G,g,l}/U^(-,var,X)_{G,g,l} when controls are included: sum across values of baseline treatment d of M^+_(d,l)* a term in brackets in companion paper. 
      if (!is.null(controls)) {
        ## Loop over values of d_sq_int_XX:sum across values of baseline treatment  
        for (l in levels_d_sq_XX) {
          df[paste0("combined",increase_XX,"_temp_",l,"_",i,"_XX")] <- 0
          ## Loop over controls: inner product of M^+_(d,l)* term in brackets in companion paper.
          for (j in 1:count_controls) {
            ## Loop over k: computation of term in brackets in companion paper.
            for (k in 1:count_controls) {
              ## Computation of all cross products between elements of jth line of Den^{-1}_d and term multiplying it. Plus, summing over k, to have jth coordinate of vector Den^{-1}_d*...
              df[[paste0("in_brackets_",l,"_",j,"_XX")]] <- df[[paste0("in_brackets_",l,"_",j,"_XX")]] + get(paste0("inv_Denom_",l,"_XX"))[j,k] * df[[paste0("in_sum_",k,"_",l,"_XX")]] *
              (df$d_sq_int_XX == l & df$F_g_XX >= 3)
            }
            ## Withdrawing theta_d
            df[[paste0("in_brackets_",l,"_",j,"_XX")]] <- df[[paste0("in_brackets_",l,"_",j,"_XX")]] - get(paste0("coefs_sq_",l,"_XX"))[j,1]
          ## Computation of all cross products between coordinates of M^+_(d,l) and those of term in brackets
          df[[paste0("combined",increase_XX,"_temp_",l,"_",i,"_XX")]] <-  df[[paste0("combined",increase_XX,"_temp_",l,"_",i,"_XX")]] + df[[paste0("M",increase_XX,"_",l,"_",j,"_",i,"_XX")]] * df[[paste0("in_brackets_",l,"_",j,"_XX")]]
          }
          #### Final sum over the status quo treatment (outer sum over d in the formula)
          df[[paste0("part2_switch",increase_XX,"_",i,"_XX")]] <- as.numeric(df[[paste0("part2_switch",increase_XX,"_",i,"_XX")]] + df[[paste0("combined",increase_XX,"_temp_",l,"_",i,"_XX")]] * (df$d_sq_int_XX == l))
        }

        #### Making the adjustement to U^(+,var)_{G,g,l} when controls are included
        if (increase_XX == 1) {
          df[[paste0("U_Gg",i,"_temp_var_XX")]] <- df[[paste0("U_Gg",i,"_temp_var_XX")]] - df[[paste0("part2_switch1_",i,"_XX")]]
        } else {
          df[[paste0("U_Gg",i,"_temp_var_XX")]] <- df[[paste0("U_Gg",i,"_temp_var_XX")]] + df[[paste0("part2_switch0_",i,"_XX")]]
        }
      }

      df[[paste0("U_Gg",i,"_temp_var_XX")]] <- as.numeric(df[[paste0("U_Gg",i,"_temp_var_XX")]])

      # Summing the U_{G,g,l} over time periods for each group
      df <- df %>% group_by(.data$group_XX) %>% 
          mutate(!!paste0("U_Gg",i,"_var_XX"):= sum(.data[[paste0("U_Gg", i,"_temp_var_XX")]], na.rm = TRUE))
    }

    ###### 4. Computing adjustements for the normalized and trends_lin options 
    #### Compute \delta^D for normalized option
    if (normalized == TRUE) {
      if (is.null(continuous)){
        df$sum_temp_XX <- ifelse(df$time_XX >= df$F_g_XX & df$time_XX <= df$F_g_XX - 1 + i & df$S_g_XX == increase_XX, df$treatment_XX - df$d_sq_XX, NA)
      } else {
	      ## Redefine this with original treatment if continuous is defined (treatment was binarized and staggerized)
        df$sum_temp_XX <- ifelse(df$time_XX >= df$F_g_XX & df$time_XX <= df$F_g_XX - 1 + i & df$S_g_XX == increase_XX, df$treatment_XX_orig - df$d_sq_XX_orig, NA)
      }
      df <- df %>% group_by(.data$group_XX) %>% 
          mutate(!!paste0("sum_treat_until_",i,"_XX") := sum(.data$sum_temp_XX, na.rm = TRUE)) %>%
          dplyr::select(-.data$sum_temp_XX)
      df[[paste0("delta_D_",i,"_cum_temp_XX")]] <- ifelse(
        df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
        (df$N_gt_XX/get(paste0("N",increase_XX,"_",i,"_XX"))) * (
          df$S_g_XX * df[[paste0("sum_treat_until_",i,"_XX")]] +
          (1 - df$S_g_XX) * (-df[[paste0("sum_treat_until_",i,"_XX")]]) 
        ), NA)
      assign(paste0("delta_norm_",i,"_XX"), sum(df[[paste0("delta_D_",i,"_cum_temp_XX")]], na.rm = TRUE))    
    }
  }
  ## At this point we have all the U_Gg_`i'_XX we need. Thus, we can sum them now for the trends_lin option.

  #### Trends_lin option 
  #### As trends_lin relies on summing up the l effects we can only estimate effect l when we can estimate all prior l-1 effects, we verify if this condition is met.
  Ntrendslin <- 1
  for (i in 1:l_u_a_XX) {
    Ntrendslin <- min(Ntrendslin, get(paste0("N", increase_XX, "_",i,"_XX")), na.rm = TRUE)
  }

  #### Compute the U_Gg_\ell for trends_lin
  if (isTRUE(trends_lin) & Ntrendslin != 0) {
	  ## Initializing at 0
    df[paste0("U_Gg",l_u_a_XX,"_TL")] <- 0
    df[paste0("U_Gg",l_u_a_XX,"_var_TL")] <- 0
	  ## summing up the U_Gg's up to the l-th (current) effect
    for (i in 1:l_u_a_XX) {
      df[[paste0("U_Gg",l_u_a_XX,"_TL")]] <- ifelse(!is.na(df[[paste0("U_Gg",i,"_XX")]]),
      df[[paste0("U_Gg",l_u_a_XX,"_TL")]] + df[[paste0("U_Gg",i,"_XX")]],
      df[[paste0("U_Gg",l_u_a_XX,"_TL")]])
      df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] <- ifelse(!is.na(df[[paste0("U_Gg",i,"_var_XX")]]),
      df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] + df[[paste0("U_Gg",i,"_var_XX")]],
      df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]])
    }

	  ## replacing the U_Gg's with the one adjusted for group specific linear trends
    df[[paste0("U_Gg",l_u_a_XX,"_XX")]] <- df[[paste0("U_Gg",l_u_a_XX,"_TL")]] 
    df[[paste0("U_Gg",l_u_a_XX,"_var_XX")]] <- df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] 
    df[[paste0("U_Gg",l_u_a_XX,"_TL")]] <- NULL
    df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] <- NULL
  }

  ###### 5. Data preparation steps to generate variables necessary for computation of placebo effects	

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
          paste0("dist_to_switch_pl_",i,"_wXX"),
          paste0("dof_cohort_pl_",i,"_ns_t_XX"), paste0("dof_cohort_pl_",i,"_ns_t_XX"),
          paste0("count_cohort_pl_",i,"_ns_t_XX"), paste0("count_cohort_pl_",i,"_ns_t_XX"),
          paste0("total_cohort_pl_",i,"_ns_t_XX"), paste0("total_cohort_pl_",i,"_ns_t_XX"),
          paste0("mean_cohort_pl_",i,"_ns_t_XX"), paste0("mean_cohort_pl_",i,"_ns_t_XX"),

          paste0("dof_cohort_pl_",i,"_ns_t_XX"), paste0("dof_cohort_pl_",i,"_s_t_XX"),
          paste0("dof_cohort_pl_",i,"_s0_t_XX"), paste0("dof_cohort_pl_",i,"_s1_t_XX"),
          paste0("dof_cohort_pl_",i,"_s2_t_XX"),
          paste0("count_cohort_pl_",i,"_ns_t_XX"), paste0("count_cohort_pl_",i,"_s_t_XX"),
          paste0("count_cohort_pl_",i,"_s0_t_XX"), paste0("count_cohort_pl_",i,"_s1_t_XX"),
          paste0("count_cohort_pl_",i,"_s2_t_XX"),
          paste0("total_cohort_pl_",i,"_ns_t_XX"), paste0("total_cohort_pl_",i,"_s_t_XX"),
          paste0("total_cohort_pl_",i,"_s0_t_XX"), paste0("total_cohort_pl_",i,"_s1_t_XX"),
          paste0("total_cohort_pl_",i,"_s2_t_XX"),
          paste0("mean_cohort_pl_",i,"_ns_t_XX"), paste0("mean_cohort_pl_",i,"_s_t_XX"),
          paste0("mean_cohort_pl_",i,"_s0_t_XX"), paste0("mean_cohort_pl_",i,"_s1_t_XX"),
          paste0("mean_cohort_pl_",i,"_s2_t_XX")
        ))) 

        # The steps to compute the placebos are:
        # 1. to place the corresponding outcome (y_{F_g-1} - y_{F_g - l - 1})) values in the same row of that (y_{F_g + l -1} - y_{F_g - 1}) of the symmetric DID_l. 
        # 2. The other variables, such as N_gt, N0_l or N1_l, remain unchanged, except that we have to check if diff_y_placebo ( = y_{F_g - 2l -2}- y_{F_g - l -1}) exists. 

        ## Computing the long differences for the placebos
        df <- df %>% 
          group_by(.data$group_XX) %>% 
          mutate(!!paste0("diff_y_pl_", i, "_XX") := lag(.data$outcome_XX, 2*i) - lag(.data$outcome_XX, i)) %>% 
          ungroup()
        
        ## Identifying the controls (g,t)s in the estimation of placebo i
        df[paste0("never_change_d_pl_", i, "_XX")] <- df[[paste0("never_change_d_", i, "_XX")]] * (!is.na(df[[paste0("diff_y_pl_", i,"_XX")]]))
        df[paste0("never_change_d_pl_", i, "_wXX")] <- df[[paste0("never_change_d_pl_", i, "_XX")]] * df$N_gt_XX

        ## number of control groups for g at t
        df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
        df <- df %>% group_by(.data$joint_trends_XX) %>%
        mutate(!!paste0("N_gt_control_placebo_", i, "_XX") := sum(.data[[paste0("never_change_d_pl_", i, "_wXX")]], na.rm = TRUE)) %>% ungroup()

        #### Creating a variable counting the number of groups \ell periods away from switch at t
        df[paste0("dist_to_switch_pl_", i, "_XX")] <- df[[paste0("distance_to_switch_", i, "_XX")]] * (!is.na(df[[paste0("diff_y_pl_",i,"_XX")]]))
        # v1.0.1.
        if (isTRUE(same_switchers_pl)) {
          df[[paste0("dist_to_switch_pl_", i, "_XX")]] <- df[[paste0("dist_to_switch_pl_", i, "_XX")]] * df$fillin_g_pl_XX
        }

        df[paste0("dist_to_switch_pl_", i, "_wXX")] <- 
        df[[paste0("dist_to_switch_pl_", i, "_XX")]] * df$N_gt_XX

        ## Creating a variable counting the number of groups \ell periods away from switch at t
        df <- df %>% group_by(.data$time_XX) %>% mutate(!!paste0("N", increase_XX,"_t_placebo_", i, "_XX") := sum(.data[[paste0("dist_to_switch_pl_", i, "_wXX")]], na.rm = TRUE))

        #### Computing N^1_\ell/N^0_\ell. for the placebos
        ## Initializing the N1_`i'_XX/N0_`i'_XX scalar at 0. 
        assign(paste0("N",increase_XX,"_placebo_",i,"_XX"), 0)
        for (t in t_min_XX:T_max_XX) {
          assign(paste0("N",increase_XX,"_placebo_",i,"_XX"), 
            get(paste0("N",increase_XX,"_placebo_",i,"_XX")) + mean(df[[paste0("N", increase_XX,"_t_placebo_", i, "_XX")]][df$time_XX == t], na.rm = TRUE))
        }
        assign(paste0("N",increase_XX,"_placebo_",i,"_XX"), get(paste0("N",increase_XX,"_placebo_",i,"_XX")))

        ## Creating N^1_{t,\ell,g}/N^0_{t,\ell,g} for the placebos: Variable counting number of groups \ell periods away from switch at t, and with same D_{g,1} and trends_nonparam.
        df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
        df <- df %>% group_by(.data$joint_trends_XX) %>%
            mutate(!!paste0("N",increase_XX,"_t_placebo_",i,"_g_XX") := sum(.data[[paste0("dist_to_switch_pl_",i,"_wXX")]], na.rm = TRUE)) %>% ungroup()

        ## Creating all the adjustment terms to compute estimators with controls, and their variances 
        if (!is.null(controls)) {

          ## Initialize intermediate Variable needed later
          df[[paste0("part2_pl_switch",increase_XX,"_",i,"_XX")]] <- 0

          ## Computing the long differences of the control variables (X_g_t - X_g_t-l)
          count_controls <- 0
          for (var in controls) {
            count_controls <- count_controls + 1
            df[paste0("diff_X",count_controls,"_placebo_",i,"_XX")]  <- lag(df[[var]],2*i)-lag(df[[var]], i)
            ## Computing N_g_t * (X_g_t - X_g_t-l)
            df[paste0("diff_X",count_controls,"_pl_",i,"_N_XX")] <- df$N_gt_XX * 
                df[[paste0("diff_X",count_controls,"_placebo_",i,"_XX")]]

            ## index l corresponds to d in the paper
            for (l in levels_d_sq_XX) {
              
              df[[paste0("m",increase_XX,"_pl_g_",l,"_",count_controls,"_",i,"_XX")]] <- 
                (i <= df$T_g_XX - 2 & df$d_sq_int_XX == l) * (G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) * ((df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * (df[[paste0("diff_X",count_controls,"_pl_",i,"_N_XX")]]))

              df <- df %>% group_by(.data$group_XX) %>% 
                  mutate(!!paste0("m_pl",increase_XX,"_",l,"_", count_controls,"_",i,"_XX") := sum(.data[[paste0("m",increase_XX,"_pl_g_",l,"_", count_controls,"_",i,"_XX")]], na.rm = TRUE))
              df[[paste0("m_pl",increase_XX,"_",l,"_",count_controls,"_",i,"_XX")]] <- ifelse(df$first_obs_by_gp_XX == 1, df[[paste0("m_pl",increase_XX,"_",l,"_", count_controls,"_",i,"_XX")]], NA)
              
              df[paste0("M_pl",increase_XX,"_",l,"_", count_controls,"_",i,"_XX")] <- 
                  sum(df[[paste0("m_pl",increase_XX,"_",l,"_", count_controls,"_",i,"_XX")]], na.rm = TRUE) / G_XX          

              if (get(paste0("useful_res_", l, "_XX")) > 1) {
                df[[paste0("diff_y_pl_", i, "_XX")]] <- ifelse(df$d_sq_int_XX == l,              
                  df[[paste0("diff_y_pl_", i, "_XX")]] - get(paste0("coefs_sq_", l, "_XX"))[count_controls, 1] * df[[paste0("diff_X", count_controls, "_placebo_", i, "_XX")]]
                  , df[[paste0("diff_y_pl_", i, "_XX")]])              
                df[[paste0("in_brackets_pl_",l,"_",count_controls,"_XX")]] <- 0

              }
            }
          }
        }
        ## Computing the variables used for the demeaning of outcome's placebo long difference diff_y_`i',
        ## which we will use to create the U_g^{var} variables, which are used to compute 
        ## the placebos' variances. 

        ## Generate placebo long-differences of outcomes time N_{g,t}, and dummy for (g,t) such that
        ## diff_y_`i' non missing and N_gt non missing.

        df[paste0("diff_y_pl_", i,"_N_gt_XX")] <- df[[paste0("diff_y_pl_", i,"_XX")]] * df$N_gt_XX
        df[paste0("dof_y_pl_", i,"_N_gt_XX")]  <- as.numeric(df$N_gt_XX != 0 & !is.na(df[[paste0("diff_y_pl_", i,"_XX")]]))

        ## For never switchers, demeaning wrt to cohorts defined by D_{g,1}, `trends_nonparam' \mathcal{D}_k in companion paper)

        df <- joint_trends(df, c("d_sq_XX", "time_XX"), trends_nonparam)

        # Denominator of the mean
        df <- df %>% group_by(.data$joint_trends_XX) %>% 
            mutate(!!paste0("count_cohort_pl_",i,"_ns_t_XX") := 
            sum(.data[["N_gt_XX"]][
              !is.na(.data[[paste0("diff_y_pl_",i,"_XX")]]) &
              .data[[paste0("never_change_d_pl_",i,"_XX")]] == 1 &
              .data[[paste0("N",increase_XX,"_t_placebo_",i,"_XX")]] > 0 &
              !is.na(.data[[paste0("N",increase_XX,"_t_placebo_",i,"_XX")]])
            ],na.rm = TRUE)) %>% ungroup()
        df[[paste0("count_cohort_pl_",i,"_ns_t_XX")]] <- ifelse(
              !is.na(df[[paste0("diff_y_pl_",i,"_XX")]]) &
              df[[paste0("never_change_d_pl_",i,"_XX")]] == 1 &
              df[[paste0("N",increase_XX,"_t_placebo_",i,"_XX")]] > 0 &
              !is.na(df[[paste0("N",increase_XX,"_t_placebo_",i,"_XX")]]),
              df[[paste0("count_cohort_pl_",i,"_ns_t_XX")]], NA)

        # Numerator of the mean
        df <- df %>% group_by(.data$joint_trends_XX) %>% 
            mutate(!!paste0("total_cohort_pl_",i,"_ns_t_XX") := 
            sum(.data[[paste0("diff_y_pl_",i,"_N_gt_XX")]][
              .data[[paste0("never_change_d_pl_",i,"_XX")]] == 1 &
              .data[[paste0("N",increase_XX,"_t_placebo_",i,"_XX")]] > 0 &
              !is.na(.data[[paste0("N",increase_XX,"_t_placebo_",i,"_XX")]])
            ],na.rm = TRUE)) %>% ungroup()
        df[[paste0("total_cohort_pl_",i,"_ns_t_XX")]] <- ifelse(
              df[[paste0("never_change_d_pl_",i,"_XX")]] == 1 &
              df[[paste0("N",increase_XX,"_t_placebo_",i,"_XX")]] > 0 &
              !is.na(df[[paste0("N",increase_XX,"_t_placebo_",i,"_XX")]]),
              df[[paste0("total_cohort_pl_",i,"_ns_t_XX")]], NA)

        # Mean
        df[[paste0("mean_cohort_pl_",i,"_ns_t_XX")]] <- df[[paste0("total_cohort_pl_",i,"_ns_t_XX")]] /
          df[[paste0("count_cohort_pl_",i,"_ns_t_XX")]]

        # Counting number of groups for DOF adjustment
        df <- df %>% group_by(.data$joint_trends_XX) %>% 
            mutate(!!paste0("dof_cohort_pl_",i,"_ns_t_XX") := 
            sum(.data[[paste0("dof_y_pl_",i,"_N_gt_XX")]][
              !is.na(.data[[paste0("diff_y_pl_",i,"_XX")]]) &
              .data[[paste0("never_change_d_pl_",i,"_XX")]] == 1 &
              .data[[paste0("N",increase_XX,"_t_placebo_",i,"_XX")]] > 0 &
              !is.na(.data[[paste0("N",increase_XX,"_t_placebo_",i,"_XX")]])
            ],na.rm = TRUE)) %>% ungroup()
        df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]] <- ifelse(
              !is.na(df[[paste0("diff_y_pl_",i,"_XX")]]) &
              df[[paste0("never_change_d_pl_",i,"_XX")]] == 1 &
              df[[paste0("N",increase_XX,"_t_placebo_",i,"_XX")]] > 0 &
              !is.na(df[[paste0("N",increase_XX,"_t_placebo_",i,"_XX")]]),
              df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]], NA)

        ## For switchers, for placebos we no longer need to distinguish depending on whether the option less_conservative_se specified or not, we always demean wrt to cohorts defined by D_{g,1}, F_g, D_{g,F_g}, `trends_nonparam' (\mathcal{C}_k in companion paper).

        df <- joint_trends(df, c("d_sq_XX", "F_g_XX", "d_fg_XX", paste0("dist_to_switch_pl_",i,"_XX")), trends_nonparam)

        # Denominator
        df <- df %>% group_by(.data$joint_trends_XX) %>% 
            mutate(!!paste0("count_cohort_pl_",i,"_s_t_XX") := 
            sum(.data[["N_gt_XX"]],na.rm = TRUE)) %>% ungroup()
        df[[paste0("count_cohort_pl_",i,"_s_t_XX")]] <- ifelse(
              df[[paste0("dist_to_switch_pl_",i,"_XX")]] == 1,
              df[[paste0("count_cohort_pl_",i,"_s_t_XX")]], NA)

        # Numerator
        df <- df %>% group_by(.data$joint_trends_XX) %>% 
            mutate(!!paste0("total_cohort_pl_",i,"_s_t_XX") := 
            sum(.data[[paste0("diff_y_pl_",i,"_N_gt_XX")]],na.rm = TRUE)) %>% ungroup()
        df[[paste0("total_cohort_pl_",i,"_s_t_XX")]] <- ifelse(
              df[[paste0("dist_to_switch_pl_",i,"_XX")]] == 1,
              df[[paste0("total_cohort_pl_",i,"_s_t_XX")]], NA)

        df[paste0("mean_cohort_pl_",i,"_s_t_XX")] <- df[[paste0("total_cohort_pl_",i,"_s_t_XX")]] /
            df[[paste0("count_cohort_pl_",i,"_s_t_XX")]]

        # Counting number of groups for DOF adjustment
        df <- df %>% group_by(.data$joint_trends_XX) %>% 
            mutate(!!paste0("dof_cohort_pl_",i,"_s_t_XX") := 
            sum(.data[[paste0("dof_y_pl_",i,"_N_gt_XX")]],na.rm = TRUE)) %>% ungroup()
        df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]] <- ifelse(
              df[[paste0("dist_to_switch_pl_",i,"_XX")]] == 1,
              df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]], NA)

        #### From those parts, generate variables for the demeaning and the DOF adjustment 
        ## E_hat_pl_(g,t), defined from parts depending on the cohort definition 
        df[[paste0("E_hat_g_t_pl_",i,"_XX")]] <- ifelse(df$F_g_XX > df$time_XX, df[[paste0("mean_cohort_pl_",i,"_ns_t_XX")]] * (df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]] >= 2), NA)
        df[[paste0("E_hat_g_t_pl_",i,"_XX")]] <- ifelse(df$time_XX == df$F_g_XX -1 + i,  df[[paste0("mean_cohort_pl_",i,"_s_t_XX")]] * (df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]] >= 2), df[[paste0("E_hat_g_t_pl_",i,"_XX")]])

        ## DOF_pl_(g,t) for DOF adjustement, defined from parts depending on the cohort definition 
        df[[paste0("DOF_g_t_pl_",i,"_XX")]] <- ifelse(df$F_g_XX-1+i==df$time_XX,1+(df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]] >= 2) * ((sqrt( df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]]/(df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]] -1)))-1), NA)
        df[[paste0("DOF_g_t_pl_",i,"_XX")]] <- ifelse(df$F_g_XX > df$time_XX, 1 + (df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]] >= 2) * ((sqrt( df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]]/(df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]] -1)))-1), df[[paste0("DOF_g_t_pl_",i,"_XX")]])

        #### 6. Computing U_Gg_\ell variables for the placebos (similar to part 4, less commented)
        df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] <- i <= df$T_g_XX - 1

        if (get(paste0("N",increase_XX,"_placebo_",i,"_XX")) != 0) {

        
          df[paste0("U_Gg_pl_",i,"_temp_XX")] <- 
          df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] *
           (G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) * df$N_gt_XX * (df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) *  df[[paste0("diff_y_pl_",i,"_XX")]]* (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) 

          df <- df %>% group_by(.data$group_XX) %>% 
              mutate(!!paste0("U_Gg_placebo_",i,"_XX") := sum(.data[[paste0("U_Gg_pl_",i,"_temp_XX")]], na.rm = TRUE))
          df[[paste0("U_Gg_placebo_",i,"_XX")]] <- df[[paste0("U_Gg_placebo_",i,"_XX")]] * df$first_obs_by_gp_XX

          df[paste0("count",i,"_pl_core_XX")] <- ifelse(!is.na(df[[paste0("U_Gg_pl_",i,"_temp_XX")]]) & df[[paste0("U_Gg_pl_",i,"_temp_XX")]] != 0 | df[[paste0("U_Gg_pl_",i,"_temp_XX")]] == 0 & df[[paste0("diff_y_pl_",i,"_XX")]]==0 & (df[[paste0("dist_to_switch_pl_",i,"_XX")]] != 0 | df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]] != 0 & df[[paste0("never_change_d_pl_",i,"_XX")]] != 0), df$N_gt_XX, 0)
          df[paste0("count",i,"_pl_core_XX")] <- as.numeric(df[[paste0("count",i,"_pl_core_XX")]])

          df[paste0("U_Gg_pl_",i,"_temp_var_XX")] <- 0

          df[paste0("U_Gg_pl_",i,"_temp_var_XX")] <- df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) * (df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df$N_gt_XX * df[[paste0("DOF_g_t_pl_",i,"_XX")]] * (df[[paste0("diff_y_pl_",i,"_XX")]] - df[[paste0("E_hat_g_t_pl_",i,"_XX")]])

          if (!is.null(controls)) {
            for (l in levels_d_sq_XX) {
              df[paste0("combined_pl",increase_XX,"_temp_",l,"_",i,"_XX")] <- 0
              for (j in 1:count_controls) {
                for (k in 1:count_controls) {
                  df[[paste0("in_brackets_pl_",l,"_",j,"_XX")]] <- df[[paste0("in_brackets_pl_",l,"_",j,"_XX")]] + get(paste0("inv_Denom_",l,"_XX"))[j,k] * df[[paste0("in_sum_",k,"_",l,"_XX")]] *(df$d_sq_int_XX == l & df$F_g_XX >= 3)
                }
                df[[paste0("in_brackets_pl_",l,"_",j,"_XX")]] <- df[[paste0("in_brackets_pl_",l,"_",j,"_XX")]] - get(paste0("coefs_sq_",l,"_XX"))[j,1]
                df[[paste0("combined_pl",increase_XX,"_temp_",l,"_",i,"_XX")]] <-  df[[paste0("combined_pl",increase_XX,"_temp_",l,"_",i,"_XX")]] + df[[paste0("M_pl",increase_XX,"_",l,"_",j,"_",i,"_XX")]] * df[[paste0("in_brackets_pl_",l,"_",j,"_XX")]]
              }
              df[[paste0("part2_pl_switch",increase_XX,"_",i,"_XX")]] <- as.numeric(df[[paste0("part2_pl_switch",increase_XX,"_",i,"_XX")]] + df[[paste0("combined_pl",increase_XX,"_temp_",l,"_",i,"_XX")]] * (df$d_sq_int_XX == l))
            }
            if (increase_XX == 1) {
              df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]] <- df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]] - 
                  df[[paste0("part2_pl_switch1_",i,"_XX")]]
            } else {
              df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]] <- df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]] + 
                  df[[paste0("part2_pl_switch0_",i,"_XX")]]
            }
          }

          df[paste0("U_Gg_pl_",i,"_temp_var_XX")] <- 
              as.numeric(df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]])

          df <- df %>% group_by(.data$group_XX) %>% 
              mutate(!!paste0("U_Gg_pl_",i,"_var_XX"):= sum(.data[[paste0("U_Gg_pl_", i,"_temp_var_XX")]], na.rm = TRUE))
        }

        ###### 7. Computing adjustements for the normalized and trends_lin options for placebos (similar to part 4, not commented) 
        if (normalized == TRUE) {
          if (is.null(continuous)) {
            df$sum_temp_pl_XX <- ifelse(df$time_XX >= df$F_g_XX & df$time_XX <= df$F_g_XX - 1 + i & !is.na(df[[paste0("diff_y_pl_",i,"_XX")]]) & df$S_g_XX == increase_XX, df$treatment_XX - df$d_sq_XX, NA)
          } else {
            df$sum_temp_pl_XX <- ifelse(df$time_XX >= df$F_g_XX & df$time_XX <= df$F_g_XX - 1 + i & !is.na(df[[paste0("diff_y_pl_",i,"_XX")]]) & df$S_g_XX == increase_XX, df$treatment_XX_orig - df$d_sq_XX_orig, NA)
          }
          df <- df %>% group_by(.data$group_XX) %>% 
              mutate(!!paste0("sum_treat_until_",i,"_pl_XX") := sum(.data$sum_temp_pl_XX, na.rm = TRUE)) %>%
              dplyr::select(-.data$sum_temp_pl_XX)
          df[[paste0("delta_D_pl_",i,"_cum_temp_XX")]] <- ifelse(
            df[[paste0("dist_to_switch_pl_",i,"_XX")]] == 1,
            (df$N_gt_XX/get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) * (
              df$S_g_XX * df[[paste0("sum_treat_until_",i,"_pl_XX")]] +
              (1 - df$S_g_XX) * (-df[[paste0("sum_treat_until_",i,"_pl_XX")]]) 
            ), NA)
          assign(paste0("delta_norm_pl_",i,"_XX"), sum(df[[paste0("delta_D_pl_",i,"_cum_temp_XX")]], na.rm = TRUE))    
        }
      }
    }

    Ntrendslin_pl <- 1
    for (i in 1:l_placebo_u_a_XX) {
      Ntrendslin_pl <- min(Ntrendslin_pl, get(paste0("N",increase_XX,"_placebo_",i,"_XX")), na.rm = TRUE)
    }

    if (isTRUE(trends_lin) & Ntrendslin_pl != 0) {
      df[[paste0("U_Gg_pl_",l_u_a_XX,"_TL")]] <- 0
      df[[paste0("U_Gg_pl_",l_u_a_XX,"_var_TL")]] <- 0
      for (i in 1:l_u_a_XX) {
        df[[paste0("U_Gg_pl_",l_u_a_XX,"_TL")]] <- ifelse(!is.na(df[[paste0("U_Gg_placebo_",i,"_XX")]]),
        df[[paste0("U_Gg_pl_",l_u_a_XX,"_TL")]] + df[[paste0("U_Gg_placebo_",i,"_XX")]],
        df[[paste0("U_Gg_pl_",l_u_a_XX,"_TL")]])
        
        df[[paste0("U_Gg_pl_",l_u_a_XX,"_var_TL")]] <- ifelse(!is.na(df[[paste0("U_Gg_pl_",i,"_var_XX")]]),
        df[[paste0("U_Gg_pl_",l_u_a_XX,"_var_TL")]] + df[[paste0("U_Gg_pl_",i,"_var_XX")]],
        df[[paste0("U_Gg_pl_",l_u_a_XX,"_var_TL")]])
      }

      df[[paste0("U_Gg_placebo_",l_u_a_XX,"_XX")]] <- df[[paste0("U_Gg_pl_",l_u_a_XX,"_TL")]] 
      df[[paste0("U_Gg_pl_",l_u_a_XX,"_var_XX")]] <- df[[paste0("U_Gg_pl_",l_u_a_XX,"_var_TL")]] 
      df[[paste0("U_Gg_pl_",l_u_a_XX,"_TL")]] <- NULL
      df[[paste0("U_Gg_pl_",l_u_a_XX,"_var_TL")]] <- NULL
    }

  }
  ###### 8. Computing Average Total Effect estimator

  ## Average total effect not estimated if trends_lin option requested:
  if (isFALSE(trends_lin)) {
    ## Computing the sum of the N1_`i'_XX for the weights w. 
    assign(paste0("sum_N",increase_XX,"_l_XX"), 0)
    for (i in 1:l_u_a_XX) {
      assign(paste0("sum_N",increase_XX,"_l_XX"), get(paste0("sum_N",increase_XX,"_l_XX"))+
      get(paste0("N",increase_XX,"_",i,"_XX")))
    }

    ## Dropping/Initializing needed variables
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

        #### Computing the weights w_(+,\ell)
        assign(paste0("w_",i,"_XX"), get(paste0("N",increase_XX,"_",i,"_XX")) / get(paste0("sum_N",increase_XX,"_l_XX")))

        #### Computing the delta^D_{+,\ell}s, which enter in the denominator of \hat{\delta}_+/\hat{\delta}_-.
        if (is.null(continuous)) {
          df[paste0("delta_D_",i,"_temp_XX")] <- df$N_gt_XX/get(paste0("N",increase_XX,"_",i,"_XX")) * ((df$treatment_XX - df$d_sq_XX) * df$S_g_XX + (1 - df$S_g_XX) * (df$d_sq_XX - df$treatment_XX))
        } else {
          ## Redefine this with original treatment if continuous option specified (treatment was binarized and staggerized)
          df[paste0("delta_D_",i,"_temp_XX")] <- df$N_gt_XX/get(paste0("N",increase_XX,"_",i,"_XX")) * ((df$treatment_XX_orig - df$d_sq_XX_orig) * df$S_g_XX + (1 - df$S_g_XX) * (df$d_sq_XX_orig - df$treatment_XX_orig))
        }

        df[[paste0("delta_D_",i,"_temp_XX")]] <- ifelse(df[[paste0("distance_to_switch_",i,"_XX")]] == 1, 
            df[[paste0("delta_D_",i,"_temp_XX")]], NA)
        df[[paste0("delta_D_",i,"_temp_XX")]][is.na(
          df[[paste0("delta_D_",i,"_temp_XX")]])] <- 0

        df[paste0("delta_D_",i,"_XX")] <- sum(df[[paste0("delta_D_",i,"_temp_XX")]], na.rm = TRUE)
        df <- df %>% dplyr::select(-.data[[paste0("delta_D_",i,"_temp_XX")]])

        ## Computing the numerator of U^+_{G,g}: summing up the U_{G,g,l}s, after weighting them
        df$U_Gg_num_XX <- df$U_Gg_num_XX + get(paste0("w_",i,"_XX")) * df[[paste0("U_Gg",i,"_XX")]]

        ## Computing the numerator of the "alternative" U_{G,g}s for the variance : summing up the U_{G,g,l}_vars, after weighting them
        df$U_Gg_num_var_XX <- df$U_Gg_num_var_XX + get(paste0("w_",i,"_XX")) * df[[paste0("U_Gg",i,"_var_XX")]]

        ## Computing the denominator of U^+_{G,g}: summing up the delta^D_{+,l}s, after weighting them
        df$U_Gg_den_XX <- df$U_Gg_den_XX + get(paste0("w_",i,"_XX")) * df[[paste0("delta_D_",i,"_XX")]]
      }
    }

    ## Computing the U^+_{G,g}s.
    df$U_Gg_XX <- df$U_Gg_num_XX/df$U_Gg_den_XX

    ## Computing the U^+_{G,g}_vars.
    df$U_Gg_var_XX <- df$U_Gg_num_var_XX/df$U_Gg_den_XX
  }

  ## Update passthrough constants
  for (e in names(const)) {
    const[[e]] <- get(e)
  }
  const[[paste0("sum_N",increase_XX,"_l_XX")]] <- get(paste0("sum_N",increase_XX,"_l_XX"))

  for (i in 1:l_u_a_XX) {
    if (isTRUE(normalized)) {
      const[[paste0("delta_norm_",i,"_XX")]] <- get(paste0("delta_norm_",i,"_XX"))
    }
  }
  if (placebo != 0) {
    if (l_placebo_u_a_XX >= 1) {
      for (i in 1:l_placebo_u_a_XX) {
          if (isTRUE(normalized)) {
            const[[paste0("delta_norm_pl_",i,"_XX")]] <- get(paste0("delta_norm_pl_",i,"_XX"))
          }
      }
    }
  }

  data <- list(
    df,
    const
  )
  data_names <- c("df", "const")

  })

  names(data) <- data_names
  return(data)
}