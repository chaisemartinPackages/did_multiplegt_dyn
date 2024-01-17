#' Internal function of did_multiplegt_dyn
#' @param df df
#' @param Y Y
#' @param G G
#' @param T T
#' @param D D
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
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang := 
#' @importFrom rlang .data
#' @noRd
did_multiplegt_dyn_core <- function(
    df, 
    Y, 
    G, 
    T, 
    D, 
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
    less_conservative_se
    ) {
  
  # Inherited Globals #
  L_u_XX <- globals$L_u_XX
  L_placebo_u_XX <- globals$L_placebo_u_XX
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
  df <- df %>% dplyr::select(-dplyr::any_of(c("num_g_paths_0_XX", "cohort_fullpath_0_XX")))

  # Estimating the DID_{+,l} or DID_{-,l}
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

      paste0("dof_cohort_",i,"_ns_t_XX"), paste0("dof_cohort_",i,"_ns_t_XX"),
      paste0("count_cohort_",i,"_ns_t_XX"), paste0("count_cohort_",i,"_ns_t_XX"),
      paste0("total_cohort_",i,"_ns_t_XX"), paste0("total_cohort_",i,"_ns_t_XX"),
      paste0("mean_cohort_",i,"_ns_t_XX"), paste0("mean_cohort_",i,"_ns_t_XX")
     ))) 

    df <- df[order(df$group_XX, df$time_XX), ]
    df <- df %>% 
      group_by(.data$group_XX) %>% 
      mutate(!!paste0("diff_y_", i, "_XX") := .data$outcome_XX - lag(.data$outcome_XX, i)) %>% 
      ungroup()

    if (isTRUE(less_conservative_se)) {
      df$d_fg0_XX <- df$d_sq_XX
      df <- df %>% group_by(.data$group_XX) %>% 
          mutate(!!paste0("d_fg",i,"_XX") := mean(
            .data$outcome_XX[.data$time_XX == .data$F_g_XX -1 + i], na.rm = TRUE))
      df[[paste0("d_fg",i,"_XX")]] <- ifelse(is.na(df[[paste0("d_fg",i,"_XX")]]),
          df[[paste0("d_fg",i-1,"_XX")]], df[[paste0("d_fg",i,"_XX")]])

      
    }

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

    # Identifying the control (g,t)s in the estimation of dynamic effect i 
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
    assign(paste0("N",increase_XX,"_",i,"_XX"), 0)
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

    ## Cohort never switchers
    # DOF
    df <- joint_trends(df, "d_sq_XX", trends_nonparam)
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

    # Denominator
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

    # Numerator
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

    # Estimator for the expectation (no need for bysort or any conditioning as the Numerator and denominator are generated along the same set of conditions)
    df[[paste0("mean_cohort_",i,"_ns_t_XX")]] <- df[[paste0("total_cohort_",i,"_ns_t_XX")]] /
       df[[paste0("count_cohort_",i,"_ns_t_XX")]]

    # if (isFALSE(less_conservative_se)) {
    df <- joint_trends(df, c("d_sq_XX", "F_g_XX", "d_fg_XX", paste0("distance_to_switch_",i,"_XX")), trends_nonparam)
    df <- df %>% group_by(.data$joint_trends_XX) %>% 
        mutate(!!paste0("dof_cohort_",i,"_s_t_XX") := 
        sum(.data[[paste0("dof_y_",i,"_N_gt_XX")]],na.rm = TRUE)) %>% ungroup()
    df[[paste0("dof_cohort_",i,"_s_t_XX")]] <- ifelse(
          df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
          df[[paste0("dof_cohort_",i,"_s_t_XX")]], NA)

    # Denominator
    df <- df %>% group_by(.data$joint_trends_XX) %>% 
        mutate(!!paste0("count_cohort_",i,"_s_t_XX") := 
        sum(.data[["N_gt_XX"]],na.rm = TRUE)) %>% ungroup()
    df[[paste0("count_cohort_",i,"_s_t_XX")]] <- ifelse(
          df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
          df[[paste0("count_cohort_",i,"_s_t_XX")]], NA)

    # Numerator
    df <- df %>% group_by(.data$joint_trends_XX) %>% 
        mutate(!!paste0("total_cohort_",i,"_s_t_XX") := 
        sum(.data[[paste0("diff_y_",i,"_N_gt_XX")]],na.rm = TRUE)) %>% ungroup()
    df[[paste0("total_cohort_",i,"_s_t_XX")]] <- ifelse(
          df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
          df[[paste0("total_cohort_",i,"_s_t_XX")]], NA)

    df[paste0("mean_cohort_",i,"_s_t_XX")] <- df[[paste0("total_cohort_",i,"_s_t_XX")]] /
         df[[paste0("count_cohort_",i,"_s_t_XX")]]
      #}

    if (get(paste0("N",increase_XX,"_",i,"_XX")) != 0) {
      df[paste0("dummy_U_Gg",i,"_XX")] <- as.numeric(i <= df$T_g_XX - 1)
    
      df[paste0("U_Gg",i,"_temp_XX")] <- df[[paste0("dummy_U_Gg",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_",i,"_XX"))) * as.numeric(df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df$N_gt_XX * (df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]])
      df[paste0("U_Gg",i,"_temp_XX")] <- df[[paste0("U_Gg",i,"_temp_XX")]] *  df[[paste0("diff_y_",i,"_XX")]]
      df <- df %>% group_by(.data$group_XX) %>% 
          mutate(!!paste0("U_Gg",i,"_XX") := sum(.data[[paste0("U_Gg",i,"_temp_XX")]], na.rm = TRUE))
      df[[paste0("U_Gg",i,"_XX")]] <- df[[paste0("U_Gg",i,"_XX")]] * df$first_obs_by_gp_XX

      # Counting the number of groups for which we can estimate U_Gg`i'_temp_XX - to help compute the "N" displayed by the command //

      df[paste0("count",i,"_core_XX")] <- ifelse(
        (!is.na(df[[paste0("U_Gg",i,"_temp_XX")]]) & df[[paste0("U_Gg",i,"_temp_XX")]] != 0) | 
        (df[[paste0("U_Gg",i,"_temp_XX")]] == 0 & df[[paste0("diff_y_",i,"_XX")]]==0 & (df[[paste0("distance_to_switch_",i,"_XX")]] != 0 | df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]] != 0 & df[[paste0("never_change_d_",i,"_XX")]] != 0)), df$N_gt_XX, 0)
      df[paste0("count",i,"_core_XX")] <- as.numeric(df[[paste0("count",i,"_core_XX")]])

      # Computing the "alternative" U_{G,g,l} which will be used for the computation of the variance only - these are like the above U_{G,g,l}s, except that the outcome differences are demeaned, and there is a DOF adjustment when possible//

      df[paste0("U_Gg",i,"_temp_var_XX")] <- 0

      # For controls
      df[[paste0("U_Gg",i,"_temp_var_XX")]] <- ifelse(
        df[[paste0("never_change_d_",i,"_XX")]] == 1 & df[[paste0("dof_cohort_",i,"_ns_t_XX")]] == 1, 
      df[[paste0("dummy_U_Gg",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_",i,"_XX"))) * (df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df[[paste0("diff_y_",i,"_N_gt_XX")]], 
      df[[paste0("U_Gg",i,"_temp_var_XX")]])

      df[[paste0("U_Gg",i,"_temp_var_XX")]] <- ifelse(
      df[[paste0("never_change_d_",i,"_XX")]] == 1 & df[[paste0("dof_cohort_",i,"_ns_t_XX")]] > 1 & 
        !is.na(df[[paste0("dof_cohort_",i,"_ns_t_XX")]]), 
      df[[paste0("dummy_U_Gg",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_",i,"_XX"))) *(df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df$N_gt_XX * (df[[paste0("diff_y_",i,"_XX")]] - (df[[paste0("mean_cohort_",i,"_ns_t_XX")]] * (sqrt(df[[paste0("dof_cohort_",i,"_ns_t_XX")]] /(df[[paste0("dof_cohort_",i,"_ns_t_XX")]]-1))) * df[[paste0("never_change_d_",i,"_XX")]])), 
      df[[paste0("U_Gg",i,"_temp_var_XX")]])

      # For switchers
      # if (isFALSE(less_conservative_se)) {
      df[[paste0("U_Gg",i,"_temp_var_XX")]] <- ifelse(
        df[[paste0("distance_to_switch_",i,"_XX")]] == 1 & df[[paste0("dof_cohort_",i,"_s_t_XX")]] == 1, 
      df[[paste0("dummy_U_Gg",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_",i,"_XX"))) *(df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df[[paste0("diff_y_",i,"_N_gt_XX")]], 
      df[[paste0("U_Gg",i,"_temp_var_XX")]])

      df[[paste0("U_Gg",i,"_temp_var_XX")]] <- ifelse(
        df[[paste0("distance_to_switch_",i,"_XX")]] == 1 & df[[paste0("dof_cohort_",i,"_s_t_XX")]] > 1 & 
            !is.na(df[[paste0("dof_cohort_",i,"_s_t_XX")]]), 
      df[[paste0("dummy_U_Gg",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_",i,"_XX"))) *(df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df$N_gt_XX * (df[[paste0("diff_y_",i,"_XX")]] - df[[paste0("mean_cohort_",i,"_s_t_XX")]] * sqrt(df[[paste0("dof_cohort_",i,"_s_t_XX")]] /(df[[paste0("dof_cohort_",i,"_s_t_XX")]]-1) * df[[paste0("distance_to_switch_",i,"_XX")]])), 
      df[[paste0("U_Gg",i,"_temp_var_XX")]])
      #}

      df[[paste0("U_Gg",i,"_temp_var_XX")]] <- as.numeric(df[[paste0("U_Gg",i,"_temp_var_XX")]])

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
      assign(paste0("delta_norm_",i,"_XX"), sum(df[[paste0("delta_D_",i,"_cum_temp_XX")]], na.rm = TRUE))    
    }
  }

  ### Trends_lin option #
  Ntrendslin <- 1
  for (i in 1:l_u_a_XX) {
    Ntrendslin <- min(Ntrendslin, get(paste0("N", increase_XX, "_",i,"_XX")), na.rm = TRUE)
  }

  if (isTRUE(trends_lin) & Ntrendslin != 0) {
    df[paste0("U_Gg",l_u_a_XX,"_TL")] <- 0
    df[paste0("U_Gg",l_u_a_XX,"_var_TL")] <- 0
    for (i in 1:l_u_a_XX) {
      df[[paste0("U_Gg",l_u_a_XX,"_TL")]] <- ifelse(!is.na(df[[paste0("U_Gg",i,"_XX")]]),
      df[[paste0("U_Gg",l_u_a_XX,"_TL")]] + df[[paste0("U_Gg",i,"_XX")]],
      df[[paste0("U_Gg",l_u_a_XX,"_TL")]])
      df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] <- ifelse(!is.na(df[[paste0("U_Gg",i,"_var_XX")]]),
      df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] + df[[paste0("U_Gg",i,"_var_XX")]],
      df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]])

    }

    df[[paste0("U_Gg",l_u_a_XX,"_XX")]] <- df[[paste0("U_Gg",l_u_a_XX,"_TL")]] 
    df[[paste0("U_Gg",l_u_a_XX,"_var_XX")]] <- df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] 
    df[[paste0("U_Gg",l_u_a_XX,"_TL")]] <- NULL
    df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] <- NULL
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
          paste0("dist_to_switch_pl_",i,"_wXX"),
          paste0("dof_cohort_pl_",i,"_ns_t_XX"), paste0("dof_cohort_pl_",i,"_ns_t_XX"),
          paste0("count_cohort_pl_",i,"_ns_t_XX"), paste0("count_cohort_pl_",i,"_ns_t_XX"),
          paste0("total_cohort_pl_",i,"_ns_t_XX"), paste0("total_cohort_pl_",i,"_ns_t_XX"),
          paste0("mean_cohort_pl_",i,"_ns_t_XX"), paste0("mean_cohort_pl_",i,"_ns_t_XX")
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
        assign(paste0("N",increase_XX,"_placebo_",i,"_XX"), 0)
        for (t in t_min_XX:T_max_XX) {
          assign(paste0("N",increase_XX,"_placebo_",i,"_XX"), 
            get(paste0("N",increase_XX,"_placebo_",i,"_XX")) + mean(df[[paste0("N", increase_XX,"_t_placebo_", i, "_XX")]][df$time_XX == t], na.rm = TRUE))
        }
        assign(paste0("N",increase_XX,"_placebo_",i,"_XX"), get(paste0("N",increase_XX,"_placebo_",i,"_XX")))

        # Computing N^0_{t,l,g} and N^1_{t,l,g}
        df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
        df <- df %>% group_by(.data$joint_trends_XX) %>%
            mutate(!!paste0("N",increase_XX,"_t_placebo_",i,"_g_XX") := sum(.data[[paste0("dist_to_switch_pl_",i,"_wXX")]], na.rm = TRUE)) %>% ungroup()

        # Computing the mean of differences of outcomes for non-treated and treated separately
        df[paste0("diff_y_pl_", i,"_N_gt_XX")] <- df[[paste0("diff_y_pl_", i,"_XX")]] * df$N_gt_XX
        df[paste0("dof_y_pl_", i,"_N_gt_XX")]  <- as.numeric(df$N_gt_XX != 0 & !is.na(df[[paste0("diff_y_pl_", i,"_XX")]]))

        ## Cohort never switchers
        # DOF
        df <- joint_trends(df, "d_sq_XX", trends_nonparam)
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

        # Denominator
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

        # Numerator
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

        # Estimator for the expectation (no need for bysort or any conditioning as the Numerator and denominator are generated along the same set of conditions)
        df[[paste0("mean_cohort_pl_",i,"_ns_t_XX")]] <- df[[paste0("total_cohort_pl_",i,"_ns_t_XX")]] /
          df[[paste0("count_cohort_pl_",i,"_ns_t_XX")]]

        # if (isFALSE(less_conservative_se)) {
        df <- joint_trends(df, c("d_sq_XX", "F_g_XX", "d_fg_XX", paste0("dist_to_switch_pl_",i,"_XX")), trends_nonparam)
        df <- df %>% group_by(.data$joint_trends_XX) %>% 
            mutate(!!paste0("dof_cohort_pl_",i,"_s_t_XX") := 
            sum(.data[[paste0("dof_y_pl_",i,"_N_gt_XX")]],na.rm = TRUE)) %>% ungroup()
        df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]] <- ifelse(
              df[[paste0("dist_to_switch_pl_",i,"_XX")]] == 1,
              df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]], NA)

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
          #}

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
            df[[paste0("never_change_d_pl_",i,"_XX")]] == 1 & df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]] == 1, 
          df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) * (df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) * (df$time_XX >= i + 2 & df$time_XX <= df$T_g_XX) * df[[paste0("diff_y_pl_",i,"_N_gt_XX")]], 
          df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]])

          df[paste0("U_Gg_pl_",i,"_temp_var_XX")] <- ifelse(
            df[[paste0("never_change_d_pl_",i,"_XX")]] == 1 & df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]] > 1 & !is.na(df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]]), 
          df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) *(df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) * (df$time_XX >= i + 2 & df$time_XX <= df$T_g_XX) * df$N_gt_XX * (df[[paste0("diff_y_pl_",i,"_XX")]] - (df[[paste0("mean_cohort_pl_",i,"_ns_t_XX")]] * sqrt(df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]] /(df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]]-1)) * df[[paste0("never_change_d_pl_",i,"_XX")]])), 
          df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]])

          # if (isFALSE(less_conservative_se)) {
          df[paste0("U_Gg_pl_",i,"_temp_var_XX")] <- ifelse(
            df[[paste0("dist_to_switch_pl_",i,"_XX")]] == 1 & df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]] == 1, 
          df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) *(df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) * (df$time_XX >= i + 2 & df$time_XX <= df$T_g_XX) * df[[paste0("diff_y_pl_",i,"_N_gt_XX")]], 
          df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]])

          df[paste0("U_Gg_pl_",i,"_temp_var_XX")] <- ifelse(
            df[[paste0("dist_to_switch_pl_",i,"_XX")]] == 1 & df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]] > 1 & !is.na(df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]]), 
          df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) *(df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) * (df$time_XX >= i + 2 & df$time_XX <= df$T_g_XX) * df$N_gt_XX * (df[[paste0("diff_y_pl_",i,"_XX")]] - (df[[paste0("mean_cohort_pl_",i,"_s_t_XX")]] * sqrt(df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]] /(df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]]-1)) * df[[paste0("dist_to_switch_pl_",i,"_XX")]])), 
          df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]])
          #}

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
          assign(paste0("delta_norm_pl_",i,"_XX"), sum(df[[paste0("delta_D_pl_",i,"_cum_temp_XX")]], na.rm = TRUE))    
        }
      }
    }

    Ntrendslin_pl <- 1
    for (i in 1:l_placebo_u_a_XX) {
      Ntrendslin_pl <- min(Ntrendslin_pl, get(paste0("N",increase_XX,"_placebo_",i,"_XX")), na.rm = TRUE)
    }

    if (isTRUE(trends_lin) & Ntrendslin_pl == 0) {
      assign(paste0("N",increase_XX,"_placebo_",l_placebo_u_a_XX), 0)
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
  # End of the Placebo computation 

  # For the estimation of \hat{\delta}

  if (isFALSE(trends_lin)) {
    assign(paste0("sum_N",increase_XX,"_l_XX"), 0)
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
  }

  # Update passthrough constants
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