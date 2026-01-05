#' Internal function of did_multiplegt_dyn that computes U_Gg_plus_XX, U_Gg_minus_XX, U_Gg_var_plus_XX, and U_Gg_var_minus_XX. These are essential variables for the computation of the DID_\ell estimators and their variances.
#' @param df df
#' @param outcome outcome
#' @param group group
#' @param time time
#' @param treatment treatment
#' @param effects effects
#' @param placebo placebo
#' @param switchers_core switchers_core
#' @param trends_nonparam trends_nonparam
#' @param controls controls
#' @param same_switchers same_switchers
#' @param same_switchers_pl same_switchers_pl
#' @param only_never_switchers only_never_switchers
#' @param normalized normalized
#' @param globals globals
#' @param const constants
#' @param trends_lin trends_lin
#' @param controls_globals controls_globals
#' @param less_conservative_se less_conservative_se
#' @param continuous continuous
#' @import data.table
#' @importFrom stats na.omit predict setNames
#' @importFrom MASS ginv
#' @importFrom fixest feols
#' @importFrom dplyr lag n_distinct
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
    cluster,
    switchers_core = NULL, 
    trends_nonparam, 
    controls, 
    same_switchers, 
    same_switchers_pl, 
    only_never_switchers,
    normalized,
    globals,
    const,
    trends_lin,
    controls_globals,
    less_conservative_se,
    continuous
    ) {

  # CRAN Compliance
  F_g_XX <- NULL
  N_gt_XX <- NULL
  T_d_XX <- NULL
  cum_fillin_XX <- NULL
  d_fg_XX_temp <- NULL
  d_sq_int_XX <- NULL
  dum_fillin_temp_XX <- NULL
  dum_fillin_temp_pl_XX <- NULL
  dummy_XX <- NULL
  fillin_g_XX<- NULL
  fillin_g_pl_XX <- NULL
  group_XX <- NULL
  num_g_paths_0_XX <- NULL
  outcome_XX <- NULL
  path_0_XX<- NULL
  relevant_y_missing_XX <- NULL
  sum_temp_XX <- NULL
  sum_temp_pl_XX <- NULL
  time_XX <- NULL
  
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
        l_placebo_u_a_XX <- min(placebo, L_placebo_u_XX)
    }
    increase_XX <- 1
  }

  if (switchers_core == "out") {
    l_u_a_XX <- min(L_a_XX, effects, na.rm = TRUE)
    if (placebo != 0) {
        l_placebo_u_a_XX <- min(placebo, L_placebo_a_XX)
    }
    increase_XX <- 0
  }

  ## Initializing values of baseline treatment
  levels_d_sq_XX <- levels(as.factor(df$d_sq_int_XX))
  df$num_g_paths_0_XX <- df$cohort_fullpath_0_XX <- NULL

  ## 2. Data preparation steps to generate variables necessary for computation of event-study effects	
  ## Loop over the number of dynamic effects

  for (i in 1:l_u_a_XX) {

    df[[paste0("distance_to_switch_",i,"_XX")]] <- NULL
    df[[paste0("never_change_d_",i,"_XX")]] <- NULL
    df[[paste0("N",increase_XX,"_t_",i,"_XX")]] <- NULL
    df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]] <- NULL
    df[[paste0("N_gt_control_",i,"_XX")]] <- NULL
    df[[paste0("diff_y_",i,"_XX")]] <- NULL
    df[[paste0("diff_y_",i,"_XX_temp")]] <- NULL
    df[[paste0("dummy_U_Gg",i,"_XX")]] <- NULL
    df[[paste0("U_Gg",i,"_temp_XX")]] <- NULL
    df[[paste0("U_Gg",i,"_XX")]] <- NULL
    df[[paste0("count",i,"_core_XX")]] <- NULL
    df[[paste0("mean_diff_y_",i,"_nd_sq_t_XX")]] <- NULL
    df[[paste0("mean_diff_y_",i,"_d_sq_t_XX")]] <- NULL
    df[[paste0("U_Gg",i,"_temp_var_XX")]] <- NULL
    df[[paste0("U_Gg",i,"_var_XX")]] <- NULL
    df[[paste0("U_Gg",i,"_var_2_XX")]] <- NULL
    df[[paste0("count_var_",i,"_ntreat_XX_temp")]] <- NULL
    df[[paste0("count_var_",i,"_ntreat_XX")]] <- NULL
    df[[paste0("count_var_",i,"_treat_XX_temp")]] <- NULL
    df[[paste0("count_var_",i,"_treat_XX")]] <- NULL
    df[[paste0("avg_diff_y_",i,"_tnp_XX")]] <- NULL
    df[[paste0("count_diff_y_",i,"_nd_sq_t_XX")]] <- NULL
    df[[paste0("count_diff_y_",i,"_d_sq_t_XX")]] <- NULL
    df[[paste0("never_change_d_",i,"_wXX")]] <- NULL
    df[[paste0("distance_to_switch_",i,"_wXX")]] <- NULL
    df[[paste0("dof_cohort_",i,"_ns_t_XX")]] <- NULL
    df[[paste0("dof_cohort_",i,"_s_t_XX")]] <- NULL
    df[[paste0("dof_cohort_",i,"_s0_t_XX")]] <- NULL
    df[[paste0("dof_cohort_",i,"_s1_t_XX")]] <- NULL
    df[[paste0("dof_cohort_",i,"_s2_t_XX")]] <- NULL
    df[[paste0("count_cohort_",i,"_ns_t_XX")]] <- NULL
    df[[paste0("count_cohort_",i,"_s_t_XX")]] <- NULL
    df[[paste0("count_cohort_",i,"_s0_t_XX")]] <- NULL
    df[[paste0("count_cohort_",i,"_s1_t_XX")]] <- NULL
    df[[paste0("count_cohort_",i,"_s2_t_XX")]] <- NULL
    df[[paste0("total_cohort_",i,"_ns_t_XX")]] <- NULL
    df[[paste0("total_cohort_",i,"_s_t_XX")]] <- NULL
    df[[paste0("total_cohort_",i,"_s0_t_XX")]] <- NULL
    df[[paste0("total_cohort_",i,"_s1_t_XX")]] <- NULL
    df[[paste0("total_cohort_",i,"_s2_t_XX")]] <- NULL
    df[[paste0("mean_cohort_",i,"_ns_t_XX")]] <- NULL
    df[[paste0("mean_cohort_",i,"_s_t_XX")]] <- NULL
    df[[paste0("mean_cohort_",i,"_s0_t_XX")]] <- NULL
    df[[paste0("mean_cohort_",i,"_s1_t_XX")]] <- NULL
    df[[paste0("mean_cohort_",i,"_s2_t_XX")]] <- NULL

    ## Creating long difference of outcome
    df <- as.data.table(df)
    df <- df[order(df$group_XX, df$time_XX), ]
    df[, lagout := data.table::shift(outcome_XX, i), by = group_XX]
    df[, paste0("diff_y_", i, "_XX") := outcome_XX - lagout]
    df[, lagout := NULL]

    
    ## Creating treatment paths if less_conservative_se option specified
    if (isTRUE(less_conservative_se)) {

      ## Creating a time-invariant, group-level variable, containing g's treatment at F_g-1+\ell
      df$d_fg_XX_temp <- ifelse(df$time_XX == df$F_g_XX +i-1,
          df$treatment_XX, NA)
      df[, paste0("d_fg",i,"_XX") := mean(d_fg_XX_temp, na.rm = TRUE), by = "group_XX"]

      ## This variable might be missing, for groups whose treatment never changes, and for groups not observed \ell periods after treatment change. Then, we impute their treatment at F_g-1+\ell-1. Inconsequential, just to avoid missing values. We also need to initialize a variable d_fg0_XX, when \ell=1.      
      if (i == 1) {
        df$d_fg0_XX <- df$d_sq_XX
        df[, path_0_XX := .GRP, by = c("d_fg0_XX", "F_g_XX")]
      }

      df[[paste0("d_fg",i,"_XX")]] <- ifelse(is.na(df[[paste0("d_fg",i,"_XX")]]),
          df[[paste0("d_fg",i-1,"_XX")]], df[[paste0("d_fg",i,"_XX")]])
      df[, paste0("path_",i,"_XX") := .GRP, by = c(paste0("path_",i-1,"_XX"), paste0("d_fg",i,"_XX"))]
      
      df$d_fg_XX_temp <- NULL

      ## For each group, define a variable counting how many groups belong to the same cohort, with cohorts defined as d_fg0_XX F_g_XX, as well as the full path.
      if (i == 1) {
        df[, num_g_paths_0_XX := n_distinct(group_XX), by = path_0_XX]
        df$cohort_fullpath_0_XX <- as.numeric(df$num_g_paths_0_XX > 1)
      }
      ## For each group, generate a dummy for whether at least two groups in their cohort.
      df[, paste0("num_g_paths_",i,"_XX") := n_distinct(group_XX), by = c(paste0("path_",i,"_XX"))]
      df[[paste0("cohort_fullpath_",i,"_XX")]] <- as.numeric(df[[paste0("num_g_paths_",i,"_XX")]] > 1)
    }
    ## Identifying the control (g,t)s in the estimation of dynamic effect i 
    df[[paste0("never_change_d_", i, "_XX")]] <- as.numeric(df$F_g_XX > df$time_XX)
    df[[paste0("never_change_d_", i, "_XX")]] <- ifelse(is.na(df[[paste0("diff_y_", i, "_XX")]]), NA,  df[[paste0("never_change_d_", i, "_XX")]]) 

    if (isTRUE(only_never_switchers)) {
     df[[paste0("never_change_d_", i, "_XX")]] <- ifelse(df$F_g_XX > df$time_XX & df$F_g_XX < T_max_XX + 1 & !is.null(df[[paste0("diff_y_",i,"_XX")]]), 0, df[[paste0("never_change_d_", i, "_XX")]])
    }

    ## Creating N^g_t:
    ## number of control groups for g at t
    df[[paste0("never_change_d_", i, "_wXX")]] <- df[[paste0("never_change_d_", i, "_XX")]] * df$N_gt_XX
    df[, paste0("N_gt_control_", i, "_XX") := sum(get(paste0("never_change_d_", i, "_wXX")), na.rm = TRUE), by = c("time_XX", "d_sq_XX", trends_nonparam)]

    ###### Creating binary variable indicating whether g is \ell periods away from switch
    ## If the same_switchers option is specified:
    if (same_switchers == TRUE) {
      df <- df[order(df$group_XX, df$time_XX), ]

      df[, N_g_control_check_XX := 0]

      for(q in 1:effects){

        setorder(df,"group_XX","time_XX")
        df[, diff_y_last_XX := outcome_XX - shift(outcome_XX, n = q, type = "lag"), by = group_XX]
        df[, never_change_d_last_XX := ifelse(!is.na(diff_y_last_XX) & F_g_XX > time_XX, 1, NA_real_)]

        if (isTRUE(only_never_switchers)) {
          df[F_g_XX > time_XX & F_g_XX < T_max_XX + 1 & !is.na(diff_y_last_XX), never_change_d_last_XX := 0]
          }

        df[, N_gt_control_last_XX := sum(never_change_d_last_XX * N_gt_XX, na.rm = TRUE), by = c("time_XX", "d_sq_XX", trends_nonparam)]

        df[, N_g_control_last_m_XX := mean(ifelse(time_XX == F_g_XX - 1 + q, N_gt_control_last_XX, NA_real_), na.rm = TRUE), by = group_XX]

        df[, diff_y_relev_XX := mean(ifelse(time_XX == F_g_XX - 1 + q, diff_y_last_XX, NA_real_), na.rm = TRUE), by = group_XX]

        df[, N_g_control_check_XX := N_g_control_check_XX + as.numeric(N_g_control_last_m_XX > 0 & !is.na(diff_y_relev_XX))]
      }

      ## If the same_switchers_pl option is specified:
      if (same_switchers_pl == TRUE) {

        df[, N_g_control_check_pl_XX := 0]

        for (q in 1:placebo) {
          
          df[, diff_y_last_XX := outcome_XX - shift(outcome_XX, n = q, type = "lead"), by = group_XX]
          
          df[, never_change_d_last_XX := ifelse(!is.na(diff_y_last_XX) & F_g_XX > time_XX, 1, NA_real_)]

          if (isTRUE(only_never_switchers)) {
            df[F_g_XX > time_XX & F_g_XX < T_max_XX + 1 & !is.na(diff_y_last_XX), never_change_d_last_XX := 0]
          }
    
          df[, N_gt_control_last_XX := sum(never_change_d_last_XX * N_gt_XX, na.rm = TRUE), by = c("time_XX", "d_sq_XX", trends_nonparam)]
                    
          df[, N_g_control_last_m_XX := mean(ifelse(time_XX == F_g_XX - 1 - q, N_gt_control_last_XX, NA_real_), na.rm = TRUE), by = group_XX]
                    
          df[, diff_y_relev_XX := mean(ifelse(time_XX == F_g_XX - 1 - q, diff_y_last_XX, NA_real_), na.rm = TRUE), by = group_XX]
          
          df[, N_g_control_check_pl_XX := N_g_control_check_pl_XX + as.numeric(N_g_control_last_m_XX > 0 & !is.na(diff_y_relev_XX))]
          }

        ## Generate a variable tagging the switchers that should be dropped
        ## Is the case if at least one of the placebos or effects we try to estimate is missing:
        df$relevant_y_missing_XX <- (is.na(df$outcome_XX) & df$time_XX >= df$F_g_XX - 1 - placebo & df$time_XX <= df$F_g_XX - 1 + effects)
        ## Or if some of the controls are missing:
        if (!is.null(controls)) {
          df$relevant_y_missing_XX <- ifelse(df$fd_X_all_non_missing_XX == 0 & df$time_XX >= df$F_g_XX - 1 - placebo & df$time_XX <= df$F_g_XX - 1 + effects, 1, df$relevant_y_missing_XX)
        }
        
        # df[, cum_fillin_XX := sum(relevant_y_missing_XX, na.rm = TRUE), by = group_XX]
        # df$dum_fillin_temp_XX <- df$cum_fillin_XX == 0 & df$time_XX == df$F_g_XX - 1 + effects
        # df[, fillin_g_XX := sum(dum_fillin_temp_XX, na.rm = TRUE), by = group_XX]
        # #v1.0.1
        # if (!is.null(placebo)) {
        #   df$dum_fillin_temp_pl_XX <- 
        #       df$cum_fillin_XX == 0 & df$time_XX == df$F_g_XX - 1 - placebo
        #   df[, fillin_g_pl_XX := sum(dum_fillin_temp_pl_XX, na.rm = TRUE), by = group_XX]
        # }
        df[, fillin_g_pl_XX := N_g_control_check_pl_XX == placebo]

        ## tag switchers who have no missings from F_g_XX-1-placebo to F_g_XX-1+effects
        df[[paste0("still_switcher_",i,"_XX")]] <- 
            df$F_g_XX - 1 + effects <= df$T_g_XX & df$N_g_control_check_XX == effects

        df[[paste0("distance_to_switch_", i, "_XX")]] <- 
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

        # df[, cum_fillin_XX := sum(relevant_y_missing_XX, na.rm = TRUE), by = group_XX]
        # df$dum_fillin_temp_XX <- df$cum_fillin_XX == 0 & df$time_XX == df$F_g_XX - 1 + effects
        # df[, fillin_g_XX := sum(dum_fillin_temp_XX, na.rm = TRUE), by = group_XX]

        ## tag switchers who have no missings from F_g_XX-1 to F_g_XX-1+effects
        df[[paste0("still_switcher_",i,"_XX")]] <- df$F_g_XX - 1 + effects <= df$T_g_XX & df$N_g_control_check_XX == effects
        df[[paste0("distance_to_switch_", i, "_XX")]] <- 
        ifelse(!is.na(df[[paste0("diff_y_", i, "_XX")]]),
        df[[paste0("still_switcher_", i, "_XX")]] == 1 &
        df$time_XX == df$F_g_XX- 1 + i & i <= df$L_g_XX & df$S_g_XX == increase_XX &
        df[[paste0("N_gt_control_", i, "_XX")]] > 0 & !is.na(df[[paste0("N_gt_control_", i, "_XX")]]), NA)
      }
    }
    else {
      ## If the same_switchers option is not specified
      i <- i  # assume i is already defined
      diff_var <- paste0("diff_y_",         i, "_XX")
      Ngt_var  <- paste0("N_gt_control_",   i, "_XX")
      dst_var  <- paste0("distance_to_switch_", i, "_XX")

      # start as missing (.) for everyone
      df[, (dst_var) := NA_integer_]

      # emulate: "if diff_y_i_XX != ."
      df[!is.na(get(diff_var)),
        (dst_var) := as.integer(
          time_XX == (F_g_XX - 1 + i) &      # time_XX == F_g_XX - 1 + i
          i <= L_g_XX &                      # `i' <= L_g_XX
          S_g_XX == increase_XX &            # S_g_XX == increase_XX
          get(Ngt_var) > 0 &                 # N_gt_control_i_XX > 0
          !is.na(get(Ngt_var))               # N_gt_control_i_XX != .
        )
      ]
    }

    #### Creating a variable counting the number of groups \ell periods away from switch at t
    df[[paste0("distance_to_switch_", i, "_wXX")]] <- 
    df[[paste0("distance_to_switch_", i, "_XX")]] * df$N_gt_XX
    df[, paste0("N", increase_XX,"_t_", i, "_XX") := sum(get(paste0("distance_to_switch_", i, "_wXX")), na.rm = TRUE), by = time_XX]
    df[, paste0("N_dw", increase_XX,"_t_", i, "_XX") := sum(get(paste0("distance_to_switch_", i, "_XX")), na.rm = TRUE), by = time_XX]

    #### Computing N^1_\ell/N^0_\ell.
    assign(paste0("N",increase_XX,"_",i,"_XX"), 0)
    assign(paste0("N",increase_XX,"_dw_",i,"_XX"), 0)
    for (t in t_min_XX:T_max_XX) {
      assign(paste0("N",increase_XX,"_",i,"_XX"), 
        get(paste0("N",increase_XX,"_",i,"_XX")) + mean(df[[paste0("N", increase_XX,"_t_", i, "_XX")]][df$time_XX == t], na.rm = TRUE))
      assign(paste0("N",increase_XX,"_dw_",i,"_XX"), 
        get(paste0("N",increase_XX,"_dw_",i,"_XX")) + mean(df[[paste0("N_dw", increase_XX,"_t_", i, "_XX")]][df$time_XX == t], na.rm = TRUE))
    }

    #### Creating N^1_{t,\ell,g}/N^0_{t,\ell,g}: Variable counting number of groups \ell periods away from switch at t, and with same D_{g,1} and trends_nonparam.
    df[, paste0("N",increase_XX,"_t_",i,"_g_XX") := sum(get(paste0("distance_to_switch_",i,"_wXX")), na.rm = TRUE), by = c("time_XX", "d_sq_XX", trends_nonparam)]

    #### Creating all the adjustment terms to compute estimators with controls, and their variances
    if (!is.null(controls)) {

      #### Initialize intermediate Variable needed later	
      df[[paste0("part2_switch",increase_XX,"_",i,"_XX")]] <- 0

      ## generation of the T_d variable = max_{g:D_g,1 = d} F_g - 1: 
      ## last period when treatment effects can still be estimated for groups with baseline treatment equal to d
      df[, T_d_XX := max(F_g_XX, na.rm = TRUE), by = d_sq_int_XX]
      df$T_d_XX <- df$T_d_XX - 1

      ## Computing the long differences of the control variables (X_g_t - X_g_t-l)
      count_controls <- 0
      for (var in controls) {
        count_controls <- count_controls + 1
        df[[paste0("diff_X",count_controls,"_", i, "_XX")]]  <- df[[var]] - lag(df[[var]], i)

        ## Computing N_g_t * (X_g_t - X_g_t-l)
        df[[paste0("diff_X",count_controls,"_",i,"_N_XX")]] <- df$N_gt_XX * 
            df[[paste0("diff_X",count_controls,"_",i,"_XX")]]

        ## index l corresponds to d in the paper
        for (l in levels_d_sq_XX) { 

          # ## intermediate variable to count the number of groups within each not yet switched cohort          
          # df$dummy_XX <- NULL
          # df$dummy_XX <- as.numeric(df$F_g_XX > df$time_XX & df$d_sq_int_XX == l)

          ## Computing coordinates of vectors m^+_{g,d,\ell} and m^-_{g,d,\ell}
          # small m
          ## Creating variable inside the summation across t in m^+_{g,d,\ell}/m^-_{g,d,\ell}
          df[[paste0("m",increase_XX,"_g_",l,"_",count_controls,"_",i,"_XX")]] <- 
            (i <= df$T_g_XX - 2 & df$d_sq_int_XX == l) * (G_XX / get(paste0("N",increase_XX,"_",i,"_XX"))) * ((df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * (df[[paste0("diff_X",count_controls,"_",i,"_N_XX")]]))

          ## Summing that variable across t, and leaving one non missing observation per g	
          df[, paste0("m",increase_XX,"_",l,"_", count_controls,"_",i,"_XX") := sum(get(paste0("m",increase_XX,"_g_",l,"_", count_controls,"_",i,"_XX")), na.rm = TRUE), by = "group_XX"]
          df[[paste0("m",increase_XX,"_",l,"_",count_controls,"_",i,"_XX")]] <- ifelse(
            df$first_obs_by_gp_XX == 1, df[[paste0("m",increase_XX,"_",l,"_", count_controls,"_",i,"_XX")]], NA)
          
          ## Computing coordinates of vectors M^+_{d,\ell} and M^-_{d,\ell}
          df[[paste0("M",increase_XX,"_",l,"_", count_controls,"_",i,"_XX")]] <- 
              sum(df[[paste0("m",increase_XX,"_",l,"_", count_controls,"_",i,"_XX")]], na.rm = TRUE) / G_XX          

          ## number of groups within each not yet switched cohort
          df[, dummy_XX := ifelse(F_g_XX > time_XX & d_sq_int_XX == l & !is.na(diff_y_XX) ,1,0)]
          df[, paste0("E_hat_denom_", count_controls,"_", l, "_XX") := sum(dummy_XX[d_sq_int_XX == l], na.rm = TRUE), by = c("time_XX", "d_sq_int_XX")]
          df[[paste0("E_hat_denom_", count_controls,"_", l, "_XX")]] <- ifelse(
            df$d_sq_int_XX == l, df[[paste0("E_hat_denom_", count_controls,"_", l, "_XX")]], NA)

          ## Add the indicator for at least two groups in the cohort to E_y_hat_gt_`l'_XX (demeaning is possible)
          df[[paste0("E_y_hat_gt_",l,"_XX")]] <- df[[paste0("E_y_hat_gt_int_",l,"_XX")]] *
              (df[[paste0("E_hat_denom_",count_controls,"_",l,"_XX")]] >= 2)

          ## Computing the summation from t=2 to F_g-1 that appears in the last term 
          ## of U^{+,var,X}_{g,\ell} and U^{-,var,X}_{g,\ell}, defined in the companion paper.          
          df[[paste0("N_c_",l,"_temp_XX")]] <- df[[paste0("N_c_",l,"_XX")]] <- df[[paste0("in_sum_temp_",count_controls,"_",l,"_XX")]]

          df[[paste0("N_c_",l,"_temp_XX")]] <- df$N_gt_XX* (df$d_sq_int_XX == l & df$time_XX >= 2 &
              df$time_XX <= df$T_d_XX & df$time_XX < df$F_g_XX & !is.na(df$diff_y_XX))
          df[[paste0("N_c_",l,"_XX")]] <- sum(df[[paste0("N_c_",l,"_temp_XX")]], na.rm = TRUE)

          df[[paste0("in_sum_temp_",count_controls,"_",l,"_XX")]] <- 
            (df[[paste0("prod_X",count_controls,"_Ngt_XX")]] * (1 + 
              (df[[paste0("E_hat_denom_", count_controls,"_",l,"_XX")]] >= 2) *
              (sqrt(df[[paste0("E_hat_denom_", count_controls,"_",l,"_XX")]] / 
              (df[[paste0("E_hat_denom_", count_controls,"_",l,"_XX")]] - 1)) -1)) *
              (df$diff_y_XX - df[[paste0("E_y_hat_gt_",l,"_XX")]]) *
              (df$time_XX >= 2 & df$time_XX <= df$F_g_XX - 1)) / df[[paste0("N_c_",l,"_XX")]]

          df[, paste0("in_sum_",count_controls,"_",l,"_XX") := sum(get(paste0("in_sum_temp_",count_controls,"_",l,"_XX")), na.rm=TRUE), by = "group_XX"]

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


  


    if (is.null(trends_nonparam)) trends_nonparam <- character(0)
    if (is.null(cluster)) cluster <- ""

    # ----- dynamic column names -----
    diff_y_i_col     <- sprintf("diff_y_%s_XX", i)
    diff_y_N_col     <- sprintf("diff_y_%s_N_gt_XX", i)
    dof_ns_col       <- sprintf("dof_ns_%s_XX", i)
    dof_s_col        <- sprintf("dof_s_%s_XX",  i)
    never_col        <- sprintf("never_change_d_%s_XX", i)
    dist_col         <- sprintf("distance_to_switch_%s_XX", i)
    Ninc_col         <- sprintf("N%s_t_%s_XX", increase_XX, i)  # N{increase_XX}_t_{i}_XX

    count_col        <- sprintf("count_cohort_%s_ns_t_XX", i)
    total_col        <- sprintf("total_cohort_%s_ns_t_XX", i)
    mean_col         <- sprintf("mean_cohort_%s_ns_t_XX",  i)
    dof_cohort_col   <- sprintf("dof_cohort_%s_ns_t_XX",   i)

    cluster_dof_col  <- sprintf("cluster_dof_%s_ns_XX",    i)

    by_cols <- c("d_sq_XX", trends_nonparam, "time_XX")

    # ----- capture drop -----
    to_drop <- intersect(c(diff_y_N_col, dof_ns_col, dof_s_col,
                            count_col, total_col, mean_col, dof_cohort_col, cluster_dof_col),
                          names(df))
    if (length(to_drop)) df[, (to_drop) := NULL]

    # ----- gen diff_y_i_N_gt_XX = N_gt_XX * diff_y_i_XX
    df[, (diff_y_N_col) := N_gt_XX * get(diff_y_i_col)]

    # ----- gen dof_ns_i_XX -----
    df[, (dof_ns_col) := as.integer(
      !is.na(N_gt_XX) & N_gt_XX != 0 &
      !is.na(get(diff_y_i_col)) &
      !is.na(get(never_col)) & get(never_col) == 1 &
      !is.na(get(Ninc_col))  & get(Ninc_col)  > 0
    )]

    # ----- gen dof_s_i_XX -----
    df[, (dof_s_col) := as.integer(
      !is.na(N_gt_XX) & N_gt_XX != 0 &
      !is.na(get(dist_col)) & get(dist_col) == 1
    )]

    # ----- grouped totals (Stata: bys ... : gegen ...) only where dof_ns_i_XX == 1 -----
    # Mean's denominator: total(N_gt_XX)
    df[get(dof_ns_col) == 1,
        (count_col) := sum(N_gt_XX, na.rm = TRUE),
        by = by_cols]

    # Mean's numerator: total(diff_y_i_N_gt_XX)
    df[get(dof_ns_col) == 1,
        (total_col) := sum(get(diff_y_N_col), na.rm = TRUE),
        by = by_cols]

    # Mean = numerator / denominator
    df[, (mean_col) := get(total_col) / get(count_col)]

    # ----- DOF counting -----
    if (identical(cluster, "") || is.na(cluster)) {
      # No cluster: sum of the indicator within group (== number of eligible rows)
      df[get(dof_ns_col) == 1,
          (dof_cohort_col) := sum(get(dof_ns_col), na.rm = TRUE),
          by = by_cols]
    } else {
      # With cluster: unique count of clusters among eligible rows
      if (cluster_dof_col %in% names(df)) df[, (cluster_dof_col) := NULL]
      df[, (cluster_dof_col) := ifelse(get(dof_ns_col) == 1, get(cluster), NA)]
      df[!is.na(get(cluster_dof_col)),
          (dof_cohort_col) := uniqueN(get(cluster_dof_col)),
          by = by_cols]
    }
    df[[paste0("diff_y_", i,"_N_gt_XX")]] <- df[[paste0("diff_y_", i,"_XX")]] * df$N_gt_XX
    df[[paste0("dof_y_", i,"_N_gt_XX")]]  <- as.numeric(df$N_gt_XX != 0 & !is.na(df[[paste0("diff_y_", i,"_XX")]]))
  
    ## For switchers, if option less_conservative_se not specified, demeaning wrt to cohorts defined by D_{g,1}, F_g, D_{g,F_g}, `trends_nonparam' (\mathcal{C}_k in companion paper).
    if (isFALSE(less_conservative_se)) {

      if (is.null(trends_nonparam)) trends_nonparam <- character(0)
      if (is.null(cluster)) cluster <- ""

      # ----- dynamic names -----
      dof_s_col      <- sprintf("dof_s_%s_XX", i)
      diff_y_N_col   <- sprintf("diff_y_%s_N_gt_XX", i)

      count_col      <- sprintf("count_cohort_%s_s_t_XX",  i)
      total_col      <- sprintf("total_cohort_%s_s_t_XX",  i)
      mean_col       <- sprintf("mean_cohort_%s_s_t_XX",   i)
      dof_cohort_col <- sprintf("dof_cohort_%s_s_t_XX",    i)
      cluster_dofcol <- sprintf("cluster_dof_%s_s_XX",     i)

      # by: d_sq_XX F_g_XX d_fg_XX `trends_nonparam'
      by_cols <- c("d_sq_XX", "F_g_XX", "d_fg_XX", trends_nonparam)

      # ----- clean any prior artifacts (capture drop) -----

      # ----- Mean's denominator: total(N_gt_XX) where dof_s_i_XX == 1 -----
      df[get(dof_s_col) == 1,
        (count_col) := sum(N_gt_XX, na.rm = TRUE),
        by = by_cols]

      # ----- Mean's numerator: total(diff_y_i_N_gt_XX) where dof_s_i_XX == 1 -----
      df[get(dof_s_col) == 1,
        (total_col) := sum(get(diff_y_N_col), na.rm = TRUE),
        by = by_cols]

      # ----- Mean -----
      df[, (mean_col) := get(total_col) / get(count_col)]

      # ----- DoF counting (Diego 16-06-25: adjust with clusters if provided) -----
      if (identical(cluster, "") || is.na(cluster)) {
        # No cluster: sum of indicator within group (i.e., number of eligible rows)
        df[get(dof_s_col) == 1,
          (dof_cohort_col) := sum(get(dof_s_col), na.rm = TRUE),
          by = by_cols]
      } else {
        # With cluster: unique number of clusters among eligible rows
        if (cluster_dofcol %in% names(df)) df[, (cluster_dofcol) := NULL]
        df[, (cluster_dofcol) := ifelse(get(dof_s_col) == 1, get(cluster), NA)]
        df[!is.na(get(cluster_dofcol)),
          (dof_cohort_col) := uniqueN(get(cluster_dofcol)),
          by = by_cols]
      }

    } else {
      ## For switchers, if option less_conservative_se specified, demeaning wrt to cohorts defined by D_{g,1} F_g, D_{g,F_g},..., D_{g,F_g+\ell}, if that cohort has at least two groups, if not: demeaning wrt to cohorts defined by D_{g,1} F_g, D_{g,F_g}, if that cohort has at least two groups, if not: demeaning wrt D_{g,1} F_g.

      # by D_{g,1}, F_g, `trends_nonparam':
      ### Denominator of the mean
      df[, paste0("count_cohort_",i,"_s0_t_XX") := sum(N_gt_XX[get(paste0("distance_to_switch_",i,"_XX")) == 1], na.rm = TRUE), by = c("path_0_XX", trends_nonparam)]

      df[[paste0("count_cohort_",i,"_s0_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
            df[[paste0("count_cohort_",i,"_s0_t_XX")]], NA)

      ### Numerator of the mean
      df[, paste0("total_cohort_",i,"_s0_t_XX") := sum(get(paste0("diff_y_",i,"_N_gt_XX"))[get(paste0("distance_to_switch_",i,"_XX")) == 1], na.rm = TRUE), by = c("path_0_XX", trends_nonparam)]

      df[[paste0("total_cohort_",i,"_s0_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
            df[[paste0("total_cohort_",i,"_s0_t_XX")]], NA)

      ### DOF
      df[, paste0("dof_cohort_",i,"_s0_t_XX") := sum(get(paste0("dof_y_",i,"_N_gt_XX"))[get(paste0("distance_to_switch_",i,"_XX")) == 1], na.rm = TRUE), by = c("path_0_XX", trends_nonparam)]

      df[[paste0("dof_cohort_",i,"_s0_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
        df[[paste0("dof_cohort_",i,"_s0_t_XX")]], NA)

      # by D_{g,1}, F_g, D_{g,F_g}, `trends_nonparam':

      ### Denominator of the mean
      df[, paste0("count_cohort_",i,"_s1_t_XX") := sum(N_gt_XX[get(paste0("distance_to_switch_",i,"_XX")) == 1], na.rm = TRUE), by = c("path_1_XX", trends_nonparam)]

      df[[paste0("count_cohort_",i,"_s1_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
            df[[paste0("count_cohort_",i,"_s1_t_XX")]], NA)

      ### Numerator of the mean
      df[, paste0("total_cohort_",i,"_s1_t_XX") := sum(get(paste0("diff_y_",i,"_N_gt_XX"))[get(paste0("distance_to_switch_",i,"_XX")) == 1], na.rm = TRUE), by = c("path_1_XX", trends_nonparam)]

      df[[paste0("total_cohort_",i,"_s1_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
            df[[paste0("total_cohort_",i,"_s1_t_XX")]], NA)

      ### DOF
      df[, paste0("dof_cohort_",i,"_s1_t_XX") := sum(get(paste0("dof_y_",i,"_N_gt_XX"))[get(paste0("distance_to_switch_",i,"_XX")) == 1], na.rm = TRUE), by = c("path_1_XX", trends_nonparam)]

      df[[paste0("dof_cohort_",i,"_s1_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
        df[[paste0("dof_cohort_",i,"_s1_t_XX")]], NA)

      # by D_{g,1}, F_g, D_{g,F_g},..., D_{g,F_g+\ell}, `trends_nonparam':

      ### Denominator of the mean
      df[, paste0("count_cohort_",i,"_s2_t_XX") := sum(N_gt_XX[get(paste0("distance_to_switch_",i,"_XX")) == 1], na.rm = TRUE), by = c(paste0("path_",i,"_XX"), trends_nonparam)]

      df[[paste0("count_cohort_",i,"_s2_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
            df[[paste0("count_cohort_",i,"_s2_t_XX")]], NA)

      ### Numerator of the mean
      df[, paste0("total_cohort_",i,"_s2_t_XX") := sum(get(paste0("diff_y_",i,"_N_gt_XX"))[get(paste0("distance_to_switch_",i,"_XX")) == 1], na.rm = TRUE), by = c(paste0("path_",i,"_XX"), trends_nonparam)]

      df[[paste0("total_cohort_",i,"_s2_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
            df[[paste0("total_cohort_",i,"_s2_t_XX")]], NA)

      ### DOF
      df[, paste0("dof_cohort_",i,"_s2_t_XX") := sum(get(paste0("dof_y_",i,"_N_gt_XX"))[get(paste0("distance_to_switch_",i,"_XX")) == 1], na.rm = TRUE), by = c(paste0("path_",i,"_XX"), trends_nonparam)]

      df[[paste0("dof_cohort_",i,"_s2_t_XX")]] <- ifelse(
            df[[paste0("distance_to_switch_",i,"_XX")]] == 1,
        df[[paste0("dof_cohort_",i,"_s2_t_XX")]], NA)

      ## Mean
      df[[paste0("mean_cohort_",i,"_s_t_XX")]] <- ifelse( df[[paste0("cohort_fullpath_",i,"_XX")]] == 1, df[[paste0("total_cohort_",i,"_s2_t_XX")]] / df[[paste0("count_cohort_",i,"_s2_t_XX")]], NA)
      df[[paste0("mean_cohort_",i,"_s_t_XX")]] <- ifelse(df[[paste0("cohort_fullpath_",i,"_XX")]] == 0 & df$cohort_fullpath_1_XX == 1, df[[paste0("total_cohort_",i,"_s1_t_XX")]] /df[[paste0("count_cohort_",i,"_s1_t_XX")]], df[[paste0("mean_cohort_",i,"_s_t_XX")]])
      df[[paste0("mean_cohort_",i,"_s_t_XX")]] <- ifelse(df$cohort_fullpath_1_XX == 0, df[[paste0("total_cohort_",i,"_s0_t_XX")]] /df[[paste0("count_cohort_",i,"_s0_t_XX")]], df[[paste0("mean_cohort_",i,"_s_t_XX")]])

      ## Counting number of groups for DOF adjustment
      df[[paste0("dof_cohort_",i,"_s_t_XX")]] <- ifelse( df[[paste0("cohort_fullpath_",i,"_XX")]] == 1, df[[paste0("dof_cohort_",i,"_s2_t_XX")]], NA)
      df[[paste0("dof_cohort_",i,"_s_t_XX")]] <- ifelse(df[[paste0("cohort_fullpath_",i,"_XX")]] == 0 & df$cohort_fullpath_1_XX == 1, df[[paste0("dof_cohort_",i,"_s1_t_XX")]], df[[paste0("dof_cohort_",i,"_s_t_XX")]])
      df[[paste0("dof_cohort_",i,"_s_t_XX")]] <- ifelse(df$cohort_fullpath_1_XX == 0, df[[paste0("dof_cohort_",i,"_s0_t_XX")]], df[[paste0("dof_cohort_",i,"_s_t_XX")]])
    }
    
    
    ## Modif Clément 30/6/2025:
    ## If a switcher is the only one in their cohort or if a not-yet-switcher is 
    ## the only one in their cohort, we demean wrt union of switchers and not-yet switchers, 
    ## provided switchers and not-yet-switchers do not all come from the same cluster.
    if (is.null(trends_nonparam)) trends_nonparam <- character(0)
    if (is.null(cluster)) cluster <- ""
    
    # ----- dynamic column names -----
    dof_ns_s_col   <- sprintf("dof_ns_s_%s_XX", i)
    count_col      <- sprintf("count_cohort_%s_ns_s_t_XX", i)
    total_col      <- sprintf("total_cohort_%s_ns_s_t_XX", i)
    mean_col       <- sprintf("mean_cohort_%s_ns_s_t_XX", i)
    dof_cohort_col <- sprintf("dof_cohort_%s_ns_s_t_XX", i)
    cluster_dofcol <- sprintf("cluster_dof_%s_ns_s_XX", i)
    
    dof_s_col      <- sprintf("dof_s_%s_XX", i)
    dof_ns_col     <- sprintf("dof_ns_%s_XX", i)
    diff_y_col     <- sprintf("diff_y_%s_N_gt_XX", i)
    
    by_cols <- c("d_sq_XX", trends_nonparam, "time_XX")
    
    # ----- drop previous artifacts (safe) -----
    drop_these <- intersect(c(dof_ns_s_col, count_col, total_col, mean_col, dof_cohort_col, cluster_dofcol),
                            names(df))
    if (length(drop_these)) df[, (drop_these) := NULL]
    
    # ----- dof_ns_s_i_XX = (dof_s_i_XX==1 | dof_ns_i_XX==1)
    df[, (dof_ns_s_col) := as.integer( (get(dof_s_col) == 1) | (get(dof_ns_col) == 1) ) ]
    
    # ----- Mean's denominator: total(N_gt_XX) within groups, only where dof_ns_s==1
    df[get(dof_ns_s_col) == 1,
       (count_col) := sum(N_gt_XX, na.rm = TRUE),
       by = by_cols]
    
    # ----- Mean's numerator: total(diff_y_i_N_gt_XX) within groups, only where dof_ns_s==1
    df[get(dof_ns_s_col) == 1,
       (total_col) := sum(get(diff_y_col), na.rm = TRUE),
       by = by_cols]
    
    # ----- Mean
    df[, (mean_col) := get(total_col) / get(count_col)]
    
    # ----- Counting number of groups for DOF adjustment
    if (identical(cluster, "") || is.na(cluster)) {
      # No cluster: total of the indicator within group (i.e., number of rows with dof_ns_s==1)
      df[get(dof_ns_s_col) == 1,
         (dof_cohort_col) := sum(get(dof_ns_s_col), na.rm = TRUE),
         by = by_cols]
    } else {
      # With cluster: unique number of clusters among rows with dof_ns_s==1
      # Create helper column restricted to those rows
      if (cluster_dofcol %in% names(df)) df[, (cluster_dofcol) := NULL]
      df[, (cluster_dofcol) := ifelse(get(dof_ns_s_col) == 1, get(cluster), NA)]
      
      df[!is.na(get(cluster_dofcol)),
         (dof_cohort_col) := uniqueN(get(cluster_dofcol)),
         by = by_cols]
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    #### From those parts, generate variables for the demeaning and the DOF adjustment 
    ## E_hat_(g,t), defined from parts depending on the cohort definition 
    # df[[paste0("E_hat_gt_",i,"_XX")]] <- ifelse(df$F_g_XX > df$time_XX, df[[paste0("mean_cohort_",i,"_ns_t_XX")]] * (df[[paste0("dof_cohort_",i,"_ns_t_XX")]] >= 2), NA)
    # df[[paste0("E_hat_gt_",i,"_XX")]] <- ifelse(df$time_XX == df$F_g_XX -1 + i,  df[[paste0("mean_cohort_",i,"_s_t_XX")]] * (df[[paste0("dof_cohort_",i,"_s_t_XX")]] >= 2), df[[paste0("E_hat_gt_",i,"_XX")]])
    
    compute_DOF_gt_with_nans <- function(DT, i, type_sect = "effect") {
      stopifnot(is.data.table(DT))

      # ----- dynamic column names -----
      if (identical(type_sect, "effect")) {
        DOF_col      <- sprintf("DOF_gt_%s_XX", i)
        dof_s_t_col  <- sprintf("dof_cohort_%s_s_t_XX", i)
        dof_ns_t_col <- sprintf("dof_cohort_%s_ns_t_XX", i)
        dof_ns_s_t_col <- sprintf("dof_cohort_%s_ns_s_t_XX", i)
      } else {
        DOF_col      <- sprintf("DOF_gt_pl_%s_XX", i)
        dof_s_t_col  <- sprintf("dof_cohort_pl_%s_s_t_XX", i)
        dof_ns_t_col <- sprintf("dof_cohort_pl_%s_ns_t_XX", i)
        dof_ns_s_t_col <- sprintf("dof_cohort_pl_%s_ns_s_t_XX", i)
      }
      
      # start like Stata: create var, default NA
      if (DOF_col %in% names(DT)) DT[, (DOF_col) := NULL]
      DT[, (DOF_col) := NA_real_]

      # gen DOF_gt_i_XX = 1 if (time_XX < F_g_XX | F_g_XX - 1 + i == time_XX)
      DT[(time_XX < F_g_XX) | (F_g_XX - 1 + as.integer(i) == time_XX),
        (DOF_col) := 1]

      # replace with sqrt(F_g_XX / (dof_cohort_i_s_t_XX - 1))
      # We need to review this since we are changing the column fro, F_g_XX to dof_cohort_`i'_s_t_XX
    #   DT[(F_g_XX - 1 + as.integer(i) == time_XX) & (get(dof_s_t_col) > 1),
    #     (DOF_col) := sqrt(F_g_XX / (get(dof_s_t_col) - 1))]
      DT[(F_g_XX - 1 + as.integer(i) == time_XX) & (get(dof_s_t_col) > 1),
        (DOF_col) := sqrt(get(dof_s_t_col) / (get(dof_s_t_col) - 1))]

      # replace with sqrt(dof_cohort_i_ns_t_XX / (dof_cohort_i_ns_t_XX - 1))
      DT[(time_XX < F_g_XX) & (get(dof_ns_t_col) > 1),
        (DOF_col) := sqrt(get(dof_ns_t_col) / (get(dof_ns_t_col) - 1))]

      # replace with sqrt(dof_cohort_i_ns_s_t_XX / (dof_cohort_i_ns_s_t_XX - 1))
      DT[(get(dof_ns_s_t_col) >= 2) &
          ( ((F_g_XX - 1 + as.integer(i) == time_XX) & (get(dof_s_t_col) == 1)) |
              ((time_XX < F_g_XX) & (get(dof_ns_t_col) == 1)) ),
        (DOF_col) := sqrt(get(dof_ns_s_t_col) / (get(dof_ns_s_t_col) - 1))]
      DT[ is.na(get(dof_s_t_col)) & is.na(get(dof_ns_t_col)) & is.na(get(dof_ns_s_t_col)),
            (DOF_col) := NA_real_]

      return(DT)
    }
    
    compute_E_hat_gt_with_nans <- function(df, i, type_sect = "effect") {
      stopifnot(is.data.table(df))
      
      if (type_sect == "effect") {
        E_hat    <- sprintf("E_hat_gt_%s_XX", i)
        mean_ns  <- sprintf("mean_cohort_%s_ns_t_XX", i)
        mean_s   <- sprintf("mean_cohort_%s_s_t_XX", i)
        mean_nss <- sprintf("mean_cohort_%s_ns_s_t_XX", i)
        
        dof_ns   <- sprintf("dof_cohort_%s_ns_t_XX", i)
        dof_s    <- sprintf("dof_cohort_%s_s_t_XX", i)
        dof_nss  <- sprintf("dof_cohort_%s_ns_s_t_XX", i)
      } else {
        E_hat    <- sprintf("E_hat_gt_pl_%s_XX", i)
        mean_ns  <- sprintf("mean_cohort_pl_%s_ns_t_XX", i)
        mean_s   <- sprintf("mean_cohort_pl_%s_s_t_XX", i)
        mean_nss <- sprintf("mean_cohort_pl_%s_ns_s_t_XX", i)
        
        dof_ns   <- sprintf("dof_cohort_pl_%s_ns_t_XX", i)
        dof_s    <- sprintf("dof_cohort_pl_%s_s_t_XX", i)
        dof_nss  <- sprintf("dof_cohort_pl_%s_ns_s_t_XX", i)
      }
      
      time_col <- "time_XX"
      Fg_col   <- "F_g_XX"
      
      ## Start missing
      df[, (E_hat) := NA_real_]
      
      ## --- 1) E_hat = 0 if cond_A: (time < Fg) | ((Fg - 1 + i) == time)
      df[
        get(time_col) < get(Fg_col) |
          (get(Fg_col) - 1L + i == get(time_col)),
        (E_hat) := 0
      ]
      
      ## --- 2) replace with mean_ns if cond_B:
      ## cond_B = (time < Fg) & (ns9999 >= 2)
      ## ns9999 = ns with NA treated as 9999 ⇒ (is.na(ns) | ns >= 2)
      df[
        get(time_col) < get(Fg_col) &
          (is.na(get(dof_ns)) | get(dof_ns) >= 2),
        (E_hat) := get(mean_ns)
      ]
      
      ## E_hat = NA if mean_ns is NA / NaN
      df[
        is.na(get(mean_ns)) | is.nan(get(mean_ns)),
        (E_hat) := NA_real_
      ]
      
      ## --- 3) replace with mean_s if cond_C:
      ## cond_C = ((Fg - 1 + i) == time) & (s9999 >= 2)
      ## s9999 = s with NA treated as 9999 ⇒ (is.na(s) | s >= 2)
      df[
        (get(Fg_col) - 1L + i == get(time_col)) &
          (is.na(get(dof_s)) | get(dof_s) >= 2),
        (E_hat) := get(mean_s)
      ]
      
      ## --- 4) replace with mean_nss if cond_D:
      ## cond_D = (nss >= 2) &
      ##          ( ((Fg - 1 + i) == time & s9999 == 1) |
      ##            (time < Fg & ns9999 == 1) )
      ## s9999==1 ⇔ !is.na(s) & s==1 ;  ns9999==1 ⇔ !is.na(ns) & ns==1
      df[
        (!is.na(get(dof_nss)) & get(dof_nss) >= 2) &
          (
            ((get(Fg_col) - 1L + i == get(time_col)) &
               !is.na(get(dof_s)) & get(dof_s) == 1) |
              (get(time_col) < get(Fg_col) &
                 !is.na(get(dof_ns)) & get(dof_ns) == 1)
          ),
        (E_hat) := get(mean_nss)
      ]
      
      ## E_hat = NA if mean_nss is NA / NaN
      df[
        is.na(get(mean_nss)) | is.nan(get(mean_nss)),
        (E_hat) := NA_real_
      ]
      
      return(df)
    }
    

    df = compute_E_hat_gt_with_nans(df, i)
    df = compute_DOF_gt_with_nans(df, i)
    
    
    
    
    
    

    ###### 3. Computing U_Gg_\ell variables
    #### If the dynamic effect can be estimated (as there are switchers), we compute the U_Gg_\ell variables etc.

    if (get(paste0("N",increase_XX,"_",i,"_XX")) != 0) {

      ## Creating a dummy variable indicating whether l<=T_g_XX-1
      df[[paste0("dummy_U_Gg",i,"_XX")]] <- as.numeric(i <= df$T_g_XX - 1)
    
      ## Computing U^+_{G,g,l}
      df[[paste0("U_Gg",i,"_temp_XX")]] <- df[[paste0("dummy_U_Gg",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_",i,"_XX"))) * as.numeric(df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df$N_gt_XX * (df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]])
      df[[paste0("U_Gg",i,"_temp_XX")]] <- df[[paste0("U_Gg",i,"_temp_XX")]] *  df[[paste0("diff_y_",i,"_XX")]]

      df[, paste0("U_Gg",i,"_XX") := sum(get(paste0("U_Gg",i,"_temp_XX")), na.rm = TRUE), by = group_XX]
      df[[paste0("U_Gg",i,"_XX")]] <- df[[paste0("U_Gg",i,"_XX")]] * df$first_obs_by_gp_XX

      # Counting the number of groups for which we can estimate U_Gg`i'_temp_XX - to help compute the "N" displayed by the command 
      df[, paste0("count", i, "_core_XX") := 0]
      
      df[
        (get(paste0("U_Gg", i, "_temp_XX")) != 0 & !is.na(get(paste0("U_Gg", i, "_temp_XX")))) |
          (
            get(paste0("U_Gg", i, "_temp_XX")) == 0 &
              get(paste0("diff_y_", i, "_XX")) == 0 &
              (
                get(paste0("distance_to_switch_", i, "_XX")) != 0 |
                  (
                    get(paste0("N", increase_XX, "_t_", i, "_g_XX")) != 0 &
                      get(paste0("never_change_d_", i, "_XX")) != 0
                  )
              )
          ),
        paste0("count", i, "_core_XX") := get("N_gt_XX")
      ]
      
      
      df[[paste0("count",i,"_core_XX")]] <- as.numeric(df[[paste0("count",i,"_core_XX")]])

      ## Computing U^(+,var)_{G,g,l}
      df[[paste0("U_Gg",i,"_temp_var_XX")]] <- 0

      df[[paste0("U_Gg",i,"_temp_var_XX")]] <- df[[paste0("dummy_U_Gg",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_",i,"_XX"))) * (df[[paste0("distance_to_switch_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_",i,"_g_XX")]]/df[[paste0("N_gt_control_",i,"_XX")]]) * df[[paste0("never_change_d_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df$N_gt_XX * df[[paste0("DOF_gt_",i,"_XX")]] * (df[[paste0("diff_y_",i,"_XX")]] - df[[paste0("E_hat_gt_",i,"_XX")]])

      ## Adding the additional part of U^(+,var,X)_{G,g,l}/U^(-,var,X)_{G,g,l} when controls are included: sum across values of baseline treatment d of M^+_(d,l)* a term in brackets in companion paper. 
      if (!is.null(controls)) {
        ## Loop over values of d_sq_int_XX:sum across values of baseline treatment  
        for (l in levels_d_sq_XX) {
          df[[paste0("combined",increase_XX,"_temp_",l,"_",i,"_XX")]] <- 0
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
          df[[paste0("part2_switch",increase_XX,"_",i,"_XX")]] <- as.numeric(df[[paste0("part2_switch",increase_XX,"_",i,"_XX")]] + df[[paste0("combined",increase_XX,"_temp_",l,"_",i,"_XX")]] )
        }
      }

      # Summing the U_{G,g,l} over time periods for each group
      df[[paste0("U_Gg",i,"_temp_var_XX")]] <- as.numeric(df[[paste0("U_Gg",i,"_temp_var_XX")]])
      df[, paste0("U_Gg",i,"_var_XX") := sum(get(paste0("U_Gg", i,"_temp_var_XX")), na.rm = TRUE), by = group_XX]

      if (!is.null(controls)) {

        #### Making the adjustement to U^(+,var)_{G,g,l} when controls are included
        if (increase_XX == 1) {
          df[[paste0("U_Gg",i,"_var_XX")]] <- df[[paste0("U_Gg",i,"_var_XX")]] - df[[paste0("part2_switch1_",i,"_XX")]]
        } else {
          df[[paste0("U_Gg",i,"_var_XX")]] <- df[[paste0("U_Gg",i,"_var_XX")]] - df[[paste0("part2_switch0_",i,"_XX")]]
        }
      }

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
      df[, paste0("sum_treat_until_",i,"_XX") := sum(sum_temp_XX, na.rm = TRUE), by = group_XX]
      df$sum_temp_XX <- NULL
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
#   if (isTRUE(trends_lin) & Ntrendslin != 0) {
# 	  ## Initializing at 0
#     df[[paste0("U_Gg",l_u_a_XX,"_TL")]] <- 0
#     df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] <- 0
# 	  ## summing up the U_Gg's up to the l-th (current) effect
#     for (i in 1:l_u_a_XX) {
#       df[[paste0("U_Gg",l_u_a_XX,"_TL")]] <- ifelse(!is.na(df[[paste0("U_Gg",i,"_XX")]]),
#       df[[paste0("U_Gg",l_u_a_XX,"_TL")]] + df[[paste0("U_Gg",i,"_XX")]],
#       df[[paste0("U_Gg",l_u_a_XX,"_TL")]])
#       df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] <- ifelse(!is.na(df[[paste0("U_Gg",i,"_var_XX")]]),
#       df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] + df[[paste0("U_Gg",i,"_var_XX")]],
#       df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]])
#     }
# 
# 	  ## replacing the U_Gg's with the one adjusted for group specific linear trends
#     df[[paste0("U_Gg",l_u_a_XX,"_XX")]] <- df[[paste0("U_Gg",l_u_a_XX,"_TL")]] 
#     df[[paste0("U_Gg",l_u_a_XX,"_var_XX")]] <- df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] 
#     df[[paste0("U_Gg",l_u_a_XX,"_TL")]] <- NULL
#     df[[paste0("U_Gg",l_u_a_XX,"_var_TL")]] <- NULL
#   }
  
  if (isTRUE(trends_lin) & Ntrendslin != 0) {
    
    lu <- as.integer(l_u_a_XX)
    
    col_TL     <- sprintf("U_Gg%d_TL",      lu)
    col_var_TL <- sprintf("U_Gg%d_var_TL",  lu)
    col_XX     <- sprintf("U_Gg%d_XX",      lu)
    col_var_XX <- sprintf("U_Gg%d_var_XX",  lu)
    
    ## Drop TL columns if they exist
    cols_to_drop <- intersect(c(col_TL, col_var_TL), names(df))
    if (length(cols_to_drop) > 0L) {
      df[, (cols_to_drop) := NULL]
    }
    
    ## Initialize TL columns to 0
    df[, c(col_TL, col_var_TL) := .(0.0, 0.0)]
    
    ## (Optional debug) list columns starting with "U_Gg"
    cols_U_Gg <- grep("^U_Gg", names(df), value = TRUE)
    # print(cols_U_Gg)
    
    ## --- Version 1: loop (closest to your polars code) ---
    for (i in seq_len(lu)) {
      df[, (col_TL)     := get(col_TL)     + get(sprintf("U_Gg%d_XX",     i))]
      df[, (col_var_TL) := get(col_var_TL) + get(sprintf("U_Gg%d_var_XX", i))]
    }
    
    ## Copy totals into the final columns
    df[, (col_XX)     := get(col_TL)]
    df[, (col_var_XX) := get(col_var_TL)]
  }
  

  ###### 5. Data preparation steps to generate variables necessary for computation of placebo effects	

  if (placebo != 0) {
    if (l_placebo_u_a_XX >= 1) {
      for (i in 1:l_placebo_u_a_XX) {

        df[[paste0("diff_y_pl_",i,"_XX")]] <- NULL
        df[[paste0("U_Gg_pl_",i,"_temp_XX")]] <- NULL
        df[[paste0("U_Gg_placebo_",i,"_XX")]] <- NULL
        df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]] <- NULL
        df[[paste0("U_Gg_pl_",i,"_var_XX")]] <- NULL
        df[[paste0("mean_diff_y_pl_",i,"_nd_sq_t_XX")]] <- NULL
        df[[paste0("mean_diff_y_pl_",i,"_d_sq_t_XX") ]] <- NULL
        df[[paste0("count_diff_y_pl_",i,"_nd_sq_t_XX")]] <- NULL
        df[[paste0("count_diff_y_pl_",i,"_d_sq_t_XX")]] <- NULL
        df[[paste0("dist_to_switch_pl_",i,"_XX")]] <- NULL
        df[[paste0("never_change_d_pl_",i,"_XX")]] <- NULL
        df[[paste0("N", increase_XX,"_t_placebo_",i,"_XX")]] <- NULL
        df[[paste0("N", increase_XX,"_t_placebo_",i,"_g_XX")]] <- NULL
        df[[paste0("N_gt_control_placebo_",i,"_XX")]] <- NULL
        df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] <- NULL
        df[[paste0("never_change_d_pl_",i,"_wXX")]] <- NULL
        df[[paste0("dist_to_switch_pl_",i,"_wXX")]] <- NULL
        df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]] <- NULL
        df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]] <- NULL
        df[[paste0("count_cohort_pl_",i,"_ns_t_XX")]] <- NULL
        df[[paste0("count_cohort_pl_",i,"_ns_t_XX")]] <- NULL
        df[[paste0("total_cohort_pl_",i,"_ns_t_XX")]] <- NULL
        df[[paste0("total_cohort_pl_",i,"_ns_t_XX")]] <- NULL
        df[[paste0("mean_cohort_pl_",i,"_ns_t_XX")]] <- NULL
        df[[paste0("mean_cohort_pl_",i,"_ns_t_XX")]] <- NULL
        df[[paste0("dof_cohort_pl_",i,"_ns_t_XX")]] <- NULL
        df[[paste0("dof_cohort_pl_",i,"_s_t_XX")]] <- NULL
        df[[paste0("dof_cohort_pl_",i,"_s0_t_XX")]] <- NULL
        df[[paste0("dof_cohort_pl_",i,"_s1_t_XX")]] <- NULL
        df[[paste0("dof_cohort_pl_",i,"_s2_t_XX")]] <- NULL
        df[[paste0("count_cohort_pl_",i,"_ns_t_XX")]] <- NULL
        df[[paste0("count_cohort_pl_",i,"_s_t_XX")]] <- NULL
        df[[paste0("count_cohort_pl_",i,"_s0_t_XX")]] <- NULL
        df[[paste0("count_cohort_pl_",i,"_s1_t_XX")]] <- NULL
        df[[paste0("count_cohort_pl_",i,"_s2_t_XX")]] <- NULL
        df[[paste0("total_cohort_pl_",i,"_ns_t_XX")]] <- NULL
        df[[paste0("total_cohort_pl_",i,"_s_t_XX")]] <- NULL
        df[[paste0("total_cohort_pl_",i,"_s0_t_XX")]] <- NULL
        df[[paste0("total_cohort_pl_",i,"_s1_t_XX")]] <- NULL
        df[[paste0("total_cohort_pl_",i,"_s2_t_XX")]] <- NULL
        df[[paste0("mean_cohort_pl_",i,"_ns_t_XX")]] <- NULL
        df[[paste0("mean_cohort_pl_",i,"_s_t_XX")]] <- NULL
        df[[paste0("mean_cohort_pl_",i,"_s0_t_XX")]] <- NULL
        df[[paste0("mean_cohort_pl_",i,"_s1_t_XX")]] <- NULL
        df[[paste0("mean_cohort_pl_",i,"_s2_t_XX")]] <- NULL

        # The steps to compute the placebos are:
        # 1. to place the corresponding outcome (y_{F_g-1} - y_{F_g - l - 1})) values in the same row of that (y_{F_g + l -1} - y_{F_g - 1}) of the symmetric DID_l. 
        # 2. The other variables, such as N_gt, N0_l or N1_l, remain unchanged, except that we have to check if diff_y_placebo ( = y_{F_g - 2l -2}- y_{F_g - l -1}) exists. 

        ## Computing the long differences for the placebos
        df[, paste0("diff_y_pl_", i, "_XX") := data.table::shift(outcome_XX, 2*i) - data.table::shift(outcome_XX, i), by = group_XX]
        
        ## Identifying the controls (g,t)s in the estimation of placebo i
        df[[paste0("never_change_d_pl_", i, "_XX")]] <- df[[paste0("never_change_d_", i, "_XX")]] * (!is.na(df[[paste0("diff_y_pl_", i,"_XX")]]))
        df[[paste0("never_change_d_pl_", i, "_wXX")]] <- df[[paste0("never_change_d_pl_", i, "_XX")]] * df$N_gt_XX

        ## number of control groups for g at t
        df[, paste0("N_gt_control_placebo_", i, "_XX") := sum(get(paste0("never_change_d_pl_", i, "_wXX")), na.rm = TRUE), by = c("time_XX", "d_sq_XX", trends_nonparam)]

        #### Creating a variable counting the number of groups \ell periods away from switch at t
        dist_var     <- paste0("distance_to_switch_", i, "_XX")
        diff_pl_var  <- paste0("diff_y_pl_", i, "_XX")
        Ngt_pl_var   <- paste0("N_gt_control_placebo_", i, "_XX")
        new_var      <- paste0("dist_to_switch_pl_", i, "_XX")

        # start as NA (Stata's .)
        df[, (new_var) := NA_real_]

        # emulate:
        # gen dist_to_switch_pl_i_XX =
        #   distance_to_switch_i_XX * (diff_y_pl_i_XX!=.) *
        #   (N_gt_control_placebo_i_XX>0 & N_gt_control_placebo_i_XX!=.)

        df[!is.na(get(dist_var)),
          (new_var) :=
            get(dist_var) *
            as.integer(!is.na(get(diff_pl_var))) *
            as.integer(get(Ngt_pl_var) > 0 & !is.na(get(Ngt_pl_var)))
        ]
        
        # v1.0.1.
        if (isTRUE(same_switchers_pl)) {
          df[[paste0("dist_to_switch_pl_", i, "_XX")]] <- df[[paste0("dist_to_switch_pl_", i, "_XX")]] * df$fillin_g_pl_XX
        }

        df[[paste0("dist_to_switch_pl_", i, "_wXX")]] <- 
        df[[paste0("dist_to_switch_pl_", i, "_XX")]] * df$N_gt_XX

        ## Creating a variable counting the number of groups \ell periods away from switch at t
        df[, paste0("N", increase_XX,"_t_placebo_", i, "_XX") := sum(get(paste0("dist_to_switch_pl_", i, "_wXX")), na.rm = TRUE), by = time_XX]
        df[, paste0("N", increase_XX,"_t_placebo_", i, "_dwXX") := sum(get(paste0("dist_to_switch_pl_", i, "_XX")), na.rm = TRUE), by = time_XX]

        #### Computing N^1_\ell/N^0_\ell. for the placebos
        ## Initializing the N1_`i'_XX/N0_`i'_XX scalar at 0. 
        assign(paste0("N",increase_XX,"_placebo_",i,"_XX"), 0)
        assign(paste0("N",increase_XX,"_dw_placebo_",i,"_XX"), 0)
        for (t in t_min_XX:T_max_XX) {
          assign(paste0("N",increase_XX,"_placebo_",i,"_XX"), 
            get(paste0("N",increase_XX,"_placebo_",i,"_XX")) + mean(df[[paste0("N", increase_XX,"_t_placebo_", i, "_XX")]][df$time_XX == t], na.rm = TRUE))
          assign(paste0("N",increase_XX,"_dw_placebo_",i,"_XX"), 
            get(paste0("N",increase_XX,"_dw_placebo_",i,"_XX")) + mean(df[[paste0("N", increase_XX,"_t_placebo_", i, "_dwXX")]][df$time_XX == t], na.rm = TRUE))
        }
        assign(paste0("N",increase_XX,"_placebo_",i,"_XX"), get(paste0("N",increase_XX,"_placebo_",i,"_XX")))
        assign(paste0("N",increase_XX,"_dw_placebo_",i,"_XX"), get(paste0("N",increase_XX,"_dw_placebo_",i,"_XX")))

        ## Creating N^1_{t,\ell,g}/N^0_{t,\ell,g} for the placebos: Variable counting number of groups \ell periods away from switch at t, and with same D_{g,1} and trends_nonparam.
        df[, paste0("N",increase_XX,"_t_placebo_",i,"_g_XX") := sum(get(paste0("dist_to_switch_pl_",i,"_wXX")), na.rm = TRUE), by =  c("time_XX", "d_sq_XX", trends_nonparam)]

        ## Creating all the adjustment terms to compute estimators with controls, and their variances 
        if (!is.null(controls)) {

          ## Initialize intermediate Variable needed later
          df[[paste0("part2_pl_switch",increase_XX,"_",i,"_XX")]] <- 0

          ## Computing the long differences of the control variables (X_g_t - X_g_t-l)
          count_controls <- 0
          for (var in controls) {
            count_controls <- count_controls + 1
            df[[paste0("diff_X",count_controls,"_placebo_",i,"_XX")]]  <- lag(df[[var]],2*i)-lag(df[[var]], i)
            ## Computing N_g_t * (X_g_t - X_g_t-l)
            df[[paste0("diff_X",count_controls,"_pl_",i,"_N_XX")]] <- df$N_gt_XX * 
                df[[paste0("diff_X",count_controls,"_placebo_",i,"_XX")]]

            ## index l corresponds to d in the paper
            for (l in levels_d_sq_XX) {
              
              df[[paste0("m",increase_XX,"_pl_g_",l,"_",count_controls,"_",i,"_XX")]] <- 
                (i <= df$T_g_XX - 2 & df$d_sq_int_XX == l) * (G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) * ((df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * (df[[paste0("diff_X",count_controls,"_pl_",i,"_N_XX")]]))

              df[, paste0("m_pl",increase_XX,"_",l,"_", count_controls,"_",i,"_XX") := sum(get(paste0("m",increase_XX,"_pl_g_",l,"_", count_controls,"_",i,"_XX")), na.rm = TRUE), by = group_XX]
              df[[paste0("m_pl",increase_XX,"_",l,"_",count_controls,"_",i,"_XX")]] <- ifelse(df$first_obs_by_gp_XX == 1, df[[paste0("m_pl",increase_XX,"_",l,"_", count_controls,"_",i,"_XX")]], NA)
              
              df[[paste0("M_pl",increase_XX,"_",l,"_", count_controls,"_",i,"_XX")]] <- 
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
        
        diff_y_pl_col     <- sprintf("diff_y_pl_%s_XX", i)
        diff_y_pl_N_col   <- sprintf("diff_y_pl_%s_N_gt_XX", i)
        dof_ns_pl         <- sprintf("dof_ns_pl_%s_XX", i)
        dof_s_pl          <- sprintf("dof_s_pl_%s_XX", i)
        never_change_col  <- sprintf("never_change_d_pl_%s_XX", i)
        N_t_placebo_col   <- sprintf("N%s_t_placebo_%s_XX", increase_XX, i)
        dist_to_switch_pl <- sprintf("dist_to_switch_pl_%s_XX", i)
        
        # diff_y_pl * N_gt_XX
        df[, (diff_y_pl_N_col) := get(diff_y_pl_col) * N_gt_XX]
        
        # dof_ns_pl: indicator (0/1)
        df[, (dof_ns_pl) := as.integer(
          N_gt_XX != 0 &
            !is.na(get(diff_y_pl_col)) &
            get(never_change_col) == 1 &
            get(N_t_placebo_col) > 0 &
            !is.na(get(N_t_placebo_col))
        )]
        
        # dof_s_pl: indicator (0/1), NA -> 0
        df[, (dof_s_pl) := as.integer(N_gt_XX != 0 & get(dist_to_switch_pl) == 1)]
        df[is.na(get(dof_s_pl)), (dof_s_pl) := 0L]
        
        ## --- 2. Groups for ns-cohort averages ---
        
        if (is.null(trends_nonparam)) {
          group_cols <- c("d_sq_XX", "time_XX")
        } else {
          group_cols <- c("d_sq_XX", "time_XX", trends_nonparam)
        }
        
        count_ns_col <- sprintf("count_cohort_pl_%s_ns_t_XX", i)
        total_ns_col <- sprintf("total_cohort_pl_%s_ns_t_XX", i)
        mean_ns_col  <- sprintf("mean_cohort_pl_%s_ns_t_XX", i)
        
        # count_cohort_pl_{i}_ns_t_XX: sum of N_gt_XX over group_cols, only for dof_ns_pl==1
        df[get(dof_ns_pl) == 1L,
           (count_ns_col) := sum(N_gt_XX),
           by = group_cols]
        
        # total_cohort_pl_{i}_ns_t_XX: sum of diff_y_pl_N over group_cols, only for dof_ns_pl==1
        df[get(dof_ns_pl) == 1L,
           (total_ns_col) := sum(get(diff_y_pl_N_col)),
           by = group_cols]
        
        # mean_cohort_pl_{i}_ns_t_XX
        df[, (mean_ns_col) := get(total_ns_col) / get(count_ns_col)]
        
        ## --- 3. DOF for ns-cohort (dof_cohort_pl_{i}_ns_t_XX) ---
        
        col_dof_ns     <- dof_ns_pl
        col_dof_coh_ns <- sprintf("dof_cohort_pl_%s_ns_t_XX", i)
        
        if (is.null(trends_nonparam)) {
          group_keys <- c("d_sq_XX", "time_XX")
        } else {
          group_keys <- c("d_sq_XX", trends_nonparam, "time_XX")
        }
        
        if (is.null(cluster) || cluster == "") {
          # sum of dof_ns_pl within group_keys (only where dof_ns_pl==1)
          df[get(col_dof_ns) == 1L,
             (col_dof_coh_ns) := as.numeric(sum(get(col_dof_ns))),
             by = group_keys]
        } else {
          # with cluster: unique clusters among rows with dof_ns_pl==1
          clust_dof <- sprintf("cluster_dof_pl_%s_ns_XX", i)
          df[, (clust_dof) := ifelse(get(col_dof_ns) == 1L, get(cluster), NA)]
          
          agg_ns <- df[!is.na(get(clust_dof)),
                       .(tmp_dof = uniqueN(get(clust_dof))),
                       by = group_keys]
          
          df[agg_ns, (col_dof_coh_ns) := as.numeric(i.tmp_dof), on = group_keys]
          df[is.na(get(clust_dof)), (col_dof_coh_ns) := NA_real_]
        }
        
        ## --- 4. Switchers cohort (C_k) demeaning ---
        
        if (is.null(trends_nonparam)) {
          group_cols_sw <- c("d_sq_XX", "F_g_XX", "d_fg_XX")
        } else {
          group_cols_sw <- c("d_sq_XX", "F_g_XX", "d_fg_XX", trends_nonparam)
        }
        
        count_sw_col <- sprintf("count_cohort_pl_%s_s_t_XX", i)
        total_sw_col <- sprintf("total_cohort_pl_%s_s_t_XX", i)
        mean_s_col   <- sprintf("mean_cohort_pl_%s_s_t_XX", i)
        
        # denominator: sum N_gt_XX within group_cols_sw for dof_s_pl==1
        df[get(dof_s_pl) == 1L,
           (count_sw_col) := sum(N_gt_XX),
           by = group_cols_sw]
        
        # numerator: sum diff_y_pl_N within group_cols_sw for dof_s_pl==1
        df[get(dof_s_pl) == 1L,
           (total_sw_col) := sum(get(diff_y_pl_N_col)),
           by = group_cols_sw]
        
        # mean
        df[, (mean_s_col) := get(total_sw_col) / get(count_sw_col)]
        
        ## --- 5. DOF for switchers (dof_cohort_pl_{i}_s_t_XX) ---
        
        if (is.null(trends_nonparam)) {
          group_keys_s <- c("d_sq_XX", "F_g_XX", "d_fg_XX")
        } else {
          group_keys_s <- c("d_sq_XX", "F_g_XX", "d_fg_XX", trends_nonparam)
        }
        
        col_dof_cohs <- sprintf("dof_cohort_pl_%s_s_t_XX", i)
        col_dof_s    <- dof_s_pl
        
        if (is.null(cluster) || cluster == "") {
          df[get(col_dof_s) == 1L,
             (col_dof_cohs) := as.numeric(sum(get(col_dof_s))),
             by = group_keys_s]
        } else {
          clust_dof_s <- sprintf("cluster_dof_pl_%s_s_XX", i)
          df[, (clust_dof_s) := ifelse(get(col_dof_s) == 1L, get(cluster), NA)]
          
          agg_s <- df[!is.na(get(clust_dof_s)),
                      .(tmp_dof = uniqueN(get(clust_dof_s))),
                      by = group_keys_s]
          
          df[agg_s, (col_dof_cohs) := as.numeric(i.tmp_dof), on = group_keys_s]
          df[is.na(get(clust_dof_s)), (col_dof_cohs) := NA_real_]
        }
        
        ## --- 6. Union of switchers and not-yet switchers (ns_s) ---
        
        if (is.null(trends_nonparam)) {
          group_keys_any <- c("d_sq_XX", "time_XX")
        } else {
          group_keys_any <- c("d_sq_XX", "time_XX", trends_nonparam)
        }
        
        col_dof_s   <- dof_s_pl
        col_dof_ns  <- dof_ns_pl
        col_dof_any <- sprintf("dof_ns_s_pl_%s_XX", i)
        
        col_N     <- "N_gt_XX"
        col_diffN <- diff_y_pl_N_col
        
        col_count <- sprintf("count_cohort_pl_%s_ns_s_t_XX", i)
        col_total <- sprintf("total_cohort_pl_%s_ns_s_t_XX", i)
        col_mean  <- sprintf("mean_cohort_pl_%s_ns_s_t_XX", i)
        
        # mask_dof = (dof_s == 1) OR (dof_ns == 1), numeric 1/0
        df[, (col_dof_any) := as.numeric((get(col_dof_s) == 1L) | (get(col_dof_ns) == 1L))]
        
        # set NA where either dof_s or dof_ns is NA
        df[is.na(get(col_dof_s)) | is.na(get(col_dof_ns)),
           (col_dof_any) := NA_real_]
        
        # group sums for union cohort (only among mask_dof rows)
        df[get(col_dof_any) == 1,
           (col_count) := sum(get(col_N)),
           by = group_keys_any]
        
        df[get(col_dof_any) == 1,
           (col_total) := sum(get(col_diffN)),
           by = group_keys_any]
        
        df[, (col_mean) := get(col_total) / get(col_count)]
        
        col_dof_coh_any <- sprintf("dof_cohort_pl_%s_ns_s_t_XX", i)
        
        if (is.null(cluster) || cluster == "") {
          df[get(col_dof_any) == 1,
             (col_dof_coh_any) := as.numeric(sum(get(col_dof_any))),
             by = group_keys_any]
        } else {
          col_dof_union <- col_dof_any
          clust_dof_any <- sprintf("cluster_dof_pl_%s_ns_s_XX", i)
          
          df[, (clust_dof_any) := ifelse(get(col_dof_union) == 1, get(cluster), NA)]
          
          # Drop previous col_dof_coh_any if exists
          if (col_dof_coh_any %chin% names(df)) {
            df[, (col_dof_coh_any) := NULL]
          }
          
          agg_any <- df[!is.na(get(clust_dof_any)),
                        .(tmp_dof = uniqueN(get(clust_dof_any))),
                        by = group_keys_any]
          
          df[agg_any, (col_dof_coh_any) := as.numeric(i.tmp_dof), on = group_keys_any]
          df[is.na(get(clust_dof_any)), (col_dof_coh_any) := NA_real_]
        }
        
    
        df = compute_E_hat_gt_with_nans(df, i, type = "placebo")
        df = compute_DOF_gt_with_nans(df, i, type = "placebo")

        #### 6. Computing U_Gg_\ell variables for the placebos (similar to part 4, less commented)
        df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] <- i <= df$T_g_XX - 1

        if (get(paste0("N",increase_XX,"_placebo_",i,"_XX")) != 0) {

        
          df[[paste0("U_Gg_pl_",i,"_temp_XX")]] <- 
          df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] *
           (G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) * df$N_gt_XX * (df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) *  df[[paste0("diff_y_pl_",i,"_XX")]]* (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) 

          df[, paste0("U_Gg_placebo_",i,"_XX") := sum(get(paste0("U_Gg_pl_",i,"_temp_XX")), na.rm = TRUE), by = group_XX]
          df[[paste0("U_Gg_placebo_",i,"_XX")]] <- df[[paste0("U_Gg_placebo_",i,"_XX")]] * df$first_obs_by_gp_XX

          df[[paste0("count",i,"_pl_core_XX")]] <- ifelse(!is.na(df[[paste0("U_Gg_pl_",i,"_temp_XX")]]) & df[[paste0("U_Gg_pl_",i,"_temp_XX")]] != 0 | df[[paste0("U_Gg_pl_",i,"_temp_XX")]] == 0 & df[[paste0("diff_y_pl_",i,"_XX")]]==0 & (df[[paste0("dist_to_switch_pl_",i,"_XX")]] != 0 | df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]] != 0 & df[[paste0("never_change_d_pl_",i,"_XX")]] != 0), df$N_gt_XX, 0)
          df[[paste0("count",i,"_pl_core_XX")]] <- as.numeric(df[[paste0("count",i,"_pl_core_XX")]])

          df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]] <- 0

          df[[paste0("U_Gg_pl_",i,"_temp_var_XX")]] <- df[[paste0("dummy_U_Gg_pl_",i,"_XX")]] * (G_XX / get(paste0("N",increase_XX,"_placebo_",i,"_XX"))) * (df[[paste0("dist_to_switch_pl_",i,"_XX")]] - (df[[paste0("N",increase_XX,"_t_placebo_",i,"_g_XX")]]/df[[paste0("N_gt_control_placebo_",i,"_XX")]]) * df[[paste0("never_change_d_pl_",i,"_XX")]]) * (df$time_XX >= i + 1 & df$time_XX <= df$T_g_XX) * df$N_gt_XX * df[[paste0("DOF_gt_pl_",i,"_XX")]] * (df[[paste0("diff_y_pl_",i,"_XX")]] - df[[paste0("E_hat_gt_pl_",i,"_XX")]])

          if (!is.null(controls)) {
            for (l in levels_d_sq_XX) {
              df[[paste0("combined_pl",increase_XX,"_temp_",l,"_",i,"_XX")]] <- 0
              for (j in 1:count_controls) {
                for (k in 1:count_controls) {
                  df[[paste0("in_brackets_pl_",l,"_",j,"_XX")]] <- df[[paste0("in_brackets_pl_",l,"_",j,"_XX")]] + get(paste0("inv_Denom_",l,"_XX"))[j,k] * df[[paste0("in_sum_",k,"_",l,"_XX")]] *(df$d_sq_int_XX == l & df$F_g_XX >= 3)
                }
                df[[paste0("in_brackets_pl_",l,"_",j,"_XX")]] <- df[[paste0("in_brackets_pl_",l,"_",j,"_XX")]] - get(paste0("coefs_sq_",l,"_XX"))[j,1]
                df[[paste0("combined_pl",increase_XX,"_temp_",l,"_",i,"_XX")]] <-  df[[paste0("combined_pl",increase_XX,"_temp_",l,"_",i,"_XX")]] + df[[paste0("M_pl",increase_XX,"_",l,"_",j,"_",i,"_XX")]] * df[[paste0("in_brackets_pl_",l,"_",j,"_XX")]]
              }
              df[[paste0("part2_pl_switch",increase_XX,"_",i,"_XX")]] <- as.numeric(df[[paste0("part2_pl_switch",increase_XX,"_",i,"_XX")]]) + df[[paste0("combined_pl",increase_XX,"_temp_",l,"_",i,"_XX")]] 
            }
          }

          df[, paste0("U_Gg_pl_",i,"_var_XX"):= sum(get(paste0("U_Gg_pl_", i,"_temp_var_XX")), na.rm = TRUE), by = group_XX]

          if(!is.null(controls)){
            if (increase_XX == 1) {
              df[[paste0("U_Gg_pl_",i,"_var_XX")]] <- df[[paste0("U_Gg_pl_",i,"_var_XX")]] - 
                  df[[paste0("part2_pl_switch1_",i,"_XX")]]
            } else {
              df[[paste0("U_Gg_pl_",i,"_var_XX")]] <- df[[paste0("U_Gg_pl_",i,"_var_XX")]] - 
                  df[[paste0("part2_pl_switch0_",i,"_XX")]]
            }
          }

        }

        ###### 7. Computing adjustements for the normalized and trends_lin options for placebos (similar to part 4, not commented) 
        if (normalized == TRUE) {
          if (is.null(continuous)) {
            df$sum_temp_pl_XX <- ifelse(df$time_XX >= df$F_g_XX & df$time_XX <= df$F_g_XX - 1 + i & df$S_g_XX == increase_XX, df$treatment_XX - df$d_sq_XX, NA) ## FELIX 18.03.2025 (delete the !is.na(diff_y_pl_`i'_XX) condition)
          } else {
            df$sum_temp_pl_XX <- ifelse(df$time_XX >= df$F_g_XX & df$time_XX <= df$F_g_XX - 1 + i & df$S_g_XX == increase_XX, df$treatment_XX_orig - df$d_sq_XX_orig, NA) ## FELIX 18.03.2025 (delete the !is.na(diff_y_pl_`i'_XX) condition)
          }
          df[, paste0("sum_treat_until_",i,"_pl_XX") := sum(sum_temp_pl_XX, na.rm = TRUE), by = group_XX]
          df$sum_temp_pl_XX <- NULL
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
    if (isTRUE(trends_lin)) {
      for (i in 1:l_placebo_u_a_XX) {
        Ntrendslin_pl <- min(Ntrendslin_pl, get(paste0("N",increase_XX,"_placebo_",i,"_XX")), na.rm = TRUE)
        }
    }
    
    if (trends_lin && Ntrendslin_pl != 0) {
      
      lp <- as.integer(l_placebo_u_a_XX)
      
      # Column names
      col_TL        <- sprintf("U_Gg_pl_%d_TL",       lp)
      col_var_TL    <- sprintf("U_Gg_pl_%d_var_TL",   lp)
      col_placebo   <- sprintf("U_Gg_placebo_%d_XX",  lp)
      col_pl_var_XX <- sprintf("U_Gg_pl_%d_var_XX",   lp)
      
      # Drop old TL columns if they exist
      cols_to_drop <- intersect(c(col_TL, col_var_TL), names(df))
      if (length(cols_to_drop) > 0L) {
        df[, (cols_to_drop) := NULL]
      }
      
      # Initialize TL columns to 0
      df[, c(col_TL, col_var_TL) := .(0.0, 0.0)]
      
      # Accumulate over i = 1..lp  (same structure as your Polars loop)
      for (i in seq_len(lp)) {
        df[, (col_TL) :=
             get(col_TL) + get(sprintf("U_Gg_placebo_%d_XX",    i))]
        df[, (col_var_TL) :=
             get(col_var_TL) + get(sprintf("U_Gg_pl_%d_var_XX", i))]
      }
      
      # Copy back into the final placebo columns
      df[, (col_placebo)   := get(col_TL)]
      df[, (col_pl_var_XX) := get(col_var_TL)]
    }

  }
  
  ###### 8. Computing Average Total Effect estimator
  if (!trends_lin) {
    ## 1) Compute sum_N{increase_XX}_l_XX  and store it with assign()
    total_key <- sprintf("sum_N%s_l_XX", increase_XX)
    
    sum_N <- sum(vapply(
      seq_len(as.integer(l_u_a_XX)),
      function(j) {
        get(sprintf("N%s_%s_XX", increase_XX, j))
      },
      numeric(1)
    ))
    
    assign(total_key, sum_N)
    
    ## 2) Initialize needed columns to 0
    init_cols <- c("U_Gg_XX", "U_Gg_num_XX", "U_Gg_den_XX",
                   "U_Gg_num_var_XX", "U_Gg_var_XX")
    df[, (init_cols) := 0]
    
    ## Loop over l = 1,...,l_u_a_XX
    for (i in seq_len(as.integer(l_u_a_XX))) {
      
      N_name          <- sprintf("N%s_%s_XX", increase_XX, i)
      N_increase      <- get(N_name)
      sum_N_increase  <- get(total_key)
      
      delta_temp      <- sprintf("delta_D_%s_temp_XX", i)
      delta           <- sprintf("delta_D_%s_XX", i)
      delta_g         <- sprintf("delta_D_g_%s_XX", i)
      dist_to_switch  <- sprintf("distance_to_switch_%s_XX", i)
      
      # Only run if N_increase != 0
      if (!is.null(N_increase) && N_increase != 0) {
        
        ## 1. Compute weight and store it with assign()
        w_i <- N_increase / sum_N_increase
        assign(sprintf("w_%s_XX", i), w_i)
        
        ## 2. Compute delta_D_temp
        if (is.null(continuous)) {
          
          df[, (delta_temp) := 0]
          
          df[get(dist_to_switch) == 1,
             (delta_temp) :=
               (N_gt_XX / N_increase) *
               ((treatment_XX - d_sq_XX) * S_g_XX +
                  (1 - S_g_XX) * (d_sq_XX - treatment_XX))
          ]
          
        } else if (continuous > 0) {
          
          den_col <- sprintf("N%s_%s_XX", increase_XX, i)
          
          df[, (delta_temp) := 0]
          
          df[get(dist_to_switch) == 1,
             (delta_temp) :=
               (N_gt_XX / get(den_col)) *
               ((treatment_XX_orig - d_sq_XX_orig) * S_g_XX +
                  (1 - S_g_XX) * (d_sq_XX_orig - treatment_XX_orig))
          ]
        }
        
        ## 3. Aggregate delta_D (sum over all obs → scalar replicated in column)
        total_delta <- df[, sum(get(delta_temp), na.rm = TRUE)]
        df[, (delta) := total_delta]
        
        ## 4. Compute delta_D_g
        df[, (delta_g) := get(delta_temp) * (N_increase / N_gt_XX)]
        
        ## 5. Drop temp
        df[, (delta_temp) := NULL]
        
        ## 6. Update U_Gg_* numerators and denominators (cumulative over i)
        U_col_i     <- sprintf("U_Gg%s_XX", i)
        U_var_col_i <- sprintf("U_Gg%s_var_XX", i)
        
        df[
          ,
          `:=`(
            U_Gg_num_XX     = U_Gg_num_XX     + w_i * get(U_col_i),
            U_Gg_num_var_XX = U_Gg_num_var_XX + w_i * get(U_var_col_i),
            U_Gg_den_XX     = U_Gg_den_XX     + w_i * get(delta)
          ),
          by = "group_XX"
        ]
      }
    }
    
    ## Final ratios
    df[, U_Gg_XX     := U_Gg_num_XX     / U_Gg_den_XX]
    df[, U_Gg_var_XX := U_Gg_num_var_XX / U_Gg_den_XX]
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
