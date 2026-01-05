#' Internal function of did_multiplegt_dyn
#' @param df df
#' @param outcome outcome
#' @param group group
#' @param time time
#' @param treatment treatment
#' @param effects effects
#' @param placebo placebo
#' @param ci_level ci_level
#' @param switchers switchers
#' @param only_never_switchers only_never_switchers
#' @param trends_nonparam trends_nonparam
#' @param weight weight
#' @param controls controls
#' @param dont_drop_larger_lower dont_drop_larger_lower
#' @param drop_if_d_miss_before_first_switch drop_if_d_miss_before_first_switch
#' @param cluster cluster
#' @param same_switchers same_switchers
#' @param same_switchers_pl same_switchers_pl
#' @param effects_equal effects_equal
#' @param save_results save_results
#' @param normalized normalized
#' @param predict_het predict_het
#' @param trends_lin trends_lin
#' @param less_conservative_se less_conservative_se
#' @param continuous continuous
#' @param data_only data_only
#' @import data.table
#' @importFrom dplyr lag
#' @importFrom matlib Ginv 
#' @importFrom plm pdata.frame make.pbalanced
#' @importFrom utils write.csv
#' @importFrom stats pchisq qnorm sd weighted.mean as.formula df.residual lm nobs qt relevel
#' @importFrom stats na.omit predict setNames
#' @importFrom MASS ginv
#' @importFrom fixest feols
#' @import lmtest
#' @import sandwich
#' @importFrom car linearHypothesis
#' @returns A list with the final estimation dataframe and other relevant matrices and scalars.
#' @noRd 
did_multiplegt_main <- function(
  df, 
  outcome, 
  group, 
  time, 
  treatment, 
  effects, 
  placebo, 
  ci_level, 
  switchers, 
  only_never_switchers,
  trends_nonparam, 
  weight, 
  controls, 
  dont_drop_larger_lower, 
  drop_if_d_miss_before_first_switch, 
  cluster, 
  same_switchers, 
  same_switchers_pl,
  effects_equal, 
  save_results, 
  normalized,
  predict_het,
  trends_lin,
  less_conservative_se,
  continuous,
  data_only = FALSE
  ) {

suppressWarnings({

  ###### 0. Pre-allocate variables that are generated via data.table (to satisfy CRAN requirements)
  gr_id <- NULL
  weight_XX <- NULL
  F_g_XX <- NULL
  F_g_trunc_XX <- NULL
  N_gt_XX <- NULL
  T_g_XX <- NULL
  U_Gg_var_global_XX <- NULL
  Yg_Fg_min1_XX <- NULL
  Yg_Fg_min2_XX <- NULL
  avg_diff_temp_XX <- NULL
  avg_post_switch_treat_XX   <- NULL
  avg_post_switch_treat_XX_temp <- NULL
  clust_U_Gg_var_global_XX <- NULL
  cluster_XX <- NULL
  cluster_var_g_XX <- NULL
  controls_time_XX <- NULL
  count_time_post_switch_XX<- NULL
  count_time_post_switch_XX_temp <- NULL
  counter <- NULL
  counter_temp <- NULL
  d_F_g_XX <- NULL
  d_F_g_temp_XX <- NULL
  d_fg_XX <- NULL
  d_sq_XX <- NULL
  d_sq_int_XX <- NULL
  d_sq_temp_XX <- NULL
  diff_y_XX <- NULL
  ever_change_d_XX <- NULL
  fd_X_all_non_missing_XX <- NULL
  first_obs_by_clust_XX <- NULL
  first_obs_by_gp_XX <- NULL
  group_XX <- NULL
  last_obs_D_bef_switch_XX <- NULL
  last_obs_D_bef_switch_t_XX <- NULL
  max_time_d_nonmiss_XX <- NULL
  mean_D <- NULL
  mean_Y <- NULL
  min_time_d_miss_aft_ynm_XX <- NULL
  min_time_d_nonmiss_XX <- NULL
  min_time_y_nonmiss_XX <- NULL
  never_change_d_XX <- NULL
  sd_het <- NULL
  sum_weights_control_XX<- NULL
  temp_F_g_XX <- NULL
  time_XX <- NULL
  time_d_miss_XX <- NULL
  time_d_nonmiss_XX<- NULL
  time_y_nonmiss_XX <- NULL
  treatment_XX_v1 <- NULL
  var_F_g_XX<- NULL


  ######## 1. Checking that syntax correctly specified
  #### Add a stop message: same_switchers_pl only works when same_switchers is specified.
  if (same_switchers == FALSE & same_switchers_pl == TRUE) {
    stop("The same_switchers_pl option only works if same_switchers is specified as well!")
  }


  #### Continous option: checking that polynomial order specified, and putting it into degree_pol scalar.
  if (!is.null(continuous)) {
    degree_pol <- continuous
  }

  ######## 2. Data preparation steps
  #### Renaming the variables in the dataset 
  original_names <- c(c(outcome, group, time, treatment), trends_nonparam, weight, controls, cluster, unlist(predict_het[1]))
  df <- subset(df, select = original_names)
  df <- data.table::setnames(df, old = c(outcome, group, time, treatment), new = c("outcome", "group", "time", "treatment"))
  df <- data.table(df)

  #### Grouping together trends_nonparam variables
  #if (!is.null(trends_nonparam)) {
  #  df$trends_nonparam_XX <- df[trends_nonparam]
  #}

  #### Patching the cluster variable: by default, the command clusters at group level. If the user specifies clustering by group, the clustering option goes to NULL.
  if (!is.null(cluster)) {
    if (paste0(cluster) == paste0(group)) {
      # df$cluster_XX <- df[[cluster]]
      cluster <- NULL
    } else{
      df$cluster_XX <- df[[cluster]]  
    }
    
  }

  #### Selecting the sample
  ## Dropping observations with missing group or time
  df <- df[ !is.na(df$group) & !is.na(df$time) ] 
  ## Dropping observations with missing controls
  if (!is.null(controls)) {
    for (var in controls) {
      df <- subset(df, !is.na(df[[var]]))
    }
  }

  #### Further sample selection steps
  ## Dropping observations with a missing clustering variable
  if (!is.null(cluster)) {
    df <- subset(df, !is.na(df$cluster_XX))
  }

  ## Dropping groups with always missing treatment or outcomes
  df[, mean_D := mean(treatment, na.rm = TRUE), by = group]
  df[, mean_Y := mean(outcome, na.rm = TRUE), by = group]
  df <- subset(df, !is.na(df$mean_Y) & !is.na(df$mean_D))
  df$mean_Y <- df$mean_D <- NULL

  #### Predict_het option for heterogeneous treatment effects analysis
  predict_het_good <- c()
  if (!is.null(predict_het)) {
    if (length(predict_het) != 2 & inherits(predict_het, "list")) {
      stop("Syntax error in predict_hat option: list with 2 elements required. Set the second element to -1 to include all the effects.")
    }
    ## Checks if predict_het and normalized are both specified
    if (isTRUE(normalized)) {
      message("The options normalized and predict_het cannot be specified together. The option predict_het will be ignored.")
    } else {
      pred_het <- unlist(predict_het[1])
      het_effects <- unlist(predict_het[2])
      ## Checks if only time-invariant variables are specified in predict_het
      for (v in pred_het) {
        df[, sd_het := fifelse(is.na(sd(get(v), na.rm = TRUE)),0,sd(get(v), na.rm = TRUE)), by = group]
        if (mean(df$sd_het) == 0) {
          predict_het_good <- c(predict_het_good, v)
        } else {
          message(sprintf("The variable %s specified in the option predict_het is time-varying, the command will therefore ignore it.", v))
        }
        df$sd_het <- NULL
      }
    }     
  }

  #### Collapse and weight
  ## Creating the weight variable 
  if (is.null(weight)) {
    df$weight_XX <- 1
  } else{
    df$weight_XX <- df[[weight]]
  }
  df$weight_XX <- ifelse(is.na(df$weight_XX), 0, df$weight_XX)

  ## Checking if the data has to be collapsed
  df$counter_temp <- 1
  df[, counter := sum(counter_temp), by = list(group, time)]
  aggregated_data <- max(df$counter) == 1
  df$counter <- df$counter_temp <- NULL 

  ## Collapsing the data if necessary
  if (aggregated_data != 1) {
    df$weight_XX <- ifelse(is.na(df$treatment), 0, df$weight_XX)
    if (is.null(cluster)) {
      df$cluster_XX <- 1
    }
    
    df_1 <- df[, lapply(.SD, weighted.mean, na.rm = TRUE, w = weight_XX), .SDcols = c("treatment", "outcome", trends_nonparam, weight, controls, predict_het_good, "cluster_XX"), by = c("group", "time")]
    df_2 <- df[, lapply(.SD, sum, na.rm = TRUE), .SDcols = "weight_XX", by = c("group", "time")]
    df <- merge(df_1, df_2, by = c("group", "time"))
    df_1 <- NULL; df_2 <- NULL;
    
    if (is.null(cluster)) {
      df$cluster_XX <- NULL
    }
  }

  ## --- Generate factorized versions of Y, G, T and D ---
  outcome <- "outcome"
  group <- "group"
  time <- "time"
  treatment <- "treatment"

  # outcome_XX = outcome
  df[, outcome_XX := outcome]

  # sort by time (like df.sort_values(time))
  setorder(df, "time")

  # group_XX and time_XX as "factorized" (1,2,3,...) in order of appearance
  df[, group_XX := .GRP, by = group]
  df[, time_XX  := .GRP, by = time ]
  df[, treatment_XX := treatment]


  # first/last date where D not missing
  df[, time_d_nonmiss_XX := ifelse(!is.na(treatment_XX), time_XX, NA_real_)]
  # first date where Y not missing
  df[, time_y_nonmiss_XX := ifelse(!is.na(outcome_XX),   time_XX, NA_real_)]

  # per-group mins & max like grp.transform('min'/'max', skipna=True)
  df[, `:=`(
    min_time_d_nonmiss_XX = {
      x <- time_d_nonmiss_XX
      if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
    },
    max_time_d_nonmiss_XX = {
      x <- time_d_nonmiss_XX
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    },
    min_time_y_nonmiss_XX = {
      x <- time_y_nonmiss_XX
      if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
    }
  ), by = group_XX]

  # first date D missing *after* Y seen
  df[, time_d_miss_XX :=
      ifelse(is.na(treatment_XX) & time_XX >= min_time_y_nonmiss_XX,
              time_XX, NA_real_)]

  # per-group: min_time_d_miss_aft_ynm_XX = grp['time_d_miss_XX'].transform('min')
  df[, min_time_d_miss_aft_ynm_XX := {
    x <- time_d_miss_XX
    if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  },
  by = group_XX]

  # drop intermediate cols (like df.drop(..., axis=1))
  df[, c("time_d_nonmiss_XX", "time_y_nonmiss_XX", "time_d_miss_XX") := NULL]


  ## --- Baseline treatment D_{g,1} ---

  # d_sq_temp_XX = treatment_XX at min_time_d_nonmiss_XX
  df[, d_sq_temp_XX :=
      ifelse(time_XX == min_time_d_nonmiss_XX, treatment_XX, NA_real_)]

  # d_sq_XX = group mean of that (only one non-NA per group, so it's the baseline)
  df[, d_sq_XX := {
    x <- d_sq_temp_XX
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  },
  by = group_XX]

  # drop temp
  df[, d_sq_temp_XX := NULL]


  ## --- Enforce "Design Restriction 2" ---

  df[, diff_from_sq_XX := treatment_XX - d_sq_XX]

  # sort by group_XX, time_XX like df.sort_values(['group_XX','time_XX'])
  setorder(df, group_XX, time_XX)

  # T_XX = int(df['time_XX'].max())
  T_XX <- as.integer(df[, max(time_XX, na.rm = TRUE)])


  if (!(dont_drop_larger_lower == TRUE)) {
    # Sort by group_XX and time_XX
    setorder(df, group_XX, time_XX)
    
    # 2. strict increase: ever_strict_increase_XX is 1 if it ever happens within group_XX
    df[, ever_strict_increase_XX :=
        as.integer(pmin(1L,
                        cumsum(diff_from_sq_XX > 0 & !is.na(treatment_XX)))),
      by = group_XX]
    
    # 3. strict decrease: ever_strict_decrease_XX is 1 if it ever happens within group_XX
    df[, ever_strict_decrease_XX :=
        as.integer(pmin(1L,
                        cumsum(diff_from_sq_XX < 0 & !is.na(treatment_XX)))),
      by = group_XX]
    
    # 4. drop rows where both == 1
    df <- df[!(ever_strict_increase_XX == 1L & ever_strict_decrease_XX == 1L)]
    
    
  }
  #### Counting number of groups
  # G_XX <- max(df$group_XX, na.rm = TRUE)

  #### Ever changed treatment 
  df$ever_change_d_XX <- abs(df$diff_from_sq_XX) > 0 & !is.na(df$treatment_XX) 
  for (i in 2:T_XX) {
    df$ever_change_d_XX[shift(df$ever_change_d_XX) == 1 & df$group_XX == shift(df$group_XX) & df$time_XX == i] <- 1
  }

  df <- data.table(df)
  #### Creating date of the first treatment change
  df <- df[order(df$group_XX, df$time_XX), ]
  df$temp_F_g_XX <- ifelse(df$ever_change_d_XX == 1 & shift(df$ever_change_d_XX) == 0, df$time_XX, 0)
  df[, F_g_XX := max(temp_F_g_XX, na.rm = TRUE), by = group_XX]
  df$temp_F_g_XX <- NULL

  #### If continuous option specified, generating polynomials of D_{g,1},
  #### storing D_{g,1} somewhere, and replacing it by 0.
  if (!is.null(continuous)) {
    for (pol_level in 1:degree_pol) {
      df[[paste0("d_sq_",pol_level,"_XX")]] <- df$d_sq_XX^pol_level
    }
    df$d_sq_XX_orig <- df$d_sq_XX
    df$d_sq_XX <- 0
  }

  ## Creating a new value with integer levels of d_sq_XX
  df[, d_sq_int_XX := .GRP, by = d_sq_XX]
  df$d_sq_int_XX <- as.numeric(as.character(df$d_sq_int_XX))

  #### Dropping values of baseline treatment such that there is no variance in F_g within
  df <- data.table(df)
  df[, var_F_g_XX := sd(F_g_XX), by = c( "d_sq_XX", trends_nonparam)]
  df <- subset(df, df$var_F_g_XX > 0)
  df$var_F_g_XX <- NULL

  #### Counting number of groups
  G_XX <- uniqueN(df$group_XX, na.rm = TRUE)


  if (nrow(df) == 0) {
    stop("No treatment effect can be estimated.\n  This is because Design Restriction 1 in de Chaisemartin & D'Haultfoeuille (2024) is not satisfied in the data, given the options requested.\n  This may be due to the fact that groups' period-one treatment is continuous, or takes a large number of values, and you have not specified the continuous option.\n  If so, you can try to specify this option.\n  If the issue persists even with this option, this means that all groups experience their first treatment change at the same date.\n  In this situation, estimators of de Chaisemartin & D'Haultfoeuille (2024) cannot be used.")
  }

  #### For each value of d_sq_XX, we drop time periods such that we do not have any control with the same baseline treatment afterwards
  #### This means the panel is no longer balanced, though it is balanced within values of the baseline treatment
  df$never_change_d_XX <- 1 - df$ever_change_d_XX 
  df[, controls_time_XX := max(never_change_d_XX), by = c("time_XX", "d_sq_XX", trends_nonparam)]
  df <- subset(df, df$controls_time_XX > 0)

  #### Computing t_min, T_max and adjusting F_g by last period pluc one for those that never change treatment
  t_min_XX <- min(df$time_XX)
  T_max_XX <- max(df$time_XX)
  df$F_g_XX[df$F_g_XX == 0] <- T_max_XX + 1



  ######## Dealing with missing treatments: most conservative option
  #### Let FMD_g denote the first date when g's treatment is missing while y has been not missing at least once, so that we know for sure that g already exists. 
  #### If that date is before the first period when g's treatment changes, we do not know when g's treatment has changed for the first time. Then, a conservative option is to drop all of g's outcomes starting at FMD_g.

  if (drop_if_d_miss_before_first_switch == TRUE) {
    df$outcome_XX <- 
      ifelse(
        ifelse(!is.na(df$min_time_d_miss_aft_ynm_XX),
              df$min_time_d_miss_aft_ynm_XX < df$F_g_XX & df$time_XX >= df$min_time_d_miss_aft_ynm_XX, FALSE), 
        NA, df$outcome_XX)
    ## RA 23/06/25 : do not filter if min_time_d_miss_aft_ynm_XX is na
  }

  ######## Dealing with missing treatments: most liberal option
  #### Let FD_g and LD_g respectively denote the first and last period where a group's treatment is non missing. Let FY_g denote the first period where a group's outcome is non missing.
  #### For groups that experience at least one treatment change, let LDBF_g denote the last date before F_g where g's treatment is non missing. We have FD_g<=LDBF_g<F_g<=LD_g, and we will deal with missing treatments depending on when they occur with respect to those four dates. 

  df$last_obs_D_bef_switch_t_XX <- ifelse(df$time_XX < df$F_g_XX & !is.na(df$treatment_XX), df$time_XX, NA)
  df[, last_obs_D_bef_switch_XX := max(last_obs_D_bef_switch_t_XX, na.rm = TRUE), by = group_XX]

  #### For groups that do not experience a treatment change, we just have FD_g<=LD_g, and we will deal with missing treatments depending on when they occur with respect to those two dates.
  #### For t<FD_g, by default we are going to consider that g joins the panel at FD_g: any non-missing outcome before FD_g replaced as missing, but all non-missing outcomes after FD_g are kept. For groups such that FY_g<FD_g, this is a "liberal" convention: those groups exist before FD_g, so one could argue that their status quo treatment is missing and they should be dropped from the analysis. We give the user the option to do that, with drop_if_d_miss_before_first_switch option

  df$outcome_XX[df$time_XX < df$min_time_d_nonmiss_XX] <- NA

  #### For groups that experience a treatment change, if D_gt missing at FD_g<t<LDBF_g, we replace their missing treatment by their status-quo treatment. Again, this is a liberal convention, so we give the user the option to not use those observations, with drop_if_d_miss_before_first_switch option.

  df$treatment_XX <- ifelse(df$F_g_XX < T_max_XX + 1 & is.na(df$treatment_XX) & df$time_XX < df$last_obs_D_bef_switch_XX & df$time_XX > df$min_time_d_nonmiss_XX, df$d_sq_XX, df$treatment_XX)

  #### For groups that experience a treatment change, if D_gt missing at LDBF_g<t<F_g (equivalent to LDBF_g<F_g-1), we cannot know the exact date when their treatment has changed, even in a binary and staggered design. Therefore, we set their outcomes at missing starting at LDBF_g+1. We also redefine their F_g as T+1 because they are effectively control groups. We also define the trunc_control_XX as LDBF_g+1 for them, because they can only be used as controls till that date.

  df$outcome_XX <- ifelse(df$F_g_XX < T_max_XX + 1 & df$time_XX > df$last_obs_D_bef_switch_XX & df$last_obs_D_bef_switch_XX < df$F_g_XX - 1, NA, df$outcome_XX)
  df$trunc_control_XX <- ifelse(df$F_g_XX < T_max_XX + 1 & df$last_obs_D_bef_switch_XX < df$F_g_XX - 1, df$last_obs_D_bef_switch_XX + 1, NA)
  df$F_g_XX[df$F_g_XX < T_max_XX + 1 & df$last_obs_D_bef_switch_XX < df$F_g_XX - 1] <- T_max_XX + 1

  #### For groups that experience a treatment change, if D_gt missing at F_g<t, we replace their missing treatment by D(g,F_g). This is again a liberal convention, but it is innocuous for the reduced-form parameters DID_l, so we do not give the user the option to overrule it (Note that overruling it could make the same_switchers option fail)

  df$d_F_g_temp_XX <- ifelse(df$time_XX == df$F_g_XX, df$treatment_XX, NA)
  df[, d_F_g_XX := mean(d_F_g_temp_XX, na.rm = TRUE), by = group_XX]
  df$treatment_XX <- ifelse(df$F_g_XX < T_max_XX + 1 & is.na(df$treatment_XX) & df$time_XX > df$F_g_XX & df$last_obs_D_bef_switch_XX == df$F_g_XX - 1, df$d_F_g_XX, df$treatment_XX)

  #### *For groups that do not experience a treatment change, if D_gt missing at FD_g<t<LD_g, we replace their missing treatment by D_g1. This is again a liberal convention, so we give the user the option to not use those observations, wi1th drop_if_d_miss_before_first_switch option.

  df$treatment_XX <- ifelse(df$F_g_XX == T_max_XX + 1 & is.na(df$treatment_XX) & df$time_XX > df$min_time_d_nonmiss_XX & df$time_XX < df$max_time_d_nonmiss_XX, df$d_sq_XX, df$treatment_XX)

  #### For groups that do not experience a treatment change, we replace all their outcomes by missing at t>LD_g. Even in a binary and staggered design, we cannot infer their treatment at t>LD_g.

  df$outcome_XX[df$F_g_XX == T_max_XX + 1 &df$time_XX > df$max_time_d_nonmiss_XX] <- NA
  df$trunc_control_XX <- ifelse(df$F_g_XX == T_max_XX + 1, df$max_time_d_nonmiss_XX + 1, df$trunc_control_XX)

  #### Store the outcome in levels, will be useful later when predict_het and trends_lin specified
  if (!is.null(predict_het)) {
    if (length(predict_het_good) > 0) {
      df$outcome_non_diff_XX <- df$outcome_XX
    }
  }

  #### When the trends_lin option is specified, drop units for which F_g_XX == 2 and redefine outcome and controls in first difference
  if (isTRUE(trends_lin)) {
    
    # make sure df is a data.table
    # if not already: data.table::setDT(df)
    
    df <- df[F_g_XX != 2]
    data.table::setorder(df, group_XX, time_XX)
    
    for (v in c("outcome_XX", controls)) {
      # first difference by group
      df[, (v) := get(v) - data.table::shift(get(v)), by = group_XX]
    }
    
    df <- df[time_XX != 1]
    t_min_XX <- min(df$time_XX)
  }


  #### Balancing the panel
  df$joint_trends_XX <- NULL
  df <- pdata.frame(df, index = c("group_XX", "time_XX")) 
  df <- make.pbalanced(df, balance.type = "fill")
  df$time_XX <- as.numeric(as.character(df$time_XX))
  df$group_XX <- as.numeric(as.character(df$group_XX))

  df <- data.table(df)

  df[, d_sq_XX := mean(d_sq_XX, na.rm = TRUE), by = group_XX]
  df[, d_sq_int_XX := mean(d_sq_int_XX, na.rm = TRUE), by = group_XX]
  df[, F_g_XX := mean(F_g_XX, na.rm = TRUE), by = group_XX]

  #### Defining N_gt, the weight of each (g,t) cell
  df$N_gt_XX <- 1
  df$N_gt_XX <- ifelse(is.na(df$outcome_XX) | is.na(df$treatment_XX), 0, df$weight_XX * df$N_gt_XX)

  #### Determining last period where g still has a control group:
  #### There is still a group with same 
  #### treatment as g's in period 1 and whose treatment has not changed since 
  #### start of panel. Definition adapted from the paper, to account for 
  #### imbalanced panel.
  df$F_g_trunc_XX <- ifelse(df$F_g_XX < df$trunc_control_XX, df$F_g_XX, df$trunc_control_XX)
  df$F_g_trunc_XX <- ifelse(is.na(df$trunc_control_XX), df$F_g_XX, df$F_g_trunc_XX)
  df$F_g_trunc_XX <- ifelse(is.na(df$F_g_XX), df$trunc_control_XX, df$F_g_trunc_XX)

  df[, T_g_XX := max(F_g_trunc_XX, na.rm = TRUE), by = c("d_sq_XX", trends_nonparam)]
  df$T_g_XX <- df$T_g_XX - 1

  #### Defining S_g: 
  #### an indicator variable for groups whose average post switch 
  #### treatment value is larger than their initial treatment D_{g,1}. 
  #### They will be considered switchers in. If S_g==0, the group is a switcher out. 
  #### For never-switchers, S_g is undefined.
  #### Definition of S_g matches that in paper, unless dont_drop_larger_lower specified.

  df$treatment_XX_v1 <- ifelse(df$time_XX >= df$F_g_XX & df$time_XX <= df$T_g_XX, df$treatment_XX, NA)

  # Assuming df is your data.table
  df[, avg_post_switch_treat_XX_temp := 
      ifelse(time_XX >= F_g_XX & time_XX <= T_g_XX, treatment_XX, NA_real_)]

  # Count of non-missing treatment observations in the post-switch period
  df[, count_time_post_switch_XX_temp := 
      ifelse(time_XX >= F_g_XX & time_XX <= T_g_XX, !is.na(treatment_XX), FALSE)]

  # Sum within group
  df[, avg_post_switch_treat_XX_temp := sum(avg_post_switch_treat_XX_temp, na.rm = TRUE), by = group_XX]
  df[, count_time_post_switch_XX := sum(count_time_post_switch_XX_temp, na.rm = TRUE), by = group_XX]

  # Divide sum by count to get the group-specific average
  df[, avg_post_switch_treat_XX_temp := avg_post_switch_treat_XX_temp / count_time_post_switch_XX]

  # Get the mean of that average across group
  df[, avg_post_switch_treat_XX := mean(avg_post_switch_treat_XX_temp, na.rm = TRUE), by = group_XX]





  df[, avg_post_switch_treat_XX_temp := sum(treatment_XX_v1, na.rm = TRUE), by = group_XX]
  df$treatment_XX_v1 <- NULL

  df$count_time_post_switch_XX_temp <- (df$time_XX >= df$F_g_XX & df$time_XX <= df$T_g_XX & !is.na(df$treatment_XX))
  df[, count_time_post_switch_XX := sum(count_time_post_switch_XX_temp, na.rm = TRUE), by = group_XX]
  df$avg_post_switch_treat_XX_temp <- df$avg_post_switch_treat_XX_temp / df$count_time_post_switch_XX
  df[, avg_post_switch_treat_XX := mean(avg_post_switch_treat_XX_temp, na.rm = TRUE), by = group_XX]
  df$avg_post_switch_treat_XX_temp <- NULL

  #### When a group is a switching group, but its average post-treatment treatment 
  #### value is exactly equal to its baseline treatment, we cannnot classify it as 
  #### a swicher in or a switcher out, but it is not a control either. 
  #### As such, we drop it from the estimation. Those groups are referred to 
  #### as no-first-stage-switchers. This issue can only arise 
  #### if dont_drop_larger_lower specified. 
  #### if continuous is specified we do this according to the original 
  #### baseline treatment and not to the one set to 0 to correctly
  #### track if a group is switcher in or switcher out.

  if (is.null(continuous)) {
    df <- subset(df, !(df$avg_post_switch_treat_XX == df$d_sq_XX & !is.na(df$avg_post_switch_treat_XX) & df$F_g_XX != df$T_g_XX + 1 & !is.na(df$F_g_XX) & !is.na(df$T_g_XX)))
    df$S_g_XX <- as.numeric(df$avg_post_switch_treat_XX > df$d_sq_XX)
    df$S_g_XX <- ifelse(df$F_g_XX != T_max_XX + 1, df$S_g_XX, NA)
  } else {
    df <- subset(df, !(df$avg_post_switch_treat_XX == df$d_sq_XX_orig  & !is.na(df$avg_post_switch_treat_XX) & df$F_g_XX != df$T_g_XX + 1 & !is.na(df$F_g_XX) & !is.na(df$T_g_XX)))
    df$S_g_XX <- as.numeric(df$avg_post_switch_treat_XX > df$d_sq_XX_orig)
    df$S_g_XX <- ifelse(df$F_g_XX != T_max_XX + 1, df$S_g_XX, NA)
  }

  #### Define another version where S_g=-1 for switchers out, which we need 
  #### when predict_het or continuous specified.
  if (length(predict_het) > 0 | !is.null(continuous)) {
    df$S_g_het_XX <- ifelse(df$S_g_XX == 0, -1, df$S_g_XX)
  }

  #### If continuous option specified: binarizing and staggerizing treatment,
  #### and adding time_FEs interacted with D_{g,1} as controls
  if (!is.null(continuous)) {
    ## Binarizing and staggerizing treatment
    df$treatment_temp_XX <- ifelse(!is.na(df$S_g_het_XX), 
                                  as.numeric((df$F_g_XX <= df$time_XX) * df$S_g_het_XX), NA)
    df$treatment_XX_orig <- df$treatment_XX
    df$treatment_XX <- df$treatment_temp_XX
    ## Enriching controls
    time_fe_XX <- levels(factor(df$time_XX))
    for (j in 2:length(time_fe_XX)) { 
      for (k in 1:degree_pol) {
        df[[paste0("time_fe_XX_",j,"_bt",k,"_XX")]] <- (df$time_XX >= j) * 
          df[[paste0("d_sq_",k,"_XX")]]
        
        ## Symmetry with Stata version ##
        #df[[paste0("time_fe_XX_",time_fe_XX[j],"_bt",k,"_XX")]] <- (df$time_XX >= time_fe_XX[j]) * 
        #df[[paste0("d_sq_",k,"_XX")]]
        controls <- c(controls, paste0("time_fe_XX_",j,"_bt",k,"_XX"))
      }
    }
  }

  #### Creating treatment at F_g: D_{g,F_g}
  df$d_fg_XX <- ifelse(df$time_XX == df$F_g_XX, df$treatment_XX, NA)
  df[, d_fg_XX := mean(d_fg_XX, na.rm = TRUE), by = group_XX]
  df$d_fg_XX <- ifelse(is.na(df$d_fg_XX) & df$F_g_XX == T_max_XX + 1, df$d_sq_XX, df$d_fg_XX)

  #### Creating the variable L_g_XX = T_g_XX - F_g_XX so that we can compute L_u or L_a afterwards
  df$L_g_XX <- df$T_g_XX - df$F_g_XX + 1

  #### Creating the equivalent variable L_g_placebo_XX for placebos
  if (placebo > 0) {
    df$L_g_placebo_XX <- ifelse(df$F_g_XX >= 3, ifelse(df$L_g_XX > df$F_g_XX - 2, df$F_g_XX - 2, df$L_g_XX), NA)
    df$L_g_placebo_XX <- ifelse(df$L_g_placebo_XX == Inf, NA, df$L_g_placebo_XX)
  }

  #### Tagging first observation of each group_XX
  df <- df[order(df$group_XX, df$time_XX), ]
  df[, first_obs_by_gp_XX := seq_len(.N) == 1, by = group_XX]
  df$first_obs_by_gp_XX <- as.numeric(df$first_obs_by_gp_XX)

  #### If cluster option if specified, flagging first obs in cluster and checking if the cluster variable is weakly coarser than the group one.
  if (!is.null(cluster)) {
    
    ## convert to integer (easier to manipulate)
    # df[, cluster_XX := as.integer(as.factor(cluster_XX))]
    
    ## complete missing clusters based on the min
    df$cluster_XX <- as.numeric( df$cluster_XX )
    df[, cluster_group_XX := min(cluster_XX, na.rm = TRUE), by = group_XX]
    df[, cluster_XX := fifelse(is.na(cluster_XX), cluster_group_XX, cluster_XX)]
    
    df[, first_obs_by_clust_XX := seq_len(.N) == 1, by = cluster_XX]
    df$first_obs_by_clust_XX <- as.numeric(df$first_obs_by_clust_XX)
    
    df[, cluster_var_g_XX := sd(cluster_XX, na.rm = TRUE), by = group_XX]
    ## Error message for clustering: non-nested case
    if (max(df$cluster_var_g_XX, na.rm = TRUE) > 0) {
      stop("The group variable should be nested within the clustering variable.")
    }
  }

  #### Declaring the data as panel after the changes above
  df <- pdata.frame(df, index = c("group_XX", "time_XX")) 
  df$time_XX <- as.numeric(as.character(df$time_XX))
  df$group_XX <- as.numeric(as.character(df$group_XX))
  df$diff_y_XX <- diff(df$outcome_XX)
  df$diff_d_XX <- diff(df$treatment_XX)
  df <- as.data.table(df)
  setorder(df, group_XX, time_XX)
  ######## 3. Necessary pre-estimation steps when the controls option is specified
  ######### CHeck this part

  if (!is.null(controls) && length(controls) > 0L) {
    
    ## 1) First differences of each control + missing flag
    count_controls <- 0L
    df[, fd_X_all_non_missing_XX := 1L]
    
    for (var in controls) {
      count_controls <- count_controls + 1L
      diff_col <- sprintf("diff_X%d_XX", count_controls)
      
      # group-wise first difference (assumes df sorted by time within group_XX)
      df[, (diff_col) := get(var) - shift(get(var)), by = group_XX]
      
      # if diff is NA, mark as missing in fd_X_all_non_missing_XX
      df[is.na(get(diff_col)), fd_X_all_non_missing_XX := 0L]
    }
    
    ## 2) Residualization prep
    count_controls <- 0L
    mycontrols_XX <- character(0L)
    
    for (var in controls) {
      count_controls <- count_controls + 1L
      diff_col <- sprintf("diff_X%d_XX", count_controls)
      
      # remove any helpers
      helper_cols <- intersect(
        c("sum_weights_control_XX", "avg_diff_temp_XX", "diff_y_wXX"),
        names(df)
      )
      if (length(helper_cols) > 0L) {
        df[, (helper_cols) := NULL]
      }
      
      # grouping keys
      grp_cols <- c("time_XX", "d_sq_XX",
                    if (!is.null(trends_nonparam)) trends_nonparam else character(0L))
      
      ## 2a) sum of N_gt for controls (within mask)
      # mask: not-yet-switched & valid diff_y & all control diffs non-missing
      df[, `_N_for_ctrl` := 0]   # or 0.0 if you want numeric
      
      ## fill only where condition holds
      df[ever_change_d_XX == 0 & !is.na(diff_y_XX) & fd_X_all_non_missing_XX == 1L,
        `_N_for_ctrl` := N_gt_XX]
      
      df[, sum_weights_control_XX := sum(`_N_for_ctrl`, na.rm = TRUE),
        by = grp_cols]
      
      df[!(ever_change_d_XX == 0 &
            !is.na(diff_y_XX) &
            fd_X_all_non_missing_XX == 1L),
        sum_weights_control_XX := NA_real_]
      
      df[, `_N_for_ctrl` := NULL]
      
      ## 2b) weighted sum of first-diffs
      df[, avg_diff_temp_XX := N_gt_XX * get(diff_col)]
      
      avg_col <- sprintf("avg_diff_X%d_XX", count_controls)
      
      df[, `_avg_diff_temp_masked` := 0 ]
      df[ ever_change_d_XX == 0 & !is.na(diff_y_XX) & fd_X_all_non_missing_XX == 1L, 
          `_avg_diff_temp_masked` := avg_diff_temp_XX ]
      
      df[, (avg_col) := sum(`_avg_diff_temp_masked`, na.rm = TRUE),
        by = grp_cols]
      
      df[!(ever_change_d_XX == 0 &
            !is.na(diff_y_XX) &
            fd_X_all_non_missing_XX == 1L),
        (avg_col) := NA_real_]
      
      df[, `_avg_diff_temp_masked` := NULL]
      
      # divide by sum_weights_control_XX
      df[, (avg_col) := get(avg_col) / sum_weights_control_XX]
      
      ## 2c) residual (sqrt(N) * (deltaX - avg deltaX))
      resid_col <- sprintf("resid_X%d_time_FE_XX", count_controls)
      
      df[, (resid_col) := sqrt(N_gt_XX) * (get(diff_col) - get(avg_col))]
      df[is.na(get(resid_col)), (resid_col) := 0]
      
      mycontrols_XX <- c(mycontrols_XX, resid_col)
      
      ## 2d) prepare product with deltaY
      df[, diff_y_wXX := sqrt(N_gt_XX) * diff_y_XX]
      
      prod_col <- sprintf("prod_X%d_Ngt_XX", count_controls)
      df[, (prod_col) := sqrt(N_gt_XX) * get(resid_col)]
      df[is.na(get(prod_col)), (prod_col) := 0]
    }
    
    ## Dictionaries / storage
    levels_d_sq_XX <- sort(unique(na.omit(df$d_sq_int_XX)))
    store_singular <- setNames(rep(FALSE, length(levels_d_sq_XX)),
                              as.character(levels_d_sq_XX))
    store_noresidualization_XX <- integer(0L)
    levels_d_sq_XX_final <- integer(0L)
    
    ## Loop over each baseline-treatment level
    for (l in levels_d_sq_XX) {
      useful <- df[d_sq_int_XX == l, uniqueN(F_g_XX)]
      assign(paste0("useful_res_", l, "_XX"),useful)
      
      if (useful > 1L) {
        mask_rows <- df$ever_change_d_XX == 0 &
          !is.na(df$diff_y_XX) &
          df$fd_X_all_non_missing_XX == 1L &
          df$d_sq_int_XX == l
        
        data_XX <- df[mask_rows]
        
        if (nrow(data_XX) == 0L) {
          store_singular[as.character(l)] <- TRUE
          store_noresidualization_XX <- c(store_noresidualization_XX, l)
          assign(paste0("useful_res_", l, "_XX"),1L)
          next
        }
        
        Y_vec <- data_XX$diff_y_wXX
        X_mat <- as.matrix(data_XX[, ..mycontrols_XX])
        YX <- cbind(Y_vec, X_mat, 1)
        
        overall <- crossprod(YX)
        val <- sum(overall)
        
        if (is.na(val)) {
          # Singular: cannot invert or accumulate
          store_singular[as.character(l)] <- TRUE
          store_noresidualization_XX <- c(store_noresidualization_XX, l)
          assign(paste0("useful_res_", l, "_XX"),1L)
        } else {
          k <- length(mycontrols_XX)
          idx_controls <- 1:k + 1L
          
          M <- overall[idx_controls, idx_controls, drop = FALSE]
          v <- overall[idx_controls, 1, drop = FALSE]
          
          # Theta
          theta_d <- ginv(M) %*% v
          assign(paste0("coefs_sq_", l, "_XX"), theta_d)
          levels_d_sq_XX_final <- c(levels_d_sq_XX_final, l)
          
          # check invertibility
          if (abs(det(M)) <= 1e-16) {
            store_singular[as.character(l)] <- TRUE
          }
          
          # rmax over all F_g_XX
          rmax <- df[, max(F_g_XX, na.rm = TRUE)]
          
          col_temp <- sprintf("N_c_%d_temp_XX", l)
          df[, (col_temp) :=
              time_XX >= 2 &
              time_XX <= (rmax - 1) &
              time_XX < F_g_XX &
              !is.na(diff_y_XX)]
          
          rsum <- df[get(col_temp) == TRUE, sum(N_gt_XX, na.rm = TRUE)]
          
          # dict_glob[[sprintf("inv_Denom_%d_XX", l)]] <- ginv(M) * rsum * G_XX
          assign(paste0("inv_Denom_",l,"_XX"),ginv(M) * rsum * G_XX)
        }
      }
    }
    
    ## Reconstruct store_singular_XX string using original d_sq_XX levels
    levels_d_sq_bis_XX <- sort(unique(na.omit(df$d_sq_XX)))
    singular_levels <- integer(0L)
    
    for (l in levels_d_sq_bis_XX) {
      key <- as.character(l)
      if (!is.null(store_singular[key]) && isTRUE(store_singular[key])) {
        singular_levels <- c(singular_levels, l)
      }
    }
    
    if (length(singular_levels) > 0L) {
      store_singular_XX <- paste(singular_levels, collapse = " ")
      warning(
        "Some control variables are not taken into account for groups with baseline treatment equal to:",
        store_singular_XX
      )
      warning(
        "1. For these groups, the regression of Y evolution X evolution and time-FE had fewer observations than regressors."
      )
      warning(
        "2. For these groups, one or more controls were perfectly collinear (no time variation)."
      )
    }
    
    ## Drop levels where residualization failed entirely
    if (length(store_noresidualization_XX) > 0L) {
      df <- df[!d_sq_int_XX %in% store_noresidualization_XX]
    }
    
    ## 3. Prepare for FE residualization regressions
    df[, time_FE_XX := as.integer(time_XX)]
    
    if (!"row_id" %in% names(df)) {
      df[, row_id := .I]
    }
    
    ## 4. Loop over each baseline-treatment level we actually residualized
    for (l in levels_d_sq_XX_final) {
      outcol <- sprintf("E_y_hat_gt_int_%d_XX", l)
      
      mask_rows <- df$d_sq_int_XX == l & df$F_g_XX > df$time_XX
      data_reg <- df[mask_rows]
      
      if (nrow(data_reg) == 0L) {
        df[, (outcol) := NA_real_]
        next
      }
      
      # Ensure time_FE_XX is factor and relevel so "2" is reference (ordering)
      data_reg[, time_FE_XX := factor(time_FE_XX)]
      levs <- levels(data_reg$time_FE_XX)
      if ("2" %in% levs) {
        new_order <- c("2", setdiff(levs, "2"))
        data_reg[, time_FE_XX := factor(time_FE_XX, levels = new_order)]
      }
      
      fe_terms <- sprintf("diff_X%d_XX", seq_len(count_controls))
      formula_str <- paste(
        "diff_y_XX ~",
        paste(c(fe_terms), collapse = " + "),
        "- 1",
        " | time_FE_XX"
      )
      form <- as.formula(formula_str)
      
      # weighted least squares
      model <- feols(form, data = data_reg, weights = weight_XX)
      data_reg[, y_hat := predict(model, newdata = data_reg)]
      
      # join predictions back by row_id
      pred_df <- data_reg[, .(row_id, y_hat)]
      setnames(pred_df, "y_hat", outcol)
      
      df <- merge(df, pred_df, by = c("row_id"), all.x = TRUE)
    }
    
    ## 5. Drop any numeric dummy columns if they exist
    for (t in 2:as.integer(T_max_XX)) {
      col <- sprintf("time_FE_XX%d", t)
      if (col %in% names(df)) {
        df[, (col) := NULL]
      }
    }
    
    ## 6. Drop temporary factor column
    if ("time_FE_XX" %in% names(df)) {
      df[, time_FE_XX := NULL]
    }
  }




  ###### 4. Performing the estimation and storing the results 
  ## Computing L_u/L_a, maximum number of event-study effects that can be computed
  ## for the switchers in/out, to compare them to number of effects requested,
  ## and finally determine the number of effects to be estimated.
  ## Same thing for the placebos.

  ## Initialize L_u_XX/L_a_XX
  L_u_XX <- NA
  L_a_XX <- NA
  L_placebo_u_XX <- NA
  L_placebo_a_XX <- NA

  ## For switchers in
  if (switchers == "" | switchers == "in") {
    L_u_XX <- max(df$L_g_XX[df$S_g_XX == 1], na.rm = TRUE)  
    if (length(df$S_g_XX[df$S_g_XX == 1 & !is.na(df$S_g_XX)]) == 0) {
      L_u_XX <- 0
    }
    ## For placebos
    if (placebo != 0) {
      L_placebo_u_XX <- max(df$L_g_placebo_XX[df$S_g_XX == 1], na.rm = TRUE)  
      L_placebo_u_XX <- ifelse(L_placebo_u_XX < 0,0, L_placebo_u_XX)
      ## If the trends_lin option was specified, L_placebo_u_XX should be decreased by 1
      ## because data starts at period 2 instead of 1.
      if (isTRUE(trends_lin)) {
        L_placebo_u_XX <- L_placebo_u_XX - 1
      }
    }
  }

  ## For switchers out
  if (switchers == "" | switchers == "out") {
    L_a_XX <- max(df$L_g_XX[df$S_g_XX == 0], na.rm = TRUE)  
    if (length(df$L_g_XX[df$S_g_XX == 0 & !is.na(df$S_g_XX)]) == 0) {
      L_a_XX <- 0
    }
    if (placebo != 0) {
      L_placebo_a_XX <- max(df$L_g_placebo_XX[df$S_g_XX == 0], na.rm = TRUE)  
      L_placebo_a_XX <- ifelse(L_placebo_a_XX < 0,0, L_placebo_a_XX)
      if (isTRUE(trends_lin)) {
        L_placebo_a_XX <- L_placebo_a_XX - 1
      }
    }
  }

  ## Error message if Design restriction 1 is not met
  if (
    (switchers == "in" & (is.na(L_u_XX) | L_u_XX == 0)) | 
    (switchers == "out" & (is.na(L_a_XX) | L_a_XX == 0)) | 
    (switchers == "" &  ((is.na(L_u_XX) | L_u_XX == 0) & (is.na(L_a_XX) | L_a_XX == 0)))
  ) {
    stop("No treatment effect can be estimated.\n  This is because Design Restriction 1 in de Chaisemartin & D'Haultfoeuille (2024) is not satisfied in the data, given the options requested.\n  This may be due to the fact that groups' period-one treatment is continuous, or takes a large number of values, and you have not specified the continuous option.\n  If so, you can try to specify this option.\n  If the issue persists even with this option, this means that all groups experience their first treatment change at the same date.\n  In this situation, estimators of de Chaisemartin & D'Haultfoeuille (2024) cannot be used.")
  }

  ## Checking that the number of dynamic and placebo effects requested by user
  ## are feasible, and correcting them if they are not. 

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

  # If the number of effects or placebos initially asked by user was too large, display error message
  if (l_XX < effects) {
    message(sprintf("The number of effects requested is too large. The number of effects which can be estimated is at most %.0f. The command will therefore try to estimante %.0f effect(s)", l_XX, l_XX))
  }

  if (placebo != 0) {
    if (l_placebo_XX < placebo & effects >= placebo) {
      message(sprintf("The number of placebos which can be estimated is at most %.0f.The command will therefore try to estimate %.0f placebo(s).", l_placebo_XX, l_placebo_XX))
    }
    if (effects < placebo) {
      message(sprintf("The number of placebo requested cannot be larger than the number of effects requested. The command cannot compute more than %.0f placebo(s).", l_placebo_XX))
    }
  }

  ## Adjustment to add more placebos (did_multiplegt_dyn_all_pl)
  max_pl_u_XX <- max_pl_a_XX <- max_pl_gap_u_XX <- max_pl_gap_a_XX <- 0
  df$pl_gap_XX <- ifelse(!is.na(df$S_g_XX), df$F_g_XX - 2 - df$L_g_XX, NA)
  if (switchers == "" | switchers == "in") {
    max_pl_u_XX <- max(df$F_g_XX[df$S_g_XX == 1], na.rm = TRUE) - 2
    max_pl_gap_u_XX <- max(df$pl_gap_XX[df$S_g_XX == 1], na.rm = TRUE)
  }
  if (switchers == "" | switchers == "out") {
    max_pl_a_XX <- max(df$F_g_XX[df$S_g_XX == 0], na.rm = TRUE) - 2
    max_pl_gap_a_XX <- max(df$pl_gap_XX[df$S_g_XX == 0], na.rm = TRUE)    
  }
  max_pl_XX <- max(max_pl_u_XX, max_pl_a_XX)
  max_pl_gap_XX <- max(max_pl_gap_u_XX, max_pl_gap_a_XX)
  max_pl_u_XX <- max_pl_a_XX <- max_pl_gap_u_XX <- max_pl_gap_a_XX <- df$pl_gap_XX <- NULL

  ## Generating default values for the variables which will be aggregated 
  ## after Program 2 below has been run for switchers in and for switchers out.

  inh_obj <- c()
  for (k in 1:l_XX) {
    df[[paste0("U_Gg", k, "_plus_XX")]] <- 0
    df[[paste0("U_Gg", k, "_minus_XX")]] <- 0
    df[[paste0("count", k, "_plus_XX")]] <- 0
    df[[paste0("count", k, "_minus_XX")]] <- 0
    df[[paste0("U_Gg_var_", k, "_in_XX")]] <- 0
    df[[paste0("U_Gg_var_", k, "_out_XX")]] <- 0
    df[[paste0("delta_D_g_", k, "_plus_XX")]] <- 0
    df[[paste0("delta_D_g_", k, "_minus_XX")]] <- 0
  }
  assign("sum_for_var_in_XX", 0)
  assign("sum_for_var_out_XX", 0)
  inh_obj <- c(inh_obj, "sum_for_var_in_XX", "sum_for_var_out_XX")
  if (placebo != 0) {
    for (k in 1:l_XX) {
      df[[paste0("U_Gg_pl_", k, "_plus_XX")]] <- 0
      df[[paste0("U_Gg_pl_", k, "_minus_XX")]] <- 0
      df[[paste0("count", k, "_pl_plus_XX")]] <- 0
      df[[paste0("count", k, "_pl_minus_XX")]] <- 0
      df[[paste0("U_Gg_var_pl_", k, "_in_XX")]] <- 0
      df[[paste0("U_Gg_var_pl_", k, "_out_XX")]] <- 0
    }
    assign("sum_for_var_placebo_in_XX", 0)
    assign("sum_for_var_placebo_out_XX", 0)
    inh_obj <- c(inh_obj, "sum_for_var_placebo_in_XX", "sum_for_var_placebo_out_XX")
  }

  for (i in 1:l_XX) {
    ### dw = de-weighted
    assign(paste0("N1_", i, "_XX"), 0) 
    assign(paste0("N1_", i, "_XX_new"), 0) 
    assign(paste0("N1_dw_", i, "_XX"), 0)
    assign(paste0("N0_", i, "_XX"), 0) 
    assign(paste0("N0_", i, "_XX_new"), 0)
    assign(paste0("N0_dw_", i, "_XX"), 0) 
    inh_obj <- c(inh_obj,paste0("N1_", i, "_XX"),paste0("N1_", i, "_XX_new"), 
                paste0("N0_", i, "_XX"), paste0("N0_", i, "_XX_new"),paste0("N1_dw_", i, "_XX"),
                paste0("N0_dw_", i, "_XX"))
    if (normalized == TRUE) {
      assign(paste0("delta_D_",i,"_in_XX"), 0)
      assign(paste0("delta_D_",i,"_out_XX"), 0)      
      inh_obj <- c(inh_obj,paste0("delta_D_",i,"_in_XX"),paste0("delta_D_",i,"_out_XX"))
    }
    if (placebo != 0) {
      assign(paste0("N1_placebo_", i, "_XX"), 0) 
      assign(paste0("N1_placebo_", i, "_XX_new"), 0) 
      assign(paste0("N1_dw_placebo_", i, "_XX"), 0) 
      assign(paste0("N0_placebo_", i, "_XX"), 0) 
      assign(paste0("N0_placebo_", i, "_XX_new"), 0) 
      assign(paste0("N0_dw_placebo_", i, "_XX"), 0) 
      inh_obj <- c(inh_obj,paste0("N1_placebo_", i, "_XX"),paste0("N1_placebo_", i, "_XX_new"), 
                  paste0("N0_placebo_", i, "_XX"), paste0("N0_placebo_", i, "_XX_new"), paste0("N1_dw_placebo_", i, "_XX"),
                  paste0("N0_dw_placebo_", i, "_XX"))
      if (normalized == TRUE) {
        assign(paste0("delta_D_pl_",i,"_in_XX"), 0)
        assign(paste0("delta_D_pl_",i,"_out_XX"), 0)      
        inh_obj <- c(inh_obj,paste0("delta_D_pl_",i,"_in_XX"),paste0("delta_D_pl_",i,"_out_XX"))
      }
    }
  }

  df$U_Gg_plus_XX <- 0
  df$U_Gg_minus_XX <- 0
  assign("U_Gg_den_plus_XX", 0)
  assign("U_Gg_den_minus_XX", 0)
  assign("sum_N1_l_XX", 0)
  assign("sum_N0_l_XX", 0)
  inh_obj <- c(inh_obj,"U_Gg_den_plus_XX", "U_Gg_den_minus_XX", "sum_N1_l_XX", "sum_N0_l_XX")
  df$U_Gg_var_plus_XX <- 0
  df$U_Gg_var_minus_XX <- 0  

  # Scalars previously passed as inherited objects
  # Their values will be changes through the next routines
  const <- NULL
  for (v in inh_obj) {
    const[[v]] <- get(v)
  }

  # Saving useful scalars to the Global Environment
  # Their values will be not changes through the next routines
  gs <- c("L_u_XX", "L_a_XX", "l_XX", "t_min_XX", "T_max_XX", "G_XX")
  if (placebo != 0) {
    gs <- c(gs, "L_placebo_u_XX", "L_placebo_a_XX")
  }
  # Add inheritance of controls #
  globals <- NULL
  for (v in gs) {
    globals[[v]] <- get(v)
  }

  controls_globals <- NULL
  if (!is.null(controls)) {
    controls_globals <- list()
    for (l in levels_d_sq_XX) {
      controls_globals <- append(controls_globals, get(paste0("useful_res_", l, "_XX")))
      names(controls_globals)[length(controls_globals)] <- paste0("useful_res_", l, "_XX")
      controls_globals <- append(controls_globals, list(get(paste0("coefs_sq_", l, "_XX"))))
      names(controls_globals)[length(controls_globals)] <- paste0("coefs_sq_", l, "_XX")
      controls_globals <- append(controls_globals, list(get(paste0("inv_Denom_",l,"_XX"))))
      names(controls_globals)[length(controls_globals)] <- paste0("inv_Denom_", l, "_XX")
    }
  }

  ## Initialize variable to earmark switchers by the number of the event-study effect
  df$switchers_tag_XX <- NA

  ## Store the data prior to estimation if requested
  if (isTRUE(data_only)) {
    data <- list(df,l_XX,T_max_XX)
    names(data) <- c("df", "l_XX", "T_max_XX")    
    return(data)
  }



  ## Perform the estimation: call the program did_multiplegt_dyn_core, 
  ## for switchers in and for switchers out, and store the results.

  if (switchers == "" | switchers == "in") {
    if (!is.na(L_u_XX) & L_u_XX != 0) {
      
      ## Perform the estimation of effects and placebos outside of the loop on 
      ## number of effects if trends_lin not specified
      if (isFALSE(trends_lin)) {
        data <- did_multiplegt_dyn_core(df, outcome = "outcome_XX", 
                                        group = "group_XX", time = "time_XX", cluster = cluster,
              treatment = "treatment_XX", effects = l_XX, placebo = l_placebo_XX, 
              switchers_core = "in", trends_nonparam = trends_nonparam, 
              controls = controls, same_switchers = same_switchers, 
              same_switchers_pl = same_switchers_pl, only_never_switchers = only_never_switchers, 
              normalized = normalized, globals = globals, const = const, 
              trends_lin = trends_lin, controls_globals = controls_globals, 
              less_conservative_se = less_conservative_se, continuous = continuous)
        
        df <- data$df
        data$df <- NULL
        for (e in names(data$const)) {
          const[[e]] <- data$const[[e]]
          assign(e, const[[e]])
        }
        
        # Store the number of the event-study effect for switchers-in
        for (k in 1:l_XX) {
          df$switchers_tag_XX[df[[paste0("distance_to_switch_",k,"_XX")]] == 1] <- k
        }
      }
      
      for (i in 1:l_XX) {
        ## Perform the estimation of effects inside of the loop on number of effects 
        ## if trends_lin is specified
        ## Note that if the option trends_lin was specified, same_switchers must also be specified.
        
        if (isTRUE(trends_lin)) {
          data <- did_multiplegt_dyn_core(df, outcome = "outcome_XX", group = "group_XX", 
                    time = "time_XX", treatment = "treatment_XX", cluster = cluster, 
                    effects = i, placebo = 0, switchers_core = "in", 
                    trends_nonparam = trends_nonparam, controls = controls, 
                    same_switchers = TRUE, same_switchers_pl = FALSE, 
                    only_never_switchers = only_never_switchers, normalized = normalized, 
                    globals = globals, const = const, trends_lin = trends_lin, 
                    controls_globals = controls_globals, 
                    less_conservative_se = less_conservative_se, continuous = continuous)
          
          df <- data$df
          data$df <- NULL
          for (e in names(data$const)) {
            const[[e]] <- data$const[[e]]
            assign(e, const[[e]])
          }          
          
          ## Store the number of the event-study effect for switchers-in
          col_dist_i <- sprintf("distance_to_switch_%d_XX", i)
          df[get(col_dist_i) == 1, switcher_tag_XX := i]
        }
        
        ## Store variables necessary for computation of effects.
        ## N.B.: in the case of unbalanced panels, it can happen that the U_Gg`i'_XX are not computed by program 2 (for example when y is missing). Consequently, for the command not to display an error message and continue running, we need to verify the variable is created, which is conditional on  N1_`i'_XX!=0.
        
        if (get(paste0("N1_",i,"_XX")) != 0) {
          df[[paste0("U_Gg",i,"_plus_XX")]] <- df[[paste0("U_Gg",i,"_XX")]]
          df[[paste0("count",i,"_plus_XX")]] <- df[[paste0("count",i,"_core_XX")]]
          df[[paste0("U_Gg_var_",i,"_in_XX")]] <- df[[paste0("U_Gg",i,"_var_XX")]]
          assign(paste0("N1_",i,"_XX_new"), get(paste0("N1_",i,"_XX")))
          const[[paste0("N1_",i,"_XX_new")]] <- get(paste0("N1_",i,"_XX_new"))
          
          if (normalized == TRUE) {
            assign(paste0("delta_D_",i,"_in_XX"), get(paste0("delta_norm_",i,"_XX")))
            const[[paste0("delta_D_",i,"_in_XX")]] <- get(paste0("delta_D_",i,"_in_XX"))
          }
          
          if (isFALSE(trends_lin)) {
            df[[paste0("delta_D_g_",i,"_plus_XX")]] <- df[[paste0("delta_D_g_",i,"_XX")]]
          }
        }
        
      }
      
      # Same as above for placebos.
      if (l_placebo_XX != 0) {
        for (i in 1:l_placebo_XX) {
          
          if (isTRUE(trends_lin)) {
            data <- did_multiplegt_dyn_core(df, 
                outcome = "outcome_XX", group = "group_XX", time = "time_XX", 
                cluster = cluster,
                treatment = "treatment_XX", effects = i, placebo = i, 
                switchers_core = "in", trends_nonparam = trends_nonparam, 
                controls = controls, same_switchers = TRUE, 
                same_switchers_pl = TRUE, only_never_switchers = only_never_switchers, 
                normalized = normalized, globals = globals, const = const, 
                trends_lin = trends_lin, controls_globals = controls_globals, 
                less_conservative_se = less_conservative_se, continuous = continuous)
            
            df <- data$df
            data$df <- NULL
            for (e in names(data$const)) {
              const[[e]] <- data$const[[e]]
              assign(e, const[[e]])
            }
            
            col_dist_i <- sprintf("distance_to_switch_%d_XX", i)
            df[get(col_dist_i) == 1, switcher_tag_XX := i]
            
          }
          
          if (get(paste0("N1_placebo_",i,"_XX")) != 0) {
            df[[paste0("U_Gg_pl_",i,"_plus_XX")]] <- df[[paste0("U_Gg_placebo_",i,"_XX")]]
            df[[paste0("count",i,"_pl_plus_XX")]] <- df[[paste0("count",i,"_pl_core_XX")]]
            df[[paste0("U_Gg_var_pl_",i,"_in_XX")]] <- df[[paste0("U_Gg_pl_",i,"_var_XX")]]
            assign(paste0("N1_placebo_",i,"_XX_new"), get(paste0("N1_placebo_",i,"_XX")))
            const[[paste0("N1_placebo_",i,"_XX_new")]] <- get(paste0("N1_placebo_",i,"_XX_new"))
            
            if (normalized == TRUE) {
              assign(paste0("delta_D_pl_",i,"_in_XX"), get(paste0("delta_norm_pl_",i,"_XX")))
              const[[paste0("delta_D_pl_",i,"_in_XX")]] <- get(paste0("delta_D_pl_",i,"_in_XX"))
            }
          }
          
        }
      }
      
      # Store variables necessary for computation of average effect.
      if (isFALSE(trends_lin)) {
        if (sum_N1_l_XX != 0) {
          df$U_Gg_plus_XX <- df$U_Gg_XX
          df$U_Gg_den_plus_XX <- df$U_Gg_den_XX
          df$U_Gg_var_plus_XX <- df$U_Gg_var_XX
        }
      }
    }
  }



  ######################## Puedes Volver aqui en cualquier momento ###############



  ## Same thing as above, for switchers out
  if (switchers == "" | switchers == "out") {
    if (!is.na(L_a_XX) & L_a_XX != 0) {

      if (isFALSE(trends_lin)) {
        data <- did_multiplegt_dyn_core(df, 
        outcome = "outcome_XX", group = "group_XX", time = "time_XX", 
        treatment = "treatment_XX", effects = l_XX, cluster = cluster,
        placebo = l_placebo_XX, switchers_core = "out", 
        trends_nonparam = trends_nonparam, controls = controls, 
        same_switchers = same_switchers, same_switchers_pl = same_switchers_pl, 
        only_never_switchers = only_never_switchers, normalized, globals = globals, 
        const = const, trends_lin = trends_lin, controls_globals = controls_globals, 
        less_conservative_se, continuous = continuous)

        df <- data$df
        data$df <- NULL
        for (e in names(data$const)) {
          const[[e]] <- data$const[[e]]
          assign(e, const[[e]])
        }

        for (k in 1:l_XX) {
          ## Store the number of the event-study effect for switchers-out
          df$switchers_tag_XX[df[[paste0("distance_to_switch_",k,"_XX")]] == 1] <- k
        }
      }

      for (i in 1:l_XX) {

        if (isTRUE(trends_lin)) {
          data <- did_multiplegt_dyn_core(df, outcome = "outcome_XX", 
              group = "group_XX", time = "time_XX", treatment = "treatment_XX", 
              effects = i, placebo = 0, switchers_core = "out", cluster = cluster,
              trends_nonparam = trends_nonparam, controls = controls, 
              same_switchers = TRUE, same_switchers_pl = FALSE, 
              only_never_switchers = only_never_switchers, normalized = normalized, 
              globals = globals, const = const, trends_lin = trends_lin, 
              controls_globals = controls_globals, 
              less_conservative_se = less_conservative_se, continuous = continuous)

          df <- data$df
          data$df <- NULL
          for (e in names(data$const)) {
            const[[e]] <- data$const[[e]]
            assign(e, const[[e]])
          }

          ## Store the number of the event-study effect for switchers-out
          df$switchers_tag_XX[df[[paste0("distance_to_switch_",i,"_XX")]] == 1] <- i
        }

        if (get(paste0("N0_",i,"_XX")) != 0) {
          df[[paste0("U_Gg",i,"_minus_XX")]] <- - df[[paste0("U_Gg",i,"_XX")]]
          df[[paste0("count",i,"_minus_XX")]] <- df[[paste0("count",i,"_core_XX")]]
          df[[paste0("U_Gg_var_",i,"_out_XX")]] <- - df[[paste0("U_Gg",i,"_var_XX")]] ## tks
          assign(paste0("N0_",i,"_XX_new"), get(paste0("N0_",i,"_XX")))
          const[[paste0("N0_",i,"_XX_new")]] <- get(paste0("N0_",i,"_XX_new"))

          if (normalized == TRUE) {
            assign(paste0("delta_D_",i,"_out_XX"), get(paste0("delta_norm_",i,"_XX")))
            const[[paste0("delta_D_",i,"_out_XX")]] <- get(paste0("delta_D_",i,"_out_XX"))
          }

          if (isFALSE(trends_lin)) {
            df[[paste0("delta_D_g_",i,"_minus_XX")]] <- df[[paste0("delta_D_g_",i,"_XX")]]
          }
        }
      }

      if (l_placebo_XX != 0) {
        for (i in 1:l_placebo_XX) {

          if (isTRUE(trends_lin)) {
            data <- did_multiplegt_dyn_core(df, outcome = "outcome_XX", 
                group = "group_XX", time = "time_XX", treatment = "treatment_XX", 
                effects = i, placebo = i, switchers_core = "out", 
                cluster = cluster,
                trends_nonparam = trends_nonparam, controls = controls, 
                same_switchers = TRUE, same_switchers_pl = TRUE, 
                only_never_switchers = only_never_switchers, normalized = normalized, 
                globals = globals, const = const, trends_lin = trends_lin, 
                controls_globals = controls_globals, 
                less_conservative_se = less_conservative_se, 
                continuous = continuous)

            df <- data$df
            data$df <- NULL
            for (e in names(data$const)) {
              const[[e]] <- data$const[[e]]
              assign(e, const[[e]])
            }
            df$switchers_tag_XX[df[[paste0("distance_to_switch_",i,"_XX")]] == 1] <- i
          }

          if (get(paste0("N0_placebo_",i,"_XX")) != 0) {
            df[[paste0("U_Gg_pl_",i,"_minus_XX")]] <- - df[[paste0("U_Gg_placebo_",i,"_XX")]]
            df[[paste0("count",i,"_pl_minus_XX")]] <- df[[paste0("count",i,"_pl_core_XX")]]
            df[[paste0("U_Gg_var_pl_",i,"_out_XX")]] <- - df[[paste0("U_Gg_pl_",i,"_var_XX")]] ## tks
            assign(paste0("N0_placebo_",i,"_XX_new"), get(paste0("N0_placebo_",i,"_XX")))
            const[[paste0("N0_placebo_",i,"_XX_new")]] <- get(paste0("N0_placebo_",i,"_XX_new"))

            if (normalized == TRUE) {
              assign(paste0("delta_D_pl_",i,"_out_XX"), get(paste0("delta_norm_pl_",i,"_XX")))
              const[[paste0("delta_D_pl_",i,"_out_XX")]] <- get(paste0("delta_D_pl_",i,"_out_XX"))
            }
          }
        }
      }

      if (isFALSE(trends_lin)) {
        if (sum_N0_l_XX != 0) {
          df$U_Gg_minus_XX <- - df$U_Gg_XX
          df$U_Gg_den_minus_XX <- df$U_Gg_den_XX
          df$U_Gg_var_minus_XX <- - df$U_Gg_var_XX ## tks
        }
      }
    }
  }
  rownames <- c()

  ###### 5. Computing the estimators and their variances

  # Creation of the matrix which stores all the estimators (DID_l, DID_pl, delta, etc.), their sd and the CIs
  mat_res_XX <- matrix(NA, nrow = l_XX + l_placebo_XX + 1, ncol = 9)

  ####  Computing DID_\ell
  ## Loop over the number of effects to be estimated
  for (i in 1:l_XX) {
    df[[paste0("U_Gg",i,"_global_XX")]] <- get(paste0("N1_",i,"_XX_new")) / (get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new"))) * df[[paste0("U_Gg",i,"_plus_XX")]] + get(paste0("N0_",i,"_XX_new")) / (get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new"))) * df[[paste0("U_Gg",i,"_minus_XX")]]
    df[[paste0("U_Gg",i,"_global_XX")]][df$first_obs_by_gp_XX == 0] <- NA

    df[[paste0("count",i,"_global_XX")]] <- ifelse(df[[paste0("count",i,"_plus_XX")]] > df[[paste0("count",i,"_minus_XX")]], df[[paste0("count",i,"_plus_XX")]], df[[paste0("count",i,"_minus_XX")]])
    df[[paste0("count",i,"_global_XX")]] <- ifelse(is.na(df[[paste0("count",i,"_plus_XX")]]), df[[paste0("count",i,"_minus_XX")]], df[[paste0("count",i,"_global_XX")]])
    df[[paste0("count",i,"_global_XX")]] <- ifelse(is.na(df[[paste0("count",i,"_minus_XX")]]), df[[paste0("count",i,"_plus_XX")]], df[[paste0("count",i,"_global_XX")]])
    df[[paste0("count",i,"_global_XX")]][df[[paste0("count",i,"_global_XX")]] == -Inf] <- NA
    df[[paste0("count",i,"_global_dwXX")]] <- as.numeric(!is.na(df[[paste0("count",i,"_global_XX")]]) & df[[paste0("count",i,"_global_XX")]] > 0)

    ## Computing aggregated delta_D (difference between treatments received wrt status quo, from F_g-1 to F_g-1+\ell), only needed for the normalized estimator
    if (normalized == TRUE) {
      assign(paste0("delta_D_",i,"_global_XX"),
            ( get(paste0("N1_",i,"_XX_new")) /( get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new")) ) ) * get(paste0("delta_D_",i,"_in_XX")) +  (get(paste0("N0_",i,"_XX_new"))/(get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new")))) * get(paste0("delta_D_",i,"_out_XX")))
    }

    ## Counting number of switchers DID_\ell applies to
    assign(paste0("N_switchers_effect_",i,"_XX"), get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new")))
    assign(paste0("N_switchers_effect_",i,"_dwXX"), get(paste0("N1_dw_",i,"_XX")) + get(paste0("N0_dw_",i,"_XX")))
    mat_res_XX[i,8] <- get(paste0("N_switchers_effect_",i,"_XX"))
    mat_res_XX[i,6] <- get(paste0("N_switchers_effect_",i,"_dwXX"))
    mat_res_XX[i,9] <- i
    assign(paste0("N_switchers_effect_",i), get(paste0("N_switchers_effect_",i,"_XX")))

    ## Counting number of observations used in the computation of DID_\ell
    df[[paste0("N_effect_",i,"_XX")]] <- sum(df[[paste0("count",i,"_global_XX")]], na.rm = TRUE)
    assign(paste0("N_effect_",i,"_XX"), mean(df[[paste0("N_effect_",i,"_XX")]]))
    assign(paste0("N_effect_",i,"_dwXX"), sum(df[[paste0("count",i,"_global_dwXX")]], na.rm = TRUE))
    assign(paste0("N_effect_",i), get(paste0("N_effect_",i,"_XX")))
    mat_res_XX[i,7] <- get(paste0("N_effect_",i,"_XX"))
    mat_res_XX[i,5] <- get(paste0("N_effect_",i,"_dwXX"))

    ## Error message if DID_\ell cannot be estimated
    if (get(paste0("N_switchers_effect_",i,"_XX")) == 0 | get(paste0("N_effect_",i,"_XX")) == 0) {
      message(paste0("Effect_",i,"cannot be estimated. There is no switcher or no control for this effect."))
    }

    ## Averaging the U_Gg\ell to compute DID_\ell
    df[[paste0("DID_",i,"_XX")]] <- sum(df[[paste0("U_Gg",i,"_global_XX")]], na.rm = TRUE)
    df[[paste0("DID_",i,"_XX")]] <- df[[paste0("DID_",i,"_XX")]] / G_XX

    ## Computing DID^n_\ell if the option was specified
    if (normalized == TRUE) {
      df[[paste0("DID_",i,"_XX")]] <- df[[paste0("DID_",i,"_XX")]] / get(paste0("delta_D_",i,"_global_XX"))
    }
    assign(paste0("DID_",i,"_XX"), mean(df[[paste0("DID_",i,"_XX")]]))

    ## Set DID_\ell missing when there is no switcher or no control
    if ((switchers == "" & get(paste0("N1_",i,"_XX_new")) == 0 & get(paste0("N0_",i,"_XX_new")) == 0) | (switchers == "out" & get(paste0("N0_",i,"_XX_new")) == 0 ) |
        (switchers == "in" & get(paste0("N1_",i,"_XX_new")) == 0 )) {
      assign(paste0("DID_",i,"_XX"), NA)
    }

    ## Store DID_\ell estimates in the ereturn and in the results matrix which will be printed out
    assign(paste0("Effect_",i), get(paste0("DID_",i,"_XX")))
    mat_res_XX[i,1] <- get(paste0("DID_",i,"_XX"))
    rownames <- append(rownames, paste0("Effect_",i, strrep(" ",(12 - nchar(paste0("Effect_",i))))))
  }

  ###### Computing the average total effect
  U_Gg_den_plus_XX <- ifelse(is.na(mean(df$U_Gg_den_plus_XX, na.rm = TRUE)), 0, mean(df$U_Gg_den_plus_XX, na.rm = TRUE))
  U_Gg_den_minus_XX <- ifelse(is.na(mean(df$U_Gg_den_minus_XX, na.rm = TRUE)), 0, mean(df$U_Gg_den_minus_XX, na.rm = TRUE))

  #### The average effect cannot be estimated when the trends_lin option is specified so the whole part will be skipped in that case
  if (isFALSE(trends_lin)) {

    ## Computing the weight w_+.
    if (switchers == "") {
      w_plus_XX <- U_Gg_den_plus_XX * sum_N1_l_XX / (U_Gg_den_plus_XX * sum_N1_l_XX + U_Gg_den_minus_XX * sum_N0_l_XX)
    }

    ## When the "switchers" option is used, full weight is put to either switchers out or in
    if (switchers == "out") {
      w_plus_XX <- 0
    }
    if (switchers == "in") {
      w_plus_XX <- 1
    }

    ## Aggregating the U_Gg for switchers in and out
    df$U_Gg_global_XX <- w_plus_XX * df$U_Gg_plus_XX + (1 - w_plus_XX) * df$U_Gg_minus_XX
    df$U_Gg_global_XX[df$first_obs_by_gp == 0] <- NA

    ## Averaging the U_Gg to compute average total effect
    df$delta_XX <- sum(df$U_Gg_global_XX, na.rm = TRUE) / G_XX
    delta_XX <- mean(df$delta_XX, na.rm = TRUE)
    assign("Av_tot_effect", delta_XX)

    ## Completing the results matrix
    # Storing the results
    mat_res_XX[l_XX+1,1] <- delta_XX
    # Number of switchers
    N_switchers_effect_XX <- 0
    N_switchers_effect_dwXX <- 0
    for (i in 1:l_XX) {
      N_switchers_effect_XX <- N_switchers_effect_XX + get(paste0("N_switchers_effect_",i,"_XX"))
      N_switchers_effect_dwXX <- N_switchers_effect_dwXX + get(paste0("N_switchers_effect_",i,"_dwXX"))
    }
    mat_res_XX[l_XX+1,8] <- N_switchers_effect_XX
    mat_res_XX[l_XX+1,6] <- N_switchers_effect_dwXX
    mat_res_XX[l_XX+1,9] <- 0
    assign("N_switchers_effect_average", N_switchers_effect_XX)
    # Number of observations used in the estimation
    df$count_global_XX <- 0
    for (i in 1:l_XX) {
      df$count_global_XX <- ifelse(!is.na(df[[paste0("count",i,"_global_XX")]]), ifelse(df$count_global_XX > df[[paste0("count",i,"_global_XX")]], df$count_global_XX, df[[paste0("count",i,"_global_XX")]]), df$count_global_XX)
    }
    df$count_global_dwXX <- as.numeric(!is.na(df$count_global_XX) & df$count_global_XX > 0)
    N_effect_XX <- sum(df$count_global_XX, na.rm = TRUE)
    N_effect_dwXX <- sum(df$count_global_dwXX, na.rm = TRUE)
    mat_res_XX[l_XX+1,7] <- N_effect_XX
    mat_res_XX[l_XX+1,5] <- N_effect_dwXX
    assign("N_avg_total_effect", N_effect_XX)
  }
  rownames <- append(rownames, paste0("Av_tot_eff", strrep(" ",(12 - nchar("Av_tot_eff")))))
  mat_res_XX[l_XX+1,9] <- 0

  #### Computing the placebo estimators (same steps as for the DID_\ell, not commented)

  if (l_placebo_XX != 0) {
    for (i in 1:l_placebo_XX) {
      df[[paste0("U_Gg_pl_",i,"_global_XX")]] <- get(paste0("N1_placebo_",i,"_XX_new")) / (get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new"))) * df[[paste0("U_Gg_pl_",i,"_plus_XX")]] + get(paste0("N0_placebo_",i,"_XX_new")) / (get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new"))) * df[[paste0("U_Gg_pl_",i,"_minus_XX")]]
      df[[paste0("U_Gg_pl_",i,"_global_XX")]][df$first_obs_by_gp_XX == 0] <- NA

      df[[paste0("count",i,"_pl_global_XX")]] <- ifelse(df[[paste0("count",i,"_pl_plus_XX")]] > df[[paste0("count",i,"_pl_minus_XX")]], df[[paste0("count",i,"_pl_plus_XX")]], df[[paste0("count",i,"_pl_minus_XX")]])
      df[[paste0("count",i,"_pl_global_XX")]] <- ifelse(is.na(df[[paste0("count",i,"_pl_plus_XX")]]), df[[paste0("count",i,"_pl_minus_XX")]], df[[paste0("count",i,"_pl_global_XX")]])
      df[[paste0("count",i,"_pl_global_XX")]] <- ifelse(is.na(df[[paste0("count",i,"_pl_minus_XX")]]), df[[paste0("count",i,"_pl_plus_XX")]], df[[paste0("count",i,"_pl_global_XX")]])
      df[[paste0("count",i,"_pl_global_XX")]][df[[paste0("count",i,"_pl_global_XX")]] == -Inf] <- NA
      df[[paste0("count",i,"_pl_global_dwXX")]] <- as.numeric(!is.na(df[[paste0("count",i,"_pl_global_XX")]]) & df[[paste0("count",i,"_pl_global_XX")]] > 0)

      if (normalized == TRUE) {
        assign(paste0("delta_D_pl_",i,"_global_XX"),
              ( get(paste0("N1_placebo_",i,"_XX_new")) /( get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new")) ) ) * get(paste0("delta_D_pl_",i,"_in_XX")) +  (get(paste0("N0_placebo_",i,"_XX_new"))/(get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new")))) * get(paste0("delta_D_pl_",i,"_out_XX")))
      }

      df[[paste0("DID_placebo_",i,"_XX")]] <- sum(df[[paste0("U_Gg_pl_",i,"_global_XX")]], na.rm = TRUE)
      df[[paste0("DID_placebo_",i,"_XX")]] <- df[[paste0("DID_placebo_",i,"_XX")]] / G_XX

      if (normalized == TRUE) {
        df[[paste0("DID_placebo_",i,"_XX")]] <- df[[paste0("DID_placebo_",i,"_XX")]] / get(paste0("delta_D_pl_",i,"_global_XX"))
      }
      assign(paste0("DID_placebo_",i,"_XX"), mean(df[[paste0("DID_placebo_",i,"_XX")]]))

      if ((switchers == "" & get(paste0("N1_placebo_",i,"_XX_new")) == 0 & get(paste0("N0_placebo_",i,"_XX_new")) == 0) | (switchers == "out" & get(paste0("N0_placebo_",i,"_XX_new")) == 0 ) | (switchers == "in" & get(paste0("N1_placebo_",i,"_XX_new")) == 0 )) {
        assign(paste0("DID_placebo_",i,"_XX"), NA)
      }

      assign(paste0("Placebo_",i), get(paste0("DID_placebo_",i,"_XX")))
      mat_res_XX[l_XX+1+i,1] <- get(paste0("DID_placebo_",i,"_XX"))
      rownames <- append(rownames, paste0("Placebo_",i, strrep(" ",(12 - nchar(paste0("Placebo_",i))))))

      assign(paste0("N_switchers_placebo_",i,"_XX"), get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new")))
      assign(paste0("N_switchers_placebo_",i,"_dwXX"), get(paste0("N1_dw_placebo_",i,"_XX")) + get(paste0("N0_dw_placebo_",i,"_XX")))
      mat_res_XX[l_XX+1+i,8] <- get(paste0("N_switchers_placebo_",i,"_XX"))
      mat_res_XX[l_XX+1+i,6] <- get(paste0("N_switchers_placebo_",i,"_dwXX"))
      mat_res_XX[l_XX+1+i,9] <- -i
      assign(paste0("N_switchers_placebo_",i), get(paste0("N_switchers_placebo_",i,"_XX")))
      df[[paste0("N_placebo_",i,"_XX")]] <- sum(df[[paste0("count",i,"_pl_global_XX")]], na.rm = TRUE)
      assign(paste0("N_placebo_",i,"_XX"), mean(df[[paste0("N_placebo_",i,"_XX")]]))
      assign(paste0("N_placebo_",i), get(paste0("N_placebo_",i,"_XX")))
      mat_res_XX[l_XX + 1 + i,7] <- get(paste0("N_placebo_",i,"_XX"))
      mat_res_XX[l_XX + 1 + i,5] <- sum(df[[paste0("count",i,"_pl_global_dwXX")]], na.rm = TRUE)

      if (get(paste0("N_switchers_placebo_",i,"_XX")) == 0 | get(paste0("N_placebo_",i,"_XX")) == 0) {
        message(paste0("Placebo_",i,"cannot be estimated. There is no switcher or no control for this placebo."))
      }

    }
  }

  ####  Computing the variance of DID_\ell

  ## Patch significance level for t statistics
  ci_level <- ci_level / 100
  z_level <- qnorm(ci_level + (1 - ci_level)/2)

  df <- data.table(df)
  ## Loop over the number of effects to be estimated
  for (i in 1:l_XX) {
    if ((switchers == "" & (get(paste0("N1_",i,"_XX_new")) != 0 | get(paste0("N0_",i,"_XX_new"))) != 0) | (switchers == "out" & get(paste0("N0_",i,"_XX_new")) != 0 ) | (switchers == "in" & get(paste0("N1_",i,"_XX_new")) != 0 )) {

      ## Aggregating the U_Gg_var_\ell for switchers in and out
      df[[paste0("U_Gg_var_glob_",i,"_XX")]] <- df[[paste0("U_Gg_var_",i,"_in_XX")]] * (get(paste0("N1_",i,"_XX_new")) / (get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new")))) + df[[paste0("U_Gg_var_",i,"_out_XX")]] * (get(paste0("N0_",i,"_XX_new")) / (get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new"))))
      
      
      compute_sigma_hat2_glob <- function(DT, i, cluster = NULL) {
        stopifnot(is.data.table(DT))
        if (!is.numeric(G_XX) || length(G_XX) != 1L) {
          stop("`G_XX` must be a numeric scalar.")
        }
        
        # ----- dynamic names -----
        U_col          <- sprintf("U_Gg_var_glob_%s_XX", i)
        eff_sq_col     <- sprintf("U_Gg_var_glob_eff%s_sqrd_XX", i)
        
        clust_sum_col  <- sprintf("clust_U_Gg_var_glob_%s_XX",  i)
        clust_sq_col   <- sprintf("clust_U_Gg_var_glob_%s_2_XX", i)
        
        # required base columns
        needed <- c(U_col, "first_obs_by_gp_XX")
        clustered <- !(is.null(cluster) || identical(cluster, ""))
        
        if (clustered) {
          needed <- c(needed, "first_obs_by_clust_XX", cluster)
        }
        
        miss <- setdiff(needed, names(DT))
        if (length(miss)) stop("Missing required columns in DT: ", paste(miss, collapse = ", "))
        
        # ----- NON-CLUSTERED CASE -----
        if (clustered==FALSE) {
          # gen U_Gg_var_glob_eff`i'_sqrd_XX = U_Gg_var_glob_`i'_XX^2 * first_obs_by_gp_XX
          DT[, (eff_sq_col) := (get(U_col)^2) * first_obs_by_gp_XX]
          
          # sum ... ; scalar sum_for_var_`i'_XX = r(sum) / G_XX^2
          num <- DT[, sum(get(eff_sq_col), na.rm = TRUE)]
          sum_for_var <- num / (G_XX^2)
          
          return(list(DT, sum_for_var))
        }
        
        # ----- CLUSTERED CASE -----
        # capture drop clust_* if present
        for (col in c(clust_sum_col, clust_sq_col)) {
          if (col %in% names(DT)) DT[, (col) := NULL]
        }
        
        # replace U_Gg_var_glob_`i'_XX = U_Gg_var_glob_`i'_XX * first_obs_by_gp_XX
        DT[, (U_col) := get(U_col) * first_obs_by_gp_XX]
        
        # by `cluster`: gegen clust_U_Gg_var_glob_`i'_XX = total(U_Gg_var_glob_`i'_XX)
        DT[, (clust_sum_col) := sum(get(U_col), na.rm = TRUE), by = c("cluster_XX")]
        
        # gen clust_U_Gg_var_glob_`i'_2_XX = clust_U_Gg_var_glob_`i'_XX^2 * first_obs_by_clust_XX
        DT[, (clust_sq_col) := (get(clust_sum_col)^2) * first_obs_by_clust_XX]
        
        # sum ... ; scalar sum_for_var_`i'_XX = r(sum) / G_XX^2
        num1 <- DT[, sum(get(clust_sq_col), na.rm = TRUE)]
        sum_for_var <- num1 / (G_XX^2)
        
        # replace U_Gg_var_glob_`i'_XX = clust_U_Gg_var_glob_`i'_XX
        DT[, (U_col) := get(clust_sum_col)]
        return(list(DT, sum_for_var))
      }
      
      list1 <- compute_sigma_hat2_glob( df, i = i, cluster = cluster )
      df <- list1[[1]]
      assign(paste0("sum_for_var_",i,"_XX"),list1[[2]])
      
      
      
      
      
      # ## Compute \hat{\sigma}^2_l without clustering
      # if (is.null(cluster)) {
      #   df[[paste0("U_Gg_var_glob_eff",i,"_sqrd_XX")]] <- df[[paste0("U_Gg_var_glob_",i,"_XX")]]^2 * df$first_obs_by_gp_XX
      #   assign(paste0("sum_for_var_",i,"_XX"), sum(df[[paste0("U_Gg_var_glob_eff",i,"_sqrd_XX")]], na.rm = TRUE) / G_XX^2)
      # } else {
      #   ## Compute \hat{\sigma}^2_l with clustering: sum U_Gg_var_\ell within a cluster, and then take average of square.
      #   df[[paste0("U_Gg_var_glob_",i,"_XX")]] <- df[[paste0("U_Gg_var_glob_",i,"_XX")]] * df$first_obs_by_gp_XX
      # 
      #   ## Sum within cluster
      #   df[, paste0("clust_U_Gg_var_glob_",i,"_XX") := sum(get(paste0("U_Gg_var_glob_",i,"_XX")), na.rm = TRUE), by = cluster_XX]
      #   ## Compute average of square
      #   df[[paste0("clust_U_Gg_var_glob_",i,"_2_XX")]] <-
      #     df[[paste0("clust_U_Gg_var_glob_", i, "_XX")]]^2 * df$first_obs_by_clust_XX
      #   assign(paste0("sum_for_var_",i,"_XX"),
      #          sum(df[[paste0("clust_U_Gg_var_glob_", i,"_2_XX")]], na.rm = TRUE) / G_XX^2)
      #   df[[paste0("U_Gg_var_glob_",i,"_XX")]] <-
      #     df[[paste0("clust_U_Gg_var_glob_",i,"_XX")]]
      # }

      ## Compute SE
      assign(paste0("se_",i,"_XX"), sqrt(get(paste0("sum_for_var_",i,"_XX"))))

      ## Normalize SE if normalized option was specified
      if (normalized == TRUE) {
        assign(paste0("se_",i,"_XX"), get(paste0("se_",i,"_XX")) / get(paste0("delta_D_",i,"_global_XX")))
      }

      ## Storing the results
      mat_res_XX[i,2] <- get(paste0("se_",i,"_XX"))
      assign(paste0("se_effect_",i),get(paste0("se_",i,"_XX")))
      assign(paste0("LB_CI_",i,"_XX"), get(paste0("DID_",i,"_XX")) - z_level * get(paste0("se_",i,"_XX")))
      mat_res_XX[i,3] <- get(paste0("LB_CI_",i,"_XX"))
      assign(paste0("UB_CI_",i,"_XX"), get(paste0("DID_",i,"_XX")) + z_level * get(paste0("se_",i,"_XX")))
      mat_res_XX[i,4] <- get(paste0("UB_CI_",i,"_XX"))
    }
  }

  ##  Computing the variances of the placebo estimators (same steps as for the DID_\ell, not commented)
  if (l_placebo_XX != 0) {
    for (i in 1:l_placebo_XX) {
      if ((switchers == "" & (get(paste0("N1_placebo_",i,"_XX_new")) != 0 | get(paste0("N0_placebo_",i,"_XX_new")) != 0)) | (switchers == "out" & get(paste0("N0_placebo_",i,"_XX_new")) != 0 ) | (switchers == "in" & get(paste0("N1_placebo_",i,"_XX_new")) != 0 )) {

        df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]] <- df[[paste0("U_Gg_var_pl_",i,"_in_XX")]] * (get(paste0("N1_placebo_",i,"_XX_new")) / (get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new")))) + df[[paste0("U_Gg_var_pl_",i,"_out_XX")]] * (get(paste0("N0_placebo_",i,"_XX_new")) / (get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new"))))

        if (is.null(cluster)) {
          df[[paste0("U_Gg_var_glob_pl_",i,"_2_XX")]] <- df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]]^2 * df$first_obs_by_gp_XX
          assign(paste0("sum_for_var_placebo_",i,"_XX"), sum(df[[paste0("U_Gg_var_glob_pl_",i,"_2_XX")]], na.rm = TRUE) / G_XX^2)
        } else {
          df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]] <- df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]] * df$first_obs_by_gp_XX
          df[, paste0("clust_U_Gg_var_glob_pl_",i,"_XX") := sum(get(paste0("U_Gg_var_glob_pl_",i,"_XX")), na.rm = TRUE), by = cluster_XX]
          df[[paste0("clust_U_Gg_var_glob_pl_",i,"_2_XX")]] <-
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
        assign(paste0("se_placebo_",i),get(paste0("se_placebo_",i,"_XX")))

        assign(paste0("LB_CI_placebo_",i,"_XX"), get(paste0("DID_placebo_",i,"_XX")) - z_level * get(paste0("se_placebo_",i,"_XX")))
        mat_res_XX[l_XX + 1 + i,3] <- get(paste0("LB_CI_placebo_",i,"_XX"))
        assign(paste0("UB_CI_placebo_",i,"_XX"), get(paste0("DID_placebo_",i,"_XX")) + z_level * get(paste0("se_placebo_",i,"_XX")))
        mat_res_XX[l_XX + 1 + i,4] <- get(paste0("UB_CI_placebo_",i,"_XX"))
      }
    }
  }

  ##  Computing the variance of the average total effect (same steps as for the DID_\ell, not commented)
  if (isFALSE(trends_lin)) {
    if ((switchers=="" & (sum_N1_l_XX!=0|sum_N0_l_XX!=0))|(switchers =="out" & sum_N0_l_XX!=0)|(switchers=="in" & sum_N1_l_XX !=0)) {

      df$U_Gg_var_global_XX <- w_plus_XX * df$U_Gg_var_plus_XX + (1 - w_plus_XX) * df$U_Gg_var_minus_XX

      if (is.null(cluster)) {
        df$U_Gg_var_global_2_XX <- df$U_Gg_var_global_XX^2 * df$first_obs_by_gp_XX
        assign("sum_for_var_XX", sum(df$U_Gg_var_global_2_XX, na.rm = TRUE) / G_XX^2)
      } else {
        df$U_Gg_var_global_XX <- df$U_Gg_var_global_XX * df$first_obs_by_gp_XX
        df[, clust_U_Gg_var_global_XX := sum(U_Gg_var_global_XX, na.rm = TRUE), by = cluster_XX]
        df$clust_U_Gg_var_global_XX <- df$clust_U_Gg_var_global_XX^2 * df$first_obs_by_clust_XX
        assign("sum_for_var_XX", sum(df$clust_U_Gg_var_global_XX)/G_XX^2)
      }

      assign("se_XX", sqrt(sum_for_var_XX))
      mat_res_XX[l_XX + 1,2] <- se_XX
      assign("se_avg_total_effect",se_XX)

      # CI level
      LB_CI_XX <- delta_XX - z_level * se_XX
      mat_res_XX[l_XX + 1,3] <- LB_CI_XX
      UB_CI_XX <- delta_XX + z_level * se_XX
      mat_res_XX[l_XX + 1,4] <- UB_CI_XX
    }
  }

  ## Average number of cumulated effects
  for (i in 1:l_XX) {
    df[[paste0("delta_D_g_",i,"_XX")]] <- NULL
  }
  df$M_g_XX <- ifelse(l_XX <= df$T_g_XX - df$F_g_XX + 1, l_XX, df$T_g_XX - df$F_g_XX + 1)

  #### Calling variables delta_D_g_`i'_XX here like that does not work because switcher in/out are run one after another!!!

  ## second sum over g: total ... if F_g_XX<=T_g_XX
  ## actually I think it can be done in one total as we sum over the periods within groups and then across groups which are all different cells
  ## generate one variable that stores all the different delta_D_g_`i'_XX

  df$delta_D_g_XX <- 0
  for (j in 1:l_XX) {
    df$delta_D_g_XX_temp <- ifelse(df[[paste0("delta_D_g_",j,"_plus_XX")]] != 0, df[[paste0("delta_D_g_",j,"_plus_XX")]],
                                  df[[paste0("delta_D_g_",j,"_minus_XX")]])
    df$delta_D_g_XX_temp <- ifelse(df$delta_D_g_XX_temp == 0, NA, df$delta_D_g_XX_temp)
    df$delta_D_g_XX <- ifelse(df$switchers_tag_XX == j, df$delta_D_g_XX + df$delta_D_g_XX_temp, df$delta_D_g_XX)
  }
  df$delta_D_g_num_XX <- df$delta_D_g_XX * (df$M_g_XX - (df$switchers_tag_XX - 1))
  delta_D_num_total <- sum(df$delta_D_g_num_XX, na.rm = TRUE)
  delta_D_denom_total <- sum(df$delta_D_g_XX, na.rm = TRUE)
  delta_D_avg_total <- delta_D_num_total / delta_D_denom_total
  ###### 6. Computing p-values from the tests

  # If the option cluster is specified, we have previously replaced U_Gg_var_glob_pl_`i'_XX by clust_U_Gg_var_glob_pl_`i'_XX, and U_Gg_var_glob_`i'_XX by clust_U_Gg_var_glob_`i'_XX.
  # Now, we must also replace first_obs_by_gp_XX by first_obs_by_clust_XX
  if (!is.null(cluster)) {
    df$first_obs_by_gp_XX <- df$first_obs_by_clust_XX
  }

  ###### Performing a test to see whether all effects are jointly equal to 0
  all_Ns_not_zero <- NA
  all_delta_not_zero <- NA
  p_jointeffects <- NULL
  ## Test can only be run when at least two effects requested:
  if (l_XX != 0 & l_XX > 1) {
    ## If test is feasible, initalize scalar at 0
    all_Ns_not_zero <- 0
    all_delta_not_zero <- 0

    ## Count the number of estimated effects included in the test
    for (i in 1:l_XX) {
      if ( (switchers == "" & (get(paste0("N1_",i,"_XX_new"))!= 0 | get(paste0("N0_",i,"_XX_new"))!= 0 )) | (switchers == "out" & get(paste0("N0_",i,"_XX_new")) != 0) | (switchers == "in" & get(paste0("N1_",i,"_XX_new")) != 0) ) {
        all_Ns_not_zero <- all_Ns_not_zero + 1
      }

      if (isTRUE(normalized)) {
        if (get(paste0("delta_D_",i,"_global_XX")) != 0 & !is.na(get(paste0("delta_D_",i,"_global_XX")))) {
          all_delta_not_zero <- all_delta_not_zero + 1
        }
      }
    }

    ## Test can only be run when all requested effects could be computed:
    if ((all_Ns_not_zero == l_XX & isFALSE(normalized)) | (isTRUE(normalized) & all_Ns_not_zero == l_XX & all_delta_not_zero == l_XX)) {

      ## Creating a vector with all dynamic effect estimates
      didmgt_Effects <- matrix(0, nrow = l_XX, ncol = 1)

      ## Creating a matrix where the variances and the covariances of the effects will be stored.
      didmgt_Var_Effects <- matrix(0, nrow = l_XX, ncol = l_XX)

      ## Fill those matrices
      for (i in 1:l_XX) {
        didmgt_Effects[i,1] <- get(paste0("DID_",i,"_XX"))
        didmgt_Var_Effects[i,i] <- get(paste0("se_",i,"_XX"))^2

        if (i < l_XX) {
          for (j in (i+1):l_XX) {
            ## Create variables necessary to compute the covariances
            if (normalized == FALSE) {
              df[[paste0("U_Gg_var_",i,"_",j,"_XX")]] <- df[[paste0("U_Gg_var_glob_",i,"_XX")]] +  df[[paste0("U_Gg_var_glob_",j,"_XX")]]
            } else {
              df[[paste0("U_Gg_var_",i,"_",j,"_XX")]] <- df[[paste0("U_Gg_var_glob_",i,"_XX")]] / get(paste0("delta_D_",i,"_global_XX")) +  df[[paste0("U_Gg_var_glob_",j,"_XX")]] / get(paste0("delta_D_",j,"_global_XX"))
            }

            ## Estimate the covariances
            df[[paste0("U_Gg_var_",i,"_",j,"_2_XX")]] <- df[[paste0("U_Gg_var_",i,"_",j,"_XX")]]^2 * df$first_obs_by_gp_XX
            assign(paste0("var_sum_",i,"_",j,"_XX"), sum( df[[paste0("U_Gg_var_",i,"_",j,"_2_XX")]], na.rm = TRUE) / G_XX^2)
            assign(paste0("cov_",i,"_",j,"_XX"), (get(paste0("var_sum_",i,"_",j,"_XX")) - get(paste0("se_",i,"_XX"))^2 - get(paste0("se_",j,"_XX"))^2)/2)

            ## Store the results
            didmgt_Var_Effects[i,j] <- get(paste0("cov_",i,"_",j,"_XX"))
            didmgt_Var_Effects[j,i] <- get(paste0("cov_",i,"_",j,"_XX"))
          }
        }
      }

      ## Compute P-value for the F-test on joint nullity of all effects
      didmgt_Var_Effects_inv <- Ginv(didmgt_Var_Effects)
      didmgt_chi2effects <- t(didmgt_Effects) %*% didmgt_Var_Effects_inv  %*% didmgt_Effects
      p_jointeffects <- 1 - pchisq(didmgt_chi2effects[1,1], df = l_XX)
    } else {
      p_jointeffects <- NA
      ## Error message if not all of the specified effects could be estimated
      message("Some effects could not be estimated. Therefore, the test of joint nullity of the effects could not be computed.")
    }
  }


  ###### Performing a test to see whether all placebos are jointly equal to 0
  all_Ns_pl_not_zero <- NA
  all_delta_pl_not_zero <- NA
  ## Test can only be run when at least two placebos requested:
  if (l_placebo_XX != 0 & l_placebo_XX > 1) {
    ## If test is feasible, initalize scalar at 0
    all_Ns_pl_not_zero <- 0
    all_delta_pl_not_zero <- 0

    ## Count the number of estimated placebos included in the test
    for (i in 1:l_placebo_XX) {
      if ( (switchers == "" & (get(paste0("N1_placebo_",i,"_XX_new"))!= 0 | get(paste0("N0_placebo_",i,"_XX_new"))!= 0 )) | (switchers == "out" & get(paste0("N0_placebo_",i,"_XX_new")) != 0) | (switchers == "in" & get(paste0("N1_placebo_",i,"_XX_new")) != 0) ) {
        all_Ns_pl_not_zero <- all_Ns_pl_not_zero + 1
      }

      if (isTRUE(normalized)) {
        if (get(paste0("delta_D_pl_",i,"_global_XX")) != 0 & !is.na(get(paste0("delta_D_pl_",i,"_global_XX")))) {
          all_delta_pl_not_zero <- all_delta_pl_not_zero + 1
        }
      }
    }

    ## Test can only be run when all requested placebos could be computed:
    if ((all_Ns_pl_not_zero == l_placebo_XX & isFALSE(normalized)) | (isTRUE(normalized) & all_Ns_pl_not_zero == l_placebo_XX & all_delta_pl_not_zero == l_placebo_XX)) {

      ## Creating a vector with all placebo estimates
      didmgt_Placebo <- matrix(0, nrow = l_placebo_XX, ncol = 1)

      ## Creating a matrix where the variances and the covariances of the placebos will be stored.
      didmgt_Var_Placebo <- matrix(0, nrow = l_placebo_XX, ncol = l_placebo_XX)

      ## Fill those matrices
      for (i in 1:l_placebo_XX) {
        didmgt_Placebo[i,1] <- get(paste0("DID_placebo_",i,"_XX"))
        didmgt_Var_Placebo[i,i] <- get(paste0("se_placebo_",i,"_XX"))^2

        if (i < l_placebo_XX) {
          for (j in (i+1):l_placebo_XX) {
            ## Create variables necessary to compute the covariances
            if (normalized == FALSE) {
              df[[paste0("U_Gg_var_pl_",i,"_",j,"_XX")]] <- df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]] +  df[[paste0("U_Gg_var_glob_pl_",j,"_XX")]]
            } else {
              df[[paste0("U_Gg_var_pl_",i,"_",j,"_XX")]] <- df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]] / get(paste0("delta_D_pl_",i,"_global_XX")) +  df[[paste0("U_Gg_var_glob_pl_",j,"_XX")]] / get(paste0("delta_D_pl_",j,"_global_XX"))
            }

            ## Estimate the covariances
            df[[paste0("U_Gg_var_pl_",i,"_",j,"_2_XX")]] <- df[[paste0("U_Gg_var_pl_",i,"_",j,"_XX")]]^2 * df$first_obs_by_gp_XX
            assign(paste0("var_sum_pl_",i,"_",j,"_XX"), sum( df[[paste0("U_Gg_var_pl_",i,"_",j,"_2_XX")]], na.rm = TRUE) / G_XX^2)
            assign(paste0("cov_pl_",i,"_",j,"_XX"), (get(paste0("var_sum_pl_",i,"_",j,"_XX")) - get(paste0("se_placebo_",i,"_XX"))^2 - get(paste0("se_placebo_",j,"_XX"))^2)/2)

            ## Store the results
            didmgt_Var_Placebo[i,j] <- get(paste0("cov_pl_",i,"_",j,"_XX"))
            didmgt_Var_Placebo[j,i] <- get(paste0("cov_pl_",i,"_",j,"_XX"))
          }
        }
      }

      ## Compute P-value for the F-test on joint nullity of all placebos
      didmgt_Var_Placebo_inv <- Ginv(didmgt_Var_Placebo)
      didmgt_chi2placebo <- t(didmgt_Placebo) %*% didmgt_Var_Placebo_inv  %*% didmgt_Placebo
      p_jointplacebo <- 1 - pchisq(didmgt_chi2placebo[1,1], df = l_placebo_XX)
    } else {
      p_jointplacebo <- NA
      ## Error message if not all of the specified placebos could be estimated
      message("Some placebos could not be estimated. Therefore, the test of joint nullity of the placebos could not be computed.")
    }
  }

  ###### Testing for effect heterogeneity
  if (!is.null(predict_het)) {
    ## Define number of effects we want to calculate
    if (length(predict_het_good) > 0) {
      if (-1 %in% het_effects) {
        het_effects <- 1:l_XX
      }
      all_effects_XX <- c(1:l_XX)[het_effects]
      if (NA %in% all_effects_XX) {
        ## error if specified effects not matching with those actually calculated
        stop("Error in predict_het second argument: please specify only numbers that are smaller or equal to the number you request in effects()")
      }

      # Preliminaries: Yg Fg1
      df$Yg_Fg_min1_XX <- ifelse(df$time_XX == df$F_g_XX - 1, df$outcome_non_diff_XX, NA)
      df[, Yg_Fg_min1_XX := mean(Yg_Fg_min1_XX, na.rm = TRUE), by = group_XX]
      df$feasible_het_XX <- !is.na(df$Yg_Fg_min1_XX)
      if (!is.null(trends_lin)) {
        df$Yg_Fg_min2_XX <- ifelse(df$time_XX == df$F_g_XX - 2, df$outcome_non_diff_XX, NA)
        df[, Yg_Fg_min2_XX := mean(Yg_Fg_min2_XX, na.rm = TRUE), by = group_XX]
        df$Yg_Fg_min2_XX <- ifelse(is.nan(df$Yg_Fg_min2_XX), NA, df$Yg_Fg_min2_XX)

        df$feasible_het_XX <- df$feasible_het_XX & !is.na(df$Yg_Fg_min2_XX)
      }
      df <- df[order(df$group_XX, df$time_XX), ]
      df[, gr_id := seq_len(.N), by = group_XX]

      lhyp <- c()
      for (v in predict_het_good) {
        lhyp <- c(lhyp, paste0(v, "=0"))
      }

      het_res <- data.frame()
      ## Loop the procedure over all requested effects for which potential heterogeneity should be predicted
      for (i in all_effects_XX) {
        # Generation of factor dummies for regression
        het_sample <- subset(df, df$F_g_XX - 1 + i <= df$T_g_XX & df$feasible_het_XX)
        het_sample <- subset(het_sample, select = c("F_g_XX", "d_sq_XX", "S_g_XX", trends_nonparam))
        het_interact <- ""
        for (v in c("F_g_XX", "d_sq_XX", "S_g_XX", trends_nonparam)) {
          if (length(levels(as.factor(het_sample[[v]]))) > 1) {
            df[[paste0(v,"_h")]] <- factor(df[[v]])
            for (l in levels(df[[paste0(v,"_h")]])) {
              df[[paste0(v,"_h",l)]] <- as.numeric(df[[v]] == l)
            }
            het_interact <- paste0(het_interact,":",v,"_h")
          }
        }
        het_interact <- substr(het_interact,2,nchar(het_interact))
        het_sample <- NULL

        # Yg,Fg-1 + l
        df[[paste0("Yg_Fg_", i, "_XX")]] <- ifelse(df$time_XX == df$F_g_XX - 1 + i, df$outcome_non_diff_XX, NA)
        df[, paste0("Yg_Fg_",i,"_XX") := mean(get(paste0("Yg_Fg_",i,"_XX")), na.rm = TRUE), by = group_XX]

        df$diff_het_XX <- df[[paste0("Yg_Fg_",i,"_XX")]] - df$Yg_Fg_min1_XX
        if (isTRUE(trends_lin)) {
          df$diff_het_XX <- df$diff_het_XX - i * (df$Yg_Fg_min1_XX - df$Yg_Fg_min2_XX)
        }

        df[[paste0("prod_het_",i,"_XX")]] <- df$S_g_het_XX * df$diff_het_XX
        df$diff_het_XX <- NULL

        # keep one observation by group to not artificially increase sample
        df[[paste0("prod_het_",i,"_XX")]] <- ifelse(df$gr_id == 1, df[[paste0("prod_het_",i,"_XX")]], NA)

        # In order to perform the test with coeftest, we need a vector of non missing regression coefficients. To avoid collinearity, we run the regression two times: the first time with the full set of regressors (F_g_XX_h#d_sq_XX_h#S_g_XX_h), then with just the non-collinear variables.
        het_reg <- paste0("prod_het_",i,"_XX ~ ")
        for (v in predict_het_good) {
          het_reg <- paste0(het_reg,v," + ")
        }
        het_reg <- paste0(het_reg, het_interact)
        het_sample <- subset(df, df$F_g_XX - 1 + i <= df$T_g_XX)
        model <- lm(as.formula(het_reg), data = het_sample, weights = het_sample$weight_XX)
        het_reg <- gsub(het_interact, "", het_reg)
        for (k in names(model$coefficients)) {
          if (!(k %in% c("(Intercept)", predict_het_good))) {
            if (!is.na(model$coefficients[[k]])) {
              het_reg <- paste0(het_reg, " + ", k)
            }
          }
        }
        model <- lm(as.formula(het_reg), data = het_sample, weights = het_sample$weight_XX)
        model_r <- matrix(coeftest(model, vcov. = vcovHC(model, type = "HC1"))[2:(length(predict_het_good)+1), 1:3], ncol = 3)
        f_stat <- linearHypothesis(model,lhyp, vcov = vcovHC(model, type = "HC1"))[["Pr(>F)"]][2]
        t_stat <- qt(0.975, df.residual(model))
        het_sample <- NULL

        ## Output Part of the predict_het option
        het_res <- rbind(het_res, data.frame(
          effect = matrix(i, nrow = length(predict_het_good)),
          covariate = predict_het_good,
          Estimate = model_r[1:nrow(model_r),1],
          SE = model_r[1:nrow(model_r),2],
          t = model_r[1:nrow(model_r),3],
          LB = model_r[1:nrow(model_r),1] - t_stat * model_r[1:nrow(model_r),2],
          UB = model_r[1:nrow(model_r),1] + t_stat * model_r[1:nrow(model_r),2],
          N = matrix(nobs(model), nrow = length(predict_het_good)),
          pF = matrix(f_stat, nrow = length(predict_het_good))
        ))
      }
      for (v in c("F_g_XX", "d_sq_XX", "S_g_XX", trends_nonparam)) {
        df[[paste0(v,"_h")]] <- NULL
      }
      het_res <- het_res[order(het_res$covariate, het_res$effect), ]
    }

    if (l_placebo_XX > 0) {
      if (-1 %in% predict_het[2]) {
        all_effects_pl_XX <- 1:l_placebo_XX
      } else {
        if (max(het_effects) > l_placebo_XX) {
          stop("You specified some numbers in predict_het that exceed the number of placebos possible to estimate! Please specify only numbers that are smaller or equal to the number of placebos you requested.")
        } else {
          all_effects_pl_XX <- het_effects
        }
      }

      for (i in all_effects_pl_XX) {
        # Generation of factor dummies for regression
        het_sample <- subset(df, df$F_g_XX - 1 + i <= df$T_g_XX & df$feasible_het_XX)
        het_sample <- subset(het_sample, select = c("F_g_XX", "d_sq_XX", "S_g_XX", trends_nonparam))
        het_interact <- ""
        for (v in c("F_g_XX", "d_sq_XX", "S_g_XX", trends_nonparam)) {
          if (length(levels(as.factor(het_sample[[v]]))) > 1) {
            df[[paste0(v,"_h")]] <- factor(df[[v]])
            for (l in levels(df[[paste0(v,"_h")]])) {
              df[[paste0(v,"_h",l)]] <- as.numeric(df[[v]] == l)
            }
            het_interact <- paste0(het_interact,":",v,"_h")
          }
        }
        het_interact <- substr(het_interact,2,nchar(het_interact))
        het_sample <- NULL

        # Yg,Fg-1 + l
        df[[paste0("Yg_Fg_pl_", i, "_XX")]] <- ifelse(df$time_XX == df$F_g_XX - 1 - i, df$outcome_non_diff_XX, NA)
        df[, paste0("Yg_Fg_pl_",i,"_XX") := mean(get(paste0("Yg_Fg_pl_",i,"_XX")), na.rm = TRUE), by = group_XX]

        df$diff_het_pl_XX <- df[[paste0("Yg_Fg_pl_",i,"_XX")]] - df$Yg_Fg_min1_XX
        if (isTRUE(trends_lin)) {
          df$diff_het_pl_XX <- df$diff_het_pl_XX - i * (df$Yg_Fg_min1_XX - df$Yg_Fg_min2_XX)
        }

        # Now we can generate
        df[[paste0("prod_het_pl_",i,"_XX")]] <- df$S_g_het_XX * df$diff_het_pl_XX
        df$diff_het_pl_XX <- NULL

        # keep one observation by group to not artificially increase sample
        df[[paste0("prod_het_pl_",i,"_XX")]] <- ifelse(df$gr_id == 1, df[[paste0("prod_het_pl_",i,"_XX")]], NA)

        # In order to perform the test with coeftest, we need a vector of non missing regression coefficients. To avoid collinearity, we run the regression two times: the first time with the full set of regressors (F_g_XX_h#d_sq_XX_h#S_g_XX_h), then with just the non-collinear variables.
        het_reg <- paste0("prod_het_pl_",i,"_XX ~ ")
        for (v in predict_het_good) {
          het_reg <- paste0(het_reg,v," + ")
        }
        het_reg <- paste0(het_reg, het_interact)
        het_sample <- subset(df, df$F_g_XX - 1 + i <= df$T_g_XX)
        model <- lm(as.formula(het_reg), data = het_sample, weights = het_sample$weight_XX)
        het_reg <- gsub(het_interact, "", het_reg)
        for (k in names(model$coefficients)) {
          if (!(k %in% c("(Intercept)", predict_het_good))) {
            if (!is.na(model$coefficients[[k]])) {
              het_reg <- paste0(het_reg, " + ", k)
            }
          }
        }
        model <- lm(as.formula(het_reg), data = het_sample, weights = het_sample$weight_XX)
        model_r <- matrix(coeftest(model, vcov. = vcovHC(model, type = "HC1"))[2:(length(predict_het_good)+1), 1:3], ncol = 3)
        f_stat <- linearHypothesis(model,lhyp, vcov = vcovHC(model, type = "HC1"))[["Pr(>F)"]][2]
        t_stat <- qt(0.975, df.residual(model))
        het_sample <- NULL

        ## Output Part of the predict_het option
        het_res <- rbind(het_res, data.frame(
          effect = matrix(-i, nrow = length(predict_het_good)),
          covariate = predict_het_good,
          Estimate = model_r[1:nrow(model_r),1],
          SE = model_r[1:nrow(model_r),2],
          t = model_r[1:nrow(model_r),3],
          LB = model_r[1:nrow(model_r),1] - t_stat * model_r[1:nrow(model_r),2],
          UB = model_r[1:nrow(model_r),1] + t_stat * model_r[1:nrow(model_r),2],
          N = matrix(nobs(model), nrow = length(predict_het_good)),
          pF = matrix(f_stat, nrow = length(predict_het_good))
        ))
      }
      for (v in c("F_g_XX", "d_sq_XX", "S_g_XX", trends_nonparam)) {
        df[[paste0(v,"_h")]] <- NULL
      }
      het_res <- het_res[order(het_res$covariate, het_res$effect), ]
    }
  }

  ###### Performing a test that all DID_\ell effects are equal (similar structure as test on placebos, not commented, except for the small differences with placebos)
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

      ## Creating a matrix of demeaned effects: null being tested = joint equality, not jointly 0
      didmgt_D <- didmgt_identity - matrix(1/l_XX, nrow = l_XX - 1, ncol = l_XX)
      didmgt_test_effects <- didmgt_D %*% didmgt_Effects
      didmgt_test_var <- didmgt_D %*% didmgt_Var_Effects %*% t(didmgt_D)
      # Enforcing symmetry
      didmgt_test_var <- (didmgt_test_var + t(didmgt_test_var)) / 2

      didmgt_chi2_equal_ef <- t(didmgt_test_effects) %*% Ginv(didmgt_test_var) %*% didmgt_test_effects
      p_equality_effects <-
        1 - pchisq(didmgt_chi2_equal_ef[1,1], df = l_XX - 1)
      assign("p_equality_effects", p_equality_effects, inherits = TRUE)
    } else {
      message("Some effects could not be estimated. Therefore, the test of equality of effects could not be computed.")
    }
  }

  ###### Storing coefficients, variances and covariances of the estimators
  l_tot_XX <- l_XX + l_placebo_XX
  didmgt_vcov <- matrix(NA, nrow = l_tot_XX, ncol = l_tot_XX)
  mat_names <-
    colnames(didmgt_vcov) <- rownames(didmgt_vcov) <- sapply(1:l_tot_XX, function(x) ifelse(x <= l_XX, paste0("Effect_",x), paste0("Placebo_",x-l_XX)))
  for (i in 1:l_XX) {
    if (isFALSE(normalized)) {
      df[[paste0("U_Gg_var_comb_",i,"_XX")]] <- ifelse(is.null(df[[paste0("U_Gg_var_glob_",i,"_XX")]]), NA, df[[paste0("U_Gg_var_glob_",i,"_XX")]])
    } else {
      df[[paste0("U_Gg_var_comb_",i,"_XX")]] <-  ifelse(is.null(df[[paste0("U_Gg_var_glob_",i,"_XX")]]), NA, df[[paste0("U_Gg_var_glob_",i,"_XX")]]/ get(paste0("delta_D_",i,"_global_XX")))
    }
  }
  if (l_placebo_XX != 0) {
    for (i in 1:l_placebo_XX) {
      if (isFALSE(normalized)) {
        df[[paste0("U_Gg_var_comb_",l_XX + i,"_XX")]] <- ifelse(is.null(df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]]), NA, df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]])
      } else {
        df[[paste0("U_Gg_var_comb_",l_XX + i,"_XX")]] <- ifelse(is.null(df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]]), NA, df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]]/get(paste0("delta_D_pl_",i,"_global_XX")))
      }
    }
  }

  for (i in 1:l_tot_XX) {
    didmgt_vcov[i,i] <- mat_res_XX[i + (i>l_XX),2]^2
    j <- 1
    while (j < i) {
      df[[paste0("U_Gg_var_comb_",i,"_",j,"_2_XX")]] <- (df[[paste0("U_Gg_var_comb_",i,"_XX")]] + df[[paste0("U_Gg_var_comb_",j,"_XX")]])^2 * df$first_obs_by_gp_XX
      var_temp <- sum(df[[paste0("U_Gg_var_comb_",i,"_",j,"_2_XX")]], na.rm = TRUE)/G_XX^2
      didmgt_vcov[i,j] <- didmgt_vcov[j,i] <- (var_temp - mat_res_XX[i + (i>l_XX),2]^2 - mat_res_XX[j + (j>l_XX),2]^2)/2
      df[[paste0("U_Gg_var_comb_",i,"_",j,"_2_XX")]] <- var_temp <- NULL
      j <- j + 1
    }
  }

  ###### Returning the results of the estimation

  ## All the results from the estimations and tests are attached to the did_multiplegt_dyn object as its "results" branch (or as the "_by_level_n$results" for n in 1:length(levels(by)) with the by option)
  ## The whole estimation dataset plus some scalars are by default stored and passed to other functions for post-estimation features.

  mat_res_XX[,1:4] <- mat_res_XX[,1:4]
  mat_res_XX[,5:8] <- mat_res_XX[,5:8]
  rownames(mat_res_XX) <- rownames
  colnames(mat_res_XX) <- c("Estimate", "SE", "LB CI", "UB CI", "N", "Switchers", "N.w", "Switchers.w", "Time")

  # Saving the results if requested
  if (!is.null(save_results)) {
    write.csv(mat_res_XX, save_results, row.names = TRUE, col.names = TRUE)
  }

  Effect_mat <- matrix(mat_res_XX[1:l_XX, 1:(ncol(mat_res_XX) -1)], ncol = ncol(mat_res_XX)-1, nrow = l_XX)
  rownames(Effect_mat) <- rownames[1:l_XX]
  colnames(Effect_mat) <- c("Estimate", "SE", "LB CI", "UB CI", "N", "Switchers", "N.w", "Switchers.w")

  ATE_mat <- matrix(mat_res_XX[l_XX + 1, 1:(ncol(mat_res_XX) -1)], ncol = ncol(mat_res_XX)-1, nrow = 1)
  rownames(ATE_mat) <- rownames[l_XX+1]
  colnames(ATE_mat) <- c("Estimate", "SE", "LB CI", "UB CI", "N", "Switchers", "N.w", "Switchers.w")

  out_names <- c("N_Effects", "N_Placebos", "Effects", "ATE", "delta_D_avg_total", "max_pl", "max_pl_gap")
  did_multiplegt_dyn <- list(
    l_XX,
    l_placebo_XX,
    Effect_mat,
    ATE_mat,
    delta_D_avg_total,
    max_pl_XX,
    max_pl_gap_XX
  )
  if (!is.null(p_jointeffects)) {
    did_multiplegt_dyn <- append(did_multiplegt_dyn, p_jointeffects)
    out_names <- c(out_names, "p_jointeffects")
  }
  if (isTRUE(effects_equal)) {
    did_multiplegt_dyn <- append(did_multiplegt_dyn, p_equality_effects)
    out_names <- c(out_names, "p_equality_effects")
  }
  if (placebo != 0) {
    Placebo_mat <- matrix(mat_res_XX[(l_XX+2):nrow(mat_res_XX), 1:(ncol(mat_res_XX) -1)], ncol = ncol(mat_res_XX) -1, nrow = l_placebo_XX)
    rownames(Placebo_mat) <- rownames[(l_XX+2):nrow(mat_res_XX)]
    colnames(Placebo_mat) <- c("Estimate", "SE", "LB CI", "UB CI", "N", "Switchers", "N.w", "Switchers.w")


    did_multiplegt_dyn <- append(did_multiplegt_dyn, list(Placebo_mat))
    out_names <- c(out_names, "Placebos")
    if (placebo > 1) {
      if (l_placebo_XX > 1) {
        did_multiplegt_dyn <- append(did_multiplegt_dyn, p_jointplacebo)
        out_names <- c(out_names, "p_jointplacebo")
      }
    }
  }
  if (!is.null(predict_het)) {
    if (length(predict_het_good) > 0) {
      did_multiplegt_dyn <- append(did_multiplegt_dyn, list(het_res))
      out_names <- c(out_names, "predict_het")
    }
  }

  # Uncomment for debugging #
  #did_multiplegt_dyn <- append(did_multiplegt_dyn, list(df))
  #out_names <- c(out_names, "debug")

  names(did_multiplegt_dyn) <- out_names

  delta <- list()
  if (isTRUE(normalized)) {
    for (i in 1:l_XX) {
      delta[[paste0("delta_D_",i,"_global_XX")]] <-
        get(paste0("delta_D_", i, "_global_XX"))
    }
  }

  coef <- list(b = mat_res_XX[-(l_XX+1), 1], vcov = didmgt_vcov)

  ret <- list(
    df,
    did_multiplegt_dyn,
    delta,
    l_XX,
    T_max_XX,
    mat_res_XX
  )
  ret_names <- c("df", "did_multiplegt_dyn", "delta", "l_XX", "T_max_XX", "mat_res_XX")
  if (placebo!= 0) {
    ret <- append(ret, l_placebo_XX)
    ret_names <- c(ret_names, "l_placebo_XX")
  }
  ret <- append(ret, list(coef))
  ret_names <- c(ret_names, "coef")

  names(ret) <- ret_names
  ret
  })
  }
