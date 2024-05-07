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
#' @import dplyr
#' @importFrom matlib Ginv 
#' @importFrom plm pdata.frame make.pbalanced
#' @importFrom data.table shift setnames
#' @importFrom magrittr %>%
#' @importFrom utils write.csv
#' @importFrom rlang :=
#' @importFrom rlang .data
#' @importFrom stats pchisq qnorm sd weighted.mean as.formula df.residual lm nobs qt relevel
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
  df <- data.frame(df)
  df <- df %>% select_at(original_names)
  df <- data.table::setnames(df, old = c(outcome, group, time, treatment), new = c("outcome", "group", "time", "treatment"))

  #### Grouping together trends_nonparam variables
  if (!is.null(trends_nonparam)) {
    df$trends_nonparam_XX <- df[trends_nonparam]
  }

  #### Patching the cluster variable: by default, the command clusters at group level. If the user specifies clustering by group, the clustering option goes to NULL.
  if (!is.null(cluster)) {
    if (paste0(cluster) == paste0(group)) {
      cluster <- NULL
    }
    df$cluster_XX <- df[[cluster]]
  }

  #### Selecting the sample
  ## Dropping observations with missing group or time
  df <- df %>% filter(!is.na(.data$group) & !is.na(.data$time)) 
  ## Dropping observations with missing controls
  if (!is.null(controls)) {
    for (var in controls) {
      df <- subset(df, !is.na(df[var]))
    }
  }

  #### Further sample selection steps
  ## Dropping observations with a missing clustering variable
  if (!is.null(cluster)) {
    df <- subset(df, !is.na(df$cluster_XX))
  }
  ## Dropping groups with always missing treatment or outcomes
  df <- df %>% group_by(.data$group) %>% mutate(mean_D = mean(.data$treatment, na.rm = TRUE)) %>% ungroup()
  df <- df %>% group_by(.data$group) %>% mutate(mean_Y = mean(.data$outcome, na.rm = TRUE)) %>% ungroup()
  df <- df %>% filter(!is.na(.data$mean_Y) & !is.na(.data$mean_D))  %>% 
      dplyr::select(-.data$mean_Y, -.data$mean_D) 

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
        df <- df %>% group_by(.data$group) %>% mutate(sd_het = sd(.data[[v]], na.rm = TRUE)) %>% ungroup()
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
  } 
  else {
    df$weight_XX <- df[[weight]]
  }
  df$weight_XX <- ifelse(is.na(df$weight_XX), 0, df$weight_XX)

  ## Checking if the data has to be collapsed
  df$counter_temp <- 1
  df <- df %>% group_by(.data$group, .data$time) %>% 
      mutate(counter = sum(.data$counter_temp)) %>% ungroup()
  aggregated_data <- max(df$counter) == 1
  df <- df %>% dplyr::select(-.data$counter, -.data$counter_temp) 

  ## Collapsing the data if necessary
  if (aggregated_data != 1) {
    df$weight_XX <- ifelse(is.na(df$treatment), 0, df$weight_XX)
    if (is.null(cluster)) {
      df$cluster_XX <- 1
    }
    . <- NULL
    df_1 <- df %>%  group_by(.data$group, .data$time) %>%
    summarise_at(vars(treatment, outcome, trends_nonparam, weight, controls, predict_het_good, .data$cluster_XX), funs(weighted.mean(., w=.data$weight_XX))) %>% ungroup()
    df_2 <- df %>% group_by(.data$group, .data$time) %>%
        summarise_at(vars(.data$weight_XX), funs(sum(., na.rm = TRUE))) %>% ungroup()
    df <- merge(df_1, df_2, by = c("group", "time"))
    df_1 <- NULL; df_2 <- NULL;
    df$trends_nonparam_XX <- df[trends_nonparam]
    if (is.null(cluster)) {
      df <- df %>% dplyr::select(-.data$cluster_XX)
    }
  }

  #### Generate factorized versions of Y, G, T and D
  df$outcome_XX <- df$outcome
  df <- df %>% group_by(.data$group) %>% mutate(group_XX = cur_group_id()) %>% ungroup()
  df <- df %>% group_by(.data$time) %>% mutate(time_XX = cur_group_id()) %>% ungroup()
  df$treatment_XX <- df$treatment

  #### Declaring that the dataset is a panel
  df <- pdata.frame(df, index = c("group_XX", "time_XX")) 
  df$time_XX <- as.numeric(as.character(df$time_XX))
  df$group_XX <- as.numeric(as.character(df$group_XX))

  #### Creating variables useful to deal with imbalanced panels
  ## G's first and last date when D not missing
  df$time_d_nonmiss_XX <- ifelse(!is.na(df$treatment_XX), df$time_XX, NA)
  df <- df %>% group_by(.data$group_XX) %>% 
      mutate(min_time_d_nonmiss_XX = min(.data$time_d_nonmiss_XX, na.rm = TRUE)) %>% ungroup()
  df <- df %>% group_by(.data$group_XX) %>% 
      mutate(max_time_d_nonmiss_XX = max(.data$time_d_nonmiss_XX, na.rm = TRUE)) %>% ungroup()
  ## G's first date when Y not missing
  df$time_y_nonmiss_XX <- ifelse(!is.na(df$outcome_XX), df$time_XX, NA)
  df <- df %>% group_by(.data$group_XX) %>%
     mutate(min_time_y_nonmiss_XX = min(.data$time_y_nonmiss_XX, na.rm = TRUE)) %>% ungroup()
  ## G's first date when D missing after Y has been not missing
  df$time_d_miss_XX <- ifelse(is.na(df$treatment_XX) & df$time_XX >= df$min_time_y_nonmiss_XX,df$time_XX, NA)
  df <- df %>% group_by(.data$group_XX) %>%
       mutate(min_time_d_miss_aft_ynm_XX = min(.data$time_d_miss_XX, na.rm = TRUE))
  df <- df %>% dplyr::select(-.data$time_d_nonmiss_XX, -.data$time_y_nonmiss_XX, 
      -.data$time_d_miss_XX)

  #### Creating baseline (status quo) treatment
  #### D_{g,1} in paper, redefined to account for imbalanced panels:
  #### g's treatment at the first period where g's treatment is not missing

  df$d_sq_XX <- ifelse(df$time_XX == df$min_time_d_nonmiss_XX,df$treatment_XX,NA)
  df <- df  %>% group_by(.data$group_XX)  %>% 
      mutate(d_sq_XX = mean(.data$d_sq_XX, na.rm = TRUE)) %>% ungroup()

  #### Enforcing Design Restriction 2 in the paper.
  #### If the option dont_drop_larger_lower was not specified, 
  #### drop (g,t) cells such that at t, g has experienced both a strictly lower 
  #### and a strictly higher treatment than its baseline treatment.

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

  #### Counting number of groups
  G_XX <- max(df$group_XX, na.rm = TRUE)

  #### Ever changed treatment 
  df$ever_change_d_XX <- abs(df$diff_from_sq_XX) > 0 & !is.na(df$treatment_XX) 
  for (i in 2:T_XX) {
    df$ever_change_d_XX[shift(df$ever_change_d_XX) == 1 & df$group_XX == shift(df$group_XX) & df$time_XX == i] <- 1
  }

  #### Creating date of the first treatment change
  df <- df[order(df$group_XX, df$time_XX), ]
  df$temp_F_g_XX <- ifelse(df$ever_change_d_XX == 1 & shift(df$ever_change_d_XX) == 0,
     df$time_XX, 0)
  df <- df  %>% group_by(.data$group_XX)  %>% mutate(F_g_XX = max(.data$temp_F_g_XX, na.rm = TRUE))  %>% dplyr::select(-.data$temp_F_g_XX) %>% ungroup()

  #### If continuous option specified, generating polynomials of D_{g,1},
  #### storing D_{g,1} somewhere, and replacing it by 0.
  if (!is.null(continuous)) {
    for (pol_level in 1:degree_pol) {
      df[paste0("d_sq_",pol_level,"_XX")] <- df$d_sq_XX^pol_level
    }
    df$d_sq_XX_orig <- df$d_sq_XX
    df$d_sq_XX <- 0
  }

  ## Creating a new value with integer levels of d_sq_XX
  df <- df %>% group_by(.data$d_sq_XX) %>% mutate(d_sq_int_XX = cur_group_id()) %>% ungroup()
  df$d_sq_int_XX <- as.numeric(as.character(df$d_sq_int_XX))

  #### Dropping values of baseline treatment such that there is no variance in F_g within
  df <- joint_trends(df, "d_sq_XX", trends_nonparam)
  df <- df %>% group_by(.data$joint_trends_XX) %>% mutate(var_F_g_XX = sd(.data$F_g_XX)) %>% ungroup()
  df <- subset(df, df$var_F_g_XX > 0)
  df <- df %>% dplyr::select(-.data$var_F_g_XX)

  if (nrow(df) == 0) {
      stop("No treatment effect can be estimated.\n  This is because Design Restriction 1 in de Chaisemartin & D'Haultfoeuille (2024) is not satisfied in the data, given the options requested.\n  This may be due to the fact that groups' period-one treatment is continuous, or takes a large number of values, and you have not specified the continuous option.\n  If so, you can try to specify this option.\n  If the issue persists even with this option, this means that all groups experience their first treatment change at the same date.\n  In this situation, estimators of de Chaisemartin & D'Haultfoeuille (2024) cannot be used.")
  }
  
  #### For each value of d_sq_XX, we drop time periods such that we do not have any control with the same baseline treatment afterwards
  #### This means the panel is no longer balanced, though it is balanced within values of the baseline treatment
  df$never_change_d_XX <- 1 - df$ever_change_d_XX 
  df <- joint_trends(df, c("time_XX", "d_sq_XX"), trends_nonparam)
  df <- df %>% group_by(.data$joint_trends_XX) %>% mutate(controls_time_XX = max(.data$never_change_d_XX)) %>% ungroup()
  df <- subset(df, df$controls_time_XX > 0)

  #### Computing t_min, T_max and adjusting F_g by last period pluc one for those that never change treatment
  t_min_XX <- min(df$time_XX)
  T_max_XX <- max(df$time_XX)
  df$F_g_XX[df$F_g_XX == 0] <- T_max_XX + 1

  ######## Dealing with missing treatments: most conservative option
  #### Let FMD_g denote the first date when g's treatment is missing while y has been not missing at least once, so that we know for sure that g already exists. 
  #### If that date is before the first period when g's treatment changes, we do not know when g's treatment has changed for the first time. Then, a conservative option is to drop all of g's outcomes starting at FMD_g.

  if (drop_if_d_miss_before_first_switch == TRUE) {
    df$outcome_XX <- ifelse(
      df$min_time_d_miss_aft_ynm_XX < df$F_g_XX & df$time_XX >= df$min_time_d_miss_aft_ynm_XX, NA, df$outcome_XX)
  }

  ######## Dealing with missing treatments: most liberal option
  #### Let FD_g and LD_g respectively denote the first and last period where a group's treatment is non missing. Let FY_g denote the first period where a group's outcome is non missing.
  #### For groups that experience at least one treatment change, let LDBF_g denote the last date before F_g where g's treatment is non missing. We have FD_g<=LDBF_g<F_g<=LD_g, and we will deal with missing treatments depending on when they occur with respect to those four dates. 

  df$last_obs_D_bef_switch_t_XX <- ifelse(df$time_XX < df$F_g_XX & !is.na(df$treatment_XX), df$time_XX, NA)
  df <- df %>% group_by(.data$group_XX) %>% 
  mutate(last_obs_D_bef_switch_XX = max(.data$last_obs_D_bef_switch_t_XX, na.rm = TRUE)) %>% ungroup()

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
  df <- df %>% group_by(.data$group_XX) %>% mutate(d_F_g_XX = mean(.data$d_F_g_temp_XX, na.rm = TRUE)) %>% ungroup()
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
    df <- subset(df, !(df$F_g_XX == 2))
    df <- df[order(df$group_XX, df$time_XX), ]
    for (v in c("outcome_XX", controls)) {
      df <- df %>% group_by(.data$group_XX) %>%
        mutate(!!paste0(v,"_L") := lag(.data[[v]], n = 1, order_by = .data$group_XX)) %>% ungroup()
      df[[v]] <- df[[v]] - df[[paste0(v,"_L")]]
      df[[paste0(v,"_L")]] <- NULL
    }
    df <- subset(df, !(df$time_XX == 1))
    t_min_XX <- min(df$time_XX)
  }

  #### Balancing the panel
  df <- df %>% select(-any_of(c("joint_trends_XX", "trends_nonparam_XX")))
  df <- pdata.frame(df, index = c("group_XX", "time_XX")) 
  df <- make.pbalanced(df, balance.type = "fill")
  df$time_XX <- as.numeric(as.character(df$time_XX))
  df$group_XX <- as.numeric(as.character(df$group_XX))
  df <- df %>% group_by(.data$group_XX) %>% mutate(d_sq_XX = mean(.data$d_sq_XX, na.rm = TRUE)) %>% ungroup()

  #### Defining N_gt, the weight of each (g,t) cell
  df$N_gt_XX <- 1
  df$N_gt_XX <- ifelse(is.na(df$outcome_XX) | is.na(df$treatment_XX), 0, df$weight_XX * df$N_gt_XX)

  #### Determining last period where g still has a control group:
  #### There is still a group with same 
  #### treatment as g's in period 1 and whose treatment has not changed since 
  #### start of panel. Definition adapted from the paper, to account for 
  #### imbalanced panel.
  df <- df %>% rowwise() %>% mutate(F_g_trunc_XX = min(.data$F_g_XX, .data$trunc_control_XX, na.rm = TRUE))
  df$F_g_trunc_XX <- ifelse(is.na(df$trunc_control_XX), df$F_g_XX, df$F_g_trunc_XX)

  df <- joint_trends(df, "d_sq_XX", trends_nonparam)
  df <- df %>% group_by(.data$joint_trends_XX) %>% 
        mutate(T_g_XX = max(.data$F_g_trunc_XX, na.rm = TRUE)) %>% ungroup()
  df$T_g_XX <- df$T_g_XX - 1

  #### Defining S_g: 
  #### an indicator variable for groups whose average post switch 
  #### treatment value is larger than their initial treatment D_{g,1}. 
  #### They will be considered switchers in. If S_g==0, the group is a switcher out. 
  #### For never-switchers, S_g is undefined.
  #### Definition of S_g matches that in paper, unless dont_drop_larger_lower specified.

  df$treatment_XX_v1 <- ifelse(df$time_XX >= df$F_g_XX & df$time_XX <= df$T_g_XX, df$treatment_XX, NA)
  df <- df %>% group_by(.data$group_XX) %>% 
      mutate(avg_post_switch_treat_XX_temp = sum(.data$treatment_XX_v1, na.rm = TRUE)) %>% ungroup()
  df <- df %>% dplyr::select(-.data$treatment_XX_v1)

  df$count_time_post_switch_XX_temp <- (df$time_XX >= df$F_g_XX & df$time_XX <= df$T_g_XX & !is.na(df$treatment_XX))
 
  df <- df %>% group_by(.data$group_XX) %>%
      mutate(count_time_post_switch_XX = sum(.data$count_time_post_switch_XX_temp, na.rm = TRUE)) %>% ungroup()
  
  df$avg_post_switch_treat_XX_temp <- df$avg_post_switch_treat_XX_temp / df$count_time_post_switch_XX

  df <- df %>% group_by(.data$group_XX) %>% 
      mutate(avg_post_switch_treat_XX = mean(.data$avg_post_switch_treat_XX_temp, na.rm = TRUE)) %>%
      dplyr::select(-.data$avg_post_switch_treat_XX_temp) %>% ungroup()

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
    df <- subset(df, !(df$avg_post_switch_treat_XX == df$d_sq_XX & df$F_g_XX != df$T_g_XX + 1))
    df$S_g_XX <- as.numeric(df$avg_post_switch_treat_XX > df$d_sq_XX)
    df$S_g_XX <- ifelse(df$F_g_XX != T_max_XX + 1, df$S_g_XX, NA)
  } else {
    df <- subset(df, !(df$avg_post_switch_treat_XX == df$d_sq_XX_orig & df$F_g_XX != df$T_g_XX + 1))
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
          #df[paste0("time_fe_XX_",time_fe_XX[j],"_bt",k,"_XX")] <- (df$time_XX >= time_fe_XX[j]) * 
              #df[[paste0("d_sq_",k,"_XX")]]
          controls <- c(controls, paste0("time_fe_XX_",j,"_bt",k,"_XX"))
      }
    }
  }

  #### Creating treatment at F_g: D_{g,F_g}
  df$d_fg_XX <- ifelse(df$time_XX == df$F_g_XX, df$treatment_XX, NA)
  df <- df %>% group_by(.data$group_XX) %>% mutate(d_fg_XX = mean(.data$d_fg_XX, na.rm = TRUE)) %>% ungroup()
  df$d_fg_XX <- ifelse(is.na(df$d_fg_XX) & df$F_g_XX == T_max_XX + 1, df$d_sq_XX, df$d_fg_XX)

  #### Creating the variable L_g_XX = T_g_XX - F_g_XX so that we can compute L_u or L_a afterwards
  df$L_g_XX <- df$T_g_XX - df$F_g_XX + 1

  #### Creating the equivalent variable L_g_placebo_XX for placebos
  if (placebo > 0) {
    df <- df %>% rowwise() %>% 
        mutate(L_g_placebo_XX = min(.data$L_g_XX[.data$F_g_XX >= 3], 
          .data$F_g_XX[.data$F_g_XX >= 3] - 2)) 
    df$L_g_placebo_XX <- ifelse(df$L_g_placebo_XX == Inf, NA, df$L_g_placebo_XX)
  }

  #### Tagging first observation of each group_XX
  df <- df[order(df$group_XX, df$time_XX), ]
  df <- df %>% group_by(.data$group_XX) %>% 
      mutate(first_obs_by_gp_XX = row_number() == 1) %>% ungroup()
  df$first_obs_by_gp_XX <- as.numeric(df$first_obs_by_gp_XX)

  #### If cluster option if specified, flagging first obs in cluster and checking if the cluster variable is weakly coarser than the group one.
  if (!is.null(cluster)) {
    df <- df %>% group_by(.data$cluster_XX) %>%
        mutate(first_obs_by_clust_XX = row_number() == 1) %>% ungroup()
    df$first_obs_by_clust_XX <- as.numeric(df$first_obs_by_clust_XX)

    df <- df %>% group_by(.data$group_XX) %>%
        mutate(cluster_var_g_XX = sd(.data$cluster_XX)) %>% ungroup()
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

  ######## 3. Necessary pre-estimation steps when the controls option is specified

  if (!is.null(controls)) {
    count_controls <- 0
    df$fd_X_all_non_missing_XX <- 1
    for (var in controls) {
      count_controls <- count_controls + 1
      #### Computing the first difference of control variables 
      df[paste0("diff_X", count_controls, "_XX")] <- diff(df[[var]])
      df$fd_X_all_non_missing_XX <- ifelse(
         is.na(df[[paste0("diff_X", count_controls, "_XX")]]), 0, df$fd_X_all_non_missing_XX)
    }

    #### Residualization and computation of term entering variance
    count_controls <- 0
    mycontrols_XX <- c()
    prod_controls_y <- ""
    for (var in controls) {
      count_controls <- count_controls + 1
      df <- df %>% select(-any_of(c("sum_weights_contol_XX","avg_diff_temp_XX", "diff_y_wXX")))

      ## Computing \Delta X_{.,t}: average of controls' first difference at time t, among groups whose treatment has not changed yet.

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

      ## Computing \Delta\Dot{X}_{g,t}, the difference between the first differences of covariates and the average of their first-difference, which gives us the residuals of a regression of covariates on time fixed effects. 
      ## Multiply by sqrt(N_gt_XX) to replicate weighted regression

      df[paste0("resid_X", count_controls, "_time_FE_XX")] <- sqrt(df$N_gt_XX) * 
          (df[[paste0("diff_X", count_controls,"_XX")]] - 
            df[[paste0("avg_diff_X",count_controls, "_XX")]])
      df[[paste0("resid_X", count_controls, "_time_FE_XX")]] <- ifelse(
        is.na(df[paste0("resid_X", count_controls, "_time_FE_XX")]), 0, 
        df[[paste0("resid_X", count_controls, "_time_FE_XX")]])

      ## Storing the obtained residuals for the computation of theta_d
      mycontrols_XX <- c(mycontrols_XX, paste0("resid_X", count_controls, "_time_FE_XX"))

      ## Generating the product between \Delta\Dot{X}_{g,t} and \Delta Y_{g,t}
      ## Multiply by sqrt(N_gt_XX) to replicate weighted regression
      df$diff_y_wXX <- sqrt(df$N_gt_XX) * df$diff_y_XX

      # Note: resid_X`count_controls'_time_FE_XX is already multiplied by sqrt(N_gt_XX)
      df[paste0("prod_X",count_controls,"_Ngt_XX")] <- df[[paste0("resid_X",count_controls,"_time_FE_XX")]] * sqrt(df$N_gt_XX)
      df[[paste0("prod_X",count_controls,"_Ngt_XX")]] <- ifelse(
        is.na(df[[paste0("prod_X",count_controls,"_Ngt_XX")]]),0,
        df[[paste0("prod_X",count_controls,"_Ngt_XX")]])
      
    }

    ## Computing the Den_d matrices and their inverts,
    ## and creating locals storing the status quos for which Den_d^{-1} not defined. 

    store_singular_XX <- ""
    store_noresidualization_XX <- c()
    levels_d_sq_XX_final <- c()
    levels_d_sq_XX <- levels(as.factor(df$d_sq_int_XX))

    for (l in levels_d_sq_XX) {
        data_XX <- df
        assign(paste0("store_singular_", l, "_XX"), 0)
        assign(paste0("useful_res_", l, "_XX"), 
            length(levels(as.factor(data_XX$F_g_XX[data_XX$d_sq_int_XX == l]))))

        ## A baseline treatment is relevant iff it is taken by at least two groups with different values of F_g_XX and non-missing diff_y_XX, otherwise we do not need to perform the residualization for this specific baseline treatment.
        if (get(paste0("useful_res_",l,"_XX")) > 1) {

          # Isolate the observations used for the computation of theta_d
          data_XX <- subset(data_XX, data_XX$ever_change_d_XX == 0 & !is.na(data_XX$diff_y_XX) & data_XX$fd_X_all_non_missing_XX == 1 & data_XX$d_sq_int_XX == l)

          # R version of the accum function
          #-- The final matrix should be order k + 1 with k n. of controls
          # Using the matrix accum function, to regress the first difference of outcome on the first differences of covariates. We will obtain the vectors of coefficients \theta_d s, where d indexes values of the baseline treatment.
          Y_vec <- as.matrix(data_XX$diff_y_wXX)
          X_vec <- as.matrix(data_XX[mycontrols_XX])
          W_vec <- as.matrix(data_XX$weight_XX)
          YX_vec <- cbind(Y_vec, X_vec, matrix(1, nrow = length(Y_vec), ncol = 1))
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

            # Computing the vectors of coefficients \theta_d for each value of the baseline tratment
            assign(paste0("coefs_sq_", l, "_XX"), Ginv(as.matrix(didmgt_XX), tol = 10^(-16)) %*% didmgt_XY)
            levels_d_sq_XX_final <- c(levels_d_sq_XX_final, l)

            # Computing the matrix Denom^{-1}
            # Check if the matrix is invertible
            det_XX <- det(as.matrix(didmgt_XX))
            if (abs(det_XX) <= 10^(-16)) {
              assign(paste0("store_singular_", l, "_XX"), 1)
            } 

            assign(paste0("inv_Denom_",l,"_XX"), Ginv(as.matrix(didmgt_XX), tol = 10^(-16)) * G_XX)
            }
            
          }
        }

    # Fill up store_singular_XX, with correct values of statu quo and not the levels
    levels_d_sq_bis_XX <- levels(as.factor(df$d_sq_XX))
    index_sing_XX <- 0 
    for (l in levels_d_sq_bis_XX) {
      index_sing_XX <- index_sing_XX + 1
      if(get(paste0("store_singular_", index_sing_XX,"_XX")) == 1) {
        store_singular_XX <- paste(store_singular_XX, l)
      }
    }

    # Display errors if one of the Den_d^{-1} is not defined
    if (store_singular_XX != "") {
      message(sprintf("Some control variables are not taken into account for groups with baseline treatment equal to: %s. This may occur in the following situations:", store_singular_XX))
      message("1. For groups with those values of the baseline treatment, the regression of the outcome first difference on the controls' first differences and time fixed effects has fewer observations than variables. Note that for each value of the baseline treatment, those regressions are estimated among (g,t)s such that g has not changed treatment yet at t.")
      message("2. For groups with those values of the baseline treatment, two or more of your control variables are perfectly collinear in the sample where the regression is run, for instance because those control variables do not vary over time.")
    }

    # Values of baseline treatment such that residualization could not be performed at all are dropped.
    for (l in store_noresidualization_XX) {
      df <- subset(df, df$d_sq_int_XX != l)
    }

    # Running residualization regression to compute predicted values
    fe_reg <- "diff_y_XX ~ "
    indep_var <- c()
    for (c in 1:count_controls) {
      fe_reg <- paste(fe_reg,"+",paste0("diff_X",c,"_XX")) 
      indep_var <- c(indep_var, paste0("diff_X",c,"_XX"))
    }
    for (t in 2:T_max_XX) {
      df[[paste0("time_FE_XX", t)]] <- as.numeric(df$time_XX == t)
      indep_var <- c(indep_var, paste0("time_FE_XX", t))
    }
    fe_reg <- paste(fe_reg,"+ time_FE_XX -1") 
    for (l in levels_d_sq_XX_final) {
      df[[paste0("E_y_hat_gt_int_",l,"_XX")]] <- 0

      data_reg <- subset(df, df$d_sq_int_XX == l &  df$F_g_XX > df$time_XX & df$time_XX != t_min_XX)  
      data_reg$time_FE_XX <- as.factor(data_reg$time_XX)
      data_reg <- within(data_reg, time_FE_XX <- relevel(time_FE_XX, ref = 2))
      model <- lm(as.formula(fe_reg),  data = data_reg, weights = data_reg$weight_XX) 

      for (v in names(model$coefficients)) {
          df$to_add <- df[[v]] * model$coefficients[[v]] 
          df$to_add[is.na(df$to_add)] <- 0
          df[[paste0("E_y_hat_gt_int_",l,"_XX")]] <- df[[paste0("E_y_hat_gt_int_",l,"_XX")]] + df$to_add
          df$to_add <- NULL
      }
      df[[paste0("E_y_hat_gt_int_",l,"_XX")]] <- ifelse( 
        df$d_sq_int_XX == l &  df$F_g_XX > df$time_XX & df$time_XX != t_min_XX, df[[paste0("E_y_hat_gt_int_",l,"_XX")]], NA)
       data_reg <- NULL
    }
    for (t in 2:T_max_XX) {
      df[[paste0("time_FE_XX", t)]] <- NULL
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

  ## Generating default values for the variables which will be aggregated 
  ## after Program 2 below has been run for switchers in and for switchers out.

  inh_obj <- c()
  df[paste0("U_Gg", 1:l_XX, "_plus_XX")] <- lapply(1:l_XX, function(i) 0)
  df[paste0("U_Gg", 1:l_XX, "_minus_XX")] <- lapply(1:l_XX, function(i) 0)
  df[paste0("count", 1:l_XX, "_plus_XX")] <- lapply(1:l_XX, function(i) 0)
  df[paste0("count", 1:l_XX, "_minus_XX")] <- lapply(1:l_XX, function(i) 0)
  df[paste0("U_Gg_var_", 1:l_XX, "_in_XX")] <- lapply(1:l_XX, function(i) 0)
  df[paste0("U_Gg_var_", 1:l_XX, "_out_XX")] <- lapply(1:l_XX, function(i) 0)
  df[paste0("delta_D_g_", 1:l_XX, "_plus_XX")] <- lapply(1:l_XX, function(i) 0)
  df[paste0("delta_D_g_", 1:l_XX, "_minus_XX")] <- lapply(1:l_XX, function(i) 0)
  assign("sum_for_var_in_XX", 0)
  assign("sum_for_var_out_XX", 0)
  inh_obj <- c(inh_obj, "sum_for_var_in_XX", "sum_for_var_out_XX")
  if (placebo != 0) {
    df[paste0("U_Gg_pl_", 1:l_XX, "_plus_XX")] <- lapply(1:l_XX, function(i) 0)
    df[paste0("U_Gg_pl_", 1:l_XX, "_minus_XX")] <- lapply(1:l_XX, function(i) 0)
    df[paste0("count", 1:l_XX, "_pl_plus_XX")] <- lapply(1:l_XX, function(i) 0)
    df[paste0("count", 1:l_XX, "_pl_minus_XX")] <- lapply(1:l_XX, function(i) 0)
    df[paste0("U_Gg_var_pl_", 1:l_XX, "_in_XX")] <- lapply(1:l_XX, function(i) 0)
    df[paste0("U_Gg_var_pl_", 1:l_XX, "_out_XX")] <- lapply(1:l_XX, function(i) 0)
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
  
  # For switchers in
  if (switchers == "" | switchers == "in") {
    if (!is.na(L_u_XX) & L_u_XX != 0) {

      ## Perform the estimation of effects and placebos outside of the loop on 
      ## number of effects if trends_lin not specified
      if (isFALSE(trends_lin)) {
        data <- did_multiplegt_dyn_core(df, outcome = "outcome_XX", group = "group_XX", time = "time_XX", treatment = "treatment_XX", effects = l_XX, placebo = l_placebo_XX, switchers_core = "in", trends_nonparam = trends_nonparam, controls = controls, same_switchers = same_switchers, same_switchers_pl = same_switchers_pl, only_never_switchers = only_never_switchers, normalized = normalized, globals = globals, const = const, trends_lin = trends_lin, controls_globals = controls_globals, less_conservative_se = less_conservative_se, continuous = continuous)

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
          data <- did_multiplegt_dyn_core(df, outcome = "outcome_XX", group = "group_XX", time = "time_XX", treatment = "treatment_XX", effects = i, placebo = 0, switchers_core = "in", trends_nonparam = trends_nonparam, controls = controls, same_switchers = TRUE, same_switchers_pl = FALSE, only_never_switchers = only_never_switchers, normalized = normalized, globals = globals, const = const, trends_lin = trends_lin, controls_globals = controls_globals, less_conservative_se = less_conservative_se, continuous = continuous)

          df <- data$df
          data$df <- NULL
          for (e in names(data$const)) {
            const[[e]] <- data$const[[e]]
            assign(e, const[[e]])
          }          

	        ## Store the number of the event-study effect for switchers-in
          df$switchers_tag_XX[df[[paste0("distance_to_switch_",i,"_XX")]] == 1] <- i
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
            data <- did_multiplegt_dyn_core(df, outcome = "outcome_XX", group = "group_XX", time = "time_XX", treatment = "treatment_XX", effects = i, placebo = i, switchers_core = "in", trends_nonparam = trends_nonparam, controls = controls, same_switchers = TRUE, same_switchers_pl = TRUE, only_never_switchers = only_never_switchers, normalized = normalized, globals = globals, const = const, trends_lin = trends_lin, controls_globals = controls_globals, less_conservative_se = less_conservative_se, continuous = continuous)

            df <- data$df
            data$df <- NULL
            for (e in names(data$const)) {
              const[[e]] <- data$const[[e]]
              assign(e, const[[e]])
            }          
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

  ## Same thing as above, for switchers out
  if (switchers == "" | switchers == "out") {
    if (!is.na(L_a_XX) & L_a_XX != 0) {

      if (isFALSE(trends_lin)) {
        data <- did_multiplegt_dyn_core(df, outcome = "outcome_XX", group = "group_XX", time = "time_XX", treatment = "treatment_XX", effects = l_XX, placebo = l_placebo_XX, switchers_core = "out", trends_nonparam = trends_nonparam, controls = controls, same_switchers = same_switchers, same_switchers_pl = same_switchers_pl, only_never_switchers = only_never_switchers, normalized, globals = globals, const = const, trends_lin = trends_lin, controls_globals = controls_globals, less_conservative_se, continuous = continuous)

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
          data <- did_multiplegt_dyn_core(df, outcome = "outcome_XX", group = "group_XX", time = "time_XX", treatment = "treatment_XX", effects = i, placebo = 0, switchers_core = "out", trends_nonparam = trends_nonparam, controls = controls, same_switchers = TRUE, same_switchers_pl = FALSE, only_never_switchers = only_never_switchers, normalized = normalized, globals = globals, const = const, trends_lin = trends_lin, controls_globals = controls_globals, less_conservative_se = less_conservative_se, continuous = continuous)

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
          df[[paste0("U_Gg_var_",i,"_out_XX")]] <- df[[paste0("U_Gg",i,"_var_XX")]]
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
            data <- did_multiplegt_dyn_core(df, outcome = "outcome_XX", group = "group_XX", time = "time_XX", treatment = "treatment_XX", effects = i, placebo = i, switchers_core = "out", trends_nonparam = trends_nonparam, controls = controls, same_switchers = TRUE, same_switchers_pl = TRUE, only_never_switchers = only_never_switchers, normalized = normalized, globals = globals, const = const, trends_lin = trends_lin, controls_globals = controls_globals, less_conservative_se = less_conservative_se, continuous = continuous)

            df <- data$df
            data$df <- NULL
            for (e in names(data$const)) {
              const[[e]] <- data$const[[e]]
              assign(e, const[[e]])
            }
          }

          if (get(paste0("N0_placebo_",i,"_XX")) != 0) {
            df[[paste0("U_Gg_pl_",i,"_minus_XX")]] <- - df[[paste0("U_Gg_placebo_",i,"_XX")]]
            df[[paste0("count",i,"_pl_minus_XX")]] <- df[[paste0("count",i,"_pl_core_XX")]]
            df[[paste0("U_Gg_var_pl_",i,"_out_XX")]] <- df[[paste0("U_Gg_pl_",i,"_var_XX")]]
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
          df$U_Gg_var_minus_XX <- df$U_Gg_var_XX
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
    df[paste0("U_Gg",i,"_global_XX")] <- get(paste0("N1_",i,"_XX_new")) / (get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new"))) * df[[paste0("U_Gg",i,"_plus_XX")]] + get(paste0("N0_",i,"_XX_new")) / (get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new"))) * df[[paste0("U_Gg",i,"_minus_XX")]]
    df[[paste0("U_Gg",i,"_global_XX")]][df$first_obs_by_gp_XX == 0] <- NA

    df <- df %>% rowwise() %>% 
        mutate(!!paste0("count",i,"_global_XX") :=  max(.data[[paste0("count",i,"_plus_XX")]], .data[[paste0("count",i,"_minus_XX")]], na.rm = TRUE)) 
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
    df[paste0("N_effect_",i,"_XX")] <- sum(df[[paste0("count",i,"_global_XX")]], na.rm = TRUE)
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
    df[paste0("DID_",i,"_XX")] <- sum(df[[paste0("U_Gg",i,"_global_XX")]], na.rm = TRUE)
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
      df <- df %>% rowwise() %>% 
      mutate(count_global_XX = max(.data$count_global_XX, .data[[paste0("count",i,"_global_XX")]], na.rm = TRUE))
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
      df[paste0("U_Gg_pl_",i,"_global_XX")] <- get(paste0("N1_placebo_",i,"_XX_new")) / (get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new"))) * df[[paste0("U_Gg_pl_",i,"_plus_XX")]] + get(paste0("N0_placebo_",i,"_XX_new")) / (get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new"))) * df[[paste0("U_Gg_pl_",i,"_minus_XX")]]
      df[[paste0("U_Gg_pl_",i,"_global_XX")]][df$first_obs_by_gp_XX == 0] <- NA

      df <- df %>% rowwise() %>% 
          mutate(!!paste0("count",i,"_pl_global_XX") :=  max(.data[[paste0("count",i,"_pl_plus_XX")]], .data[[paste0("count",i,"_pl_minus_XX")]], na.rm = TRUE)) 
      df[[paste0("count",i,"_pl_global_XX")]][df[[paste0("count",i,"_pl_global_XX")]] == -Inf] <- NA
      df[[paste0("count",i,"_pl_global_dwXX")]] <- as.numeric(!is.na(df[[paste0("count",i,"_pl_global_XX")]]) & df[[paste0("count",i,"_pl_global_XX")]] > 0)

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

      assign(paste0("Placebo_",i), get(paste0("DID_placebo_",i,"_XX")))
      mat_res_XX[l_XX+1+i,1] <- get(paste0("DID_placebo_",i,"_XX")) 
      rownames <- append(rownames, paste0("Placebo_",i, strrep(" ",(12 - nchar(paste0("Placebo_",i))))))

      assign(paste0("N_switchers_placebo_",i,"_XX"), get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new")))
      assign(paste0("N_switchers_placebo_",i,"_dwXX"), get(paste0("N1_dw_placebo_",i,"_XX")) + get(paste0("N0_dw_placebo_",i,"_XX")))
      mat_res_XX[l_XX+1+i,8] <- get(paste0("N_switchers_placebo_",i,"_XX"))
      mat_res_XX[l_XX+1+i,6] <- get(paste0("N_switchers_placebo_",i,"_dwXX"))
      mat_res_XX[l_XX+1+i,9] <- -i
      assign(paste0("N_switchers_placebo_",i), get(paste0("N_switchers_placebo_",i,"_XX")))
      df[paste0("N_placebo_",i,"_XX")] <- sum(df[[paste0("count",i,"_pl_global_XX")]], na.rm = TRUE)
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

  ## Loop over the number of effects to be estimated
  for (i in 1:l_XX) {
    if ((switchers == "" & (get(paste0("N1_",i,"_XX_new")) != 0 | get(paste0("N0_",i,"_XX_new"))) != 0) | (switchers == "out" & get(paste0("N0_",i,"_XX_new")) != 0 ) | (switchers == "in" & get(paste0("N1_",i,"_XX_new")) != 0 )) {

        ## Aggregating the U_Gg_var_\ell for switchers in and out
        df[paste0("U_Gg_var_glob_",i,"_XX")] <- df[[paste0("U_Gg_var_",i,"_in_XX")]] * (get(paste0("N1_",i,"_XX_new")) / (get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new")))) + df[[paste0("U_Gg_var_",i,"_out_XX")]] * (get(paste0("N0_",i,"_XX_new")) / (get(paste0("N1_",i,"_XX_new")) + get(paste0("N0_",i,"_XX_new"))))

        ## Compute \hat{\sigma}^2_l without clustering
        if (is.null(cluster)) {
        df[paste0("U_Gg_var_glob_eff",i,"_sqrd_XX")] <- df[[paste0("U_Gg_var_glob_",i,"_XX")]]^2 * df$first_obs_by_gp_XX
        assign(paste0("sum_for_var_",i,"_XX"), sum(df[[paste0("U_Gg_var_glob_eff",i,"_sqrd_XX")]], na.rm = TRUE) / G_XX^2) 
        } else {
        ## Compute \hat{\sigma}^2_l with clustering: sum U_Gg_var_\ell within a cluster, and then take average of square. 
          df[[paste0("U_Gg_var_glob_",i,"_XX")]] <- df[[paste0("U_Gg_var_glob_",i,"_XX")]] * df$first_obs_by_gp_XX

	        ## Sum within cluster
          df <- df %>% group_by(.data$cluster_XX) %>%
              mutate(!!paste0("clust_U_Gg_var_glob_",i,"_XX") 
                  := sum(.data[[paste0("U_Gg_var_glob_",i,"_XX")]], na.rm = TRUE)) %>% ungroup()
	        ## Compute average of square
          df[paste0("clust_U_Gg_var_glob_",i,"_2_XX")] <-
              df[[paste0("clust_U_Gg_var_glob_", i, "_XX")]]^2 * df$first_obs_by_clust_XX
          assign(paste0("sum_for_var_",i,"_XX"), 
              sum(df[[paste0("clust_U_Gg_var_glob_", i,"_2_XX")]], na.rm = TRUE) / G_XX^2)
          df[[paste0("U_Gg_var_glob_",i,"_XX")]] <- 
              df[[paste0("clust_U_Gg_var_glob_",i,"_XX")]]
        }

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
          
          df[paste0("U_Gg_var_glob_pl_",i,"_XX")] <- df[[paste0("U_Gg_var_pl_",i,"_in_XX")]] * (get(paste0("N1_placebo_",i,"_XX_new")) / (get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new")))) + df[[paste0("U_Gg_var_pl_",i,"_out_XX")]] * (get(paste0("N0_placebo_",i,"_XX_new")) / (get(paste0("N1_placebo_",i,"_XX_new")) + get(paste0("N0_placebo_",i,"_XX_new"))))

          if (is.null(cluster)) {
          df[paste0("U_Gg_var_glob_pl_",i,"_2_XX")] <- df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]]^2 * df$first_obs_by_gp_XX
          assign(paste0("sum_for_var_placebo_",i,"_XX"), sum(df[[paste0("U_Gg_var_glob_pl_",i,"_2_XX")]], na.rm = TRUE) / G_XX^2) 
          } else {
            df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]] <- df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]] * df$first_obs_by_gp_XX
            df <- df %>% group_by(.data$cluster_XX) %>%
                mutate(!!paste0("clust_U_Gg_var_glob_pl_",i,"_XX") 
                    := sum(.data[[paste0("U_Gg_var_glob_pl_",i,"_XX")]], na.rm = TRUE)) %>% ungroup()
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
        df <- df %>% group_by(.data$cluster_XX) %>%
            mutate(clust_U_Gg_var_global_XX = sum(.data$U_Gg_var_global_XX, na.rm = TRUE)) %>% ungroup()
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

    # Preliminaries: Yg,Fg−1
    df$Yg_Fg_min1_XX <- ifelse(df$time_XX == df$F_g_XX - 1, df$outcome_non_diff_XX, NA)
    df <- df %>% group_by(.data$group_XX) %>% 
        mutate(Yg_Fg_min1_XX = mean(.data$Yg_Fg_min1_XX, na.rm = TRUE)) %>% ungroup()
    df$feasible_het_XX <- !is.na(df$Yg_Fg_min1_XX)
    if (!is.null(trends_lin)) {
      df$Yg_Fg_min2_XX <- ifelse(df$time_XX == df$F_g_XX - 2, df$outcome_non_diff_XX, NA)
      df <- df %>% group_by(.data$group_XX) %>% 
          mutate(Yg_Fg_min2_XX = mean(.data$Yg_Fg_min2_XX, na.rm = TRUE)) %>% ungroup()
      df$Yg_Fg_min2_XX <- ifelse(is.nan(df$Yg_Fg_min2_XX), NA, df$Yg_Fg_min2_XX)

      df$feasible_het_XX <- df$feasible_het_XX & !is.na(df$Yg_Fg_min2_XX)
    }
    df <- df[order(df$group_XX, df$time_XX), ]
    df <- df %>% group_by(.data$group_XX) %>% mutate(gr_id = row_number()) %>% ungroup()

    lhyp <- c()
    for (v in predict_het_good) {
      lhyp <- c(lhyp, paste0(v, "=0"))
    }

    het_res <- data.frame()
    ## Loop the procedure over all requested effects for which potential heterogeneity should be predicted
    for (i in all_effects_XX) {
      # Generation of factor dummies for regression
      het_sample <- subset(df, df$F_g_XX - 1 + i <= df$T_g_XX & df$feasible_het_XX)[c("F_g_XX", "d_sq_XX", "S_g_XX", trends_nonparam)]
      het_interact <- ""
      for (v in c("F_g_XX", "d_sq_XX", "S_g_XX", trends_nonparam)) {
        if (length(levels(as.factor(het_sample[[v]]))) > 1) {
          df[paste0(v,"_h")] <- factor(df[[v]])
          for (l in levels(df[[paste0(v,"_h")]])) {
            df[[paste0(v,"_h",l)]] <- as.numeric(df[[v]] == l)
          }
          het_interact <- paste0(het_interact,":",v,"_h")
        }
      }
      het_interact <- substr(het_interact,2,nchar(het_interact))
      het_sample <- NULL

      # Yg,Fg-1 + l
      df[paste0("Yg_Fg_", i, "_XX")] <- ifelse(df$time_XX == df$F_g_XX - 1 + i, df$outcome_non_diff_XX, NA)
      df <- df %>% group_by(.data$group_XX) %>% 
          mutate(!!paste0("Yg_Fg_",i,"_XX") := mean(.data[[paste0("Yg_Fg_",i,"_XX")]], na.rm = TRUE)) %>% ungroup()

      df$diff_het_XX <- df[[paste0("Yg_Fg_",i,"_XX")]] - df$Yg_Fg_min1_XX
      if (isTRUE(trends_lin)) {
        df$diff_het_XX <- df$diff_het_XX - i * (df$Yg_Fg_min1_XX - df$Yg_Fg_min2_XX)        
      }

      # Now we can generate Sg*(Yg,Fg−1+l − Yg,Fg−1)
      df[paste0("prod_het_",i,"_XX")] <- df$S_g_het_XX * df$diff_het_XX
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
    df[[paste0("U_Gg_var_comb_",i,"_XX")]] <- df[[paste0("U_Gg_var_glob_",i,"_XX")]]
  } else {
    df[[paste0("U_Gg_var_comb_",i,"_XX")]] <- df[[paste0("U_Gg_var_glob_",i,"_XX")]]/ get(paste0("delta_D_",i,"_global_XX"))
  }
}
if (l_placebo_XX != 0) {
  for (i in 1:l_placebo_XX) {
    if (isFALSE(normalized)) {
      df[[paste0("U_Gg_var_comb_",l_XX + i,"_XX")]] <- df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]]
    } else {
      df[[paste0("U_Gg_var_comb_",l_XX + i,"_XX")]] <- df[[paste0("U_Gg_var_glob_pl_",i,"_XX")]]/get(paste0("delta_D_pl_",i,"_global_XX"))
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

mat_res_XX[,1:4] <- round(mat_res_XX[,1:4],5)
mat_res_XX[,5:8] <- round(mat_res_XX[,5:8],0)
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

out_names <- c("N_Effects", "N_Placebos", "Effects", "ATE", "delta_D_avg_total")
did_multiplegt_dyn <- list(
  l_XX,
  l_placebo_XX,
  Effect_mat,
  ATE_mat,
  delta_D_avg_total
)
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
    did_multiplegt_dyn <- append(did_multiplegt_dyn, p_jointplacebo)
    out_names <- c(out_names, "p_jointplacebo")
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