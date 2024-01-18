#' Core function for did_multiplegt_dyn
#' @importFrom haven read_dta
#' @md 
#' @description Estimation of event-study Difference-in-Difference (DID) estimators in designs with multiple groups and periods,and with a potentially non-binary treatment that may increase or decrease multiple times.
#' @param df (dataframe or matrix) the starting data set.
#' @param Y (str) is the outcome variable. 
#' @param G (str) is the group variable. 
#' @param T (str) is the time period variable. The command assumes that the time variable is evenly spaced (e.g.: the panel is at the yearly level, and no year is missing for all groups). When it is not (e.g.: the panel is at the yearly level, but three consecutive years are missing for all groups), the command can still be used, though it requires a bit of tweaking, see FAQ section below.
#' @param D (str) is the treatment variable.
#' @param effects (int) gives the number of event-study effects to be estimated. With a balanced panel of groups, the maximum number of dynamic effects one can estimate can be determined as follows. For each value of the baseline treatment d, start by computing the difference between the last period at which at least one group has had treatment d since period 1, and the first period at which a group with treatment d at period 1 changed its treatment. Add one to this difference. Then, the maximum number of dynamic effects is equal to the maximum of the obtained values, across all values of the baseline treatment. With an unbalanced panel of groups (e.g.: counties appear or disappear over time if the data is a county-level panel), this method can still be used to derive an upper bound of the maximum number of dynamic effects one can estimate.
#' @param placebo (int) gives the number of placebo estimators to be computed. Placebos compare the outcome evolution of switchers and of their controls, before switchers' treatment changes for the first time. Under the parallel trends and no-anticipation assumptions underlying the event-study estimators computed by did_multiplegt_dyn, the expectation of the placebos is equal to zero. Thus, placebos can be used to test those assumptions, by testing the null that all placebos are equal to zero. If the user requests that at least two placebos be estimated, the command computes the p-value of a joint test of that null hypothesis. The number of placebos requested can be at most equal to the number of time periods in the data minus 2, though most often only a smaller number of placebos can be computed. Also, the number of placebos requested cannot be larger than the number of effects requested.
#' @param ci_level (int) with this option you can change the level of the confidence intervals displayed in the output tables and the graphs. The default value is fixed at 95, accoring to a 95% coverage, and can be adjusted to any desired integer value. The adjusted confidence intervals then automatically apply to the treatment effect estimates, the average effect estimate, the placebo estimates, and the heterogeneity estimates. 
#' @param switchers (str in "", "in", "out") one may be interested in estimating separately the treatment effect of switchers-in, whose average treatment after they switch is larger than their baseline treatment, and of switchers-out, whose average treatment after they switch is lower than their baseline treatment. In that case, one should run the command first with the switchers(in) option, and then with the switchers(out) option.
#' @param trends_nonparam (atomic str or vector of str) when this option is specified, the DID estimators computed by the command only compare switchers to controls whose treatment has not changed yet, with the same baseline treatment, and with the same value of the varlist. Estimators with the trends_nonparam() option are unbiased even if groups experience differential trends, provided all groups with the same value of the varlist experience parallel trends. The varlist can only include time-invariant variables, and the interaction of those variables has to be coarser than the group variable. For instance, if one works with a county x year data set and one wants to allow for state-specific trends, then one should write trends_nonparam(state), where state is the state identifier.
#' @param weight (str) gives the name of a variable to be used to weight the data. For instance, if one works with a district x year data set and one wants to weight the estimation by each district x year's population, one should write weight(population), where population is the population of each district x year. If the data set is at a more disaggregated level than group x time, the command aggregates it at the group x time level internally, and weights the estimation by the number of observations in each group x time cell if the weight option is not specified, or by the sum of the weights of the observations in each group x time cell if the weight option is specified.
#' @param controls (atomic str or vector of str) gives the names of the control variables to be included in the estimation. Estimators with controls are similar to those without controls, except that the first-difference of the outcome is replaced by residuals from regressions of the first-difference of the outcome on the first-differences of the controls and time fixed effects. Those regressions are estimated in the sample of control (g,t)s: (g,t)s such that group g's treatment has not changed yet at t. Those regressions are also estimated separately for each value of the baseline treatment. Estimators with controls are unbiased even if groups experience differential trends, provided such differential trends can be fully explained by a linear model in covariates changes. To control for time-invariant covariates, one needs to interact them with the time variable T, or with time fixed effects. See Section 1.2 of the Web Appendix of de Chaisemartin and D'Haultfoeuille (2020a) for further details.
#' @param dont_drop_larger_lower (logical) by default, the command drops all the (g,t) cells such that at t, group g has experienced both a strictly larger and a strictly lower treatment than its baseline treatment. de Chaisemartin and D'Haultfoeuille (2020a) recommend this procedure, if you are interested in more details you can see their Section 3.1. The option dont_drop_larger_lower allows to overwrite this procedure and keeps (g,t) cells such that at t, group g has experienced both a strictly larger and a strictly lower treatment than its baseline treatment in the estimation sample.
#' @param save_sample (logical) if this option is specified, the command will generate a (numeric) group level variable did_sample. This variable may take on three non-missing values ("Never-switcher" (0), "Switcher-on" (1), "Switcher-off" (-1)) and is missing for all (g,t) cells that are dropped from the estimation.
#' @param drop_if_d_miss_before_first_switch (logical) This option is relevant when the treatment of some groups is missing at some time periods. Then, the command imputes some of those missing treatments. Those imputations are detailed in de Chaisemartin et al (2023a). In designs where groups' treatments can change at most once, all those imputations are justified by the design. In other designs, some of those imputations may be liberal. drop_if_d_miss_before_first_switch can be used to overrule the potentially liberal imputations that are not innocuous for the non-normalized event-study estimators. See de Chaisemartin et al (2023a) for further details.
#' @param cluster (str) can be used to cluster the estimators' standard errors. Only one clustering variable is allowed. A common practice in DID analysis is to cluster standard errors at the group level. Such clustering is implemented by default by the command. Standard errors can be clustered at a more aggregated level than the group level, but they cannot be clustered at a more disaggregated level.
#' @param same_switchers (logical) if this option is specified and the user requests that at least two effects be estimated, the command will restrict the estimation of the event-study effects to switchers for which all effects can be estimated, to avoid compositional changes.
#' @param same_switchers_pl (logical) requires same_switchers. This option applies the same restriction as the same_switchers option to placebos.
#' @param effects_equal (logical) when this option is specified and the user requests that at least two effects be estimated, the command performs an F-test that all effects are equal. When the normalized option is specified, this test can be useful to assess if the current and lagged treatments all have the same effect on the outcome or if their effects differ, see Lemma 3 of de Chaisemartin and D'Haultfoeuille (2020a).
#' @param save_results (str) if this option is specified, the command saves the estimators requested, their standard error, their 95% confidence interval, and the number of observations used in the estimation in a separate data set, at the location specified in {it:path}.
#' @param normalized (logical) when this option is not specified, the command estimates non-normalized event-study effects. Non-normalized event-study effects are average effects of having been exposed to a weakly higher treatment dose for l periods, where the magnitude and timing of the incremental treatment doses can vary across groups. When this option is specified, the command estimates normalized event-study effects, that are equal to a weighted average of the effects of the current treatment and of its l-1 first lags on the outcome. See Sections 3.1 and 3.2 of de Chaisemartin and D'Haultfoeuille (2020a) for further details.
#' @param design (2 args: float, str path) detects the treatment paths common to at least (float x 100)%  of the switchers. This option retrieves the number of periods  after the first switch from the effects() argument. Results can be printed in the Stata console specifying console as the second option.  For example, design = c(0.5, console) with effects = 5 retrieves the treatment paths experienced by at least 50% of the  switchers up to five periods after the first switch and prints the output in the Stata console. Alternatively, the output can be stored in a Excel file providing a valid file path in the string argument.
#' @param date_first_switch (2 args: str in ("", "by_baseline_treat"), str path) detects first switch dates and how many groups experienced a treatment switch on each date. The reference population is the sample of switchers having non missing treatment status from their first switch for a number of periods equal to the effects() argument. If by_baseline_treat is specified as the first argument, separate tables are displayed for each level of the status quo treatment. Results can be printed in the Stata console specifying console as the second option. Alternatively, the output can be stored in a Excel file providing a valid file path in the string argument.
#' @param normalized_weights (str in ("by_k", "by_calendar")) requires normalized = TRUE. Computes the weights attached to the effect of the kth lag (with 0<= k <= l) for each of the l effects. See Section 3.2 of de Chaisemartin and D'Haultfoeuille (2020a) for the analytical definition of the weights. The option allows only two arguments: "by_k" and "by_calendar". With the first specification, the row variable of the output table is k. Since k is the lag index, weights on the same row will not be referred to the same lag period in calendar time. For instance, with 2 effects the first row (k = 0) includes: (1) for l=1 (first column), the weight on the Fg-1+1-0 = Fg lag, (2) for l = 2 (second column), the weight on the Fg-1+2-0 = Fg+1 lag. The argument "by_calendar" solves this issue sorting the rows by calendar time.    
#' @param graph_off (logical) when this option is specified, the command does not return a graph. Regardless, a ggplot object will be still generated and stored in the did_multiplegt_dyn class object.
#' @param by (str) when this option is specified, the command estimates all the effects by the different levels of the specified variable. If the variable is a binary variable for example, then the estimation is carried out once for the sample of groups with var=0 and once for the sample of groups with var=1. The major output of interest then is a graph which shows the treatment effect evolution by the different values of the variable which allows to assess if there is heterogeneity by {it:varname}. Please note that this type of analysis is only valid if the variable is constant over time within groups. 
#' @param predict_het (list with 2 args: str or vector of str, -1 or vector of positive integers)  when this option is specified, the command outputs tables with estimators and tests on the heterogeneity of treatment effects following Appendix Section 1.5 of de Chaisemartin and D'Haultfoeuille (2020a). With the second argument set to -1, with this option you receive l (the number of effects specified in effects()) tables, each displaying the estimate of the effect of your variables defined in the first argument on the treatment effects computed by did_multiplegt_dyn. The p-value of a test on the null hypothesis that all heterogeneity estimates are equal to zero is shown below each table. If you are only interested in the impact of your first-argument variables on a subset of the l estimated treatment effects, you can specify those inside an integer vector and then the estimators and tests on the heterogeneity of treatment effects will only be carried out for those.Please note that the variables you specify in the varlist should be constant over time within groups.  
#' @param trends_lin trends_lin
#' @param less_conservative_se less_conservative_se
#' @param continuous continuous
#' @section Overview:
#' did_multiplegt_dyn estimates the effect of a treatment on an outcome, using group-(e.g. county- or state-) level panel data with multiple groups and periods. It computes the DID event-study estimators introduced in de Chaisemartin and D'Haultfoeuille (2020a).did_multiplegt_dyn can be used with a binary and absorbing (staggered) treatment but it can also be used with a non-binary treatment (discrete or continuous) that can increase or decrease multiple times, even if lagged treatments affect the outcome, and if the current and lagged treatments have heterogeneous effects, across space and/or over time. The event-study estimators computed by the command rely on a no-anticipation and parallel trends assumptions.
#' 
#' The command can be used in sharp designs, where the treatment is assigned at the group x period level, and in some fuzzy designs where the treatment varies within group x time cell, see de Chaisemartin and D'Haultfoeuille (2020a) for further details. The panel of groups may be unbalanced: not all groups have to be observed at every period. The data may also be at a more disaggregated level than the group level (e.g. individual-level wage data to measure the effect of a regional-level minimum-wage on individuals' wages).  
#' 
#' For all "switchers", namely groups that experience a change of their treatment over the study period, let F_g denote the first time period when g's treatment changes. The command computes the non-normalized event-study estimators DID_l. DID_1 is the average, across all switchers, of DID estimators comparing the F_g-1 to F_g outcome evolution of g to that of groups with the same baseline (period-one) treatment as g but whose treatment has not changed yet at F_g. More generally, DID_l is the average, across all switchers, of DID estimators comparing the F_g-1 to F_g-1+l outcome evolution of g to that of groups with the same baseline treatment as g but whose treatment has not changed yet at F_g-1+l. Non-normalized event-study effects are average effects of having been exposed to a weakly higher treatment dose for l periods, where the magnitude and timing of the incremental treatment doses can vary across groups. The command also computes the normalized event-study estimators DID^n_l, that normalize DID_l by the average total incremental treatment dose received by switchers from F_g-1 to F_g-1+l with respect to their baseline treatment. This normalization ensures that DID^n_l estimates a weighted average of the effects of the current treatment and of its l-1 first lags on the outcome. The command also computes an estimated average total effect per unit of treatment, where “total effect” refers to the sum of the effects of a treatment increment, at the time when it takes place and at later periods, see Section 3.3 of de Chaisemartin and D'Haultfoeuille (2020a) for further details.Finally, the command also computes placebo estimators, that average DIDs comparing the outcome evolution of switcher gand of its control groups, from F_g-1 to F_g-1-l, namely before g's treatment changes for the first time. Those placebos can be used to test the parallel trends and no-anticipation assumptions under which the estimators computed by {cmd:did_multiplegt_dyn} are unbiased.
#' @section Option compatibility:
#' Here are some highlights that one should be aware of when combining some options in the command:
#' 1. The option _by_ and the option _predict_het_ are not compatible unless they receive different inputs (varname). In such case (i.e, two different inputs), the command carries out the heterogeneity prediction,according to the variable specified in _predict_het_, conditional on the different values taken by the variable specified in _by_.
#' 2. If the option _by_ is specified, and ones requests the data to be saved using the option _save_results_, the command will save the estimation results as usual except that the names of the columns are indexed by the level of the variable inputed in _by_. E.g., if the variable (let's call it by_var) has 4 levels: in the saved dataset, one will have point_estimate1 (for by_var==1), point_estimate2 (for by_var==2) etc. as the estimates of effects estimated conditional on the sample such that by_var==1, by_var==2, etc. respectively.
#' 3. Options _by_ and _design_: If one requests the design to be displayed in the console, the command displays succesively the design for each level of the variable inputed in _by_. Otherwise, if one requests the design to be stored in an Excel file, the command stores each design in a specific sheet. Exactly the same reasoning applies when specifying the _by_ option together with the _date_first_switch_ option. 
#' 4. The option _normalized_ should not be specified if one wants to use _predict_het_. For more details see Lemma 6 of [de Chaisemartin, C, D'Haultfoeuille, X (2020a)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3731856).
#' @export 
did_multiplegt_dyn <- function(
    df, 
    Y, 
    G, 
    T, 
    D, 
    effects = 1, 
    placebo = 0, 
    ci_level = 95, 
    switchers = "", 
    trends_nonparam = NULL, 
    weight = NULL, 
    controls = NULL, 
    dont_drop_larger_lower = FALSE, 
    save_sample = FALSE, 
    drop_if_d_miss_before_first_switch = FALSE, 
    cluster = NULL, 
    same_switchers = FALSE, 
    same_switchers_pl = FALSE, 
    effects_equal = FALSE, 
    save_results = NULL, 
    normalized = FALSE, 
    design = NULL, 
    date_first_switch = NULL, 
    normalized_weights = NULL, 
    graph_off = FALSE,
    by = NULL,
    predict_het = NULL,
    trends_lin = FALSE,
    less_conservative_se = FALSE,
    continuous = NULL
    ) { 
  
  args <- list()
  for (v in names(formals(did_multiplegt_dyn))) {
    if (v != "df") {
      args[[v]] <- get(v)
    }
  }
  did_multiplegt_dyn <- list(args)
  f_names <- c("args")

  by_levels <- c("_no_by")
  if (!is.null(by)) {
    if (!did_multiplegt_dyn_by_check(df, G, by)) {
      cat(sprintf("The variable %s is time-variant. The command will ignore the by option", by));cat("\n");
    } else {
      by_levels <- levels(factor(df[[by]]))
      did_multiplegt_dyn <- append(did_multiplegt_dyn, list(by_levels))
      f_names <- c(f_names, "by_levels")
    }
  }

  append_design <- FALSE
  append_dfs <- FALSE
  for (b in 1:length(by_levels)) {

    if (by_levels[b] != "_no_by") {
        df_main <- subset(df, df[[by]] == by_levels[b])
    } else {
      df_main <- df
    }

    df_est <- did_multiplegt_main(df_main, Y, G, T, D, 
    effects, placebo, ci_level, switchers, trends_nonparam, 
    weight, controls, dont_drop_larger_lower, 
    drop_if_d_miss_before_first_switch, cluster, 
    same_switchers, same_switchers_pl, effects_equal, 
    save_results, normalized, predict_het, trends_lin, 
    less_conservative_se, continuous)

    temp_obj <- list(df_est$did_multiplegt_dyn)
    names(temp_obj)[length(temp_obj)] <- "results"

    if (!is.null(design)) {
      temp_obj <- append(temp_obj, list(did_multiplegt_dyn_design(df_est, design, weight, by, by_levels[b], append_design)))
      names(temp_obj)[length(temp_obj)] <- "design"
      append_design <- TRUE
    }

    if (!is.null(date_first_switch)) {
      temp_obj <- append(temp_obj, list(did_multiplegt_dyn_dfs(df_est, date_first_switch, by, by_levels[b], append_dfs)))
      names(temp_obj)[length(temp_obj)] <- "date_first_switch"
      append_dfs <- TRUE
    }

    if (!is.null(normalized_weights)) {
      temp_obj <- append(temp_obj, 
          list(did_multiplegt_dyn_normweights(df_est, normalized, normalized_weights, same_switchers, continuous)))
      names(temp_obj)[length(temp_obj)] <- "normalized_weights"
    }

    temp_obj <- append(temp_obj, list(did_multiplegt_dyn_graph(df_est)))
    names(temp_obj)[length(temp_obj)] <- "plot"

    if (isTRUE(save_sample)) {
      df_save_XX <- did_save_sample(df_est, G, T)
      df_m <- merge(df, df_save_XX, by = c(G, T))
      temp_obj <- append(temp_obj, list(df_m))
      names(temp_obj)[length(temp_obj)] <- "save_sample"
    }

    if (by_levels[b] == "_no_by") {
      did_multiplegt_dyn <- c(did_multiplegt_dyn, temp_obj)
      f_names <- c(f_names, names(temp_obj))
    } else {
      temp_obj <- append(temp_obj, list(by_levels[b]))
      names(temp_obj)[length(temp_obj)] <- "level"
      did_multiplegt_dyn <- append(did_multiplegt_dyn, list(temp_obj))
      f_names <- c(f_names, paste0("by_level_",b))
    }
  }

  names(did_multiplegt_dyn) <- f_names

  if (!is.null(by)) {
    if (isTRUE(save_sample)) {
      did_multiplegt_dyn <- adj_save_sample(did_multiplegt_dyn)
      names(did_multiplegt_dyn)[length(did_multiplegt_dyn)] <- "save_sample"
    }
    did_multiplegt_dyn <- append(did_multiplegt_dyn, list(combine_plot(did_multiplegt_dyn)))
    names(did_multiplegt_dyn)[length(did_multiplegt_dyn)] <- "plot"
  }

  if (isFALSE(graph_off)) {
    print(did_multiplegt_dyn$plot)
  }
  
  class(did_multiplegt_dyn) <- c(class(did_multiplegt_dyn), "did_multiplegt_dyn")
  return(did_multiplegt_dyn)
}
