# did_multiplegt_dyn
Estimation in Difference-in-Difference (DID) designs with multiple groups and periods.

[Short description](#Short-description) | [Setup](#Setup) | [Syntax](#Syntax) | [Description](#Description)

[Options](#Options) | [Example](#Example) | [FAQ](#FAQ) | [References](#References ) | [Authors](#Authors)

## Short description

Estimation of event-study Difference-in-Difference (DID) estimators in designs with multiple groups and periods, and with a potentially non-binary treatment that may increase 
or decrease multiple times.  This is a beta version of the command. New options will be added soon, and some of the options already provided are not fully stabilized yet.

## Setup

### Stata
```s
ssc install did_multiplegt_dyn, replace
```

### R
```s
```
[![R-CMD-check](https://github.com/DiegoCiccia/did_multiplegt_dyn/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DiegoCiccia/did_multiplegt_dyn/actions/workflows/R-CMD-check.yaml)

## Syntax

### Stata

**did_multiplegt_dyn Y G T D** [if] [in] [, **effects**(#) **normalized effects_equal placebo**(#) **ci_level**(#) **controls**(*varlist*) **trends_nonparam**(*varlist*) **weight**(*varlist*) **switchers**(*string*) **same_switchers drop_larger_lower drop_if_d_miss_before_first_switch cluster**(*varname*) **by**(*varname*) **predict_het**( *varlist,numlist*) **trends_lin normalized_weightS**(*string*) **date_first_switch**( [*by_baseline_treat*],*string*) **graphoptions**(*string*) **graph_off save_results**(*path*) **save_sample**]


## Description

**did_multiplegt_dyn** estimates the effect of a treatment on an outcome, using group-(e.g. county-
    or state-) level panel data with multiple groups and periods.  It computes the DID event-study
    estimators introduced in de Chaisemartin and D'Haultfoeuille (2020a).  **did_multiplegt_dyn** can be
    used with a binary and absorbing (staggered) treatment but it can also be used with a non-binary
    treatment (discrete or continuous) that can increase or decrease multiple times, even if lagged
    treatments affect the outcome, and if the current and lagged treatments have heterogeneous
    effects, across space and/or over time.  The event-study estimators computed by the command rely
    on a no-anticipation and parallel trends assumptions.

The command can be used in sharp designs, where the treatment is assigned at the group*period
    level, and in some fuzzy designs where the treatment varies within group*time cell, see de
    Chaisemartin and D'Haultfoeuille (2020a) for further details.  The panel of groups may be
    unbalanced:  not all groups have to be observed at every period.  The data may also be at a more
    disaggregated level than the group level (e.g. individual-level wage data to measure the effect
    of a regional-level minimum-wage on individuals' wages).

For all "switchers", namely groups that experience a change of their treatment over the study
    period, let F_g denote the first time period when g's treatment changes.  The command computes
    the non-normalized event-study estimators DID_l.  DID_1 is the average, across all switchers, of
    DID estimators comparing the F_g-1 to F_g outcome evolution of g to that of groups with the same
    baseline (period-one) treatment as g but whose treatment has not changed yet at F_g.  More
    generally, DID_l is the average, across all switchers, of DID estimators comparing the F_g-1 to
    F_g-1+l outcome evolution of g to that of groups with the same baseline treatment as g but whose
    treatment has not changed yet at F_g-1+l.  Non-normalized event-study effects are average effects
    of having been exposed to a weakly higher treatment dose for l periods, where the magnitude and
    timing of the incremental treatment doses can vary across groups.  The command also computes the
    normalized event-study estimators DID^n_l, that normalize DID_l by the average total incremental
    treatment dose received by switchers from F_g-1 to F_g-1+l with respect to their baseline
    treatment.  This normalization ensures that DID^n_l estimates a weighted average of the effects
    of the current treatment and of its l-1 first lags on the outcome.  The command also computes an
    estimated average total effect per unit of treatment, where “total effect” refers to the sum of
    the effects of a treatment increment, at the time when it takes place and at later periods, see
    Section 3.3 of de Chaisemartin and D'Haultfoeuille (2020a) for further details.  Finally, the
    command also computes placebo estimators, that average DIDs comparing the outcome evolution of
    switcher g and of its control groups, from F_g-1 to F_g-1-l, namely before g's treatment changes
    for the first time.  Those placebos can be used to test the parallel trends and no-anticipation
    assumptions under which the estimators computed by **did_multiplegt_dyn** are unbiased.

**Y** is the outcome variable.

**G** is the group variable.

**T** is the time period variable.  The command assumes that the time variable is evenly spaced
    (e.g.: the panel is at the yearly level, and no year is missing for all groups).  When it is not
    (e.g.: the panel is at the yearly level, but three consecutive years are missing for all groups),
    the command can still be used, though it requires a bit of tweaking, see FAQ section below.

**D** is the treatment variable.
    

## Options   

**effects(#)** gives the number of event-study effects to be estimated.  With a balanced panel of
    groups, the maximum number of dynamic effects one can estimate can be determined as follows.  For
    each value of the baseline treatment d, start by computing the difference between the last period
    at which at least one group has had treatment d since period 1, and the first period at which a
    group with treatment d at period 1 changed its treatment.  Add one to this difference.  Then, the
    maximum number of dynamic effects is equal to the maximum of the obtained values, across all
    values of the baseline treatment.  With an unbalanced panel of groups (e.g.: counties appear or
    disappear over time if the data is a county-level panel), this method can still be used to derive
    an upper bound of the maximum number of dynamic effects one can estimate.

**normalized:** when this option is not specified, the command estimates non-normalized event-study
    effects.  Non-normalized event-study effects are average effects of having been exposed to a
    weakly higher treatment dose for l periods, where the magnitude and timing of the incremental
    treatment doses can vary across groups.  When this option is specified, the command estimates
    normalized event-study effects, that are equal to a weighted average of the effects of the
    current treatment and of its l-1 first lags on the outcome.  See Sections 3.1 and 3.2 of de
    Chaisemartin and D'Haultfoeuille (2020a) for further details.

**effects_equal:** when this option is specified and the user requests that at least two effects be
    estimated, the command performs an F-test that all effects are equal. When the **normalized** option
    is specified, this test can be useful to assess if the current and lagged treatments all have the
    same effect on the outcome or if their effects differ, see Lemma 3 of de Chaisemartin and
    D'Haultfoeuille (2020a).

**placebo(#)** gives the number of placebo estimators to be computed.  Placebos compare the outcome
    evolution of switchers and of their controls, before switchers' treatment changes for the first
    time.  Under the parallel trends and no-anticipation assumptions underlying the event-study
    estimators computed by **did_multiplegt_dyn**, the expectation of the placebos is equal to zero.
    Thus, placebos can be used to test those assumptions, by testing the null that all placebos are
    equal to zero.  If the user requests that at least two placebos be estimated, the command
    computes the p-value of a joint test of that null hypothesis.  The number of placebos requested
    can be at most equal to the number of time periods in the data minus 2, though most often only a
    smaller number of placebos can be computed.  Also, the number of placebos requested cannot be
    larger than the number of effects requested.

**controls**(varlist) gives the names of the control variables to be included in the estimation.
    Estimators with controls are similar to those without controls, except that the first-difference
    of the outcome is replaced by residuals from regressions of the first-difference of the outcome
    on the first-differences of the controls and time fixed effects.  Those regressions are estimated
    in the sample of control (g,t)s:  (g,t)s such that group g's treatment has not changed yet at t.
    Those regressions are also estimated separately for each value of the baseline treatment.
    Estimators with controls are unbiased even if groups experience differential trends, provided
    such differential trends can be fully explained by a linear model in covariates changes.  To
    control for time-invariant covariates, one needs to interact them with the time variable **T**, or
    with time fixed effects.  See Section 1.2 of the Web Appendix of de Chaisemartin and
    D'Haultfoeuille (2020a) for further details.

**trends_nonparam**(*varlist*)**:** when this option is specified, the DID estimators computed by the
    command only compare switchers to controls whose treatment has not changed yet, with the same
    baseline treatment, and with the same value of varlist.  Estimators with the
    **trends_nonparam**(*varlist*) option are unbiased even if groups experience differential trends,
    provided all groups with the same value of varlist experience parallel trends.  varlist can only
    include time-invariant variables, and the interaction of those variables has to be coarser than
    the group variable.  For instance, if one works with a county*year data set and one wants to
    allow for state-specific trends, then one should write **trends_nonparam**(*varlist*), where state is the
    state identifier.

**weight**(varlist) gives the name of a variable to be used to weight the data.  For instance, if one
    works with a district*year data set and one wants to weight the estimation by each
    district*year's population, one should write **weight**(*population*), where population is the
    population of each district*year.  If the data set is at a more disaggregated level than
    group*time, the command aggregates it at the group*time level internally, and weights the
    estimation by the number of observations in each group*time cell if the weight option is not
    specified, or by the sum of the weights of the observations in each group*time cell if the weight
    option is specified.

**switchers**(*string*)**:** one may be interested in estimating separately the treatment effect of
    switchers-in, whose average treatment after they switch is larger than their baseline treatment,
    and of switchers-out, whose average treatment after they switch is lower than their baseline
    treatment.  In that case, one should run the command first with the **switchers**(*in*) option, and
    then with the **switchers**(*out*) option.

**same_switchers:** if this option is specified and the user requests that at least two effects be
    estimated, the command will restrict the estimation of the event-study effects to switchers for
    which all effects can be estimated, to avoid compositional changes.

**dont_drop_larger_lower:** by default, the command drops all the (g,t) cells such that at t, group g
    has experienced both a strictly larger and a strictly lower treatment than its baseline
    treatment.  de Chaisemartin and D'Haultfoeuille (2020a) recommend this procedure, if you are
    interested in more details you can see their Section 3.1.  The option **dont_drop_larger_lower**
    allows to overwrite this procedure and keeps (g,t) cells such that at t, group g has experienced
    both a strictly larger and a strictly lower treatment than its baseline treatment in the
    estimation sample.

**drop_if_d_miss_before_first_switch:** This option is relevant when the treatment of some groups is
    missing at some time periods.  Then, the command imputes some of those missing treatments.  Those
    imputations are detailed in de Chaisemartin et al (2023a).  In designs where groups' treatments
    can change at most once, all those imputations are justified by the design.  In other designs,
    some of those imputations may be liberal.  **drop_if_d_miss_before_first_switch** can be used to
    overrule the potentially liberal imputations that are not innocuous for the non-normalized
    event-study estimators.  See de Chaisemartin et al (2023a) for further details.

**cluster**(varname) can be used to cluster the estimators' standard errors. Only one clustering
    variable is allowed.  A common practice in DID analysis is to cluster standard errors at the
    group level. Such clustering is implemented by default by the command.  Standard errors can be
    clustered at a more aggregated level than the group level, but they cannot be clustered at a more
    disaggregated level.

**graphoptions**(*string*)**:**  one can use the **graphoptions**(*string*) option to modify the appearance of
    the graph produced by the command.  Options requested have to follow the syntax of Stata
    **twoway_options**.  Do not use quotation marks for text passed into the arguments of **twoway_options**.
    For instance, if you want the title of your graph to be "Graph to convince skeptical referee",
    you should type **graphoptions**(*title(Graph to convince skeptical referee)*).

**graph_off:** when this option is specified, the command does not return a graph.

**by**(*varname*)**:** when this option is specified, the command estimates all the effects by the
    different levels of varname.  If varname is a binary variable for example, then the estimation is
    carried out once for the sample of groups with varname=0 and once for the sample of groups with
    varname=1. The major output of interest then is a graph which shows the treatment effect
    evolution by the different values of varname which allows to assess if there is heterogeneity by
    varname. Please note that this type of analysis is only valid if varname is constant over time
    within groups. Further, keep in mind that there will be as many treatment effect pathways as the
    number of values varname can take, so when you use **graphoptions** you need to account for that and
    specify your individual options for the corresponding number of lines and CI's in your graph.

**predict_het**(*varlist,numlist*)**:** when this option is specified, the command outputs tables with
    estimators and tests on the heterogeneity of treatment effects following Appendix Section 1.5 of
    de Chaisemartin and D'Haultfoeuille (2020a). By default, with this option you receive **l** (the
    number of effects specified in **effects(#)**) tables, each displaying the estimate of the effect of
    your variables defined in varlist on the treatment effects computed by **did_multiplegt_dyn**.  The
    p-value of a test on the null hypothesis that all heterogeneity estimates are equal to zero is
    shown below each table. If you are only interested in the impact of your varlist variables on a
    subset of the l estimated treatment effects, you can specify those inside numlist and then the
    estimators and tests on the heterogeneity of treatment effects will only be carried out for
    those.  Please note that the variables you specify in varlist should be constant over time within
    groups and that you always have to put the comma after varlist, even if you are not including a
    specific numlist.

**ci_level(#):** with this option you can change the level of the confidence intervals displayed in
    the output tables and the graphs. The default value is fixed at 95, accoring to a 95% coverage,
    and can be adjusted to any desired integer value.  The adjusted confidence intervals then
    automatically apply to the treatment effect estimates, the average effect estimate, the placebo
    estimates, and the heterogeneity estimates.

**trends_lin:** when this option is specified, the estimation of the treatment effects allows for
    group-specific linear trends.  This is achieved by implementing a verion of the DID estimator
    that compares the evolution of first-differenced outcomes between groups that experience a change
    of their treatment over the study period to those of groups with the same baseline treatment
    whose treatment has not changed yet. For a more detailed explainaition of this procedure see
    Appendix Section 1.3 of de Chaisemartin and D'Haultfoeuille (2020a).

**normalized_weights**(*string*)**:** *requires the* **normalized** *option*. Computes the weights attached to the
    effect of the kth lag (with 0<=k<=ℓ) for each of the ℓ effects. See Section 3.2 of de
    Chaisemartin and D'Haultfoeuille (2020a) for the analytical definition of the weights. The option
    allows only two arguments:  **by_k** and **by_calendar**. With the first specification, the row variable
    of the output table is k. Since k is the lag index, weights on the same row will not be referred
    to the same lag period in calendar time. For instance, with 2 effects the first row (k = 0)
    includes: (1) for ℓ=1 (first column), the weight on the Fg-1+1-0 = Fg lag, (2) for ℓ = 2 (second
    column), the weight on the Fg-1+2-0 = Fg+1 lag.  The argument **by_calendar** solves this issue
    sorting the rows by calendar time.

**design**(*[float], string*)**:**  detects the treatment paths common to at least (float*100)% of the
    switchers. This option retrieves the number of periods after the first switch from the effects()
    argument.  Results can be printed in the Stata console specifying console as the second option.
    For example, **did_multiplegt_dyn Y G T D, effects(5) design(0.5, console)** retrieves the treatment
    paths experienced by at least 50% of the switchers up to five periods after the first switch and
    prints the output in the Stata console. Alternatively, the output can be stored in a Excel file
    providing a valid file path in the string argument.

**date_first_switch**(*[by_baseline_treat],string*)**:**  detects first switch dates and how many groups
    experienced a treatment switch on each date. The reference population is the sample of switchers
    having non missing treatment status from their first switch for a number of periods equal to the
    effects() argument.  If by_baseline_treat is specified as the first argument, separate tables are
    displayed for each level of the status quo treatment.  Results can be printed in the Stata
    console specifying console as the second option. Alternatively, the output can be stored in a
    Excel file providing a valid file path in the string argument.

**save_results**(*path*)**:** if this option is specified, the command saves the estimators requested,
    their standard error, their 95% confidence interval, and the number of observations used in the
    estimation in a separate data set, at the location specified in path.

**save_sample:** if this option is specified, the command will generate a (numeric) group level
    variable _did_sample. Differently from the traditional e(sample), this variable may take on three
    non-missing values ("Never-switcher" [0], "Switcher-on" [1], "Switcher-off" [-1]) and is missing
    for all (g,t) cells that are dropped from the estimation.

**<ins>Options compatibility and interaction:</ins>** Here are some highlights that one should be aware of when combining some options in the command:

 **i.** The option **by**(*varname*) and the option **predict_het**(*varname*) are not compatible unless they
    receive different inputs (varname). In such case (i.e, two different inputs), the command carries
    out the heterogeneity prediction, according to the variable specified in **predict_het**(),
    conditional on the different values taken by the variable specified in **by**().

 **ii.** If the option **by**() is specified, and ones requests the data to be saved using the option
    **save**(), the command will save the estimation results as usual except that the names of the
    columns are indexed by the level of the variable inputed in **by**().  E.g., if the variable (let's
    call it by_var) has 4 levels:  in the saved dataset, one will have point_estimate1 (for
    by_var==1), point_estimate2 (for by_var==2) etc.  as the estimates of effects estimated
    conditional on the sample such that by_var==1, by_var==2, etc. respectively.

**iii.** Option **by**() and Option **design**(): If one requests the design to be displayed in the console,
    the command displays succesively the design for each level of the variable inputed in **by**().
    Otherwise, if one requests the design to be stored in an Excel file, the command stores each
    design in a specific sheet.  Exactly the same reasoning applies when specifying the **by**() option
    together with the **date_first_switch**() option.

**iv.** The option **normalized** should not be specified if one wants to use **predict_het**().  For more
    details see Lemma 6 of [de Chaisemartin, C, D'Haultfoeuille, X (2020a)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3731856).
    

## Example

This example is estimating the effect of banking deregulations on loans volume, using the data of Favara and Imbs (2015)

```applescript
ssc install did_multiplegt_dyn
net get did_multiplegt_dyn
use favara_imbs_did_multiplegt_dyn.dta, clear
```

Estimating eight non-normalized event-study effects and three placebo effects of banking
deregulations on loans volume:

```applescript
did_multiplegt_dyn Dl_vloans_b county year inter_bra, effects(8) placebo(3) cluster(state_n)
```

Estimating eight normalized event-study effects and three placebo effects of banking
deregulations on loans volume, restricting the estimation to switchers for which all effects can be estimated, and testing that effects are equal:

```applescript
 did_multiplegt_dyn Dl_vloans_b county year inter_bra, effects(8) cluster(state_n) normalized same_switchers effects_equal
```

    
## :question: FAQ

 > :question: *did_multiplegt_dyn does not output exactly the same results as did_multiplegt, is there something
    wrong?*

No, the two commands can sometimes output different results.  This is mostly due to different
    conventions in the way the two commands deal with missing values.  See de Chaisemartin et al
    (2023b) for further details.

> :question: *Do I have to include group and time fixed effects as controls when using the did_multiplegt_dyn
    package?*

No, you do not have to.  Group and time fixed effects are automatically controlled for.

> :question: *My group-level panel is unbalanced: some groups (e.g. counties) are not observed in every year.
    Can I still use the command?*

You can. A frequent case of unbalancedness is when some groups are not observed over the full
    duration of the panel.  For instance, your data may be a yearly county-level panel from 1990 to
    2000, where some counties appear after 1990 while some exit before 2000.  Then, the command just
    redefines group's baseline treatment as their treatment at the first period when they are
    observed.

It may also be that some groups enter and exit the data multiple times.  For instance, you
    observe a county in 1990, 1991, 1994, 1996, and 2000. Then, the command may impute some of that
    county's missing treatments.  Those imputations are detailed in de Chaisemartin et al (2023a).
    In designs where groups' treatments can change at most once, all those imputations are justified
    by the design.  In other designs, some of those imputations may be liberal.
    **drop_if_d_miss_before_first_switch** can be used to overrule the potentially liberal imputations
    that are not innocuous for the non-normalized event-study estimators.  See de Chaisemartin et al
    (2023a) for further details.

Finally, it may also be the case that the data is fully missing at one or several time periods.
    For instance, you have data for 1990, 1991, and 1993, but 1992 is missing for every group.  Then,
    it is important to fill the gap in the data, as otherwise the estimation will assume that 1991
    and 1993 are as far apart as 1990 and 1991.  There are two ways of doing so.  First, you can
    append to your data a data set identical to your 1991 data, but with the year equal to 1992, and
    the outcome missing for every observation.  This is a conservative solution, where no first
    treatment change occurring between 1991 and 1993 will be used in the estimation, which may be
    reasonable because the year in which the change occurred is effectively unknown.  Second, you can
    append to your data a data set identical to your 1993 data, with the year equal to 1992, and the
    outcome missing for every observation.  Then, treatment changes occurring between 1991 and 1993
    will be used in the estimation, assuming they all took place between 1991 and 1992.

> :question: *Related to imbalanced panels, my outcomes (and potentially the control variables) are measured
    less frequently than the treatment.  For instance, the outcome is measured every two years, but I
    know the treatment of every group in every year.  How should I proceed?*

To fix ideas, let us first assume that the outcome is measured every two years, but you know the
    treatment of every group in every year.

You should split the sample into two subsamples, and run the command twice, one time on each of
    the subsamples.  In the first estimation, you should include all group*time cells (g,t) such that
    at t, g's treatment has never changed since the start of the panel, and all (g,t)s such that g's
    treatment has changed at least once at t, and the changed occurred at a period where the outcome
    is observed.  Since the outcome is measured every two years, in that subsample the first
    event-study effect (denoted effect_1) is the effect of being exposed to a higher treatment for
    one period, the second effect (effect_2) is the effect of being exposed to a higher treatment for
    three periods, etc.  In the second estimation, you should include all group*time cells (g,t) such
    that at t, g's treatment has never changed since the start of the panel, and all (g,t)s such that
    g's treatment has changed at least once at t, and the change occurred at a period where the
    outcome is not observed.  In that subsample, the first event-study effect (denoted effect_1) is
    the effect of being exposed to a higher treatment for two periods, the second effect (effect_2)
    is the effect of being exposed to a higher treatment for four periods, etc.  You may then combine
    the two sets of estimated effects into one event-study graph, with the only caveat that the "odd"
    and "even" effects are estimated on different subsamples.  Importantly, the two estimations have
    to be run on a dataset at the same bi-yearly level as the outcome variable: the yearly level
    treatment information should only be used to select the relevant subsamples.

If the treatment is observed three times more often than the treatment, you can follow the same
    logic, splitting the sample into three subsamples and running the command three times, etc.

A short do file with a simple example where the treatment status is observed in each period while
    the outcome is only observed every second period can be found [here](https://drive.google.com/uc?export=download&id=1NBwfsFeNltU3XSOsORdthUW49LIezm1z).

> :question: *How many control variables can I include in the estimation?*

Estimators with control variables are similar to those without controls, except that the
    first-difference of the outcome is replaced by residuals from regressions of the first-difference
    of the outcome on the first-differences of the controls and time fixed effects.  Those
    regressions are estimated in the sample of control (g,t)s:  (g,t)s such that group g's treatment
    has not changed yet at period t.  Those regressions are also estimated separately for each value
    of the baseline treatment.  If the treatment takes values 0, 1, 2, 3, and 4, one regression is
    estimated for control (g,t)s with a treatment equal to 0, one regression is estimated for control
    (g,t)s with a treatment equal to 1, etc.  The number of control variables needs to be
    significantly smaller than the number of control (g,t)s in each of those regressions.  Otherwise,
    those regressions will overfit and produce noisy estimates.  If the number of observations is
    lower than the number of variables in one of those regressions, the command will run but will not
    take into account all the controls for all values of the baseline treatment.  An error message
    will let the user know that they are encountering this situation, and may thus want to reduce
    their number of control variables.

> :question: *In my application, groups' baseline treatment is a continuous variable, meaning that all groups
    have a different period-one treatment.  Therefore, Assumption 1 in de Chaisemartin and
    D'Haultfoeuille (2020a) fails. Can I still use the command?*

 Yes you can estimate non-normalized event-study effects.  Essentially, you just need to define a
    new treatment variable equal to 0 if g's treatment has never changed at t, to 1 if g's treatment
    has changed at t and g's period-t treatment is larger than its baseline treatment, and to -1 if
    g's treatment has changed at t and g's period-t treatment is lower than its baseline treatment.
    Then, you run the command with this new treatment, including interactions of period fixed effects
    and a polynomial in the baseline treatment as control variables.  For instance, if one wants to
    model the relationship between the counterfactual outcome trend and the baseline treatment as
    quadratic, and the data has 12 periods, one needs to include 22 variables as controls: the
    baseline treatment interacted with the period 2 to 12 fixed effects, and the baseline treatment
    squared interacted with the period 2 to 12 fixed effects.  See Section 1.5 of de Chaisemartin et
    al (2022) for further details. If groups' baseline treatment is not continuous but takes many
    values, pursuing this strategy may yield more precise estimators, applying to a larger number of
    switchers, than just running the command with the original treatment, at the expense of incurring
    a potential bias if the model for the counterfactual outcome trend is misspecified.

> :question: *My design is such that treatment is binary, and groups can enter the treatment, and then leave it
    once.  Can I use the command to separately estimate the effect of joining and leaving the
    treatment?*

 Yes you can. See Section 1.5 of the Web Appendix of de Chaisemartin and D'Haultfoeuille (2020a)
    for further details.

> :question: *My design has several treatments that may all have dynamic effects.  Can I use the command to
    estimate the effect of a treatment controlling for other treatments?*

 Yes, if those treatments follow binary and staggered designs.  See Section 3.2 of the Web
    Appendix of de Chaisemartin and D'Haultfoeuille (2020b) for further details.

> :question: *Can I perform triple difference-in-differences with the command?*

Yes. Suppose for instance your third difference is across men and women in the same (g,t) cell.
    Then, for each (g,t) cell, you just need to compute the difference between the average outcome of
    men and women in cell (g,t).  Then, you simply run the command with this new outcome.

> :question: *Is it possible to compute the percentage increase in the outcome of the treated relative to its
    counterfactual?*

Yes. The command already outputs the average total treatment effect, that is, the numerator in
    the percentage increase.  To compute the denominator, you need to define a new outcome variable
    \tilde{Y} = - Y * 1{t < F_g}, where F_g is the first date at which g's treatment has changed.
    Essentially, you replace the outcome by 0 after the treatment change, and by -Y before the
    treatment change.  Finally, you can just compute un-normalized event-study estimators with
    \tilde{Y} as the outcome.

## :bookmark_tabs: References

de Chaisemartin, C, D'Haultfoeuille, X (2020a).  [Difference-in-Differences Estimators of
Intertemporal Treatment Effects](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3731856).
    
de Chaisemartin, C, D'Haultfoeuille, X (2020b).  [Two-way fixed effects regressions with several
treatments](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3751060).

de Chaisemartin, C, D'Haultfoeuille, X, Pasquier, F, Vazquez-Bare, G (2022).
[Difference-in-Differences Estimators for Treatments Continuously Distributed at Every Period](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4011782).
    
de Chaisemartin, C, D'Haultfoeuille, X, Malézieux, M, Sow, D (2023a).  [Conventions used by the
did_multiplegt_dyn Stata command with missing treatments](https://drive.google.com/file/d/1NGgScujLCCS4RrwdN-PC1SnVigfBa32h/view).
    
de Chaisemartin, C, D'Haultfoeuille, X, Malézieux, M, Sow, D (2023b).  [Comparing the
did_multiplegt and did_multiplegt_dyn Stata commands](https://drive.google.com/file/d/1d23jtOT8tiHG3mpjD_cuJy3h4I2sYFEK/view).

## Authors

Clément de Chaisemartin, Economics Department, Sciences Po, France.  
Xavier D'Haultfoeuille, CREST-ENSAE, France.  
Diego Ciccia, Economics Department, Sciences Po, France.  
Felix Knau, Economics Department, Sciences Po, France.  
Mélitine Malézieux, Economics Department, Sciences Po, France.  
Doulo Sow, CREST-ENSAE, France.  

**<ins>Contact:</ins>**  
[chaisemartin.packages@gmail.com](mailto:chaisemartin.packages@gmail.com)






