{smcl}
{* *! version 1  2023-06-29}{...}
{viewerjumpto "Syntax" "did_multiplegt_dyn##syntax"}{...}
{viewerjumpto "Description" "did_multiplegt_dyn##description"}{...}
{viewerjumpto "Options" "did_multiplegt_dyn##options"}{...}
{viewerjumpto "Examples" "did_multiplegt_dyn##examples"}{...}
{viewerjumpto "Saved results" "did_multiplegt_dyn##saved_results"}{...}

{title:Title}

{p 4 4}
{cmd:did_multiplegt_dyn} {hline 2} Estimation of event-study Difference-in-Difference (DID) estimators
in designs with multiple groups and periods,
and with a potentially non-binary treatment that may increase or decrease multiple times.
This is a beta version of the command. New options will be added soon, and some of the options
already provided are not fully stabilized yet.
{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 4}
{cmd:did_multiplegt_dyn Y G T D} {ifin}
[{cmd:,}
{cmd:effects(}{it:#}{cmd:)}
{cmd:{ul:norm}alized}
{cmd:effects_equal}
{cmd:placebo(}{it:#}{cmd:)}
{cmd:controls(}{it:varlist}{cmd:)}
{cmd:trends_nonparam(}{it:varlist}{cmd:)}
{cmd:weight(}{it:varlist}{cmd:)}
{cmd:switchers(}{it:string}{cmd:)}
{cmd:same_switchers}
{cmd:drop_larger_lower}
{cmd:drop_if_d_miss_before_first_switch}
{cmd:cluster(}{it:varname}{cmd:)}
{cmd:graphoptions(}{it:string}{cmd:)}
{cmd:graph_off}
{cmd:{ul:sav}e_results(}{it:path}{cmd:)}]
{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 4}
{cmd:did_multiplegt_dyn} estimates the effect of a treatment on an outcome,
using group-(e.g. county- or
state-) level panel data with multiple groups and periods.
It computes the DID event-study estimators introduced
in de Chaisemartin and D'Haultfoeuille (2020a).
{cmd:did_multiplegt_dyn} can be used with a binary
and absorbing (staggered) treatment
but it can also be used with a non-binary treatment
that can increase or decrease multiple times,
even if lagged treatments affect the outcome,
and if the current and lagged treatments have heterogeneous
effects, across space and/or over time.
The event-study estimators computed by the command rely on a no-anticipation
and parallel trends assumptions.
{p_end}

{p 4 4}
The command can be used in sharp designs,
where the treatment is assigned at the group*period level,
and in some fuzzy designs where the treatment varies within group*time cell,
see de Chaisemartin and D'Haultfoeuille (2020a) for further details.
The panel of groups may be unbalanced:
not all groups have to be observed at every period.
The data may also be at a more disaggregated level than
the group level (e.g. individual-level wage data
to measure the effect of a regional-level minimum-wage on individuals' wages).
{p_end}

{p 4 4}
For all "switchers", namely groups that experience a change of their treatment
over the study period, let F_g denote the first time period when g's treatment changes.
The command computes the non-normalized event-study estimators DID_l.
DID_1 is the average, across all switchers,
of DID estimators comparing the F_g-1 to F_g outcome evolution of g to that of groups
with the same baseline (period-one) treatment
as g but whose treatment has not changed yet at F_g.
More generally,
DID_l is the average, across all switchers,
of DID estimators comparing the F_g-1 to F_g-1+l
outcome evolution of g to that of groups
with the same baseline treatment as g
but whose treatment has not changed yet at F_g-1+l.
Non-normalized event-study effects are average effects of having been
exposed to a weakly higher treatment dose for l periods,
where the magnitude and timing of
the incremental treatment doses can vary across groups.
The command also computes the normalized event-study estimators DID^n_l,
that normalize DID_l by the average total incremental treatment dose
received by switchers from F_g-1 to F_g-1+l
with respect to their baseline treatment.
This normalization ensures that DID^n_l estimates
a weighted average of the effects of the current treatment and
of its l-1 first lags on the outcome.
The command also computes an estimated average total
effect per unit of treatment, where “total effect”
refers to the sum of the effects of a treatment
increment, at the time when it takes place and at later periods,
see Section 3.3 of de Chaisemartin and D'Haultfoeuille (2020a) for further details.
Finally, the command also computes placebo estimators,
that average DIDs comparing the outcome evolution of switcher g
and of its control groups, from F_g-1 to F_g-1-l,
namely before g's treatment changes for the first time.
Those placebos can be used to test
the parallel trends and no-anticipation assumptions
under which the estimators computed by {cmd:did_multiplegt_dyn} are unbiased.
{p_end}

{p 4 4}
{cmd:Y} is the outcome variable.
{p_end}

{p 4 4}
{cmd:G} is the group variable.
{p_end}

{p 4 4}
{cmd:T} is the time period variable.
The command assumes that the time variable is evenly spaced
(e.g.: the panel is at the yearly level,
and no year is missing for all groups).
When it is not (e.g.: the panel is at the yearly level,
but three consecutive years are missing for all
groups), the command can still be used,
though it requires a bit of tweaking, see FAQ section below.
{p_end}

{p 4 4}
{cmd:D} is the treatment variable.
{p_end}

{marker options}{...}
{title:Options}

{p 4 4}
{cmd:effects(}{it:#}{cmd:)} gives the number of event-study effects to be estimated.
With a balanced panel of groups,
the maximum number of dynamic effects one can estimate
can be determined as follows.
For each value of the baseline treatment d,
start by computing the difference between the last period
at which at least one group has had treatment d since period 1,
and the first period at which a group with treatment d at period 1
changed its treatment.
Add one to this difference.
Then, the maximum number of dynamic effects is equal to
the maximum of the obtained values,
across all values of the baseline treatment.
With an unbalanced panel of groups
(e.g.: counties appear or disappear over time
if the data is a county-level panel),
this method can still be used to derive an upper bound
of the maximum number of dynamic effects one can estimate.
{p_end}

{p 4 4}
{cmd:normalized}: when this option is not specified,
the command estimates non-normalized event-study effects.
Non-normalized event-study effects are average effects of having been
exposed to a weakly higher treatment dose for l periods,
where the magnitude and timing of the incremental treatment
doses can vary across groups.
When this option is specified,
the command estimates normalized event-study effects,
that are equal to a weighted average of the effects
of the current treatment and of its l-1 first lags on the outcome.
See Sections 3.1 and 3.2 of de Chaisemartin and D'Haultfoeuille (2020a)
for further details.
{p_end}

{p 4 4}
{cmd:effects_equal}: when this option is specified and the user requests
that at least two effects be estimated,
the command performs an F-test that all effects are
equal. When the {cmd:normalized} option is specified,
this test can be useful to assess if the current and
lagged treatments all have the same effect on
the outcome or if their effects differ,
see Lemma 3 of de Chaisemartin and D'Haultfoeuille (2020a).
{p_end}

{p 4 4}
{cmd:placebo(}{it:#}{cmd:)} gives the number of placebo estimators to be computed.
Placebos compare the outcome evolution of switchers and of their controls,
before switchers' treatment changes for the first time.
Under the parallel trends and no-anticipation assumptions
underlying the event-study estimators computed by {cmd:did_multiplegt_dyn}, 
the expectation of the placebos is equal to zero.
Thus, placebos can be used to test those assumptions, by testing the null that
all placebos are equal to zero.
If the user requests that at least two placebos be estimated,
the command computes the p-value of a joint test of that null hypothesis.
The number of placebos requested can be at most
equal to the number of time periods in the data minus 2,
though most often only a smaller number of placebos can be computed.
Also, the number of placebos requested cannot be larger
than the number of effects requested.
{p_end}

{p 4 4}
{cmd:controls(}{it:varlist}{cmd:)} gives the names of
the control variables to be included in the estimation.
Estimators with controls are similar
to those without controls,
except that the first-difference of the outcome is
replaced by residuals from regressions
of the first-difference of the outcome
on the first-differences of the controls and time fixed effects.
Those regressions are estimated in the sample of control (g,t)s:
(g,t)s such that group g's treatment has not changed yet at t.
Those regressions are also estimated separately
for each value of the baseline treatment.
Estimators with controls are unbiased
even if groups experience differential trends,
provided such differential trends can be
fully explained by a linear model in covariates changes.
To control for time-invariant covariates,
one needs to interact them with the time variable {cmd:T},
or with time fixed effects.
See Section 1.2 of the Web Appendix
of de Chaisemartin and D'Haultfoeuille (2020a) for further details.
{p_end}

{p 4 4}
{cmd:trends_nonparam(}{it:varlist}{cmd:)}: when this option is specified,
the DID estimators computed by the command only compare switchers
to controls whose treatment has not changed yet,
with the same baseline treatment, and with the same value of {it:varlist}.
Estimators with the {cmd:trends_nonparam(}{it:varlist}{cmd:)}
option are unbiased even if groups experience differential trends,
provided all groups with the same value of {it:varlist}
experience parallel trends.
{it:varlist} can only include time-invariant variables,
and the interaction of those variables
has to be coarser than the group variable.
For instance, if one works with a county*year data set
and one wants to allow for state-specific trends,
then one should write {cmd:trends_nonparam(}state{cmd:)},
where state is the state identifier.
{p_end}

{p 4 4}
{cmd:weight(}{it:varlist}{cmd:)} gives the name of
a variable to be used to weight the data.
For instance,
if one works with a district*year data set
and one wants to weight the estimation
by each district*year's population,
one should write {cmd:weight(}population{cmd:)},
where population is the population of each district*year.
If the data set is at a more disaggregated level than group*time,
the command aggregates it at the group*time
level internally, and weights the estimation
by the number of observations in each group*time cell
if the weight option is not specified,
or by the sum of the weights of the observations
in each group*time cell if the weight option is specified.
{p_end}

{p 4 4}
{cmd:switchers(}{it:string}{cmd:)}: one may be interested in
estimating separately the treatment effect
of switchers-in, whose average treatment after
they switch is larger than their baseline treatment,
and of switchers-out, whose average treatment after
they switch is lower than their baseline treatment.
In that case, one should run the command first with
the {cmd:switchers(}{it:in}{cmd:)} option,
and then with the {cmd:switchers(}{it:out}{cmd:)} option.
{p_end}

{p 4 4}
{cmd:same_switchers}: if this option is specified
and the user requests that at least two effects be estimated,
the command will restrict the estimation
of the event-study effects to switchers
for which all effects can be estimated,
to avoid compositional changes.
{p_end}

{p 4 4}
{cmd:drop_larger_lower}: This option is relevant when the
treatment is non-binary.
It drops all the (g,t) cells such that at t,
group g has experienced both a strictly larger and a strictly
lower treatment than its baseline treatment.
Then, the event-study effects of those (g,t) cells can be written as weighted sums
of the effects of different treatment lags,
where the effects of some lags are weighted negatively.
de Chaisemartin and D'Haultfoeuille (2020a)
recommend dropping those (g,t) cells,
see their Section 3.1 for further details.
{p_end}

{p 4 4}
{cmd:drop_if_d_miss_before_first_switch}: This option is relevant
when the treatment of some groups is missing at some time periods.
Then,
the command imputes some of those missing treatments.
Those imputations are detailed in de Chaisemartin et al (2023a).
In designs where
groups' treatments can change at most once,
all those imputations are justified by the design.
In other designs, some of those
imputations may be liberal.
{cmd:drop_if_d_miss_before_first_switch} can be used to overrule
the potentially liberal imputations that are not innocuous
for the non-normalized event-study estimators.
See de Chaisemartin et al (2023a) for further details.
{p_end}

{p 4 4}
{cmd:cluster(}{it:varname}{cmd:)} can be used to cluster the estimators' standard
errors. Only one clustering variable is allowed.
A common practice in DID analysis is to cluster standard
errors at the group level. Such clustering is implemented by default by the command.
Standard errors can be clustered at a more aggregated level than the group level,
but they cannot be clustered at a more disaggregated level.
{p_end}

{p 4 4}
{cmd:graphoptions(}{it:string}{cmd:)}:
one can use the {cmd:graphoptions(}{it:string}{cmd:)}
option to modify the appearance of the graph produced by the command.
Options requested have to follow the syntax of Stata {cmd:twoway_options}.
Do not use quotation marks for text passed into the arguments of {cmd:twoway_options}.
For instance, if you want the title of your graph to be "Graph to convince
skeptical referee", you should type
{cmd:graphoptions(}title(Graph to convince skeptical referee){cmd:)}.
{p_end}

{p 4 4}
{cmd:graph_off}: when this option is specified,
the command does not return a graph.
{p_end}

{p 4 4}
{cmd:save_results(}{it:path}{cmd:)}: if this option is specified,
the command saves the estimators requested,
their standard error,
their 95% confidence interval,
and the number of observations used in the estimation in a separate data set,
at the location specified in {it:path}.
{p_end}

{marker Example}{...}
{title:Example: estimating the effect of banking deregulations on loans volume, using the data of Favara and Imbs (2015)}

{p 4 4}
The data for this example can be downloaded by running:
{p_end}

{p 4 4}
ssc install did_multiplegt_dyn
{p_end}

{p 4 4}
net get did_multiplegt_dyn
{p_end}

{p 4 4}
use favara_imbs_did_multiplegt_dyn.dta, clear
{p_end}

{p 4 4}
Estimating eight non-normalized event-study effects
and three placebo effects of banking deregulations on loans volume:
{p_end}

{p 4 4}
did_multiplegt_dyn Dl_vloans_b county year inter_bra, effects(8) placebo(3) cluster(state_n)
{p_end}

{p 4 4}
Estimating eight normalized event-study effects
and three placebo effects of banking deregulations on loans volume,
restricting the estimation to switchers
for which all effects can be estimated, and testing that effects are equal:
{p_end}

{p 4 4}
did_multiplegt_dyn Dl_vloans_b county year inter_bra,
effects(8) cluster(state_n) normalized same_switchers effects_equal
{p_end}

{marker FAQ}{...}
{title:FAQ}

{p 4 4}
{it:did_multiplegt_dyn does not output exactly the same results as did_multiplegt,}
{it:is there something wrong?}
{p_end}

{p 4 4}
No, the two commands can sometimes output different results.
This is mostly due to different conventions
in the way the two commands deal with missing values.
See de Chaisemartin et al (2023b) for further details.
{p_end}

{p 4 4}
{it: Do I have to include group and time fixed effects as controls when using the did_multiplegt_dyn package?}
{p_end}

{p 4 4}
No, you do not have to.
Group and time fixed effects are automatically controlled for.
{p_end}

{p 4 4}
{it:My group-level panel is unbalanced: some groups (e.g. counties)}
{it:are not observed in every year. Can I still use the command?}
{p_end}

{p 4 4}
You can. A frequent case of unbalancedness
is when some groups are not observed over the full duration of the panel.
For instance, your data may be a yearly county-level panel from 1990 to 2000,
where some counties appear after 1990 while some exit before 2000.
Then, the command just redefines group's baseline treatment
as their treatment at the first period when they are observed.
{p_end}

{p 4 4}
It may also be that some groups enter and exit the data multiple times.
For instance, you observe a
county in 1990, 1991, 1994, 1996, and 2000. Then,
the command may impute some of that county's missing treatments.
Those imputations are detailed in de Chaisemartin et al (2023a).
In designs where
groups' treatments can change at most once,
all those imputations are justified by the design.
In other designs, some of those
imputations may be liberal.
{cmd:drop_if_d_miss_before_first_switch} can be used to overrule
the potentially liberal imputations
that are not innocuous for the non-normalized event-study estimators. 
See de Chaisemartin et al (2023a) for further details.
{p_end}

{p 4 4}
Finally, it may also be the case that the data
is fully missing at one or several time periods.
For instance, you have data for 1990, 1991,
and 1993, but 1992 is missing for every group.
Then, it is important to fill the gap in the data,
as otherwise the estimation will assume that 1991 and 1993 are as far apart as 1990 and 1991.
There are two ways of doing so.
First, you can append to your data a data set identical to your 1991 data,
but with the year equal to 1992,
and the outcome missing for every observation.
This is a conservative solution,
where no first treatment change occurring between 1991
and 1993 will be used in the estimation,
which may be reasonable because the year in which the change occurred is effectively unknown.
Second, you can append to your data a data set identical to your 1993 data,
with the year equal to 1992,
and the outcome missing for every observation.
Then, treatment changes occurring between 1991
and 1993 will be used in the estimation,
assuming they all took place between 1991 and 1992.
{p_end}

{p 4 4}
{it: Related to imbalanced panels,}
{it:my outcomes (and potentially the control variables) are measured}
{it:less frequently than the treatment.}
{it:For instance, the outcome is measured every two years,}
{it:but I know the treatment of every group in every year.}
{it:How should I proceed?}
{p_end}

{p 4 4}
To fix ideas,
let us first assume
that the outcome is measured every two years,
but you know the treatment of every group in every year.
{p_end}

{p 4 4}
You should split the sample into two subsamples,
and run the command twice,
one time on each of the subsamples.
In the first estimation,
you should include all group*time cells (g,t)
such that at t, g's treatment has never changed
since the start of the panel, and all (g,t)s such that g's
treatment has changed at
least once at t, and the changed occurred at a period
where the outcome is observed.
Since the outcome is measured every two years,
in that subsample the first event-study effect (denoted effect_1)
is the effect of being exposed to a higher treatment for one period,
the second effect (effect_2)
is the effect of being exposed to a higher treatment for three periods, etc.
In the second estimation,
you should include all group*time cells (g,t)
such that at t, g's treatment has never changed
since the start of the panel, and all (g,t)s such that g's
treatment has changed at
least once at t, and the change occurred at a period
where the outcome is not observed.
In that subsample, the first event-study effect (denoted effect_1)
is the effect of being exposed to a higher treatment for two periods,
the second effect (effect_2)
is the effect of being exposed to a higher treatment for four periods, etc.
You may then combine the two sets of estimated effects
into one event-study graph, with the only caveat that
the "odd" and "even"
effects are estimated on different subsamples.
Importantly, the two estimations have to be run
on a dataset at the same bi-yearly level as the outcome
variable: the yearly level treatment information
should only be used to select the relevant subsamples.
{p_end}

{p 4 4}
If the treatment is observed three times more often than the
treatment, you can follow the same logic,
splitting the sample into three subsamples
and running the command three times, etc.
{p_end}

{p 4 4}
{it:How many control variables can I include in the estimation?}
{p_end}

{p 4 4}
Estimators with control variables are similar to those without controls,
except that the first-difference of the outcome
is replaced by residuals from regressions
of the first-difference of the outcome on
the first-differences of the controls and time fixed effects.
Those regressions are estimated in the sample of control (g,t)s:
(g,t)s such that group g's treatment has not changed yet at period t.
Those regressions are also estimated separately
for each value of the baseline treatment.
If the treatment takes values 0, 1, 2, 3, and 4,
one regression is estimated
for control (g,t)s with a treatment equal to 0,
one regression is estimated for control
(g,t)s with a treatment equal to 1, etc.
The number of control variables
needs to be significantly smaller than
the number of control (g,t)s in each of those regressions.
Otherwise, those regressions will overfit and produce noisy estimates.
If the number of observations is lower than the number
of variables in one of those regressions,
the command will run but will not take into account all
the controls for all values of the baseline treatment.
An error message will let the user know that
they are encountering this situation, and may thus want to reduce their
number of control variables.
{p_end}

{p 4 4}
{it:In my application, groups' baseline treatment is a continuous variable,}
{it: meaning that all groups have a different period-one treatment.}
{it: Therefore, Assumption 1 in de Chaisemartin and D'Haultfoeuille (2020a)}
{it: fails. Can I still use the command?}
{p_end}

{p 4 4}
Yes you can estimate non-normalized event-study effects. 
Essentially, you just need to define
a new treatment variable equal to 0 if g's treatment has never changed at t, 
to 1 if g's treatment has changed at t and g's period-t treatment is larger
than its baseline treatment, and to -1 if g's treatment has changed at t and g's period-t 
treatment is lower than its baseline treatment.
Then, you run the command with this new treatment, 
including interactions of period fixed effects and a polynomial in
the baseline treatment as control variables. 
For instance, if one wants to model the relationship between the counterfactual 
outcome trend and the baseline treatment as quadratic, and the data has 12 periods,
one needs to include 22 variables as controls: the baseline treatment interacted with the
period 2 to 12 fixed effects, 
and the baseline treatment squared interacted with the period
2 to 12 fixed effects.
See Section 1.5 of 
de Chaisemartin et al (2022) for further details. If groups' baseline
treatment is not continuous but takes many values, pursuing this strategy
may yield more precise estimators, applying to a larger number of switchers, than just
running the command with the original treatment, at the expense of incurring a 
potential bias if the model for the counterfactual outcome trend is misspecified.
{p_end}

{p 4 4}
{it:My design is such that treatment is binary, and groups can enter the treatment, and then leave it once.}
{it:Can I use the command to separately estimate the effect of joining and leaving the treatment?}
{p_end}

{p 4 4}
Yes you can. See Section 1.5 of the Web Appendix
of de Chaisemartin and D'Haultfoeuille (2020a) for further details.
{p_end}

{p 4 4}
{it:My design has several treatments that may all have dynamic effects.}
{it:Can I use the command to estimate the effect of a treatment controlling for other treatments?}
{p_end}

{p 4 4}
Yes, if those treatments follow binary and staggered designs.
See Section 3.2 of the Web Appendix of de Chaisemartin and D'Haultfoeuille (2020b)
for further details.
{p_end}

{p 4 4}
{it:Can I perform triple difference-in-differences with the command?}
{p_end}

{p 4 4}
Yes. Suppose for instance your third difference is across men
and women in the same (g,t) cell. Then, 
for each (g,t) cell, you just need to compute the difference between the average 
outcome of men and women in cell (g,t).
Then, you simply run the command with this new outcome.
{p_end}

{title:References}

{p 4 4}
de Chaisemartin, C, D'Haultfoeuille, X (2020a).
{browse "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3731856":Difference-in-Differences Estimators of Intertemporal Treatment Effects}.
{p_end}
{p 4 4}
de Chaisemartin, C, D'Haultfoeuille, X (2020b).
{browse "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3751060":Two-way fixed effects regressions with several treatments}.
{p_end}
{p 4 4}
de Chaisemartin, C, D'Haultfoeuille, X, Pasquier, F, Vazquez-Bare, G (2022).
{browse "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4011782":Difference-in-Differences Estimators for Treatments Continuously Distributed at Every Period}.
{p_end}
{p 4 4}
de Chaisemartin, C, D'Haultfoeuille, X, Malézieux, M, Sow, D (2023a).
{browse "https://drive.google.com/file/d/1NGgScujLCCS4RrwdN-PC1SnVigfBa32h/view?usp=drive_link":Conventions used by the did_multiplegt_dyn Stata command with missing treatments}.
{p_end}
{p 4 4}
de Chaisemartin, C, D'Haultfoeuille, X, Malézieux, M, Sow, D (2023b).
{browse "https://drive.google.com/file/d/1d23jtOT8tiHG3mpjD_cuJy3h4I2sYFEK/view?usp=drive_link":Comparing the did_multiplegt and did_multiplegt_dyn Stata commands}.
{p_end}

{title:Authors}

{p 4 4}
Clément de Chaisemartin, Economics Department, Sciences Po, France.
{p_end}
{p 4 4}
Xavier D'Haultfoeuille, CREST-ENSAE, France.
{p_end}
{p 4 4}
Mélitine Malézieux, Economics Department, Sciences Po, France.
{p_end}
{p 4 4}
Doulo Sow, CREST-ENSAE, France.
{p_end}

{title:Contact}

{p 4 4}
{browse "mailto:chaisemartin.packages@gmail.com":chaisemartin.packages@gmail.com}
{p_end}

{marker saved_results}{...}
{title:Saved results}

{p 4 4}
{cmd:e(effect_l)}: estimated event-study effect l.
{p_end}

{p 4 4}
{cmd:e(N_effect_l)}: number of observations used in the estimation of {cmd:e(effect_l)}.
{p_end}

{p 4 4}
{cmd:e(N_switchers_effect_1)}: number of switchers {cmd:e(effect_l)} applies to.
{p_end}

{p 4 4}
{cmd:e(se_effect_l)}: estimated standard error of {cmd:e(effect_l)}.
{p_end}

{p 4 4}
{cmd:e(p_equality_effects)}: p-value of a joint test
that all effects are equal, when the option {cmd:effects_equal} is specified.
{p_end}

{p 4 4}
{cmd:e(placebo_l)}: estimated placebo l.
{p_end}

{p 4 4}
{cmd:e(N_placebo_l)}: number of observations used
in the estimation of {cmd:e(placebo_l)}.
{p_end}

{p 4 4}
{cmd:e(N_switchers_placebo_l)}: number of switchers {cmd:e(placebo_l)} applies to.
{p_end}

{p 4 4}
{cmd:e(se_placebo_l)}: estimated standard error of {cmd:e(placebo_l)}.
{p_end}

{p 4 4}
{cmd:e(p_jointplacebo)}: p-value of the joint test that all the placebos are equal to 0,
if two or more placebo estimators were requested.
{p_end}

{p 4 4}
{cmd:e(effect_average)}: estimated average total
effect per unit of treatment, where “total effect”
refers to the sum of the effects of a treatment
increment, at the time when it takes place and at later periods.
{p_end}

{p 4 4}
{cmd:e(N_effect_average)}: number of observations used
in the estimation of {cmd:e(effect_average)}.
{p_end}

{p 4 4}
{cmd:e(N_switchers_effect_average)}: number of switchers*periods
{cmd:e(effect_average)} applies to.
{p_end}

{p 4 4}
{cmd:e(se_effect_average)}: estimated standard error of {cmd:e(effect_average)}.
{p_end}

{p 4 4}
{cmd:e(cmd)}: macro equal to "did_multiplegt_dyn", the name of the command.
{p_end}
