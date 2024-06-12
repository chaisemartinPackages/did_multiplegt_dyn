# Recovering more placebos in did_multiplegt_dyn

Section 2.1 of `did_multiplegt_dyn` [companion paper](https://drive.google.com/file/d/1NGgScujLCCS4RrwdN-PC1SnVigfBa32h/view) describes the conventions used in the computation of placebo estimators. Specifically, the most important convention can be summarised as follows:

> if the $\ell$-th dynamic effect of a switching group cannot be computed, its $\ell$-th placebo will not be computed as well.

The reason for this constraint is that placebos should be informative of the correspondent dynamic effect estimators' bias coming from violations of parallel trends [(Rambachan and Roth, 2023)](https://academic.oup.com/restud/article/90/5/2555/7039335). We will give an intuition of this rationale in the next section.
A direct consequence of this approach is that users cannot require more placebos than effects, since, by definition, placebos and dynamic effects are built symmetrically.
The goal of this vignette is to show how to tweak `did_multiplegt_dyn` to compute more placebos than effects. This accrues to loosening the link between placebo and dynamic effects as described above. However, it could be useful to retrieve more placebos in application with very few post treatment periods, hence limited room for pre-trend testing.

## A simple example

Let's take the case described in the table below. There are two groups and four periods, with treatment $D_{g,t}$ reported in column. Group $1$ switches at the last period, while group $2$ never switches. 

|T|$D_{1,t}$| $D_{2,t}$ |
|---:|---:|---:|
|1|0|0|
|2|0|0|
|3|0|0|
|4|1|0|

Henceforth, we will use the notation and results from [de Chaisemartin and D'Haultfeuille (2024)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3731856). Furthermore, let $0_t$ and $1_t$ denote vectors of 0s and 1s of length $t$. 

### The link between placebos and dynamic effects
Due to the limitations of the data, the only dynamic effect we can compute is $DID_{1}$, which can be disaggregated as follows:

$$
\begin{align*}
DID_{1} = DID_{1,1} &= Y_{1,4}(0_3,1) - Y_{1,3}(0_3) - (Y_{2,4}(0_4) -  Y_{2,3}(0_3)) \\
 &= (Y_{1,4}(0_3,1) - Y_{1,4}(0_4)) + (Y_{1,4}(0_4) - Y_{1,3}(0_3)) - (Y_{2,4}(0_4) -  Y_{2,3}(0_3)).
\end{align*}
$$

The first term is the treatment effect for group $1$, i.e. the difference between the actual and the counterfactual outcome one period after the first switch. The second and third terms are the differences between the untreated outcomes of the two groups at periods $3$ and $4$. If the parallel trends assumption holds, then the last two terms should cancel out each other, since the untreated potential outcome evolution is equal for all groups with the same status quo treatment. Specifically, one could use $DID^{pl}_{1}$ to test whether the assumption that group 1 and 2 experience parallel evolutions of their untreated outcomes holds for 1 period. In that case, we would be comparing the period-2-to-3 outcome evolution of group 1 against 2. In this sense, $DID^{pl}_1$ is _informative_ of the potential violation of parallel trends assumption that could bias $DID_1$.

As noted above, $DID_2$ is not defined in this dataset, but it is possible to estimate the second placebo by comparing the period-1-to-3 outcome evolution of group 1 against 2. Still, the second placebo conveys information about potential violations of the parallel trends assumption that cannot be used to support the unbiasedness of another dynamic effect. As a result, the second placebo is quite uninformative under this rationale.

### A simple tweak to retrieve more placebos

The paragraph above explains why `did_multiplegt_dyn` does not allow to request more placebos than effects. However, it could be argued that tests for pre-trends under this rationale do not exploit in full the available data. To this end, it is possible to tweak the data and retrieve the placebos that are, by default, ignored by the program.

Suppose that we want to force `did_multiplegt_dyn` to compute $DID^{pl}_2$ using the treatment rollout in the table above. As remarked above, this boils down to allowing for the estimation of $DID_2$. A way to do that is to create an artificial period $5$ where the treatment (and outcome) of all groups is equal to 0 (henceforth, we will call them _null periods_):

|T|$D_{1,t}$| $D_{2,t}$ |
|---:|---:|---:|
|1|0|0|
|2|0|0|
|3|0|0|
|4|1|0|
|5|0|0|

Now, group 1 has two periods after the first switch. As a result, $DID_2$ is defined and $DID^{pl}_2$ is estimated by default. 
In general, users should:
1. run `did_multiplegt_dyn` to retrieve the dynamic effects estimates;
2. run `did_multiplegt_dyn` again with the augmented dataset to retrieve the placebo estimates.

Notice that the dynamic effects estimates from the second run should be disregarded, since the outcome of all groups that do not have a second period after their first switch is set to zero. 

## Retrieving more placebos - Code Example (Stata)

>[!WARNING]
>
> The following method can be used if there is <ins>at least one never-switcher group per baseline treatment</ins>.
> The companion program to run compute the remaining placebos does not support any option from `did_multiplegt_dyn`, except **effects**, **placebo**, **switchers()** and **only_never_switchers**. 
> Integration with other options will depend on future demand for this feature.

Let's generate a toy dataset with 25 groups and 4 periods. Groups can either be never-switchers or switch from 0 to 1 at the last period. 
The setting is equivalent to the example above, even though we allow for more than 2 groups to draw inference on the point estimates.
As a result, we cannot compute more than 1 dynamic effect nor placebo.

```stata
clear
set seed 0
set obs 100
gen G = mod(_n-1, 25) + 1
bys G: gen T = _n
gen D = uniform() > 0.5 & T == 4
gen Y = uniform() * (1 + D)

did_multiplegt_dyn Y G T D, effects(1) placebo(1) graph_off
```

Starting from June 2024, the Stata version of `did_multiplegt_dyn` includes a subcommand called `did_multiplegt_dyn_all_pl` to retrieve all the feasible placebo estimates. The syntax of this subcommand is the same as the main program, with a few exceptions:

+ the user can request more placebos than effects;
+ the user can request only **switchers()** and **only_never_switchers** as additional options.

Even though the user can request more placebos than effects, the number of placebo estimates that the program can compute can never exceed 
$$
\max_{g:F_g \neq T+1} F_g - 2
$$
that is, the maximum, across all switchers, of their pre-first-switch periods. As a result, an error message will be prompted if the limit above is exceeded.
If we wish to retrieve the last placebo, we can run the following line:

```
did_multiplegt_dyn_all_pl Y G T D, effects(1) placebo(2) 
```

Notice that the output is divided into two steps

1. First, the program retrieves the effects as in vanilla `did_multiplegt_dyn`;
2. Then, the dataset is filled with enough null periods to retrieve the requested placebos.

As a result, the new output will include one additional placebo compared with the previous case.
Notice that we are requesting more placebos than effects, yet, since groups switch only at period 4, the number of placebos requested within the boundaries of the maximum number of placebos we can compute, i.e. 4 - 2 = 2.
The combined output from `did_multiplegt_dyn` and `did_multiplegt_dyn_all_pl` will look like this:

<table border="0" width="*">
<tr><td>            </td><td>         (1)              </td><td>         (2)              </td></tr>
<tr><td>Effect_1    </td><td>       0.343              </td><td>       0.343              </td></tr>
<tr><td>            </td><td>     (0.219)              </td><td>     (0.219)              </td></tr>
<tr><td>Placebo_1   </td><td>     -0.0587              </td><td>     -0.0587              </td></tr>
<tr><td>            </td><td>     (0.158)              </td><td>     (0.158)              </td></tr>
<tr><td>Placebo_2   </td><td>                          </td><td>     -0.0799              </td></tr>
<tr><td>            </td><td>                          </td><td>     (0.136)              </td></tr>
<tr><td>Av_tot_eff  </td><td>       0.343              </td><td>       0.343              </td></tr>
<tr><td>            </td><td>     (0.219)              </td><td>     (0.219)              </td></tr>
<tr><td colspan=3><hr></td></tr>
Standard errors in parentheses
<br /><sup>*</sup> <i>p</i> < 0.05, <sup>**</sup> <i>p</i> < 0.01, <sup>***</sup> <i>p</i> < 0.001
</td></tr>
</table>