# Recovering more placebos in did_multiplegt_dyn

Section 2.1 of `did_multiplegt_dyn` [companion paper](https://drive.google.com/file/d/1NGgScujLCCS4RrwdN-PC1SnVigfBa32h/view) describes the conventions used in the computation of placebo estimators. Specifically, the most important convention can be summarised as follows:

> if $\ell$-th dynamic effect of a switching group cannot be computed, its $\ell$-th placebo will not be computed as well.

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

### The link between placebos and effects
Due to the limitations of the data, the only dynamic effect we can compute is $DID_{1}$, which can be disaggregated as follows:

$$
\begin{align*}
DID_{1} = DID_{1,1} &= Y_{1,4}(0_3,1) - Y_{1,3}(0_3) - (Y_{2,4}(0_4) -  Y_{2,3}(0_3)) \\
 &= (Y_{1,4}(0_3,1) - Y_{1,4}(0_4)) + (Y_{1,4}(0_4) - Y_{1,3}(0_3)) - (Y_{2,4}(0_4) -  Y_{2,3}(0_3)).
\end{align*}
$$

The first term is the treatment effect for group $1$, i.e. the difference between the actual and the counterfactual outcome one period after the first switch. The second and third terms are the differences between the untreated outcomes of the two groups at periods $3$ and $4$. If the parallel trends assumption holds for at least one period, then the last two terms should cancel out each other, since the untreated potential outcome evolution is equal for all groups with the same status quo treatment. Specifically, one could use $DID^{pl}_{1}$ to assess the assumption that group 1 and 2 experience parallel evolutions of their untreated outcomes holds for 1 period.




## Retrieving more placebos
