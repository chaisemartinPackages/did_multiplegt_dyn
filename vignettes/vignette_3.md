# Recovering more placebos in did_multiplegt_dyn

Section 2.1 of `did_multiplegt_dyn` [companion paper](https://drive.google.com/file/d/1NGgScujLCCS4RrwdN-PC1SnVigfBa32h/view) describes the conventions used in the computation of placebo estimators. Specifically, the most important convention can be summarised as follows:

> if $\ell$-th dynamic effect of a switching group cannot be computed, its $\ell$-th placebo will not be computed as well.

The reason for this constraint is that placebos should be informative of the correspondent dynamic effect estimators' bias coming from violations of parallel trends [(Rambachan and Roth, 2023)](https://academic.oup.com/restud/article/90/5/2555/7039335). We will give an intuition of this rationale in the next section.

A direct consequence of this approach is that users cannot require more placebos than effects, since, by definition, placebos and dynamic effects are built symmetrically.

The goal of this vignette is to show how to tweak `did_multiplegt_dyn` to compute more placebos than effects. This accrues to loosening the link between placebo and dynamic effects as described above. However, it could be useful to retrieve more placebos in application with very few post treatment periods, hence limited room for pre-trend testing.

## A simple example
|T|$D_{1,t}$| $D_{2,t}$ |
|---:|---:|---:|
|1|0|0|
|2|0|0|
|3|0|0|
|4|1|0|

$$ DID_{1,1} = Y_{1,4}(\bm{0}_3,1) - Y_{1,3}(0_3) - (Y_{2,4}(0_3, 0) -  Y_{2,3}(0_3)) $$

## Retrieving more placebos
