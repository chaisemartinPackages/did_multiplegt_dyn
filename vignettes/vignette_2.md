# did_multiplegt_dyn and esttab

**esttab** is one of the most used Stata packages for retrieving estimation results as formatted tables. In the light of many users' requests, we have made **`did_multiplegt_dyn compatible with esttab`**. This vignette showcases how to use **esttab** in combination with **did_multiplegt_dyn** to save your results in an external file. **esttab** allows for many formats (TeX, html, ...) which can also be accessed with a **did_multiplegt_dyn** output. In this tutorial, we will only focus on LaTeX tabulars.

+ [Setup](#setup)
+ [Integration with esttab](#integration-with-esttab)
  - [General use](#general-use)
  - [Formatting](#formatting)

## Setup

We test **esttab** via a DGP with 1000 groups and 20 periods. The treatment $D_{g,t}$ is randomly drawn from $\lbrace 0,1\rbrace$ at each $(g,t)$ cell. The outcome $Y_{g,t}$ is a function of $D_{g,t}$ and a time-dependent covariate $X_{g,t}$. We also generate two group-specific variables, $H^1_g$ and $H^2_g$, to test the output of the integration when **did_multiplegt_dyn** is run with the predict_het option.

```applescript
clear
set seed 123
scalar TT = 20
scalar GG = 1000
set obs `= TT * GG'
gen G = mod(_n-1,GG) + 1
gen T = floor((_n-1)/GG)
sort G T
gen D = uniform() > 0.5
gen X = uniform() * T
gen Y = uniform() * (1 + D + X)
gen H1 = G/GG
gen H2 = mod(G, TT)
```

## Integration with esttab

### General use
The normal use of **esttab** requires saving the estimation model (plus eventual scalars/locals) with **estimates store** and then runnig **esttab** with the models' names as arguments. Let's test these basic features with a **did_multiplegt_dyn** model. We will run two specifications:
1. a model with 5 dynamic effects, where we also test for the equality of the event study estimates;
2. the same model, but we also estimate 3 placebos, control for $X_{g,t}$ and test for the joint significance of the placebos.

This setting allows us to showcase how to retrieve the estimation results and add scalars/locals from **did_multiplegt_dyn** to **esttab**. 

First, we run **did_multiplegt_dyn** with the specifications above and save the results with **est sto** and the scalars/locals with **estadd**.

```applescript
est clear
did_multiplegt_dyn Y G T D, effects(5) graph_off effects_equal
estadd scalar p_joint = e(p_equality_effects) 
estadd local controls = "No"
est sto model_1

did_multiplegt_dyn Y G T D, effects(5) placebo(3) effects_equal controls(X) graph_off
estadd scalar p_joint = e(p_equality_effects) 
estadd scalar p_placebo = e(p_jointplacebo)
estadd local controls = "Yes"
est sto model_2
```

Then, we use **esttab** to save the results in a TeX tabular:
```applescript
esttab model_* using "filename.tex", replace booktabs se s(p_joint p_placebo controls)
```

The resulting table should look like this:
<p>
  <image src="https://github.com/chaisemartinPackages/did_multiplegt_dyn/blob/main/vignettes/assets/reg1.png" alt>
</p>

The code above is sufficient to save **did_multiplegt_dyn** output virtually with any option set, except with the **by()** or **by_path()** options, since, by design, the program will return an e(V) and e(b) only for the last level of the *by* variable. 

A special case occurs with the **predict_het()** option, since the program outputs also the results from regressing the group-level estimates of the event study effects on the variables specified as the option argument. For instance, we can run the following model specification and save the results with **esttab**:

```applescript
est clear
did_multiplegt_dyn Y G T D, predict_het(H1 H2, all) graph_off effects(3)
est sto model_1
esttab model_* using "filename.tex", replace booktabs se noobs
```

In this case, the resulting **esttab** table will be partitioned as follows: 
<p>
  <image src="https://github.com/chaisemartinPackages/did_multiplegt_dyn/blob/main/vignettes/assets/reg2.png" alt>
</p>

The first equation box (labeled after the outcome variable) contains the main results from **did_multiplegt_dyn**, while the equation boxes below show the output from **predict_het()**.

### Formatting

**esttab** allows for a vast range of customisation. In this last subsection, we show how to improve the formatting of the first table in this tutorial. Specifically, we will change
+ the scalar/local labels, via `s(..., label(...))`;
+ the format of the floats, via `b(fmt)`;
+ the coefficient labels, as to mirror the notation from de Chaisemartin and D'Haultfoeuille (2024), via `coeflabels(coef "label")`;
+ Stata to LaTeX syntax, via `substitute(\_ _)`;
+ the space between the column number and the first equation line, via `mlabels(,none) collabels(,none)`.

Run again the first code block from the [previous subsection](#integration-with-esttab). Then, run the following:

```
esttab model_* using "filename.tex", replace booktabs se s(control p_joint p_placebo, label("Controls" "Joint Eq. Effects" "Joint Sig. Placebo"))  b(%9.5fc) coeflabels(Effect_1 "$\hat{\delta}_1$" Effect_2 "$\hat{\delta}_2$" Effect_3 "$\hat{\delta}_3$" Effect_4 "$\hat{\delta}_4$" Effect_5 "$\hat{\delta}_5$" Avg_Tot_Effect "$\hat{\delta}$" Placebo_1 "$\hat{\delta}_1^{pl}$" Placebo_2 "$\hat{\delta}_2^{pl}$" Placebo_3 "$\hat{\delta}_3^{pl}$") substitute(\_ _) mlabels(,none) collabels(,none)
```

The code above yields this table:
<p>
  <image src="https://github.com/chaisemartinPackages/did_multiplegt_dyn/blob/main/vignettes/assets/reg3.png" alt>
</p>

---
