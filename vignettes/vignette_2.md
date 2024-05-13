# did_multiplegt_dyn and esttab

**esttab** is one of the most used Stata packages for retrieving estimation results as formatted tables. In the light of many users' requests, we have made **`did_multiplegt_dyn compatible with esttab`**. This vignette showcases how to use **esttab** in combination with **did_multiplegt_dyn** to save your results in an external file. **esttab** allows for many formats (TeX, html, ...) which can also be accessed with a **did_multiplegt_dyn** output. In this tutorial, we will only focus on LaTeX tabulars.

+ [Setup](#setup)
+ [Integration with esttab](#integration-with-esttab)
  - [General use](#general-use)
  - Formatting

## Setup

We test **esttab** via a DGP with 1000 groups and 20 periods. The treatment $D_{g,t}$ is randomly drawn from ${0,1}$ at each $(g,t)$ cell. The outcome $Y_{g,t}$ is a function of $D_{g,t}$ and time-dependent covariate $X_{g,t}$. We also generate two group-specific variables, $H1$ and $H2$, to test the output of the integration when **did_multiplegt_dyn** is run with the predict_het option.

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
The normal use of **esttab** requires to save estimation model (and eventual scalars/locals) with **estimates store** and then run **esttab** with the models' names as arguments. Let's test these basic features with a **did_multiplegt_dyn** model. We will run two specifications:
1. a model with 5 dynamic effects, where we also test for the equality of the event study estimates;
2. the same model where we also estimate 3 placebos, control for $X_{g,t}$ and test for the joint significance of the placebos.
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
esttab model_* using "filename.tex", replace booktabs se
```

The resulting table should look like this:
<p>
  <image src="https://github.com/DiegoCiccia/did_multiplegt_dyn/blob/main/vignettes/assets/reg1.png" alt>
</p>

The code above is sufficient to save **did_multiplegt_dyn** output virtually with any option set, except with the **by()** or **by_path()** options, since, by design, the program will return an e(V) and e(b) only for the last level of the *by* variable. A special case occurs with the **predict_het()** option, since the program outputs also the results from regressing the group-level estimates of the event study effects on the variables specified as the option argument. For instance, we can run the following model specification and save the results with **esttab**:

```applescript
est clear
did_multiplegt_dyn Y G T D, predict_het(H1 H2, all) graph_off effects(3)
est sto model_1
esttab model_* using "filename.tex", replace booktabs se noobs
```

In this case, the resulting **esttab** table will be partitioned as follows: 
<p>
  <image src="https://github.com/DiegoCiccia/did_multiplegt_dyn/blob/main/vignettes/assets/reg2.png" alt>
</p>

The first equation box (labeled by the outcome variable) contains the main results from **did_multiplegt_dyn**, while the equation boxes below show the output from **predict_het()**.
