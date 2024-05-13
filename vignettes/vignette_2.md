# did_multiplegt_dyn and esttab

`esttab` is one of the most used Stata packages for retrieving estimation results as formatted tables. In the light of many users' requests, we have made **`did_multiplegt_dyn` compatible with `esttab`**. This vignette showcases how to use `esttab` in combination with `did_multiplegt_dyn` to save your results in an external file. `esttab` allows for many formats (TeX, html, ...) which can also be accessed with a `did_multiplegt_dyn` output. In this tutorial, we will only focus on LaTeX tabulars.

+ [Setup](#setup)
+ Integration with esttab
  - General use
  - Formatting

## Setup

We test `esttab` via a DGP with 1000 groups and 20 periods. The treatment $D_{g,t}$ is randomly drawn from ${0,1}$ at each $(g,t)$ cell. The outcome $Y_{g,t}$ is a function of $D_{g,t}$ and time-dependent covariate $X_{g,t}$. We also generate two group-specific variables, H1 and H2, to test the output of the integration when `did_multiplegt_dyn` is run with the `predict_het` option.

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

<p>
  <image src="https://github.com/DiegoCiccia/did_multiplegt_dyn/blob/main/vignettes/assets/reg1.png" alt>
</p>