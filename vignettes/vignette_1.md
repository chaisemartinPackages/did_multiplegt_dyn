# Outcome observed less frequently than the treatment

In many circumstances, the outcome is observed less frequently than the treatment. For instance, electoral outcomes are recorded only during election periods, while treatment may be observed every period. In this tutorial, we show how to use of the estimator from de Chaisemartin & D'Haultfoeuille (2024) with outcomes observed less frequently than the treatment. 

#### Sections
+ [A toy example](#a-toy-example)
+ [General Case with Stata and R code](#general-case-with-stata-and-r-code)
  - [Part I: Data Generation](#part-i-data-generation)
  - [Part II: Data Adjustment](#part-ii-data-adjustment)
  - [Part III: Estimation](#part-iii-estimation)
  - [Part IV: Comparison with naive did_multiplegt_dyn](#part-iv-comparison-with-naive-did_multiplegt_dyn)
  - [Part V: Graph output](#part-v-graph-output)
+ [Requesting placebos](#requesting-placebos)


## A toy example

Let's take the case of a researcher wanting to estimate the effect of the induction of a new party leader on electoral performance. Party $A$ has changed its main candidate in 2004. Assume that we have biannual election data as follows:

| Party |period |Share |LeadChangeYr|Treatment|
| ----: |----:|----: |----:       |----:    |
| A     |2003 |0.40  |2004        |0        |
| A     |2005 |0.35  |2004        |1        |
| A     |2007 |0.25  |2004        |1        |
| B     |2003 |0.40  |.           |0        |
| B     |2005 |0.45  |.           |0        |
| B     |2007 |0.50  |.           |0        |

Under a parallel-trends assumption, for $\ell=2$ and $\ell=4$, the effect of leadership change after $\ell$ periods can be estimated by a difference-in-difference comparing the 2003-to-(2003+ $\ell$) change in the vote share of treated party $A$ to the 2003-to-(2003+ $\ell$) change in the vote share of control party $B.$ However, for $\ell=1,3$, this DID cannot be computed, because elections do not happen in 2003+ $\ell$.

Assume now that the dataset includes another party, $C$, that has experienced a leadership change in 2005:

| Party |period |Share |LeadChangeYr|Treatment|
| ----: |----:|----: |----:       |----:    |
| A     |2003 |0.40  |2004        |0        |
| A     |2005 |0.35  |2004        |1        |
| A     |2007 |0.30  |2004        |1        |
| B     |2003 |0.40  |.           |0        |
| B     |2005 |0.45  |.           |0        |
| B     |2007 |0.50  |.           |0        |
| C     |2003 |0.15  |2005        |0        |
| C     |2005 |0.10  |2005        |1        |
| C     |2007 |0.05  |2005        |1        |

In this example, naively applying the event-study estimators from de Chaisemartin & D'Haultfoeuille (2024) would yield a first event-study estimator that averages the effect of one period of exposure to treatment for party C with the effect of two periods of exposure for party A. Similarly, the second event-study estimator averages the effect of three periods of exposure to treatment for party C with the effect of four periods of exposure for party A. Instead, we can estimate the effect of $\ell$ periods of exposure to leadership change for $\ell=1$ and $\ell=3$, using a DID estimator comparing the 2003-to-(2003+ $1+\ell$) change in the vote share of treated party $C$ to the 2003-to-(2003+ $1+\ell$) change in the vote share of control party $B.$ This DID estimator uses the 2003 voting rate, the most recent non-missing outcome before C gets treated, as the baseline outcome. Combining these DID estimators with those comparing parties A and B, we are now able to estimate four event-study effects, without averaging effects of different exposure lengths, with effects $\ell=1$ and $\ell=3$ applying to Party C, while effects $\ell=2$ and $\ell=4$ apply to Party A. The figure below shows the combined event-study plot from our toy example.

<p>
  <img src="https://github.com/chaisemartinPackages/did_multiplegt_dyn/blob/main/vignettes/assets/vignette_1_Stata_fig1.jpg" alt>
</p>

In the next section, we show how to use `did_multiplegt_dyn` to produce an event-study graph without averaging effects of different exposure lengths whenever the outcome is observed less frequently than the treatment. For each step, we also report the corresponding block of [Stata](https://github.com/chaisemartinPackages/did_multiplegt_dyn/blob/main/vignettes/src/vignette_1_Stata.do) and [R](https://github.com/chaisemartinPackages/did_multiplegt_dyn/blob/main/vignettes/src/vignette_1_R.R) code.

## General Case with Stata and R code

### Part I: Data Generation 

> [!IMPORTANT]
> If you already have data to which you would like to apply this tutorial, you can skip to [Part II](#part-ii-data-adjustment). Before doing so, just make sure to rename your variables of interest as follows:
> + Outcome variable $\to$ "Y"
> + Group variable $\to$ "G"
> + Time variable $\to$ "T"
> + Treatment variable $\to$ "D"
>   
> Also, your group and time variables should be positive integers (in Stata, `egen group`, in R, `cur_group_id()`).

We use a DGP with five groups and 20 periods. We can observe the outcome every fourth period. If the outcome is observed at a different frequency in your application, you just need to replace 4 by the corresponding integer in the code below. Group 1's treatment never switches, while groups 2 to 5 switch at periods 4 to 7 (group id + 2). The untreated outcome is uniformly distributed and increases by 100 after the first treatment switch.

<table>
  <tr>
    <th>Stata</th>
    <th>R</th>
  </tr>
  <tr>
    <td>
    <pre><code>
clear
set seed 123
scalar TT = 20
scalar GG = 5
set obs `= TT * GG'
gen G = mod(_n-1,GG) + 1
gen T = floor((_n-1)/GG)
sort G T
gen D = 0
forv j=2/5 {
    replace D = 1 if G == `j' & T == `j' + 2
}
bys G: gen D_stag = sum(D)
gen Y = uniform() * (1 + 100*D_stag) if mod(T,4) == 0
drop D_stag
browse
    </pre></code>
    </td>
    <td>
    <pre><code>
library(dplyr)
set.seed(123)
TT <- 20; GG <- 5
df <- data.frame(id = 1:(GG*TT))
df$G <- ((df$id-1) %% GG)+1
df$T <- floor((df$id-1)/GG)
df$id <- NULL
df <- df[order(df$G, df$T), ]
df$D <- 0
for (v in c(2,3,4,5)) {
    df$D <- ifelse(df$G == v & df$T == v+2, 1, df$D)
}
df <- df %>% group_by(.data$G) %>% 
  mutate(D_stag = cumsum(.data$D)) %>% ungroup()
df$Y <- ifelse(df$T %% 4 == 0, 
  runif(n = nrow(df)) * (1 + 100*df$D_stag), NA)
df$D_stag <- NULL
View(df)
    </pre></code>
    </td>
  </tr>
</table>

### Part II: Data Adjustment 
For switchers (i.e. groups whose treatment changes before the end of the panel), we refer to $(g,t)$ cells such that $t$ is lower than (equal to or higher than) the first period where their treatment changes as *pre-switch cells* (resp. *post-switch cells*). We need to generate a variable through which we partition our population of $(g,t)$ cells into 5 subsamples:
+ A) All $(g,t)$ cells of groups whose treatment does not change (never-switchers) plus all pre-switch cells;
+ B) Post-switch cells of groups whose treatment changes for the first time at a time period where the outcome is non-missing;
+ C) Post-switch cells of groups whose treatment changes for the first time one period before a period where the outcome is non-missing;
+ D) Post-switch cells of groups whose treatment changes for the first time two periods before a period where the outcome is non-missing;
+ E) Post-switch cells of groups whose treatment changes for the first time three periods before a period where the outcome is non-missing.
  
Then, we will use the values of this variable to run `did_multiplegt_dyn` on subsamples A) and B), then on subsamples A) and C), etc.

#### a) Identify whether treatment has changed within each group

<table>
  <tr>
    <th>Stata</th>
    <th>R</th>
  </tr>
  <tr>
    <td>
    <pre><code>
bys G: gen D0 = D[1]
gen D_change = abs(D - D0) != 0
bys G: gen at_least_one_D_change = sum(D_change)
    </pre></code>
    </td>
    <td>
    <pre><code>
df$D0 <- df$D[(df$G-1)*length(levels(factor(df$T)))+1]
df$D_change <- as.numeric(abs(df$D - df$D0) != 0)
df <- df %>% group_by(.data$G) %>% 
  mutate(at_least_one_D_change = cumsum(.data$D_change)) %>% ungroup()
    </pre></code>
    </td>
  </tr>
</table>

#### b) Identify when treatment changed for the first time

We generate a variable F_g equal to the period when group g's treatment changes for the first time. The modulus of F_g divided by 4 is equal to zero is g's treatment changed for the first time in a period when the outcome is non missing, it is equal to 1 (resp. 2, 3) if g's treatment changed for the first time 3 (resp. 2, 1) periods before a period when the outcome is non missing. The aforementioned partition variable is 0 for never-switchers and pre-switch cells of switchers. For post-switch cells of groups switching on a period divisible by 4, it is equal to 1. Lastly, for post-switch cells of groups switching 1, 2 or 3 periods prior to a non-missing outcome period, it is equal to 4 minus their modulus plus 1.

<table>
  <tr>
    <th>Stata</th>
    <th>R</th>
  </tr>
  <tr>
    <td>
    <pre><code>
bys G: egen never_treated = max(at_least_one_D_change)
replace never_treated = 1 - never_treated
bys G: egen F_g_temp = min(T * D_change) if D_change != 0
bys G: egen F_g = mean(F_g_temp)
sum T
replace F_g = r(max) + 1 if missing(F_g)
gen subsample = (4 - mod(F_g, 4)) * (mod(F_g, 4) != 0) + 1
replace subsample = 0 if at_least_one_D_change == 0
    </pre></code>
    </td>
    <td>
    <pre><code>
df <- df %>% group_by(.data$G) %>% 
mutate(never_treated = as.numeric(sum(.data$D_change, na.rm = TRUE) == 0)) %>%
mutate(F_g = ifelse(.data$never_treated == 1, max(df$T, na.rm = TRUE) +1, 
  min(ifelse(.data$D_change == 0, NA, .data$T * .data$D_change), na.rm = TRUE))) %>% 
    ungroup()
df$subsample <- (4 - (df$F_g %% 4)) * (df$F_g %% 4 != 0) + 1
df$subsample <- df$subsample * df$at_least_one_D_change
    </pre></code>
    </td>
  </tr>
</table>

#### c) Keep if Y is not missing
<table>
  <tr>
    <th>Stata</th>
    <th>R</th>
  </tr>
  <tr>
    <td>
    <pre><code>
keep if !missing(Y)
    </pre></code>
    </td>
    <td>
    <pre><code>
df <- subset(df, !is.na(df$Y))
    </pre></code>
    </td>
  </tr>
</table>

### Part III: Estimation
Once the dataset has only non missing outcomes, we run `did_multiplegt_dyn` 4 times: 
+ First, we keep groups in subsamples A) and B) (subsample=0 or subsample=1). In that subsample, the first event-study effect is an effect of one period of exposure to treatment, the second event-study effect is an effect of five periods of exposure, etc.
+ Then, we keep groups in subsamples A) and C) (subsample=0 or subsample=2). In that subsample, the first event-study effect is an effect of two periods of exposure to treatment, the second event-study effect is an effect of six periods of exposure, etc.
+ etc. 

<table>
  <tr>
    <th>Stata</th>
    <th>R</th>
  </tr>
  <tr>
    <td>
    <pre><code>
local effects = 2
mat define res = J(4*`effects', 6, .)
local r_effects ""
forv j=1/4 {
    did_multiplegt_dyn Y G T at_least_one_D_change if inlist(subsample, 0, `j'), effects(`effects') graph_off
    forv i = 1/`effects'{
        mat adj = mat_res_XX[`i',1..6]
        forv c =1/6 {
            mat res[`j'+(`i'-1)*4,`c'] = adj[1, `c']
        }
    }
}
mat li res
    </pre></code>
    </td>
    <td>
    <pre><code>
library(DIDmultiplegtDYN)
effects <- 2
table <- NULL
for (j in 1:4) {
    temp <- did_multiplegt_dyn(
      subset(df, df$subsample %in% c(0, j)), "Y", "G", "T", "at_least_one_D_change", 
      graph_off = TRUE, effects = effects)
    rownames(temp$results$Effects) <- 
      sapply(1:temp$results$N_Effects, function(x) paste0("Effect_",  j + (x-1) * 4))
    table <- rbind(table, temp$results$Effects)
}
rown <- unlist(strsplit(rownames(table), "_")) 
table <- cbind(table, as.numeric(rown[rown != "Effect"]))
print(table[order(table[,ncol(table)]),1:(ncol(table)-1)])
    </pre></code>
    </td>
  </tr>
</table>

### Part IV: Comparison with naive did_multiplegt_dyn
Running `did_multiplegt_dyn` on the collapsed data only yields biased point estimates for the event-study coefficients. Let us take the case of a slightly different DGP, where $$Y_{g,t} = U (1 + \sum_{t'\leq t} D_{g,t'})$$ with $U \sim U(0,1)$. The outcome is now increasing in the cumulative treatment received over time.  We increase the number of groups to 1000, while keeping 20 periods, and $Y_{g,t}$ is still non missing only at every fourth period. Groups whose index is divisible by 5 as never-switchers. The following code blocks generate a random sample from this DGP.

<table>
  <tr>
    <th>Stata</th>
    <th>R</th>
  </tr>
  <tr>
    <td>
    <pre><code>
clear
set seed 123
scalar TT = 20
scalar GG = 1000
set obs `= TT * GG'
gen G = mod(_n-1,GG) + 1
gen T = floor((_n-1)/GG)
sort G T
gen D = 0
forv j= 0/3 {
    replace D = 1 if mod(G, 4) == `j' & T == `j' + 2 & mod(G, 5) != 0
}
bys G: gen D_stag_temp = sum(D)
bys G: gen D_stag = sum(D_stag_temp)
gen Y = uniform() * (1 + D_stag)
replace Y = . if mod(T, 4) != 0
drop D_stag*
bys G: gen D0 = D[1]
gen D_change = abs(D - D0) != 0
bys G: gen at_least_one_D_change = sum(D_change)
bys G: egen never_treated = max(at_least_one_D_change)
replace never_treated = 1 - never_treated
bys G: egen F_g_temp = min(T * D_change) if D_change != 0
bys G: egen F_g = mean(F_g_temp)
sum T
replace F_g = r(max) + 1 if missing(F_g)
gen subsample = (4 - mod(F_g, 4)) * (mod(F_g,4) != 0) + 1
replace subsample = 0 if at_least_one_D_change == 0
sort G T
keep if !missing(Y)
    </pre></code>
    </td>
    <td>
    <pre><code>
set.seed(123)
library(dplyr)
library(DIDmultiplegtDYN)
TT <- 20; GG <- 1000
df <- data.frame(id = 1:(GG*TT))
df$G <- ((df$id-1) %% GG)+1
df$T <- floor((df$id-1)/GG)
df$id <- NULL
df <- df[order(df$G, df$T), ]
df$D <- 0
for (v in 0:3) {
    df$D <- ifelse(df$G %% 4 == v & df$T == v+2 & df$G %% 5 != 0, 1, df$D)
}
df <- df %>% group_by(.data$G) %>% mutate(D_stag_temp = cumsum(.data$D)) %>% ungroup()
df <- df %>% group_by(.data$G) %>% mutate(D_stag = cumsum(.data$D_stag_temp)) %>% ungroup()
df$Y <- ifelse(df$T %% 4 == 0, runif(n = nrow(df)) * (1 + df$D_stag), NA)
df$D_stag <- NULL
df$D0 <- df$D[(df$G-1)*length(levels(factor(df$T)))+1]
df$D_change <- as.numeric(abs(df$D - df$D0) != 0)
df <- df %>% group_by(.data$G) %>% mutate(at_least_one_D_change = cumsum(.data$D_change)) %>% ungroup()
df <- df %>% group_by(.data$G) %>% 
mutate(never_treated = as.numeric(sum(.data$D_change, na.rm = TRUE) == 0)) %>%
mutate(F_g = ifelse(.data$never_treated == 1, max(df$T, na.rm = TRUE) +1, min(ifelse(.data$D_change == 0, NA, .data$T * .data$D_change), na.rm = TRUE))) %>% ungroup()
df$subsample <- (4 - (df$F_g %% 4)) * (df$F_g != 4) + 1
df$subsample <- df$subsample * df$at_least_one_D_change
df <- subset(df, !is.na(df$Y))
    </pre></code>
    </td>
  </tr>
</table>

From the definition of the DGP, the average treatment effect of being treated for $\ell$ periods should be equal to $\ell/2$, i.e. the average of $U$ times $\ell$. However, we get a way higher point estimate for $\ell = 1$ when we run `did_multiplegt_dyn` on the full sample:

<table>
  <tr>
    <th>Stata</th>
    <th>R</th>
  </tr>
  <tr>
    <td>
    <pre><code>
did_multiplegt_dyn Y G T at_least_one_D_change, graph_off effects(1)
    </pre></code>
    </td>
    <td>
    <pre><code>
did_multiplegt_dyn(df, "Y", "G", "T", "at_least_one_D_change", graph_off = TRUE, effects = 1)
    </pre></code>
    </td>
  </tr>
</table>

This is due to the fact that effect one averages the treatment effects of groups that have switched up to 3 periods before the first switch detected in collapsed data. Instead, running the loop above with the partitioned data yields the correct point estimates

<table>
  <tr>
    <th>Stata</th>
    <th>R</th>
  </tr>
  <tr>
    <td>
    <pre><code>
local effects = 2
mat define res = J(4*`effects', 6, .)
local r_effects ""
forv j=1/4 {
    did_multiplegt_dyn Y G T at_least_one_D_change if inlist(subsample, 0, `j'), effects(`effects') graph_off
    forv i = 1/`effects'{
        mat adj = mat_res_XX[`i',1..6]
        forv c =1/6 {
            mat res[`j'+(`i'-1)*4,`c'] = adj[1, `c']
        }
    }
}
mat li res
    </pre></code>
    </td>
    <td>
    <pre><code>
library(DIDmultiplegtDYN)
effects <- 2
table <- NULL
for (j in 1:4) {
    temp <- did_multiplegt_dyn(
      subset(df, df$subsample %in% c(0, j)), "Y", "G", "T", "at_least_one_D_change", 
      graph_off = TRUE, effects = effects)
    rownames(temp$results$Effects) <- 
      sapply(1:temp$results$N_Effects, function(x) paste0("Effect_",  j + (x-1) * 4))
    table <- rbind(table, temp$results$Effects)
}
rown <- unlist(strsplit(rownames(table), "_")) 
table <- cbind(table, as.numeric(rown[rown != "Effect"]))
print(table[order(table[,ncol(table)]),1:(ncol(table)-1)])
    </pre></code>
    </td>
  </tr>
</table>

### Part V: Graph output

The output of the loop above is a matrix, which can be turned into an event-study plot very easily. In R, we use the ggplot2 library, while in Stata we use the built-in gr commands.
<table>
  <tr>
    <th>Stata</th>
    <th>R</th>
  </tr>
  <tr>
    <td>
    <pre><code>
mat res = (0,0,0,0,0,0) \ res
svmat res
gen rel_time = _n-1 if !missing(res1)
local xtitle "Relative time to last period before treatment changes (t=0)"
local title "DID, from last period before treatment changes (t=0) to t"
tw rcap res3 res4 rel_time, lc(blue) || connected res1 rel_time, mc(blue) lc(blue) ||, xtitle("`xtitle'") title("`title'") leg(off) 
    </pre></code>
    </td>
    <td>
    <pre><code>
library(ggplot2)
table <- table[order(table[,ncol(table)]), ]
table <- rbind(rep(0, ncol(table)), table)
colnames(table)[ncol(table)] <- "Time"
table <- as.data.frame(table)
out_plot <- ggplot(table, aes(x = .data$Time, y = .data$Estimate, group = 1)) + 
geom_line(colour = "blue") +
geom_errorbar(data = ~dplyr::filter(.x, table$Estimate != 0), aes(ymin = .data[["LB CI"]], ymax = .data[["UB CI"]]), 
position=position_dodge(0.05), width = 0.2, colour = "red") + 
geom_point(colour = "blue") + 
ggtitle("DID, from last period before treatment changes (t=0) to t") + 
xlab("Relative time to last period before treatment changes (t=0)") +
theme(plot.title = element_text(hjust = 0.5))
print(out_plot)
    </pre></code>
    </td>
  </tr>
</table>

The resulting graph should look like this:
<p>
  <img src="https://github.com/chaisemartinPackages/did_multiplegt_dyn/blob/main/vignettes/assets/vignette_1_Stata_fig2.jpg" alt>
</p>

---

## Requesting Placebos

Testing for no anticipation and parallel trends assumptions becomes trickier when the outcome variable is missing at regular intervals.
The econometrics behind this case is explained in the [companion paper](https://github.com/chaisemartinPackages/did_multiplegt_dyn/blob/main/vignettes/assets/main.pdf) of this vignette.
In general, we show that in a dataset with an outcome variable missing every $k$ periods
+ the mapping of the placebos from collapsed data differs from the mapping that we used to retrieve the dynamic effects estimators:
+ the earliest placebo estimator retrievable from the current data is the $k$-th one;
+ (almost) all placebo estimators retrieved in this way can be disaggregated into the _difference_ of two normal placebo estimators.

As a result, the interpretation of the estimates changes, in that finding that the difference of two placebos is not statistically different from 0 is for sure a weaker test of the parallel trends assumption.

In any case, it is possible to tweak `did_multiplegt_dyn` to retrieve placebo estimators under the caveats discussed above. Let's change the DGP of the example to allow the computation. Run the DGP generation of the last subsection replacing

```stata
replace D = 1 if T == mod(G - 1, 5) + 3 & mod(G-1, 5) != 0
```

with 


```stata
replace D = 1 if T == mod(G - 1, 5) + 7 & mod(G-1, 5) != 0
```

Now, all the treated groups will have at least two periods before their first switch in the collapsed data, meaning that we can compute at least one placebo in all the subsamples.

Using the index adjustments that are described in detail in the [companion paper](https://github.com/chaisemartinPackages/did_multiplegt_dyn/blob/main/vignettes/assets/main.pdf), we can retrieve placebo estimators comparing the outcome evolution of switchers and not-yet-switchers from the period before the first treatment switch to 4, 5, 6 and 7 periods prior as follows:

```stata
scalar effects = 2
scalar placebo = 1
mat define res = J(4*placebo + 4 * effects, 7, .)
forv j=1/4 {
    did_multiplegt_dyn Y G T at_least_one_D_change ///
     if inlist(subsample, 0, `j'), effects(`=effects') placebo(`=placebo') graph_off
    forv i = 1/`=effects' {
        mat adj = mat_res_XX[`i',1..6]
        forv c =1/6 {
            mat res[`j'+(`i'-1)*4,`c'] = adj[1, `c']
        }
        mat res[`j'+(`i'-1)*4,7] =`j'+(`i'-1)*4
    }
    forv i = 1/`=placebo' {
        mat adj = mat_res_XX[effects + 1 + `i',1..6]
        forv c =1/6 {
            mat res[(`i'+1)*4 - `j' + (4*(effects - 1) +1),`c'] = adj[1, `c']
        }
        mat res[(`i'+1)*4 - `j' + (4*(effects - 1) +1),7] = - ((`i'+1)*4 - `j')
    }
}
mat li res
```

Notice that time stamps are now stored inside the **res** matrix since the relative time indexes are not consecutive anymore.
The graph can be produced with a similar procedure as the one described in previous section:

```stata
mat res = (0,0,0,0,0,0,0) \ res
svmat res
sort res7
tw rcap res3 res4 res7, lc(blue) || ///
connected res1 res7, mc(blue) lc(blue) || ///
, xtitle("Relative time to last period before treatment changes (t=0)") ///
title("DID, from last period before treatment changes (t=0) to t") ///
ytitle(" ") leg(off) xlabel(-7(1)8)
```

The output should look like this:
<p>
  <img src="https://github.com/chaisemartinPackages/did_multiplegt_dyn/blob/main/vignettes/assets/vignette_1_Stata_fig3.jpg" alt>
</p>
---


