# Outcome observed less frequently than the treatment

In many circumstances, the outcome is observed less frequently than the treatment. For instance, electoral outcomes are recorded only during election periods, while treatment may be observed every period. In this tutorial, we show how to use of the estimator from de Chaisemartin & D'Haultfoeuille (2024) with outcomes observed less frequently than the treatment. 

#### Sections
+ [A toy example](#a-toy-example)
+ [General Case with Stata and R code](#general-case-with-stata-and-r-code)

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
  <img src="https://github.com/DiegoCiccia/did_multiplegt_dyn/blob/main/vignettes/assets/vignette_1_Stata_fig1.jpg" alt>
</p>

In the next section, we show how to use `did_multiplegt_dyn` to produce an event-study graph without averaging effects of different exposure lengths whenever the outcome is observed less frequently than the treatment. For each step, we also report the corresponding block of [Stata](https://github.com/DiegoCiccia/did_multiplegt_dyn/blob/main/vignettes/src/vignette_1_Stata.do) and [R](https://github.com/DiegoCiccia/did_multiplegt_dyn/blob/main/vignettes/src/vignette_1_R.R) code.

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
We need to generate a variable through which we partition our population into 5 subsamples: 
+ A) Groups whose treatment never changes (never switchers)
+ B) Groups whose treatment changes for the first time at a time period where the outcome is non-missing
+ C) Groups whose treatment changes for the first time one period before a period where the outcome is non-missing
+ D) Groups whose treatment changes for the first time two periods before a period where the outcome is non-missing
+ E) Groups whose treatment changes for the first time three periods before a period where the outcome is non-missing
  
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

We generate a variable F_g equal to the period when group g's treatment changes for the first time. The modulus of F_g divided by 4 is equal to zero is g's treatment changed for the first time in a period when the outcome is non missing, it is equal to 1 (resp. 2, 3) if g's treatmet changed for the first time 3 (resp. 2, 1) periods before a period when the outcome is non missing. We use this modulus to create the aforementioned partition variable. 

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
gen subsample_temp = mod(F_g, 4)
gen subsample=4-subsample_temp
replace subsample=subsample+1
replace subsample=0 if never_treated ==1 
    </pre></code>
    </td>
    <td>
    <pre><code>
df <- df %>% group_by(.data$G) %>% 
mutate(never_treated = as.numeric(sum(.data$D_change, na.rm = TRUE) == 0)) %>%
mutate(F_g = ifelse(.data$never_treated == 1, max(df$T, na.rm = TRUE) +1, 
  min(ifelse(.data$D_change == 0, NA, .data$T * .data$D_change), na.rm = TRUE))) %>% 
    ungroup()
df$subsample <- (df$F_g %% 4) + 1
df$model_subset <- df$subsample * df$at_least_one_D_change
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
    did_multiplegt_dyn Y G T at_least_one_D_change if inlist(model_subset, 0, `j'), effects(`effects') graph_off
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
      subset(df, df$model_subset %in% c(0, j)), "Y", "G", "T", "at_least_one_D_change", 
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



