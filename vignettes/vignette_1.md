# Outcome observed less frequently than the treatment

In several circumstances, data about certain outcomes can only be collected discountinously. For instance, electoral outcomes are recorded only during election years. This may often contrast with the data availability of the treatments. Instantaneous treatments can be staggerized and imputed to subsequent periods. Continuous treatments can be monitored at very disaggregated units of time. In general, it is possible that the treatment is observed more frequently than the outcome. In this tutorial, we show how to use of the estimator from de Chaisemartin & D'Haultfoeuille (2024) with outcomes periodically missing. 

#### Sections
+ [A toy example](#a-toy-example)
+ [General Case with Stata and R code](#general-case-with-stata-and-r-code)

## A toy example

Let's take the case of a researcher wanting to estimate the effect of the insediation of a new party leader on electoral performance. Party $A$ has changed its main candidate in 2004. The natural design for this application is an event-study at the year level. Assume that we have biannual election data as follows:

| Party |Year |Share |LeadChangeYr|Treatment|
| ----: |----:|----: |----:       |----:    |
| A     |2003 |0.40  |2004        |0        |
| A     |2005 |0.35  |2004        |1        |
| A     |2007 |0.25  |2004        |1        |
| B     |2003 |0.40  |.           |0        |
| B     |2005 |0.45  |.           |0        |
| B     |2007 |0.50  |.           |0        |

The treatment variable can be defined as $1\lbrace\text{Year} \geq \text{LeadChangeYr}, \text{Party} = A\rbrace$. The effect of leadership change after $\ell$ years can be estimated by comparing the 2003-to-(2003+ $\ell$) change in the shares of treated party $A$ with the 2003-to-(2003+ $\ell$) change in the shares of control party $B.$ Notice that the year of the first treatment change (2004) is not in the data. Hence, we cannot compute the first dynamic effect. The same holds true for 2006 and 2008. With the dataset above, it is possible to estimate only the second and fourth (2007) event-study effect. 

Assume now that the dataset includes another party, $C$, that has experienced a leadership change in 2005:

| Party |Year |Share |LeadChangeYr|Treatment|
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

The treatment variable is now defined $1\lbrace\text{Year} \geq \text{LeadChangeYr}, \text{Party} \in (A, C)\rbrace$. As before, we can make two actual-versus-status quo comparisons for party $C$ against $B$, i.e. 2004-to-2005 and 2004-to-2007. However, we do not observe the outcome in 2004. This is where the parallel trends assumption comes handy. Let $Y^i_t$ be the outcome of party $i \in \lbrace A,B,C\rbrace$ in year $t$. For further clarity, we denote unobservable outcomes as $\mathring{Y}^i_t.$ The expectation of the 2004-to-2005 comparison can be decomposed as follows:

$$
\begin{align*}
E[(Y^C_{2005} - \mathring{Y}^C_{2004}) - (Y^B_{2005} - \mathring{Y}^B_{2004})] &= \\
E[(Y^C_{2005} - Y^C_{2003} +  Y^C_{2003} - \mathring{Y}^C_{2004}) - (Y^B_{2005}  - Y^B_{2003} +  Y^B_{2003} - \mathring{Y}^B_{2004})]  &= \\
E[(Y^C_{2005} - Y^C_{2003}) -  (Y^B_{2005}  - Y^B_{2003})] - E[(\mathring{Y}^C_{2004} - Y^C_{2003})  - (\mathring{Y}^B_{2004} - Y^B_{2003})] &= \\
E[(Y^C_{2005} - Y^C_{2003}) -  (Y^B_{2005}  - Y^B_{2003})] \\\
\end{align*}
$$

where in the last equality we use the assumption that all the groups experience the same untreated outcome evolution (party $C$ gets treated in 2005, while party $B$ remains always untreated). As a result, we can estimate any dynamic effect using the most recent non-missing untreated outcome as the status quo outcome as long as the corresponding actual outcome is non-missing. This allows us to use group $C$ to estimate all the other missing dynamic effects from the previous dataset. Specifically, the 2004-to-2005 effect will correspond to the first dynamic effect for party $C$ and the 2004-to-2007 effect to the third dynamic effect.

Following this method, we can estimate up to 4 dynamic effects with our data. The figure below shows the combined event-study plot from our toy example.

<p>
  <img src="https://github.com/DiegoCiccia/did_multiplegt_dyn/blob/main/vignettes/assets/vignette_1_Stata_fig1.jpg" alt>
</p>

In the next section, we show step-by-step how to adapt the method described above to any case whereby the outcome is missing at regular intervals. For each step, we also report the correspondent block of the full length [Stata](https://github.com/DiegoCiccia/did_multiplegt_dyn/blob/main/vignettes/src/vignette_1_Stata.do) and [R](https://github.com/DiegoCiccia/did_multiplegt_dyn/blob/main/vignettes/src/vignette_1_R.R) code.

## General Case with Stata and R code

### Part I: Data Generation 

> [!IMPORTANT]
> If you wish to use this tutorial on some data, make sure to rename your variables of interest as follows:
> + Outcome variable $\to$ "Y"
> + Group variable $\to$ "G"
> + Time variable $\to$ "T"
> + Treatment variable $\to$ "D"
> Also, your group and time variables should be mapped to positive integers (in Stata, `egen group`, in R, `cur_group_id()`).
>
> Then, you can skip to [Part II](#part-ii-data-adjustment).

We use a DGP with five groups and 20 periods. We can observe the outcome every fourth period. It is enough to replace 4 with any other interval to generalize the code below. Group 1 never switches, while groups 2 to 5 switch at periods 4 to 7 (group id + 2). The untreated outcome is uniformly distributed and increases by a factor of 100 after the first treatment switch.

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
We need to generate partition the $(g,t)$ cells in our data according to whether and when they switch treatment. Keep in mind that the outcome is observable every 4 periods. We need to generate a variable at the (g,t) level through which we can earmark not-yet-switchers (here including never-switchers) and switchers in either a non-missing year or 1, 2, or 3 periods before a non-missing year. We will use the values of this variable to run `did_multiplegt_dyn` on subsets of our data. To define this partition, we can move along the following steps.

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

As in the DGP description, the group 2 has at least one treatment change starting from period 4 (the period where its first treatment change occurred), group 3 from period 5 and so on.

#### b) Identify when treatment changed for the first time

We generate a variable (F_g) equal to the earliest period when there has been a treatment change. For never switchers we use the same convention as de Chaisemartin & D'Haultfoeuille (2024) and set F_g = number of periods + 1. Let's use the modulus operator to check whether the F_g falls on a 4th period or 1/2/3 years after. Each year where the outcome is not missing takes value 1. Since the outcome is non-missing every 4 years, we just check whether the period is divisible by 4. If yes, the modulus operator by 4 yields 0. If not, it yield the remainder from dividing by 4. This remainder is always between 1 and 3. So, it is enough to increase the remainder by 1 to obtain 4 distinct values such that the corresponding observations with non missing outcome (every fourth period) have value 1, observations after them take value 2 and so on. Lastly, we multiply the subset variable by the indicator for at least one change in D to obtain the partition variable.

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
    gen subsample = mod(F_g, 4) + 1
    gen model_subset = subsample * at_least_one_D_change
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
Once the dataset has only non missing outcomes, we can run `did_multiplegt_dyn` a number of times equal to the interval in the original data between non-missing outcomes (4 in our case). For the $n$-th run, we keep a subset of our data correspoding to observations having model_subset equal to 0 (not-yet-switchers) or $n$ (switchers in the $n$-th partition). In this way, we are able to redefine the indices of the effects. For instance, in groups with $n = 1$ (i.e. switching on a year where the outcome can be observed), the first effect computed by `did_multiplegt_dyn` matches their first dynamic effect, but the second matches their fifth dynamic effect (consider the case of party $C$ in the [toy example](#a-toy-example) above).

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
        did_multiplegt_dyn Y G T at_least_one_D_change ///
         if inlist(model_subset, 0, `j'), effects(`effects') graph_off
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



