# Outcome observed less frequently than the treatment

In several circumstances, data about certain outcomes can only be collected discountinously. For instance, electoral outcomes are recorded only during election years.This may often contrast with the data availability of the treatments. Instantaneous treatments can be staggerized and imputed to subsequent periods. Continuous treatments can be monitored at very disaggregated units of time. In general, it is possible that the treatment is observed more frequently than the outcome. As a result, outcomes will be missing on regular intervals. 

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

The treatment variable can be defined as $1\lbrace\text{Year} \geq \text{LeadChangeYr}, \text{Party} = A\rbrace$. The effect of leadership change after $\ell$ years can be estimated by comparing the 2003-to-(2003+ $\ell$) change in the shares of treated party $A$ with the 2003-to-(2003+ $\ell$) change in the shares of control party $B$. Notice that the year of the first treatment change (2004) is not in the data. Hence, we cannot compute the first dynamic effect. The same holds true for 2006 and 2008. With the dataset above, it is possible to estimate only the second and fourth (2007) event-study effect. 

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

The treatment variable is now defined $1\lbrace\text{Year} \geq \text{LeadChangeYr}, \text{Party} \in (A, C)\rbrace$. As before, we can make two actual-versus-status quo comparisons for party $C$ against $B$, i.e. 2004-to-2005 and 2004-to-2007. However, we do not observe the outcome in 2004. This is where the parallel trends assumption comes handy. Let $Y^i_t$ be the outcome of party $i$ in year $t$. For further clarity, we denote unobservable outcomes as $\mathring{Y}^i_t$. The expectation of the 2004-to-2005 comparison can be decomposed as follows:
$$
\begin{array}{ll}
& E[(Y^C_{2005} - \mathring{Y}^C_{2004}) - (Y^B_{2005} - \mathring{Y}^B_{2004})] \\
= & E[(Y^C_{2005} - Y^C_{2003} +  Y^C_{2003} - \mathring{Y}^C_{2004}) - (Y^B_{2005}  - Y^B_{2003} +  Y^B_{2003} - \mathring{Y}^B_{2004})]  \\
= & E[(Y^C_{2005} - Y^C_{2003}) -  (Y^B_{2005}  - Y^B_{2003})] - E[(\mathring{Y}^C_{2004} - Y^C_{2003})  - (\mathring{Y}^B_{2004} - Y^B_{2003})] \\
= & E[(Y^C_{2005} - Y^C_{2003}) -  (Y^B_{2005}  - Y^B_{2003})] \\
\end{array}
$$
where in the last equality we use the assumption that all the groups experience the same untreated outcome evolution (party $C$ gets treated only in 2005, while group $B$ remains always untreated). As a result, we can estimate any dynamic effect using the most recent non-missing untreated outcome as the status quo outcome as long as the corresponding actual outcome is non-missing. This allows us to use group $C$ to estimate all the other missing dynamic effects from the previous dataset. Specifically, the 2004-to-2005 effect will correspond to the first dynamic effect for party $C$ and the 2004-to-2007 effect to the third dynamic effect.

Following this method, we can estimate up to 4 dynamic effects with our data. Figure 1 displays the combined event-study plot from our toy example.
![vignette_1_Stata_fig1.jpg]()



## General Case with Stata and R code

### Part I: Data Generation 

<table>
  <tr>
    <th>Stata</th>
    <th>R</th>
  </tr>
  <tr>
    <td>
    <pre><code>
    clear
    set seed 0
    local TT = 20
    local GG = 5
    set obs `=`TT' * `GG''
    gen G = mod(_n-1,`GG') + 1
    gen T = floor((_n-1)/`GG')
    sort G T
    </pre></code>
    </td>
    <td>
    <pre><code>
      set.seed(0)
      TT <- 20; GG <- 5
      df <- data.frame(id = 1:(GG*TT))
      df$G <- ((df$id-1) %% GG)+1
      df$T <- floor((df$id-1)/GG)
      df$id <- NULL
      df <- df[order(df$G, df$T), ]
    </pre></code>
    </td>
  </tr>
</table>

