# Outcome observed less frequently than the treatment

In several circumstances, data about certain outcomes can only be collected discountinously. For instance, electoral outcomes are recorded only during election years.This may often contrast with the data availability of the treatments. Instantaneous treatments can be staggerized and imputed to subsequent periods. Continuous treatments can be monitored at very disaggregated units of time. In general, it is possible that the treatment is observed more frequently than the outcome. As a result, outcomes will be missing on regular intervals. 

Let's take the case of a researcher wanting to estimate the effect of the insediation of a new party leader on electoral performance. Party $A$ has changed its main candidate in 2004. The natural design for this application is an event-study at the year level. The treatment variable can be defined as $1 \{\text{Year} \geq 2004, \text{Party} = A\}$. The effect of leadership change after $\ell$ years can be stimated by comparing the 2003-to-(2003+$\ell$) change in the shares of treated party $A$ with the 2003-to-(2003+$\ell$) change in the shares of control party $B$. Assume that we have biannual election data as follows:

| Party |Year |Share |Treatment|
| ----: |----:|----: |----:    |
| A     |2003 |0.40  |0        |
| A     |2005 |0.35  |1        |
| A     |2007 |0.25  |1        |
| B     |2003 |0.40  |0        |
| B     |2005 |0.45  |0        |
| B     |2007 |0.50  |0        |

Notice that the year of the treatment change (2004) is not in the data. Hence, we cannot compute the first dynamic effect. The same holds true for 2006 and 2008. With the dataset above, it is possible to estimate only the second and fourth (2007) event-study effect. 



## Part I: Data Generation 

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

