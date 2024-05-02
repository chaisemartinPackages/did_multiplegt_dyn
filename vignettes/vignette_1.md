# Outcome observed less frequently than the treatment

In several circumstances, data about certain outcomes cannot be collected at every period where the related treatment is recorded. For instance, electoral outcomes are recorded only during election years. Let's take the case of a researcher wanting to estimate the dynamic effect of some variable or event on the electoral standing of some party or candidate using a DiD design. Treatments can be recorded more frequently than electoral data. As a result, even with yearly data, outcomes will be missing on regular intervals. 

Let's take the case of two parties, $A$ and $B$. Assume that the share of each party at some biannual election is the following:
| Party |Year |Share |
| A     |2003 |0.50  |
| A     |2005 |0.35  |
| B     |2003 |0.40  |
| B     |2005 |0.45  |

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

