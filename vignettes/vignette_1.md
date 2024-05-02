# Outcome observed less frequently than the treatment

## Part I: Data Generation 

<table style = "table-layoutfixed;width:800px;">
  <tr>
    <th>Stata</th>
    <th>R</th>
  </tr>
  <tr>
    <td>
    <pre><code style ="width:500px;">
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
    <pre><code style ="width:500px;">
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

