# Outcome observed less frequently than the treatment

## Part I: Data Generation 

<div style="background-color: rgb(50, 50, 50);">
```r
set.seed(123)
TT <- 20; GG <- 5
df <- data.frame(id = 1:(GG*TT))
df$G <- ((df$id-1) %% GG)+1
df$T <- floor((df$id-1)/GG)
df$id <- NULL
df <- df[order(df$G, df$T), ]
```
</div>
