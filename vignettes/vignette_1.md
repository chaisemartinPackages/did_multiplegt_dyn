# Outcome observed less frequently than the treatment

## Part I: Data Generation 

## Example 1

### R

<pre><code class="language-r">set.seed(123)
TT <- 20; GG <- 5
df <- data.frame(id = 1:(GG*TT))
df$G <- ((df$id-1) %% GG)+1
df$T <- floor((df$id-1)/GG)
df$id <- NULL
df <- df[order(df$G, df$T), ]
</code></pre>

<div class="md-fenced-code-tabs" id="tab-tab-group-4">
    <input name="tab-group-4" type="radio" id="tab-group-4-0_c" checked="checked" class="code-tab" data-lang="c" aria-controls="tab-group-4-0_c-panel" role="tab">
    <label for="tab-group-4-0_c" class="code-tab-label" data-lang="c" id="tab-group-4-0_c-label">C</label>
    <div class="code-tabpanel" role="tabpanel" data-lang="c" id="tab-group-4-0_c-panel" aria-labelledby="tab-group-4-0_c-label">
        ... the highlighted syntax ...
    </div>
    <input name="tab-group-4" type="radio" id="tab-group-4-1_java" class="code-tab" data-lang="java" aria-controls="tab-group-4-1_java-panel" role="tab">
    <label for="tab-group-4-1_java" class="code-tab-label" data-lang="java" id="tab-group-4-1_java-label">Java</label>
    <div class="code-tabpanel" role="tabpanel" data-lang="java" id="tab-group-4-1_java-panel" aria-labelledby="tab-group-4-1_java-label">
        ... the highlighted syntax ...
    </div>
</div>

### Stata
```{.applescript}

```
