# Mathematical Justification for the `hdfe()` Option in `did_multiplegt_dyn`

## 1. Setup and Data Generating Process

Consider a panel data setting with:
- Groups $g \in \{1, \ldots, G\}$ (e.g., firms)
- Time periods $t \in \{1, \ldots, T\}$
- High-dimensional fixed effect categories $h(g) \in \{1, \ldots, H\}$ (e.g., industries)

The outcome is generated as:

$$Y_{g,t} = \alpha_g + \gamma_t + \beta_{h(g)} \cdot t + \tau \cdot D_{g,t} + \varepsilon_{g,t}$$

Where:
- $\alpha_g$ = group fixed effect
- $\gamma_t$ = time fixed effect
- $\beta_{h(g)} \cdot t$ = category-specific linear trend (slope $\beta_h$ for category $h$)
- $\tau$ = treatment effect (ATT)
- $D_{g,t}$ = treatment indicator

## 2. The Problem: Differential Trends Violate Parallel Trends

The parallel trends assumption requires that in the absence of treatment:

$$E[Y_{g,t}(0) - Y_{g,t-1}(0) | G=g] = E[Y_{g',t}(0) - Y_{g',t-1}(0) | G=g']$$

for all groups $g, g'$ being compared.

However, with category-specific trends:

$$\Delta Y_{g,t} = Y_{g,t} - Y_{g,t-1} = \gamma_t - \gamma_{t-1} + \beta_{h(g)} + \tau \cdot \Delta D_{g,t} + \Delta \varepsilon_{g,t}$$

The first difference contains $\beta_{h(g)}$, which varies across categories. If switchers and non-switchers belong to different categories (industries), this violates parallel trends.

## 3. The `hdfe()` Solution: Absorbing Category-Specific Trends

The `hdfe()` option follows the **exact same approach** as `controls()` to ensure numerical equivalence.

### Step 1: Run Weighted Regression Among Pre-Switch Observations

For each baseline treatment level $d$, run the regression:

$$\Delta Y_{g,t} = \sum_{s} \gamma_s \cdot \mathbf{1}[t=s] + \sum_{h} \theta_h \cdot \mathbf{1}[h(g)=h] + \epsilon_{g,t}$$

**Matching `controls()` exactly:**
- Sample: Pre-switch observations only ($t < F_g$)
- Weights: $N_{g,t}$ (cell weights)
- Includes: Time fixed effects ($\gamma_s$)
- Absorbs: HDFE category effects ($\theta_h$)

The absorbed coefficients $\hat{\theta}_h$ represent the category-specific means of first differences, controlling for time effects.

### Step 2: Extract HDFE Coefficients

The absorbed FE coefficients are computed as:

$$\hat{\theta}_h = \frac{\sum_{g: h(g)=h, t<F_g} N_{g,t} \cdot (\Delta Y_{g,t} - \hat{\gamma}_t)}{\sum_{g: h(g)=h, t<F_g} N_{g,t}}$$

This is the weighted mean of time-demeaned first differences for each category.

### Step 3: Adjust Outcome Levels for Long Differences

The core `did_multiplegt_dyn` program computes long differences:

$$\Delta^{(\ell)} Y_{g,t} = Y_{g,t} - Y_{g,t-\ell}$$

To ensure these also have category-specific trends removed, we adjust the outcome levels:

$$\widetilde{Y}_{g,t} = Y_{g,t} - \hat{\theta}_{h(g)} \cdot t$$

Then:

$$\widetilde{Y}_{g,t} - \widetilde{Y}_{g,t-\ell} = (Y_{g,t} - \hat{\theta}_h \cdot t) - (Y_{g,t-\ell} - \hat{\theta}_h \cdot (t-\ell))$$
$$= Y_{g,t} - Y_{g,t-\ell} - \hat{\theta}_h \cdot \ell$$

This removes the category-specific trend contribution $\beta_h \cdot \ell$ from the $\ell$-period difference.

## 4. Equivalence to `controls()` with Category×Time Dummies

### The `controls()` Approach

When using `controls(industry_1×time, industry_2×time, ..., industry_H×time)`, the command:

1. First-differences the control variables: $\Delta X^h_{g,t} = X^h_{g,t} - X^h_{g,t-1}$
2. Since $X^h_{g,t} = \mathbf{1}[h(g)=h] \cdot t$, we have $\Delta X^h_{g,t} = \mathbf{1}[h(g)=h]$
3. Regresses $\Delta Y_{g,t}$ on $\{\Delta X^h_{g,t}\}_h$ among controls
4. Residualizes $\Delta Y_{g,t}$ by subtracting fitted values

The coefficient on $\Delta X^h$ is exactly $\bar{\Delta Y}_h$ (the category mean).

### Mathematical Equivalence

Both approaches compute:

$$\widetilde{\Delta Y}_{g,t} = \Delta Y_{g,t} - \sum_h \hat{\beta}_h \cdot \mathbf{1}[h(g)=h]$$

Where $\hat{\beta}_h = \bar{\Delta Y}_h$ is the OLS coefficient, which equals the category mean when the only regressors are category indicators.

**Therefore, `hdfe(industry)` and `controls(industry×time)` are mathematically equivalent for removing category-specific linear trends.**

## 5. Advantages of `hdfe()` over `controls()`

| Aspect | `controls()` | `hdfe()` |
|--------|--------------|----------|
| **Implementation** | Create $H$ dummy variables × time | Use `reghdfe`/`areg` to absorb |
| **Memory** | $O(N \times H)$ | $O(N + H)$ |
| **Speed** | Slower with many categories | Faster via sparse matrix methods |
| **Variance** | Accounts for estimation error | Simpler variance (absorbed FE) |
| **Flexibility** | Multiple controls possible | Dedicated to FE absorption |

## 6. Key Assumptions

For the `hdfe()` option to produce unbiased estimates:

1. **Time-invariance**: The hdfe variable $h(g)$ must be constant within groups over time
2. **Linear trends**: The category-specific effect is linear in time ($\beta_h \cdot t$)
3. **Parallel trends within categories**: After removing category trends, parallel trends holds

## 7. Empirical Validation

From our test with True ATT = 3.0:

| Method | Effect_1 | Effect_2 | Effect_3 | Placebo_1 | Placebo_2 |
|--------|----------|----------|----------|-----------|-----------|
| Baseline | 3.264 | 3.124 | 3.348 | 0.103 | 0.289 |
| **controls(ind×t)** | **3.241929** | **3.125341** | **3.027178** | **0.128168** | **0.265299** |
| **hdfe(industry)** | **3.241929** | **3.125341** | **3.027179** | **0.128168** | **0.265299** |
| trends_nonparam | 3.237 | 3.031 | 3.039 | 0.115 | 0.284 |

### Numerical Equivalence Confirmed:
- **Effect_1**: `controls()` = `hdfe()` = 3.241929 ✓
- **Effect_2**: `controls()` = `hdfe()` = 3.125341 ✓
- **Effect_3**: Difference of 0.000001 (rounding) ✓
- **Placebos**: Identical to 6 decimal places ✓

The **exact numerical agreement** between `hdfe()` and `controls()` validates the mathematical equivalence.

## 8. References

- de Chaisemartin, C., & D'Haultfœuille, X. (2024). "Difference-in-Differences Estimators of Intertemporal Treatment Effects." *Review of Economics and Statistics*.
- Correia, S. (2016). "Linear Models with High-Dimensional Fixed Effects: An Efficient and Feasible Estimator." Working Paper.
- Frisch, R., & Waugh, F. V. (1933). "Partial Time Regressions as Compared with Individual Trends." *Econometrica*.
