# Performance of the `hdfe_controls()` patch

## Short version

**The patched ado is NOT faster than the OLD manual-workaround approach.** What it buys you is **convenience** (you no longer type the dummies by hand) and **exact match** with the existing estimator. It does **not** speed up the regression itself.

If you measured the cells test (48 dummies) and found it slow, that's because the original `did_multiplegt_dyn` is slow when `controls()` holds many variables — regardless of how those variables got there.

## Why

The patch's only role is the auto-expansion in lines L282–L370 of the ado:

```
hdfe_controls(i.ac_uq_id#i.monthyear)
  -> internally generates 48 byte dummies and appends them to `controls'
```

After that point the original code runs unchanged with `count_controls = 50` (`av_wind_speed + wind_direction + 48 dummies`). Three downstream bottlenecks scale with `count_controls`:

| Stage | Operation | Cost in `count_controls = K` |
|---|---|---|
| Stage 1, L1018 | `matrix accum overall_XX = ...` builds a $(K+1)\times(K+1)$ cross-product matrix | $O(K \cdot N)$ |
| Stage 1, L1034 | `invsym(didmgt_XX)` inverts the $K \times K$ controls block | $O(K^3)$ — for $K=50$, ~125k flops, trivial |
| Stage 2, L4679 | Long-difference residualization: `forvalues k = 1/K { replace diff_y_l -= theta * diff_X_k_l }` | $O(K \cdot N \cdot \text{horizons})$ |
| **Stage 3, L4921–4946** | **Nested `forvalues j=1/K { forvalues k=1/K { gen / replace } }`** | $O(K^2 \cdot \text{baselines} \cdot \text{horizons})$, each Stata command touching $N$ rows |

Stage 3 is the dominant cost. With $K = 50$ controls, 5 baseline cells, and 4 effect/placebo horizons:

$$50^2 \times 5 \times 4 = 50{,}000 \text{ Stata commands}, \text{ each running on } \sim 13{,}000 \text{ rows}.$$

That is several hundred million row-operations dispatched one Stata command at a time — *that* is what makes the cells test slow. The auto-expansion has nothing to do with it. The OLD manual-workaround path goes through exactly the same bottleneck.

## What would actually be faster

Two genuinely different paths to a speedup, each with a real cost:

1. **Vectorize Stage 3 in Mata** — rewrite the L4921–4946 nested loop as a single Mata matrix expression. Same math, no Stata `gen`/`replace` overhead. Estimated speedup for $K = 50$: 20–50×. Cost: a non-trivial Mata implementation that has to match the existing logic bit-for-bit (and play nicely with the cluster, less_conservative_se, etc. options). I'd estimate 1–2 days of careful work.

2. **True FE absorption** (the strategy in the first iteration we tried, before reverting) — use `reghdfe` for residualization in Stage 1 and *do not* materialize the dummies. Then $K$ stays at 2 (the genuine numeric controls), and Stage 3 only loops $2^2 = 4$ times. Estimated speedup: orders of magnitude. **But** this requires also extending Stage 2 (L4679) to subtract the absorbed-FE part of the long difference, otherwise NEW will diverge from OLD as it did in the earlier test (max |diff| ≈ 0.29 on Av_tot_eff). Doing that subtraction analytically without materializing the dummies is non-trivial — it involves recovering the per-cell fixed-effect coefficients from `reghdfe` and constructing their long-difference for each event horizon. Estimated cost: a longer derivation + careful code. 3–5 days.

## What to do for now

- Use the patched ado as is. It gives correct estimates and saves typing.
- If a real workload starts costing minutes, that is the signal to commit to path (1) or (2).
- Path (1) is the lower-risk option: same math, just faster Stage 3. Path (2) is more powerful but requires opening the variance derivation. We can pick the one that fits the workload when the time comes.
