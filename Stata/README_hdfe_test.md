# Testing the `hdfe_controls()` patch

## Files

- **`did_multiplegt_dyn.ado`** — patched ado with the new `hdfe_controls()` option.
- **`did_multiplegt_dyn.ado.backup`** — pristine copy of the original. Keep this; you can restore by `cp did_multiplegt_dyn.ado.backup did_multiplegt_dyn.ado`.
- **`test_hdfe_controls.do`** — Stata script that runs OLD vs NEW on a 5-cluster x 24-month slice of `data_test.csv` and prints an element-wise comparison.

## What was patched in the ado

1. **L53 (syntax)** — added `hdfe_controls(string)` to the option list.
2. **After L73 (dependencies)** — soft check for `reghdfe` and `hdfe` when the option is requested.
3. **L919–961 (Stage 1 demeaning loop)** — wrapped in `if "`hdfe_controls'" == ""` (original gegen path) `else` (one `hdfe` call per baseline cell `d`, demeaning ΔY-style controls by `time_XX` + `trends_nonparam`×time + user's `hdfe_controls`).
4. **L984–996 (Stage 1 fitted-value regression)** — added an `else` branch that uses `reghdfe ..., absorb(...) residuals(_res_DY_l_XX)`. The downstream `predict` block at L1092 was made conditional too: it uses the OLD `predict` when `hdfe_controls` is empty, and rebuilds the fitted value as `ΔY − _res_DY_l_XX` otherwise.

Stages 2 and 3 (L4582+ and L4905+) are unchanged. They consume the rebuilt Stage-1 outputs.

## How to run

On your machine:

1. Place the **patched** `did_multiplegt_dyn.ado` in a directory (e.g., the same folder as your data).
2. In `test_hdfe_controls.do`, set:
   - `local data_path` — absolute path to `data_test.csv`.
   - `local ado_path`  — directory holding the patched ado.
   - `local Y_var`     — name of your outcome variable. **Required**, no default.
3. From Stata: `do test_hdfe_controls.do`.

## What the test does

1. Imports `data_test.csv`, keeps the 5 smallest `ac_uq_id` values, keeps the first 24 distinct `monthyear` values, drops rows missing any key variable.
2. Builds 4 manual workaround controls: dummies for the first 4 `ac_uq_id` values, each multiplied by `monthyear` — the "include `Z × t`" trick from the help file. With FD they become group-level cluster dummies.
3. Runs the OLD path:
   ```stata
   did_multiplegt_dyn Y_var unique_small_grid_id monthyear downup_ac,
       effects(4) placebo(2)
       controls(av_wind_speed wind_direction <4 cluster*monthyear dummies>)
       cluster(ac_uq_id) graph_off _no_updates
   ```
4. Runs the NEW path:
   ```stata
   did_multiplegt_dyn Y_var unique_small_grid_id monthyear downup_ac,
       effects(4) placebo(2)
       controls(av_wind_speed wind_direction)
       hdfe_controls(i.ac_uq_id#c.monthyear)
       cluster(ac_uq_id) graph_off _no_updates
   ```
5. Stores `e(b)` and `e(V)` from each, prints element-wise differences and `mreldif(b_old, b_new)`.

## What to expect

- **Stage 1 θ_d on `av_wind_speed` and `wind_direction`** is identical between OLD and NEW by Frisch–Waugh–Lovell. Adding redundant dummy columns vs. partialling them out give the same slope on the genuine numeric controls.

- **The event-study coefficients in `e(b)`** may match exactly, or may show a small discrepancy. The discrepancy, if any, comes from Stage 2's L4679: the OLD path subtracts `θ_a · ℓ · D_a` from each long-difference for every cluster dummy `a` (because they are numeric controls and L4679 loops over all of them), while the NEW path subtracts only the genuine-control contribution (because the cluster trends were absorbed in Stage 1's projector and Stage 2 only iterates over numeric controls).

  In data where the cluster distribution is balanced across switchers and not-yet-switchers within each baseline cell, the cluster trend cancels in the event-study comparison and the two paths match. In imbalanced data, a residual cluster-trend term survives and the two paths differ.

- If they differ and you want them to match exactly, the next iteration is to extend Stage 2 so that the long-difference of the absorbed FE part is also subtracted. That is doable but I want to see the test output first before committing to it.

## Restoring the original

```bash
cp /workspace/did_multiplegt_dyn.ado.backup /workspace/did_multiplegt_dyn.ado
```
