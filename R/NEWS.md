# DIDmultiplegtDYN 2.4.0

## New options

* `reset` — port of the Stata `reset(integer)` option. Balances the panel to
  groups with the maximum number of non-missing treatment observations, then
  splits each group into sub-groups whenever its treatment has remained
  constant for `reset` consecutive periods after a change. When no `cluster`
  variable is supplied, the original group identifier is used as the
  clustering level.

* `avg_time_periods` — port of the Stata `avg_time_periods` option. Reports
  the average number of time periods over which the effect of a treatment
  dose is accumulated. The returned object now exposes top-level scalars
  `avg_cumul`, `N_switchers_effect_k`, and `N_switch_avg_k` mirroring the
  Stata `e()` return.

## Performance

* Ported the Mata routines `_compute_Tg` and `_avg_cumul` to C++ via Rcpp
  (`compute_Tg_cpp`, `avg_cumul_cpp`). Inner control-search lookups use
  precomputed `(group, time) -> row_index` and `cls -> [group_idx]` maps,
  so the cost drops from `O(G^2 * ell * T)` to `O(G * ell)` with `O(1)`
  candidate-control lookups.

## Argument-check fixes

* All `inherits(x, "numeric")` checks replaced with `is.numeric(x)` so
  integer literals (e.g. `effects = 5L`, `bootstrap = 100L`) no longer
  fail the type check.
* Bootstrap reps are now validated as integers.
* `bootstrap` accepts a numeric scalar (`bootstrap = 100`), a numeric
  vector (`bootstrap = c(100, 42)`), or a named vector
  (`bootstrap = c(reps = 100, seed = 42)`). The legacy
  `list(reps, seed)` form is still accepted for backward compatibility.
* `by_path` now accepts `"all"` in addition to a positive integer or `-1`.

## Internal

* Dropped the unused `matlib::Ginv` import, which removed the transitive
  `rgl` / OpenGL dependency.
* Removed leftover unused variables in `src/did_loops.cpp` that triggered
  `-Wunused-variable` warnings.

# DIDmultiplegtDYN 2.3.4

* Previous release.
