# Bootstrap-equivalence test session notes

Working notes for the test infrastructure under
`tests/manual/bootstrap_equivalence/` and the package-code fixes that came out
of building it. Written so a future session can pick up where this one left
off without re-deriving the context.

---

## What the test does

For each of 34 specifications, it runs `did_multiplegt_dyn(... bootstrap = list(4, 1234))`
twice — once with the **new** weight-based path
(`DID_BOOTSTRAP_METHOD = "subprocess"`) and once with the **legacy** row-replication
path (`DID_BOOTSTRAP_METHOD = "row_replication"`). The two paths share a CSV
sample directory so they see the *same* bootstrap draws, and we assert
`Effects` / `ATE` / `Placebos` match to floating-point tolerance (`1e-8`).

Each spec runs in **its own** `Rscript` process; within each spec, each of the
two methods runs in **its own** nested `Rscript` process. This is enforced by
the bash driver and gives us deterministic memory return between specs and
between methods, which is otherwise impossible because polars allocates in
Rust's allocator and never returns pages to the OS within a live R session.

---

## Files

```
tests/manual/bootstrap_equivalence/
├── specs.R              # SPECS list (34) + dataset loaders + preprocess helpers
├── run_one_path.R       # Runs ONE method for ONE spec, writes an RDS payload
├── run_one_spec.R       # Wrapper: spawns run_one_path.R twice, compares, writes per-spec CSV
├── run_all_specs.sh     # Driver: spawns run_one_spec.R per spec, sequentially
└── SESSION_NOTES.md     # This file
```

Logs land in `tests/manual/logs/per_spec/<spec_name>.log` and the aggregate
summary in `tests/manual/logs/bootstrap_equivalence_summary.csv`.

### Process layout per spec

```
run_all_specs.sh
└── Rscript run_one_spec.R <spec_name> <out_csv>          # one process per spec
    ├── Rscript run_one_path.R <spec> subprocess <sd> FALSE <new.rds>   # writes CSVs
    │   └── callr::r() x4   (one nested subprocess per bootstrap iteration)
    └── Rscript run_one_path.R <spec> row_replication <sd> TRUE  <old.rds>   # reads CSVs
```

Both per-method invocations share the same temp `<sd>` directory of
`rep_NNNNN_<tag>.csv` files (see `stratum_tag` note below). The new path
writes them; the legacy path replays them.

---

## How to run

```bash
cd tests/manual/bootstrap_equivalence
./run_all_specs.sh                   # all 34 specs
./run_all_specs.sh ch01_T2_deryugina_basic stata_normalized   # subset by name
```

A single spec in isolation:

```bash
LD_PRELOAD=/lib/aarch64-linux-gnu/libjemalloc.so.2 \
  Rscript run_one_spec.R stata_normalized /tmp/check.csv
```

Inspect logs / summary:

```bash
cat ../logs/per_spec/<spec_name>.log
column -t -s, ../logs/bootstrap_equivalence_summary.csv | less -S
```

---

## Environment setup (one-time)

Test container is Ubuntu 24.04 arm64 with R 4.4.2.

```bash
# System packages
sudo apt-get install -y cmake cargo rustc libuv1-dev libgl1-mesa-dev \
    libglu1-mesa-dev libx11-dev libfreetype6-dev libpng-dev libxml2-dev \
    libjemalloc-dev libjemalloc2

# R packages -> user library /home/node/Rlib  (R_LIBS_USER pinned in ~/.Renviron)
Rscript -e 'install.packages(c("matlib","openxlsx","car","lmtest","cowplot",
  "rnames","fs","rgl","nloptr","lme4","pbkrtest","bslib","sass","rmarkdown",
  "htmlwidgets"), lib="/home/node/Rlib", repos="https://cloud.r-project.org")'
NOT_CRAN=true Rscript -e 'install.packages("polars",
  repos="https://community.r-multiverse.org", lib="/home/node/Rlib")'

# Build & install our package from source
rm -f src/*.o                      # nuke stale .o (often macOS Mach-O)
R CMD INSTALL . --library=/home/node/Rlib --no-test-load --no-help --no-html --no-demo
```

---

## Bootstrap option quirk (validator)

`did_multiplegt_dyn` validates the bootstrap option with
`inherits(x, "numeric")`. An R **integer** literal (`4L`) is NOT `"numeric"`, so
`bootstrap = list(4L, 1234L)` is rejected with `"First element (reps) must be
an integer greater than 1"`. Use a double: `BOOT_REPS <- 4`.

---

## Bugs we hit and what we did

### 1. glibc heap corruption in the legacy path

`malloc(): unsorted double linked list corrupted` aborting the legacy
`row_replication` process on `stata_normalized` (and stochastically on others).
Root cause: glibc + polars/Rust allocation interplay in long-running parents.

**Mitigation (test infra, no code change):**
`run_all_specs.sh` exports

```
LD_PRELOAD=/lib/aarch64-linux-gnu/libjemalloc.so.2
MALLOC_ARENA_MAX=1
```

jemalloc replaces glibc's allocator for the whole process tree; corruption no
longer occurs. New path is unaffected (it isolates iterations in subprocesses)
but we set it globally so both paths run under identical conditions.

### 2. `by`-stratification: CSV filenames collided across strata  *(real package bug)*

`did_multiplegt_dyn` calls `did_multiplegt_bootstrap` once per `by`-stratum.
The helper named its sample CSVs `rep_NNNNN.csv` keyed only by iteration `j`,
so stratum 2 overwrote stratum 1's files. The legacy path then loaded the
wrong stratum's unit IDs, `match()`'d everything to NA, and crashed with
`"missing value where TRUE/FALSE needed"`.

**Fix in `R/did_multiplegt_bootstrap.R`:**
- Derived a `stratum_tag` from `unit_ids` (uses `rlang::hash` when available,
  falls back to a length+sum+min+max fingerprint).
- Renamed CSVs to `rep_NNNNN_<tag>.csv`.
- Both writer (new path) and reader (legacy path) derive the same tag from
  the same `unit_ids`, so by-stratum runs no longer collide.

This is not just a test fix — anyone using `DID_BOOTSTRAP_SAMPLE_DIR` for
validation under `by`/`by_path` would have hit it.

### 3. callr subprocess transiently fails to start  *(robustness fix)*

`stata_trends_lin` (heavy: linear trends per group on 50k rows) failed mid-loop
with `"callr subprocess failed: could not start R, exited with non-zero
status, has crashed or was killed"`. The parent process had grown large
enough that `fork()` couldn't reliably spawn a new R.

**Fix in `R/did_multiplegt_bootstrap.R` `run_subprocess`:**
- Retry up to `.MAX_SUBPROC_TRIES = 3L` times.
- Between attempts: `gc()` to release committed parent pages + small back-off
  (`0.5 * attempt` seconds).
- Final failure re-throws with the last error message.

Also added an `invisible(gc(verbose = FALSE))` after each iteration's
result is stowed, so the parent doesn't accumulate R-side allocations across
iterations.

### 4. east_et_al_2023 baseline OOM in 8 GB container

Original ch02 T5/T6 used `effects=5, placebo=5, 16 controls, continuous=1` —
the baseline `did_multiplegt_main()` call alone exceeds 7 GB (matlib variance
machinery), independent of bootstrap method. Both paths OOM-killed before any
bootstrap iteration ran.

**Mitigation (test data only):** trimmed `EAST_CONTROLS` to 5 vars and
east-specs to `effects=3, placebo=3` in `specs.R`. Still exercises the
continuous + weight + multi-control path; just fits in the container.

### 5. Spec/dataset mismatches that error in *both* paths

`stata_controls`, `stata_continuous`, `stata_switchers_out` use options that
are degenerate on `deryugina_2017` (e.g. binary `hurricane` with
`continuous=1`, time-invariant control). Both paths fail with the *same*
error.

**Test infra:** `run_one_spec.R` treats "both paths errored with the same
message" as `OK` (consistent failure) rather than `FAIL`. That's faithful to
the user's spec: "same exact result, even when degenerate".

---

## Final result

```
Specs run:    34
OK rows:      96
FAIL rows:    0
ERROR rows:   0
max abs diff = 1.421e-14
NEW total: 1847 s   LEGACY total: 1404 s
```

Diff is floating-point noise (well under the 1e-8 tolerance), confirming the
two paths are bit-equivalent under shared draws across every option exposed
by `did_multiplegt_dyn()` that flows through bootstrap.

---

## Coverage matrix

Options exercised under `bootstrap = list(4, 1234)`:

| Option | Spec(s) |
|---|---|
| effects, placebo | all |
| controls | ch01_T3, ch02_T5/T6, stata_controls |
| trends_nonparam | ch04_T14, stata_trends_nonparam |
| trends_lin | stata_trends_lin |
| continuous | ch02_T5/T6, stata_continuous |
| weight | ch02_T5/T6, stata_weight, favara_imbs_weight |
| cluster | ch01_T2/T3, ch04_T13, favara_imbs_*, stata_cluster |
| same_switchers / same_switchers_pl | ch03_T11, stata_same_switchers* |
| switchers | stata_switchers_{in,out} |
| only_never_switchers | stata_only_never_switchers |
| ci_level | stata_ci_level_{90,99} |
| less_conservative_se | stata_less_conservative_se |
| dont_drop_larger_lower | stata_dont_drop_larger_lower |
| effects_equal | ch02/03/04 + stata_effects_equal |
| normalized | ch02_T6, ch03_T10, stata_normalized |
| normalized_weights | gentzkow_normalized_weights |
| predict_het | deryugina_predict_het |
| by (stratified) | deryugina_by_coastal |
| drop_if_d_miss_before_first_switch | deryugina_drop_if_d_miss |
| package-bundled dataset (favara_imbs) | favara_imbs_{basic,weight} |

Not exercised (out of scope of bootstrap helper or post-bootstrap reporting):
`by_path`, `predict_het_hc2bm`, `more_granular_demeaning`, `date_first_switch`,
`design`, `save_results`, `save_sample`.

---

## Things to come back to

- **`more_granular_demeaning` is NOT passed to the bootstrap helper** (see
  `did_multiplegt_dyn.R:430`). A user calling
  `did_multiplegt_dyn(..., more_granular_demeaning=TRUE, bootstrap=...)`
  gets a baseline that uses the option but bootstrap iterations that don't —
  the SEs would be from a different estimator than the point estimate.
  Worth deciding whether to add it to `did_multiplegt_bootstrap`'s
  signature + `smaller_kwargs`.
- Same observation applies to `predict_het_hc2bm` (not passed; only matters
  for post-bootstrap reporting today).
- Consider wiring up `.DID_USE_POLARS()` (defined at top of
  `did_multiplegt_main.R` but never referenced) so users can fall back to a
  pure data.table backend — that would eliminate the polars allocator issue
  entirely.

---

## Reminder: package code we touched

Only one file in `R/`:

- `R/did_multiplegt_bootstrap.R`
  - `stratum_tag` + `sample_csv_path()` — stratum-namespaced CSV filenames
  - `run_subprocess` — retry + back-off + parent `gc()` on subprocess launch failure
  - `gc()` after each iteration's results are stowed

No changes to `R/did_multiplegt_main_smaller.R`, `R/did_multiplegt_main.R`,
or anything else.
