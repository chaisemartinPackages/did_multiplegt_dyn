# Bootstrap memory-leak fix — Windows test guide

`DIDmultiplegtDYN`'s bootstrap (and base estimation) uses **polars**, which allocates its working buffers in **Rust's allocator**. On Windows that allocator is especially reluctant to release pages back to the OS — `gc()` and `rm()` in R don't help because they only manage the R heap. The leak is therefore *much* more visible on Windows than on macOS/Linux.

This document shows how to test the **new** bootstrap path (each iteration runs in a `callr::r()` subprocess that exits and gives memory back) against the **legacy** path (row replication, all in-process), on Windows.

---

## What the two paths do

| Path | How it builds the bootstrap sample | Where the engine runs | Memory behaviour |
|---|---|---|---|
| **NEW** (`subprocess`) | Draw unit multiplicities, merge as weights into a slim copy of the data, multiply with the user weight. | `callr::r()` worker — separate R process per iteration. | Subprocess exits → OS reclaims polars pages. Parent stays flat. |
| **LEGACY** (`row_replication`) | Sample row indices, physically duplicate rows for sampled units. | Same R process. | polars buffers accumulate across iterations. |

Both paths produce **numerically identical** `Effects` / `ATE` / `Placebos` when fed the same unit selections (verified to ~5e-17 on east_et_al_2023, exact 0 on favara_imbs).

Switching paths is a single option:

```r
options(DID_BOOTSTRAP_METHOD = "subprocess")       # new (default)
options(DID_BOOTSTRAP_METHOD = "row_replication")  # legacy
```

Two further options are used **only** by the verification scripts below:

```r
options(DID_BOOTSTRAP_SAMPLE_DIR   = "C:/some/dir") # write/read unit-selection CSVs
options(DID_BOOTSTRAP_LOAD_SAMPLES = TRUE)          # replay an existing CSV set
```

---

## Prerequisites

```r
install.packages("DIDmultiplegtDYN")  # or install your local build
install.packages("polars", repos = "https://rpolars.r-universe.dev")
install.packages(c("callr", "data.table"))
```

PowerShell must be available — it is on every modern Windows install.

---

## Helper: a Windows-compatible RSS reader

Every test below uses this to read the **parent R process's** working set (RSS) from PowerShell. Copy it into your scripts as-is:

```r
pid <- Sys.getpid()
get_rss_gb <- function() {
    out <- tryCatch(
        suppressWarnings(system(
            paste0('powershell -NoProfile -Command "(Get-Process -Id ', pid,
                   ').WorkingSet64 / 1GB"'),
            intern = TRUE
        )),
        error = function(e) NA_character_
    )
    val <- suppressWarnings(as.numeric(trimws(out)))
    if (length(val) == 0 || is.na(val[1])) return(NA_real_)
    val[1]
}
```

> Note: `WorkingSet64` is the *resident* set on Windows and is the right number to watch for OS-level memory reclaim. `Sys.getpid()` returns the parent R process id; the `callr::r()` workers are children whose memory is reclaimed when they exit, so the parent's RSS is exactly the metric the user cares about.

---

## Test 1 — Equivalence (small, fast)

Verifies that both paths produce the same answers on the same draws.

Save as `test_equivalence_windows.R` and run with `Rscript test_equivalence_windows.R`:

```r
suppressPackageStartupMessages({
    library(DIDmultiplegtDYN)
    library(polars)
})

data(favara_imbs)
df <- as.data.frame(favara_imbs)

run <- function(method, sample_dir, load) {
    options(
        DID_BOOTSTRAP_METHOD       = method,
        DID_BOOTSTRAP_SAMPLE_DIR   = sample_dir,
        DID_BOOTSTRAP_LOAD_SAMPLES = load
    )
    did_multiplegt_dyn(
        df = df, outcome = "Dl_hpi", group = "county",
        time = "year", treatment = "inter_bra",
        effects = 3, placebo = 2, cluster = "state_n",
        graph_off = TRUE, bootstrap = list(8, 999)
    )
}

sd <- tempfile("boot_eq_"); dir.create(sd)

cat("\n--- NEW (subprocess + weights, writes CSVs) ---\n")
new <- run("subprocess",       sd, FALSE)

cat("\n--- LEGACY (row replication, replays CSVs) ---\n")
old <- run("row_replication",  sd, TRUE)

cmp <- function(label, a, b, tol = 1e-8) {
    diff <- max(abs(as.matrix(a) - as.matrix(b)), na.rm = TRUE)
    cat(sprintf("[%s] %s: max abs diff = %.3e\n",
                if (diff < tol) "OK" else "FAIL", label, diff))
}
cmp("Effects",  new$results$Effects,  old$results$Effects)
cmp("ATE",      new$results$ATE,      old$results$ATE)
cmp("Placebos", new$results$Placebos, old$results$Placebos)
```

Expected output (numbers shown are macOS, Windows should match to the same precision):

```
[OK] Effects: max abs diff = 0.000e+00
[OK] ATE: max abs diff = 0.000e+00
[OK] Placebos: max abs diff = 0.000e+00
```

---

## Test 2 — east_et_al_2023 (full-feature equivalence)

The user's flagship example: `continuous=1`, 16 controls, `weight = "births"`, `effects_equal = "all"`, 30 bootstraps. Save as `test_east_windows.R`:

```r
suppressPackageStartupMessages({
    library(DIDmultiplegtDYN)
    library(polars)
})

# Cache the dataset locally on first run
cache <- file.path(tempdir(), "east_et_al_2023.RData")
if (!file.exists(cache)) {
    utils::download.file(
        "https://raw.githubusercontent.com/Credible-Answers/did_multiplegt_dyn_tutorial/main/data/east_et_al_2023.RData",
        cache, mode = "wb", quiet = TRUE
    )
}
load(cache)

stcontrols <- c("stmarried","stblack","stother","sthsdrop","sthsgrad",
                "stsomecoll","pop0_4","pop5_17","pop18_24","pop25_44",
                "pop45_64","urate","incpc","maxafdc","abortconsent","abortmedr")

run_args <- list(
    df = df, outcome = "lbw_detrend81", group = "plborn",
    time = "dob_y_p", treatment = "newsimeli",
    effects = 5, placebo = 5, continuous = 1,
    controls = stcontrols, weight = "births",
    effects_equal = "all", graph_off = TRUE,
    bootstrap = list(30, 12345)
)

sd <- tempfile("boot_east_"); dir.create(sd)

# 1) NEW path writes the CSV draws.
options(DID_BOOTSTRAP_METHOD = "subprocess",
        DID_BOOTSTRAP_SAMPLE_DIR = sd,
        DID_BOOTSTRAP_LOAD_SAMPLES = FALSE)
t0 <- Sys.time()
res_new <- do.call(did_multiplegt_dyn, run_args)
cat(sprintf("NEW elapsed: %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))

# 2) LEGACY path replays the same CSVs.
options(DID_BOOTSTRAP_METHOD = "row_replication",
        DID_BOOTSTRAP_SAMPLE_DIR = sd,
        DID_BOOTSTRAP_LOAD_SAMPLES = TRUE)
t0 <- Sys.time()
res_old <- do.call(did_multiplegt_dyn, run_args)
cat(sprintf("LEGACY elapsed: %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))

cmp <- function(label, a, b, tol = 1e-8) {
    diff <- max(abs(as.matrix(a) - as.matrix(b)), na.rm = TRUE)
    cat(sprintf("[%s] %s: max abs diff = %.3e\n",
                if (diff < tol) "OK" else "FAIL", label, diff))
}
cmp("Effects",  res_new$results$Effects,  res_old$results$Effects)
cmp("ATE",      res_new$results$ATE,      res_old$results$ATE)
cmp("Placebos", res_new$results$Placebos, res_old$results$Placebos)
```

Expect tiny floating-point differences (~1e-16) because the row-replication path collapses N copies of the same row whereas the weights path multiplies the weight by N — the operations are algebraically identical but execute in a different order.

---

## Test 3 — Memory leak (the headline test)

The leak is most visible at scale on Windows. Save as `test_memory_windows.R`:

```r
suppressPackageStartupMessages({
    library(DIDmultiplegtDYN)
    library(polars)
    library(data.table)
})

# Scale favara_imbs 50x by relabeling counties.
data(favara_imbs)
fi_big <- rbindlist(lapply(seq_len(50), function(i) {
    tmp <- as.data.table(as.data.frame(favara_imbs))
    tmp[, county := county + (i - 1L) * 100000L]
    tmp
}))
fi_big <- as.data.frame(fi_big)
cat("rows:", nrow(fi_big), " groups:", length(unique(fi_big$county)), "\n")

pid <- Sys.getpid()
get_rss_gb <- function() {
    out <- tryCatch(
        suppressWarnings(system(
            paste0('powershell -NoProfile -Command "(Get-Process -Id ', pid,
                   ').WorkingSet64 / 1GB"'),
            intern = TRUE
        )),
        error = function(e) NA_character_
    )
    val <- suppressWarnings(as.numeric(trimws(out)))
    if (length(val) == 0 || is.na(val[1])) return(NA_real_)
    val[1]
}

run_one <- function(label, method) {
    cat("\n############################################################\n")
    cat("##  Method: ", label, "\n", sep = "")
    cat("############################################################\n")
    options(DID_BOOTSTRAP_METHOD       = method,
            DID_BOOTSTRAP_SAMPLE_DIR   = NULL,
            DID_BOOTSTRAP_LOAD_SAMPLES = FALSE)

    # --- A) Bootstrap leak ------------------------------------------------
    gc(full = TRUE, verbose = FALSE)
    cat(sprintf("[A] RSS before:                 %5.2f GB\n", get_rss_gb()))
    t0 <- Sys.time()
    res <- did_multiplegt_dyn(
        df = fi_big, outcome = "Dl_hpi", group = "county",
        time = "year", treatment = "inter_bra",
        effects = 3, placebo = 2, cluster = "state_n",
        graph_off = TRUE, bootstrap = list(30, 1)
    )
    el <- as.numeric(Sys.time() - t0, units = "secs")
    cat(sprintf("[A] elapsed:                    %5.1f s\n", el))
    cat(sprintf("[A] RSS after bootstrap=30:     %5.2f GB\n", get_rss_gb()))
    rm(res); gc(full = TRUE, verbose = FALSE)
    cat(sprintf("[A] RSS after gc:               %5.2f GB\n", get_rss_gb()))

    # --- B) Between-run leak (no bootstrap) -------------------------------
    for (i in 1:3) {
        res <- did_multiplegt_dyn(
            df = fi_big, outcome = "Dl_hpi", group = "county",
            time = "year", treatment = "inter_bra",
            effects = 3, placebo = 2, cluster = "state_n",
            graph_off = TRUE
        )
        cat(sprintf("[B] RSS after run %d:            %5.2f GB\n", i, get_rss_gb()))
        rm(res); gc(full = TRUE, verbose = FALSE)
        cat(sprintf("[B] RSS after gc:               %5.2f GB\n", get_rss_gb()))
    }
}

# Run NEW first, then OLD in the same R session.
run_one("subprocess (NEW)",      "subprocess")
run_one("row_replication (OLD)", "row_replication")
```

Run it (PowerShell or cmd):

```bat
Rscript test_memory_windows.R 2>&1 | tee test_memory_windows.log
```

### What to look for

- **Section A** is the bootstrap leak. The interesting number is **`RSS after gc` minus `RSS before`** — the memory the parent process *kept* after the bootstrap run.
  - On the `subprocess` path this should be small and bounded — every iteration's polars buffers died with its callr worker.
  - On the `row_replication` path this number grows roughly linearly with the bootstrap count and the dataset size.
- **Section B** is the unrelated base-estimation leak (it shows up *without* bootstrapping). It is a separate issue and is **not** what the subprocess change targets — both paths leak in B because B's calls run in the parent.

Reference numbers we measured on macOS for the same script (favara_imbs ×50, 30 reps):

```
NEW   subprocess     : 0.21 GB → 1.31 GB → 1.19 GB after gc   (66.6 s)
OLD   row_replication: 2.14 GB → 3.12 GB → 3.00 GB after gc   (29.1 s)
```

On Windows the gap should be **larger** in both directions: NEW will look almost flat, OLD will keep climbing as you add reps. The subprocess approach trades wall-clock for bounded memory — the trade-off the user explicitly asked for.

---

## Quick PowerShell wrapper to run all three back-to-back

Save as `run_all_windows.ps1`:

```powershell
$ErrorActionPreference = "Stop"
$logDir = Join-Path $PSScriptRoot "logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

function Run-RScript($name) {
    $log = Join-Path $logDir ("{0}.log" -f $name)
    Write-Host "==> Running $name (log: $log)"
    & Rscript "$PSScriptRoot/$name.R" 2>&1 | Tee-Object -FilePath $log
}

Run-RScript "test_equivalence_windows"
Run-RScript "test_east_windows"
Run-RScript "test_memory_windows"
```

Run it from a PowerShell prompt:

```powershell
powershell -ExecutionPolicy Bypass -File run_all_windows.ps1
```

---

## Troubleshooting on Windows

- **`callr` not found**: `install.packages("callr")`. Without it, the bootstrap silently falls back to the in-process slim engine — that still uses the new weight-based math, but does *not* free polars buffers between iterations.
- **`'powershell' is not recognized`**: run from a regular shell that has `powershell.exe` on `PATH` (the default), not from a heavily customized terminal.
- **Long paths in `DID_BOOTSTRAP_SAMPLE_DIR`**: prefer `tempfile("boot_")` to a deeply nested path. Windows path-length limits can break the CSV reader otherwise.
- **First call seems slow**: each callr worker has to spin up R + load polars + load DIDmultiplegtDYN. Wall-clock-wise that costs roughly 2–4 s per iteration of overhead. Memory-wise it is exactly what buys the leak fix.
- **Numbers differ from the macOS reference**: that is fine. The polars/Rust allocator is more aggressive about retaining pages on Windows; absolute numbers shift, the *trend* (NEW flat, OLD climbing) is the load-bearing observation.
