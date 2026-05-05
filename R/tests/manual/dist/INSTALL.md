# Installing the patched `DIDmultiplegtDYN` in RStudio

This folder contains everything you need:

| File | What it is |
|---|---|
| `DIDmultiplegtDYN_2.3.4.tar.gz` | Source bundle of the patched package. Contains the new `did_multiplegt_main_smaller()` engine and the rewritten `did_multiplegt_bootstrap()`. |
| `..\test_equivalence_windows.R` | 8-rep equivalence test on `favara_imbs` (~30 s). |
| `..\test_east_windows.R` | 30-rep equivalence test on `east_et_al_2023` (slow but the headline test). |
| `..\test_memory_windows.R` | Memory-leak comparison on `favara_imbs` scaled 50×. |
| `..\run_all_windows.ps1` | PowerShell wrapper that runs all three. |

## Prerequisite: Rtools (Windows only)

The package compiles a small amount of C++ via Rcpp. Installing from source therefore needs **Rtools** matching your R version (Rtools43 for R 4.3.x, Rtools44 for R 4.4.x, Rtools45 for R 4.5.x). Get the matching installer from <https://cran.r-project.org/bin/windows/Rtools/> and run it once. RStudio picks it up automatically afterward.

To check Rtools is wired up:

```r
pkgbuild::find_rtools()  # should return TRUE and a non-empty path
```

## Install option 1 — RStudio menu

1. **Tools → Install Packages…**
2. Set **Install from:** to **Package Archive File (.zip; .tar.gz)**.
3. Click **Browse…** and pick `DIDmultiplegtDYN_2.3.4.tar.gz` from this folder.
4. Click **Install**.

## Install option 2 — R console

```r
install.packages(
    "C:/full/path/to/DIDmultiplegtDYN_2.3.4.tar.gz",
    repos = NULL,
    type  = "source"
)
```

## Install option 3 — `pak` (handles dependencies cleanly)

```r
# install.packages("pak")
pak::pkg_install("local::C:/full/path/to/DIDmultiplegtDYN_2.3.4.tar.gz")
```

`pak` will pull `polars`, `callr`, `data.table`, etc. for you.

## Required CRAN deps (if any are missing)

```r
install.packages(c("callr", "data.table"))
install.packages("polars", repos = "https://rpolars.r-universe.dev")
```

## Verifying the install

```r
library(DIDmultiplegtDYN)
packageVersion("DIDmultiplegtDYN")           # 2.3.4
exists("did_multiplegt_main_smaller",
       envir = asNamespace("DIDmultiplegtDYN"))   # TRUE
```

## Running the tests

From a PowerShell prompt **in the parent folder** (`...\tests\manual`):

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all_windows.ps1
```

Or run them one by one:

```powershell
Rscript .\test_equivalence_windows.R
Rscript .\test_east_windows.R
Rscript .\test_memory_windows.R
```

## What "works" looks like

* **Equivalence tests** — every `[OK]` line should report `max abs diff ≤ 1e-15`. (Exact 0 on `favara_imbs`, ~5e-17 on `east_et_al_2023`.)
* **Memory test** — after the bootstrap, `RSS after gc` should be **noticeably smaller** under `subprocess (NEW)` than under `row_replication (OLD)`. The bigger the dataset / replicate count, the bigger the gap.

## If something breaks

* `cannot find function "did_multiplegt_main_smaller"` → you installed an older copy of the package. Re-install from `DIDmultiplegtDYN_2.3.4.tar.gz` and restart R.
* `there is no package called 'callr'` → `install.packages("callr")`. Without `callr` the bootstrap silently falls back to the in-process slim engine and the memory leak fix doesn't apply.
* `polars` install errors → see <https://rpolars.r-universe.dev>; on R 4.3+ Windows you usually need the r-universe binary (the line above already targets it).
