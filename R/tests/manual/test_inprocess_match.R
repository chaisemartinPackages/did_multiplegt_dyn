# Quick equivalence check on a small problem (favara_imbs, 5 reps).
# Uses subprocess path *without* callr (so each iteration runs in the same
# process). This is faster than the full callr test and a good sanity
# check that the math matches before we do the slow subprocess+east test.
suppressPackageStartupMessages({
    library(DIDmultiplegtDYN)
    library(polars)
})

# Force in-process subprocess (skip callr) by temporarily masking it
# inside this script via DID_BOOTSTRAP_METHOD = "row_replication" or
# uninstalling callr. Easier: we just run NEW (with callr) once, OLD once,
# both replaying the same draws.
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

sd <- tempfile("boot_inproc_")
dir.create(sd)

cat("\n--- NEW (subprocess + weights, writes CSVs) ---\n")
new <- run("subprocess", sd, FALSE)

cat("\n--- LEGACY (row replication, replays CSVs) ---\n")
old <- run("row_replication", sd, TRUE)

cmp <- function(label, a, b, tol = 1e-8) {
    if (is.null(a) && is.null(b)) {
        cat(sprintf("[skip] %s\n", label)); return(invisible(NULL))
    }
    diff <- max(abs(as.matrix(a) - as.matrix(b)), na.rm = TRUE)
    cat(sprintf("[%s] %s: max abs diff = %.3e\n",
                if (diff < tol) "OK" else "FAIL", label, diff))
}

cat("\n--- COMPARISON ---\n")
cmp("Effects",  new$results$Effects,  old$results$Effects)
cmp("ATE",      new$results$ATE,      old$results$ATE)
cmp("Placebos", new$results$Placebos, old$results$Placebos)

cat("\nNEW Effects:\n");   print(new$results$Effects)
cat("LEGACY Effects:\n"); print(old$results$Effects)
cat("\nNEW ATE:\n");       print(new$results$ATE)
cat("LEGACY ATE:\n");      print(old$results$ATE)
cat("\nNEW Placebos:\n");  print(new$results$Placebos)
cat("LEGACY Placebos:\n"); print(old$results$Placebos)
