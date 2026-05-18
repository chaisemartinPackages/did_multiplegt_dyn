#!/usr/bin/env bash
# =============================================================================
# Bootstrap equivalence driver -- runs each spec sequentially in its own
# Rscript subprocess so polars/Rust memory is reclaimed between specs.
#
# Usage:  ./run_all_specs.sh [spec_name ...]
#   With no args, runs every spec in specs.R.
#   With names, runs only those specs.
# =============================================================================
set -u

cd "$(dirname "$0")" || exit 2

# --------------------------------------------------------------------------
# Allocator stabilisation for the LEGACY path.
#
# The legacy `row_replication` bootstrap runs several full did_multiplegt_main()
# calls in a single R process. Each call routes through polars (Rust), and the
# repeated alloc/free pattern across glibc + Rust can corrupt glibc's heap
# bookkeeping on some specs ("malloc(): unsorted double linked list corrupted",
# stochastic). The new subprocess+weights path is unaffected because every
# iteration is a fresh R process.
#
# We LD_PRELOAD jemalloc so both paths use a single, consistent allocator. This
# is a TEST-RUNNER stability fix only; the package code is untouched.
# MALLOC_ARENA_MAX=1 is kept as a defence-in-depth knob.
export MALLOC_ARENA_MAX=1
JEMALLOC_LIB=""
for cand in \
    /lib/aarch64-linux-gnu/libjemalloc.so.2 \
    /usr/lib/aarch64-linux-gnu/libjemalloc.so.2 \
    /usr/lib/x86_64-linux-gnu/libjemalloc.so.2 \
    /lib/x86_64-linux-gnu/libjemalloc.so.2; do
    if [[ -f "$cand" ]]; then JEMALLOC_LIB="$cand"; break; fi
done
if [[ -n "$JEMALLOC_LIB" ]]; then
    export LD_PRELOAD="$JEMALLOC_LIB${LD_PRELOAD:+:$LD_PRELOAD}"
    echo "Using jemalloc: $JEMALLOC_LIB"
fi

LOG_DIR="../logs/per_spec"
SUMMARY_CSV="../logs/bootstrap_equivalence_summary.csv"
mkdir -p "$LOG_DIR"
rm -f "$SUMMARY_CSV"

# ---- enumerate spec names --------------------------------------------------
if [[ $# -gt 0 ]]; then
    SPEC_NAMES=("$@")
else
    mapfile -t SPEC_NAMES < <(
        Rscript -e 'source("specs.R"); cat(names(SPECS_BY_NAME), sep="\n")'
    )
fi

N_SPECS=${#SPEC_NAMES[@]}
echo "Running $N_SPECS specs sequentially (each in its own R process)."
echo "Summary: $SUMMARY_CSV"
echo

# ---- header for summary ----------------------------------------------------
echo "spec,status,quantity,diff,t_new,t_old" > "$SUMMARY_CSV"

ok=0
fail=0
err=0
i=0
for name in "${SPEC_NAMES[@]}"; do
    i=$((i+1))
    log_file="$LOG_DIR/${name}.log"
    csv_file="$LOG_DIR/${name}.csv"
    printf "[%2d/%d] %-55s " "$i" "$N_SPECS" "$name"

    Rscript run_one_spec.R "$name" "$csv_file" \
        > "$log_file" 2>&1
    rc=$?

    if [[ -f "$csv_file" ]]; then
        # skip header, append to summary
        tail -n +2 "$csv_file" >> "$SUMMARY_CSV"
        # status of each row
        # write.csv() quotes strings: $2 is "OK"/"FAIL"/"ERROR" with quotes
        n_fail=$(awk -F, 'NR>1 && $2 ~ /FAIL/' "$csv_file" | wc -l | tr -d ' ')
        n_err=$( awk -F, 'NR>1 && $2 ~ /ERROR/' "$csv_file" | wc -l | tr -d ' ')
        n_ok=$(  awk -F, 'NR>1 && $2 ~ /OK/'    "$csv_file" | wc -l | tr -d ' ')
        if [[ "$n_fail" -gt 0 ]]; then
            echo "FAIL ($n_fail failing, $n_ok ok)"
            fail=$((fail+n_fail))
        elif [[ "$n_err" -gt 0 ]]; then
            echo "ERROR ($n_err errored)"
            err=$((err+n_err))
        else
            echo "OK"
            ok=$((ok+n_ok))
        fi
    else
        echo "ERROR (no CSV; rc=$rc) -- see $log_file"
        err=$((err+1))
        echo "$name,ERROR,NA,NA,NA,NA" >> "$SUMMARY_CSV"
    fi
done

echo
echo "===================== SUMMARY ====================="
echo "Specs run:    $N_SPECS"
echo "OK rows:      $ok"
echo "FAIL rows:    $fail"
echo "ERROR rows:   $err"
echo "Summary CSV:  $SUMMARY_CSV"
echo "Per-spec log: $LOG_DIR/<spec_name>.log"

if [[ "$fail" -gt 0 || "$err" -gt 0 ]]; then
    echo
    echo "Failing/erroring rows:"
    awk -F, 'NR==1 || ($2 !~ /OK/ && $2 != "status")' "$SUMMARY_CSV"
    exit 1
fi
exit 0
