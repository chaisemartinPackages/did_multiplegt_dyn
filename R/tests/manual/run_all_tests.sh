#!/usr/bin/env bash
# Runs the equivalence test (east_et_al_2023) and the memory leak test
# (favara_imbs scaled 50x), capturing logs side by side. Run from the
# package root.
set -u

LOG_DIR="$(cd "$(dirname "$0")/logs" && pwd)"
mkdir -p "$LOG_DIR"

stamp() { date +"%Y-%m-%d %H:%M:%S"; }

echo "[$(stamp)] === east_et_al_2023 equivalence test ==="
Rscript tests/manual/test_east_match.R 2>&1 | tee "$LOG_DIR/test_east_match.log"

echo
echo "[$(stamp)] === favara_imbs memory test (NEW vs OLD) ==="
Rscript tests/manual/test_favara_memory.R 2>&1 | tee "$LOG_DIR/test_favara_memory.log"

echo
echo "[$(stamp)] All tests done. Logs in $LOG_DIR"
