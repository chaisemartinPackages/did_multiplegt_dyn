# Runs the three Windows test scripts in order, capturing logs to .\logs\.
# Usage (from a PowerShell prompt, in this directory):
#   powershell -ExecutionPolicy Bypass -File .\run_all_windows.ps1
$ErrorActionPreference = "Stop"

$here   = Split-Path -Parent $MyInvocation.MyCommand.Path
$logDir = Join-Path $here "logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

function Run-RScript($name) {
    $log = Join-Path $logDir ("{0}.log" -f $name)
    $rs  = Join-Path $here ("{0}.R" -f $name)
    if (-not (Test-Path $rs)) {
        Write-Host "Missing $rs - skipping"; return
    }
    Write-Host ""
    Write-Host "==> Running $name (log: $log)" -ForegroundColor Cyan
    & Rscript $rs 2>&1 | Tee-Object -FilePath $log
}

Run-RScript "test_equivalence_windows"
Run-RScript "test_east_windows"
Run-RScript "test_memory_windows"

Write-Host ""
Write-Host "All tests done. Logs in $logDir" -ForegroundColor Green
