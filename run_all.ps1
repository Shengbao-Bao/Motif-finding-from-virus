# Viral RNA motif discovery - one-command reproduction
# Usage:  powershell -ExecutionPolicy Bypass -File .\run_all.ps1
#         powershell -ExecutionPolicy Bypass -File .\run_all.ps1 -Refresh   (re-download)
param(
    [switch]$Refresh,
    [int]$Seed = 42,
    [string]$Python = "python"
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $MyInvocation.MyCommand.Definition
Set-Location $root

Write-Host "[run_all] python: $((& $Python --version) -join ' ')"
Write-Host "[run_all] root:   $root"

$arguments = @("src/cli.py", "all", "--seed", $Seed)
if ($Refresh) { $arguments += "--refresh" } else { $arguments += "--offline" }

& $Python @arguments
if ($LASTEXITCODE -ne 0) { throw "pipeline failed with exit code $LASTEXITCODE" }

Write-Host ""
Write-Host "[run_all] self-checks"
& $Python "tests/test_stats.py"
if ($LASTEXITCODE -ne 0) { throw "self-checks failed" }

Write-Host ""
Write-Host "[run_all] artefacts:"
Get-ChildItem -Recurse -File results | ForEach-Object {
    "{0,-58} {1,10} B" -f $_.FullName.Substring($root.Length + 1), $_.Length
} | Sort-Object
