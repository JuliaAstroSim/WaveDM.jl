# run_julia_pipeline.ps1 — end-to-end WaveDM.jl constraints Julia pipeline
#
# Runs from a USER terminal (NOT inside Trae sandbox — that can't write to E:\.julia).
# Assumes you have already done, in your default v1.13 env:
#
#     Pkg.develop(path="E:/JuliaAstroSim/ParameterSpaceVis.jl")
#     Pkg.add("YAML")
#
# Usage:
#     powershell -ExecutionPolicy Bypass -File run_julia_pipeline.ps1
#
# Produces:
#     output/ma_constraints.pspace.yaml
#     output/ma_constraints_filtered.yaml
#     output/ma_constraints_julia.png  / .pdf  / .svg

$ErrorActionPreference = "Stop"
$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
$JULIA = "D:\Programs\julia-1.13.0\bin\julia.exe"

if (-not (Test-Path $JULIA)) {
    Write-Error "Julia not found at $JULIA; please edit this script and set `$JULIA."
}

$env:JULIA_DEPOT_PATH = "E:\.julia"
Write-Host "[1/3] Building ParameterSpace DataSet from YAML..."
& $JULIA --startup-file=no --project=@. build_constraints_dataset.jl
if ($LASTEXITCODE -ne 0) { throw "build_constraints_dataset.jl failed" }

Write-Host ""
Write-Host "[2/3] Plotting via ParameterSpaceVis (CairoMakie backend)..."
& $JULIA --startup-file=no --project=@. plot_constraints_julia.jl
if ($LASTEXITCODE -ne 0) { throw "plot_constraints_julia.jl failed" }

Write-Host ""
Write-Host "[3/3] Done. Outputs:"
Get-ChildItem -Path (Join-Path $Here "output") -Filter "ma_constraints*" |
    Select-Object Name, Length | Format-Table | Out-String | Write-Host