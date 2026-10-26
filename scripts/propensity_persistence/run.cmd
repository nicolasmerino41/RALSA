@echo off
cd /d "%~dp0\..\.."
julia --startup-file=no --project=. scripts/propensity_persistence/run_analysis.jl
