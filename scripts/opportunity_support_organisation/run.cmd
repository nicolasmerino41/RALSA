@echo off
julia --startup-file=no --project="%~dp0..\.." "%~dp0run_analysis.jl"
exit /b %errorlevel%
