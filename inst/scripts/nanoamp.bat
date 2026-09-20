@echo off
REM Run the nanoamp CLI on Windows.
Rscript --vanilla -e "library(nanoamp); nanoamp_cli()" %*
