@echo off
REM Launch the nanoamp Shiny GUI on Windows.
Rscript --vanilla -e "library(nanoamp); nanoamp_gui()"
if errorlevel 1 pause
