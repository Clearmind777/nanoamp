#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(testthat)
  library(nanoamp)
})
test_check("nanoamp")
