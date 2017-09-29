Sys.setenv("R_TESTS" = "")
library(testthat)
library(Rcpp)
library(rex2arma)

test_check("rex2arma")
