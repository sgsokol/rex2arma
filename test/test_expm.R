# test rex2arma() on a real word example, a function calculating
# matrix exponential: expm.Higham08()

library(rbenchmark)
library(RcppKalman)
library(expm)
#source("../rex2arma.inc.R"); source("../rex2arma.R")
#require(Rcpp, lib="~/.R/library")
library(Rcpp)
library(rex2arma, lib="~/.R/library")
library(rbenchmark)



# The original R code for expm(A) by a method of Higham is accessible
# from the package expm, function expm.Higham08

# The original R code was modified in 4 places before submission to rex2arma():
# 1. "I <- if (is(A, "Matrix")) Diagonal(n) else diag(n)"
#    was replaced by "I <- diag(n)" as two differnet types (here sparse and dense
#    matrices) cannot be assigned to the same variable (here I)
# 2. "d <- baS$scale" was replaced by "dd <- baS$scale" as the name "d" was already used
#    in a previous code for an integer dimension vector ("d <- dim(A)").
#    Changing variable type (here from ivec to vec) is not allowed in C++ at runtime.
#    Naturally, the subsequent use of name "d" was replaced by "dd".
# 3. "X <- X * (dd * rep(1/dd, each = n))" was replaced by "X <- X * (rep(dd, n) * rep(1/dd, each = n))"
#    More narrowly, "dd *" in the above expression was replaced by "rep(dd, n) *" as
#    RcppArmadillo does not short vector recylcling in binary operations.
#    So we have to ensure in the R code that vectors have an equal length in
#    binary operations.
# The original code can be seen in expm.Higham08.R, while the modifed
# version with the three above corrections is in expm.higham.R

# get R code

source("expm.higham.R") # create expm.higham() function

# translate it in C++ expm_cpp() function
A=diag(2) # just examples of two input parameters
balancing=TRUE
rex_code=rex2arma(expm.higham, fname="expm_cpp", exec=1, rebuild=TRUE, copy=FALSE, verbose=TRUE)

# make tests and benchmarks
# Let do it on three matrix sizes: small (s, n=10), medium (m, n=100) and
# large (l, n=1000)

set.seed(7) # make the tests reproducible

# small
n=5
As=matrix(rnorm(n*n), n)
cat(sprintf("small matrix (%d x %d) benchmark\n", n, n))
print(benchmark(
   expm.higham(As, balancing=TRUE),
   expm.higham(As, balancing=FALSE),
   expm_cpp(As, balancing=TRUE),
   expm_cpp(As, balancing=FALSE),
   RcppKalman::expm(As),
   replications=500,
   order="relative"
)[,1:4])
stopifnot(diff(range(expm.higham(As, balancing=TRUE)-expm_cpp(As, balancing=TRUE))) < 1.e-13)
stopifnot(diff(range(expm.higham(As, balancing=FALSE)-expm_cpp(As, balancing=FALSE))) < 1.e-13)
stopifnot(diff(range(expm.higham(As, balancing=FALSE)-RcppKalman::expm(As))) < 1.e-13)

# medium
n=50
As=matrix(rnorm(n*n), n)
cat(sprintf("medium matrix (%d x %d) benchmark\n", n, n))
print(benchmark(
   expm.higham(As, balancing=TRUE),
   expm.higham(As, balancing=FALSE),
   expm_cpp(As, balancing=TRUE),
   expm_cpp(As, balancing=FALSE),
   RcppKalman::expm(As),
   order="relative"
)[,1:4])
stopifnot(diff(range(expm.higham(As, balancing=TRUE)-expm_cpp(As, balancing=TRUE))) < 1.e-11)
stopifnot(diff(range(expm.higham(As, balancing=FALSE)-expm_cpp(As, balancing=FALSE))) < 1.e-11)
stopifnot(diff(range(expm.higham(As, balancing=FALSE)-RcppKalman::expm(As))) < 1.e-9)

# large
n=500
As=matrix(rnorm(n*n), n)
cat(sprintf("large matrix (%d x %d) benchmark\n", n, n))
print(benchmark(
   expm.higham(As, balancing=TRUE),
   expm.higham(As, balancing=FALSE),
   expm_cpp(As, balancing=TRUE),
   expm_cpp(As, balancing=FALSE),
   RcppKalman::expm(As),
   replications=5,
   order="relative"
)[,1:4])
stopifnot(diff(range(expm.higham(As, balancing=TRUE)-expm_cpp(As, balancing=TRUE))) < 1.e-14)
stopifnot(diff(range(expm.higham(As, balancing=FALSE)-expm_cpp(As, balancing=FALSE))) < 1.e-14)
stopifnot(diff(range(expm.higham(As, balancing=FALSE)-RcppKalman::expm(As))) < 2.e-4)
