# unit tests of rex2arma()
library("rex2arma")
library("Rcpp")
rebuild=FALSE

# prepare inputs
n=10
set.seed(7)

# scalars
s1=runif(1)
s2=runif(1)

# vectors
v1=rnorm(n)
v2=rnorm(n)

# matrices
m1=matrix(rnorm(n*n), n)
m2=matrix(rnorm(n*n), n)

# arithmetic homogeneous operations
op=c("+", "-", "*", "/", "%*%", "%o%")
vartype=c("s", "v", "m")

for (o in op) {
   for (a in vartype) {
      # exclude some meaningless cases
      #if (o == "%*%" && a == "s") {
      #   next
      #}
      src=paste(a, "1", o, a, "2", sep="")
      # pass text var
      if (!all.equal(c(rex2arma(src, r=rebuild)), c(eval(parse(t=src))))) {
         stop(sprintf("test '%s' failed", src))
      }
      # pass one shot function
      eval(parse(t=sprintf("f=function(%s1, %s2) %s", a, a, src)))
      if (!all.equal(c(rex2arma(f, r=rebuild)), c(eval(parse(t=sprintf("f(%s1, %s2)", a, a)))))) {
         stop(sprintf("test 'f(%s1, %s2)%s' failed", a, a, src))
      }
      # pass body of f in braces {...}
      fsrc=sprintf("f=function(%s1, %s2) {%s}", a, a, src)
      eval(parse(t=fsrc))
      if (!all.equal(c(rex2arma(f, r=rebuild)), c(eval(parse(t=sprintf("f(%s1, %s2)", a, a)))))) {
         stop(sprintf("test '%s' failed", a, a, fsrc))
      }
      # pass plain expression
      e=parse(t=src)
      if (!all.equal(c(rex2arma(e, r=rebuild)), c(eval(e)))) {
         stop(sprintf("test 'e'", e))
      }
      # pass braced expression
      e=parse(t=sprintf("{%s}", src))
      if (!all.equal(c(rex2arma(e, r=rebuild)), c(eval(e)))) {
         stop(sprintf("test 'e'", e))
      }
   }
}

# function specific test
f="head"
# one argument
for (a in c("v", "m")) {
   src=sprintf("%s(%s1)", f, a)
   if (!all.equal(c(rex2arma(src, r=rebuild)), c(eval(parse(t=src))))) {
      stop(sprintf("test '%s' failed", src))
   }
}
# second argument is small positive, big positive, small negative, big negative, zero
for (nh in c(3, n+3, -3, -3-n, 0)) {
   for (a in c("v", "m")) {
      src=sprintf("%s(%s1, %s)", f, a, nh)
      if (!all.equal(c(rex2arma(src, r=rebuild)), c(eval(parse(t=src))))) {
         stop(sprintf("test '%s' failed", src))
      }
   }
}
