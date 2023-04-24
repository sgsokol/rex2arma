#browser()
rebuild <- FALSE
exec <- 2 # convert to C++, execute it and return its result
verbose <- FALSE
library(rex2arma)
library(RUnit)
test_simple=function() {
   n <- 1L
   x <- 2.5
   l <- 2>1
   vn <- 1L:2L
   vx <- 2.5+1:2
   vl <- (1:2)>1
   mn <- matrix(1L:6L, 3L)
   mx <- 2.5+mn

   cases=alist(
      sinteger=n+1L,
      sdouble=x+1.1,
      slogical=l && TRUE,
      vinteger=vn+1L,
      vdouble=vx+1.1,
      #vlogical=vl & TRUE # doesn't have an equivalent in Armadillo
      minteger=mn+1L,
      mdouble=mx+1.1
   )
   for (nm in names(cases)) {
      e=cases[[nm]]
      te=deparse(e)
      checkEqualsNumeric(c(rex2arma(te, r=rebuild, exec=exec, verbose=verbose, envir=parent.frame())), c(eval(e)), msg=paste0("translating '", te, "'", collapse=""))
   }
}


# arithmetic homogeneous operations
test_arith=function() {
   # binary ops
   # prepare inputs
   n <- 10
   set.seed(7)

   # scalars
   s1 <- runif(1)
   s2 <- runif(1)

   # vectors
   v1 <- rnorm(n)
   v2 <- rnorm(n)

   # matrices
   m1 <- matrix(rnorm(n*n), n)
   m2 <- matrix(rnorm(n*n), n)
   op <- c("+", "-", "*", "/", "%*%", "%o%")
   vartype <- c("s", "v", "m")
   for (o in op) {
      for (a in vartype) {
         # exclude some meaningless cases
         #if (o == "%*%" && a == "s") {
         #   next
         #}
         src=paste(a, "1", o, a, "2", sep="")
         # pass text var
         checkEquals(c(rex2arma(src, r=rebuild, exec=exec)), c(eval(parse(text=src))), msg=src)
         
         # pass one shot function
         eval(parse(t=sprintf("f=function(%s1, %s2) %s", a, a, src)), envir=.GlobalEnv)
         checkEquals(c(rex2arma(f, r=rebuild, exec=exec, envir=parent.frame())), c(eval(parse(text=sprintf("f(%s1, %s2)", a, a)))), msg="function()")
         # pass body of f in braces {...}
         fsrc=sprintf("f=function(%s1, %s2) {%s}", a, a, src)
         eval(parse(t=fsrc), envir=.GlobalEnv)
         checkEqualsNumeric(c({rex2arma(f, fname="f_arma_", r=rebuild, exec=exec); do.call(f_arma_, lapply(paste(a, 1:2, sep=""), as.symbol))}), c(eval(parse(t=sprintf("f(%s1, %s2)", a, a)))), msg="function(){}")
         # pass plain expression
         e=parse(t=src)
         checkEquals(c(rex2arma(e, r=rebuild, exec=exec)), c(eval(e)), msg="expr")
         # pass braced expression
         e=parse(t=sprintf("{%s}", src))
         checkEquals(c(rex2arma(e, r=rebuild, exec=exec)), c(eval(e)), msg="{expr}")
      }
   }
}
