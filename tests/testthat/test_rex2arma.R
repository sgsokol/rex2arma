#browser()
rebuild <<- FALSE
exec <<- 2 # convert to C++, execute it and return its result
verbose <<- FALSE
library(rex2arma)
context("simple expressions")
n <<- 1L
x <<- 2.5
l <<- 2>1
vn <<- 1L:2L
vx <<- 2.5+1:2
vl <<- (1:2)>1
mn <<- matrix(1L:6L, 3L)
mx <<- 2.5+mn

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
   expect_equal(c(rex2arma(te, r=rebuild, exec=exec, verbose=verbose)), c(eval(e)), info=paste0("translating '", te, "'", collapse=""))
}

# binary ops
context("test text, func(), func(){expr], expr, {expr}")
# prepare inputs
n <<- 10
set.seed(7)

# scalars
s1 <<- runif(1)
s2 <<- runif(1)

# vectors
v1 <<- rnorm(n)
v2 <<- rnorm(n)

# matrices
m1 <<- matrix(rnorm(n*n), n)
m2 <<- matrix(rnorm(n*n), n)

# arithmetic homogeneous operations
op <<- c("+", "-", "*", "/", "%*%", "%o%")
vartype <<- c("s", "v", "m")
for (o in op) {
   for (a in vartype) {
      # exclude some meaningless cases
      #if (o == "%*%" && a == "s") {
      #   next
      #}
      src=paste(a, "1", o, a, "2", sep="")
      # pass text var
      test_that(src, expect_equal(c(rex2arma(src, r=rebuild, exec=exec)), c(eval(parse(text=src))), label=src))
      
      # pass one shot function
      eval(parse(t=sprintf("f=function(%s1, %s2) %s", a, a, src)), envir=.GlobalEnv)
      test_that(src, expect_equal(c(rex2arma(f, r=rebuild, exec=exec)), c(eval(parse(text=sprintf("f(%s1, %s2)", a, a)))), label="function()"))
      # pass body of f in braces {...}
      fsrc=sprintf("f=function(%s1, %s2) {%s}", a, a, src)
      eval(parse(t=fsrc), envir=.GlobalEnv)
      test_that(src, expect_equal(c({rex2arma(f, fname="f_arma_", r=rebuild, exec=exec); do.call(f_arma_, lapply(paste(a, 1:2, sep=""), as.symbol))}), c(eval(parse(t=sprintf("f(%s1, %s2)", a, a)))), label="function(){}"))
      # pass plain expression
      e=parse(t=src)
      test_that(src, expect_equal(c(rex2arma(e, r=rebuild, exec=exec)), c(eval(e)), label="expr"))
      # pass braced expression
      e=parse(t=sprintf("{%s}", src))
      test_that(src, expect_equal(c(rex2arma(e, r=rebuild, exec=exec)), c(eval(e)), label="{expr}"))
   }
}
