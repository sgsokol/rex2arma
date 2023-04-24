\name{rex2arma-package}
\alias{rex2arma-package}
\docType{package}
\title{
    Convert R Expressions to RcppArmadillo Code
}
\description{
    Automatic conversion of basic R code (straightforward mathematical formulas and alike) into RcppArmadillo code.
}
\details{
    \tabular{ll}{
        Package: \tab rex2arma\cr
        Type: \tab Package\cr
        Version: \tab 0.7.1\cr
        Date: \tab 2023-04-20\cr
        License: \tab GPL (==2)\cr
    }
Principal function is rex2arma().
    Few examples are provided in \code{Sys.glob(file.path(system.file("examples", package="rex2arma"), "ex_*.R"))}
}
\author{
    Author: Serguei Sokol

    Maintainer: <sokol@insa-toulouse.fr>
}
\keyword{ automatic code translation }
\seealso{
    \code{\link[Rcpp:Rcpp-package]{Rcpp}}
}
\examples{
    x=rnorm(3) # make x a short vector of reals for prototyping
    vadd_one=function(x) x+1
    code=rex2arma(vadd_one, exec=0)
    cat(code)
    
    # or directly
    #cat(rex2arma("x+1", fname="vadd_one_bis", exec=1)) # commented for time saving during tests
    #vadd_one_bis(1:10)
}
