\name{rex2arma-package}
\alias{rex2arma-package}
\alias{rex2arma}
\docType{package}
\title{
What the package does (short line)
~~ package title ~~
}
\description{
More about what it does (maybe more than one line)
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab rex2arma\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-27\cr
License: \tab What license is it under?\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
    x=double(3) # make x a short vector of reals for prototyping
    vadd_one=function(x) x+1
    code=rex2arma(vadd_one, exec=0)
    cat(code)
}
