% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rex2arma.inc.R
\name{st2arma}
\alias{st2arma}
\title{Translate an elementary R statement to RcppArmadillo}
\usage{
st2arma(
  st,
  call2arma = call2a,
  indent = "",
  iftern = FALSE,
  env = parent.frame(),
  ...
)
}
\arguments{
\item{st}{A statement to transalte}

\item{call2arma}{An environment used as modifiable dictionary and giving
Armadillo equivalence
for some R functions. Entries are R function names to be translated.}

\item{indent}{A character string used as indentation. It is incremented
by "   " (3 spaces) inside a \code{{...}} block.}

\item{iftern}{A logical indicating whether \code{if/else} must be
considered as a ternary operator (TRUE) or as a classical \code{if/else}
operator (FALSE).}

\item{env}{An environment where statement may be executed.}

\item{...}{Parameters passed through to \code{symdim()} calls.}
}
\value{
A string with RcppArmadillo code.
}
\description{
Translate an elementary R statement to RcppArmadillo
}
\details{
Parameter \code{env} is also passed to \code{symdim()} with \code{...}.
translate statement to arma code
}
