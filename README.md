rex2arma
========
This package is a converter of R code to RcppArmadillo code.

It is not a finished project and probably will never be. However, even in its current state it can be useful in some (rather simple) situations. It can work for straightforward translation of mathematical formulas written in R to C++ (RcppArmadillo flavor).

Non exhaustive list of caveats can be found in inst/doc/caveats.txt

Installation
============
`devtools::install_github("sgsokol", "rex2arma")`

Unit Testing
============
Requires R package `RUnit` installed.

`source(system.file("tests/RUnit.R", package="rex2arma"))`

It can take several minutes, as many C++ compilations are done during testing.

Examples
========

```
x=rnorm(5) # we have to init R variable involved in R code
code=rex2arma("x+1", fname="vadd_one")
# if you wish to see the C++ code
cat(code)

# call C++ function
vadd_one(1:10)

```

Few more complete and more complex examples of usage are available in `inst/examples/ex_*.R` files

Legal Information
=================
Author: Serguei Sokol (INRAE, France)

Licence: GPL2

Copyright: 2023 INRAE/INSA/CNRS, France
