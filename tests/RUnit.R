Sys.unsetenv("R_TESTS")
if (requireNamespace("RUnit", quietly=TRUE) && requireNamespace("rex2arma", quietly=TRUE) &&
    requireNamespace("Rcpp", quietly=TRUE) && requireNamespace("RcppArmadillo", quietly=TRUE)) {
  testSuite <- RUnit::defineTestSuite(
    name = "rex2arma unit tests",
    dirs = system.file("unitTests", package = "rex2arma"),
    testFuncRegexp = "^[Tt]est.+",
    rngKind = "Mersenne-Twister"
  )
  tests <- RUnit::runTestSuite(testSuite, verbose=FALSE,
    useOwnErrorHandler=TRUE)

  RUnit::printTextProtocol(tests)

  if (RUnit::getErrors(tests)$nFail > 0) stop("RUnit test failure")
  if (RUnit::getErrors(tests)$nErr > 0) stop("Errors in RUnit tests")
}
