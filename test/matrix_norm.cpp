//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nmat(mat A) {
   Function Matrix_norm_r_=Environment("package:Matrix")["norm"];
   //SEXP s=PROTECT(wrap("1"));
   //SEXP s=wrap("1");
   //Rf_PrintValue(s);
   //double res=as<double>(Matrix_norm_r_(wrap(A), CharacterVector("1")));
   //double res=as<double>(Matrix_norm_r_(NumericMatrix(wrap(A)), CharacterVector::create("1")));
   double res=as<double>(Matrix_norm_r_(A, "1"));
   //UNPROTECT(1);
   return res;
}
// [[Rcpp::export]]
SEXP s2(mat a) {
   Function f=Environment("package:Matrix")["norm"];
   //if (1) {
   //   f=Environment("package:Matrix")["norm"];
   //}
   
   SEXP s=wrap("1");
   return s;
}
