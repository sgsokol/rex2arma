#ifdef RETURN_COLVEC_AS_VECTOR
namespace Rcpp {
    inline SEXP wrap(const arma::vec& obj) {
        Vector<REALSXP> x(obj.begin(), obj.end());
        return x;
    }
}
#endif
