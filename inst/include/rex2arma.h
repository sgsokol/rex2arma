#ifndef rex2arma_h
#define rex2arma_h


#define def_pref4(__NAME__, __pref__)                                                                                        \
template <typename T>                                                                                                        \
T __pref__##__NAME__(                                                                                                        \
    const T& x, double mu=0, double sigma=1., bool log = false                                                               \
) {                                                                                                                          \
    T res(x);                                                                                                                \
    std::transform(x.begin(), x.end(), res.begin(),                                                                          \
        [mu, sigma, log](double xx) -> double {return ::Rf_##__pref__##__NAME__##4(xx, mu, sigma, log);});                   \
    return res;                                                                                                              \
}                                                                                                                            \
inline double __pref__##__NAME__(                                                                                            \
    const double& x, double mu=0, double sigma=1., bool log = false                                                          \
) {                                                                                                                          \
    return ::Rf_##__pref__##__NAME__##4(x, mu, sigma, log);                                                                  \
}

#define def_pref5(__NAME__, __pref__)                                                                                        \
template <typename T>                                                                                                        \
T __pref__##__NAME__(                                                                                                        \
    const T& x, double mu=0, double sigma=1., bool lower = true, bool log = false                                            \
) {                                                                                                                          \
    T res(x);                                                                                                                \
    std::transform(x.begin(), x.end(), res.begin(),                                                                          \
        [mu, sigma, lower, log](double xx) -> double {return ::Rf_##__pref__##__NAME__##5(xx, mu, sigma, lower, log);});     \
    return res;                                                                                                              \
}                                                                                                                            \
inline double __pref__##__NAME__(                                                                                            \
    const double& x, double mu=0, double sigma=1., bool lower=true, bool log = false                                         \
) {                                                                                                                          \
    return ::Rf_##__pref__##__NAME__##5(x, mu, sigma, lower, log);                                                           \
}

#define RCPP_DPQ(__NAME__)                                                                                                   \
namespace Rcpp {                                                                                                             \
    def_pref4(__NAME__, d)                                                                                                   \
    def_pref5(__NAME__, p)                                                                                                   \
    def_pref5(__NAME__, q)                                                                                                   \
}

RCPP_DPQ(norm)

template <typename T>
inline uword which_max(T v) {
   return v.index_max()+1;
}

template<typename T>
inline uword which_min(T v) {
   return v.index_min()+1;
}

template <typename T>
Col<T> _rex_arma_rep(T x, int size, int each=1) {
   Col<T> v(size*each);
   v.fill(x);
   return v;
}
template <typename T>
Col<T> _rex_arma_rep(Col<T> x, int size, int each=1) {
/*
   Rcout << "x.size()=" << x.size() << std::endl;
   Rcout << "size=" << size << std::endl;
   Rcout << "each=" << each << std::endl;
   Mat<T> m=repmat(x.t(), each, size);
   m.print("m=");
*/
   return vectorise(repmat(x.t(), each, size));
}

template<typename T, typename T2>
inline T at_rc(T m, T2 idx) {
   // index of matrix m elements by a two column matrix idx (0-based)
   return m.elem(conv_to<uvec>::from(idx.col(0)+idx.col(1)*m.n_rows));
}
// for random generator
RNGScope scope;
// auxiliary functions
Environment base_env_r_=Environment::base_env();
Function rep_r_=base_env_r_["rep"];
Function c_r_=base_env_r_["c"];

#endif
