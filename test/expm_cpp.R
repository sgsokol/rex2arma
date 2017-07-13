
cppFunction(depends='RcppArmadillo', rebuild=TRUE, includes='
template <typename T>
inline unsigned which_max(T v) {
   unsigned i;
   v.max(i);
   return i+1;
}

template<typename T>
inline unsigned which_min(T v) {
   unsigned i;
   v.min(i);
   return i+1;
}
', 
"

using namespace arma;
using namespace Rcpp;

SEXP expm_cpp(
NumericMatrix A_in_,
bool balancing) {

   // auxiliary functions
   Environment base_env_r_=Environment::base_env();
   Function rep_r_=base_env_r_[\"rep\"];
   Function c_r_=base_env_r_[\"c\"];
   // External R function declarations
   Function expm_balance_r_=Environment(\"package:expm\")[\"balance\"];
   Function Matrix_norm_r_=Environment(\"package:Matrix\")[\"norm\"];

   // Input variable declarations and conversion
   mat A(A_in_.begin(), A_in_.nrow(), A_in_.ncol(), false);

   // Output and intermediate variable declarations
   mat A2;
   mat B;
   mat B2;
   mat B4;
   mat B6;
   List baP;
   List baS;
   mat C;
   vec c_;
   ivec d;
   vec dd;
   int i;
   mat I;
   int k;
   int l;
   int n;
   double nA;
   mat P;
   ivec pp;
   double s;
   vec t;
   mat tt;
   mat U;
   mat V;
   mat X;

   // Translated code starts here
   d=ivec(IntegerVector::create(A.n_rows, A.n_cols));
   if (d.size() != 2 || d.at(0) != d.at(1)) stop(\"'A' must be a square matrix\");
   n=d.at(0);
   if (n <= 1) return wrap(exp(A));
   if (balancing) {
      baP=as<List>(expm_balance_r_(A, \"P\"));
      baS=as<List>(expm_balance_r_(as<mat>(baP[\"z\"]), \"S\"));
      A=as<mat>(baS[\"z\"]);
   }
   nA=as<double>(Matrix_norm_r_(A, \"1\"));
   I=eye<mat>(n, n);
   if (nA <= 2.1) {
      t=vec({0.015, 0.25, 0.95, 2.1});
      l=which_max(nA <= t);
      C=join_vert(join_vert(join_vert(vec({120, 60, 12, 1, 0, 0, 0, 0, 0, 0}).st(), vec({30240, 15120, 3360, 420, 30, 1, 0, 0, 0, 0}).st()), vec({17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1, 0, 0}).st()), vec({17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1}).st());
      A2=A*A;
      P=I;
      U=C.at((l)-1, 1) * I;
      V=C.at((l)-1, 0) * I;
      for (int incr_arma_=(1 <= l ? 1 : -1), k=1; k != l+incr_arma_; k+=incr_arma_) {
         P=P*A2;
         U=U + C((l)-1, ((2 * k) + 2)-1) * P;
         V=V + C((l)-1, ((2 * k) + 1)-1) * P;
      }
      U=A*U;
      X=solve(V - U, V + U);
   } else {
      s=log2(nA / 5.4);
      B=A;
      if (s > 0) {
         s=ceil(s);
         B=B / (pow(2, s));
      }
      c_=vec({64764752532480000, 32382376266240000, 7771770303897600, 1187353796428800, 129060195264000, 10559470521600, 670442572800, 33522128640, 1323241920, 40840800, 960960, 16380, 182, 1});
      B2=B*B;
      B4=B2*B2;
      B6=B2*B4;
      U=B*(B6*(c_.at(13) * B6 + c_.at(11) * B4 + c_.at(9) * B2) + c_.at(7) * B6 + c_.at(5) * B4 + c_.at(3) * B2 + c_.at(1) * I);
      V=B6*(c_.at(12) * B6 + c_.at(10) * B4 + c_.at(8) * B2) + c_.at(6) * B6 + c_.at(4) * B4 + c_.at(2) * B2 + c_.at(0) * I;
      X=solve(V - U, V + U);
      if (s > 0)       for (int incr_arma_=(1 <= s ? 1 : -1), t=1; t != s+incr_arma_; t+=incr_arma_)       X=X*X;
   }
   if (balancing) {
      dd=as<vec>(baS[\"scale\"]);
      X=X % mat(vec((as<vec>(rep_r_(dd, n)) % as<vec>(rep_r_(1 / dd, _[\"each\"]=n)))).begin(), X.n_rows, X.n_cols, false);
      pp=conv_to<ivec>::from(as<vec>(baP[\"scale\"]));
      if (as<int>(baP[\"i1\"]) > 1) {
         for (int incr_arma_=((as<int>(baP[\"i1\"]) - 1) <= 1 ? 1 : -1), i=(as<int>(baP[\"i1\"]) - 1); i != 1+incr_arma_; i+=incr_arma_) {
            tt=X(span(), (i)-1);
            X(span(), (i)-1)=X(span(), (pp.at((i)-1))-1);
            X(span(), (pp.at((i)-1))-1)=tt;
            tt=X((i)-1, span());
            X((i)-1, span())=X((pp.at((i)-1))-1, span());
            X((pp.at((i)-1))-1, span())=tt;
         }
      }
      if (as<int>(baP[\"i2\"]) < n) {
         for (int incr_arma_=((as<int>(baP[\"i2\"]) + 1) <= n ? 1 : -1), i=(as<int>(baP[\"i2\"]) + 1); i != n+incr_arma_; i+=incr_arma_) {
            tt=X(span(), (i)-1);
            X(span(), (i)-1)=X(span(), (pp.at((i)-1))-1);
            X(span(), (pp.at((i)-1))-1)=tt;
            tt=X((i)-1, span());
            X((i)-1, span())=X((pp.at((i)-1))-1, span());
            X((pp.at((i)-1))-1, span())=tt;
         }
      }
   }

   return wrap(X);

}"
)
