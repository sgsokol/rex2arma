require(Rcpp)
cppFunction(depends='RcppArmadillo',
'List f(
NumericVector res,
NumericVector x) {
   using namespace arma;
   // Variable declarations
   colvec res_arma_(res.begin(), res.size(), false);
   colvec x_arma_(x.begin(), x.size(), false);
   mat A=randu<mat>(2,3);
   mat B = repmat(A, 4, 5);
   // Translated code starts here
   res_arma_=x_arma_.t()*x_arma_;
   
   return Rcpp::List::create(
      Rcpp::Named("A")=A,
      Rcpp::Named("B")=B
   );
}'
)

require(Rcpp)
cppFunction(depends='RcppArmadillo', rebuild=TRUE,
   'List fastLm_rex(
NumericVector coef,
NumericVector df,
NumericMatrix res,
NumericMatrix s2,
NumericVector std_err,
NumericMatrix X,
NumericVector y) {
   using namespace arma;
   // Variable declarations
   colvec coef_arma_(coef.begin(), coef.size(), true);
   colvec df_arma_(df.begin(), df.size(), true);
   mat res_arma_(res.begin(), res.nrow(), res.ncol(), true);
   mat s2_arma_(s2.begin(), s2.nrow(), s2.ncol(), true);
   colvec std_err_arma_(std_err.begin(), std_err.size(), true);
   mat X_arma_(X.begin(), X.nrow(), X.ncol(), true);
   colvec y_arma_(y.begin(), y.size(), true);
   // Translated code starts here
   df_arma_=(X_arma_).n_rows-(X_arma_).n_cols;
   coef_arma_=solve(X_arma_,y_arma_);
   res_arma_=y_arma_-X_arma_*coef_arma_;
   s2_arma_=(res_arma_).t()*res_arma_/df_arma_;
   std_err_arma_=sqrt(as_scalar(s2_arma_)*diagvec(pinv((X_arma_).t()*X_arma_)));
   
   return Rcpp::List::create(
      Rcpp::Named("df")=df_arma_,
      Rcpp::Named("coef")=coef_arma_,
      Rcpp::Named("res")=res_arma_,
      Rcpp::Named("s2")=s2_arma_,
      Rcpp::Named("std_err")=std_err_arma_
   );

   }'
)
#-------
rm(fastLm_arma)
cppFunction(depends='RcppArmadillo', rebuild=T,
'List fastLm_arma(arma::colvec y, arma::mat X) {
//arma::colvec y=Rcpp::as<arma::colvec>(ys);// to arma
//arma::mat X=Rcpp::as<arma::mat>(Xs);

int df = X.n_rows - X.n_cols ;

arma::colvec coef = arma::solve ( X , y ); // fit y ~ X
arma::colvec res = y - X * coef ; // residuals

double s2 = std::inner_product ( res.begin (), res.end (), res.begin (), 0.0 )/ df ; // std.err coefs
arma::colvec std_err = arma::sqrt ( s2 * arma::diagvec ( arma::pinv ( arma::trans ( X ) * X )));

return Rcpp::List::create ( Rcpp::Named ( "df" ) = df,
Rcpp::Named ( "stderr" ) = std_err ,
Rcpp::Named ( "coefficients" ) = coef );
}')
y <- log(trees$Volume)
X <- cbind(1, log(trees$Girth))
fastLm_arma(y, X)
#--------
rm(retvec)
cppFunction(depends='RcppArmadillo', rebuild=T,
'NumericVector retvec(
arma::vec x) {
   using namespace arma;
   // Variable declarations
   //vec x(x_r_.begin(), x_r_.size(), false);
   // Translated code starts here
   double d=2.;
   x*=d;
   //std::vector<double> y(x.size()); // => OK
   // y=conv_to< Rcpp:NumericVector >::from(x); // => error
   //return Rcpp::as<NumericVector>(wrap(x)); // => 1-col matrix :(
   //return wrap(x); // => 1-col matrix :(
   return NumericVector(x.begin(), x.end()); // => OK
}'
)
retvec(1:5)
#---------
cppFunction(depends='RcppArmadillo', rebuild=TRUE,
'List fastLm_rex(
NumericMatrix X_r_,
NumericVector y_r_) {
   using namespace arma;
   // Variable declarations
   mat X(X_r_.begin(), X_r_.nrow(), X_r_.ncol(), false);
   colvec y(y_r_.begin(), y_r_.size(), false);
   // Translated code starts here
   double df=(X).n_rows-(X).n_cols;
   vec coef=solve(X,y);
   vec res=y-X*coef;
   double s2 = std::inner_product(res.begin(), res.end(),
                                   res.begin(), 0.0)/df;
   //double s2=as_scalar((res).t()*res)/df;
   vec std_err=sqrt(s2*diagvec(pinv((X).t()*X)));
   
   return Rcpp::List::create(
      Rcpp::Named("coefficients")=coef,
      Rcpp::Named("stderr")=std_err,
      Rcpp::Named("df")=df
   );

}'
)
#---------------environements
cppFunction(depends='RcppArmadillo', rebuild=TRUE,
'SEXP renv() {
   Environment glob(1);
   arma::vec foo=NumericVector(glob.get("foo"));
   Rcout << "foo+2=" << foo+2 << std::endl;
   glob.assign("foo", 4);
   return R_NilValue;
}'
)
#---------------resize a vec
cppFunction(depends='RcppArmadillo', rebuild=TRUE,
'NumericVector vecresize(arma::vec x, int n) {
   x.resize(n);
   return wrap(x);
}'
)
#---------------create a vec from repeated scalar
cppFunction(depends='RcppArmadillo', rebuild=TRUE,
'NumericVector crevec(arma::vec x, double n) {
   //arma::vec y=arma::vec(n);
   return wrap(arma::vec(n).fill(arma::as_scalar(x)));
}'
)
#---------------
cppFunction(depends='RcppArmadillo', rebuild=TRUE,
'NumericVector cvec(arma::vec x, double y) {
   //x = arma::join_cols(x, y, z);
   x << y;
   return wrap(x);
}'
)
#---------------return types
cppFunction(depends='RcppArmadillo', rebuild=TRUE,
'SEXP ret(arma::Mat<char> x) {
   return wrap(x);
}'
)
