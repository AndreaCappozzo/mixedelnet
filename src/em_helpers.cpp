#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


const double log2pi = std::log(2.0 * M_PI);

double dmvnrm_arma(arma::mat x,
                   arma::rowvec mean,
                   arma::mat sigma,
                   bool logd = false) {
  //int n = x.n_rows;
  int xdim = x.n_cols;
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  arma::vec z = rooti * arma::trans( x - mean) ;
  out      = constants - 0.5 * arma::sum(z%z) + rootisum;

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
Rcpp::List estep_lmm_cpp(arma::vec res_fixed, arma::mat Z,
                      arma::vec group_indicator, arma::mat inv_Omega,
                      double sigma2, int J)
{
  int N = Z.n_rows;
  int q = Z.n_cols;

  // Containers
  arma::mat est_second_moment(q,q);
  double est_second_moment_error = 0.0;
  arma::vec raneff_i(N);
  arma::mat mu_raneff(J,q);

  //  Fill containers
  raneff_i.zeros();
  est_second_moment.zeros();

  arma::mat Omega=inv_Omega.i();
  for ( int j = 1; j < (J+1); j++ ) {
    arma::uvec rows_j = find(group_indicator==j);
    arma::mat Z_j=Z.rows(rows_j);
    arma::vec res_fixed_j=res_fixed(rows_j);
    arma::mat first_piece = (Z_j.t()*Z_j)/ sigma2 + inv_Omega;
    arma::mat Gamma_j = first_piece.i();
    arma::vec mu_j = (Gamma_j*Z_j.t()*res_fixed_j)/sigma2;
    mu_raneff.row((j-1))=mu_j.t();
    raneff_i(rows_j)=Z_j * mu_j;
    est_second_moment += Gamma_j + mu_j * mu_j.t();
    // Rcout << est_second_moment;
    est_second_moment_error += as_scalar(arma::trace(Z_j*Gamma_j*Z_j.t()));
  }

  return Rcpp::List::create( Named("est_second_moment") = est_second_moment,
                             Named("est_second_moment_error") = est_second_moment_error,
                             Named("mu_raneff") = mu_raneff,
                             Named("raneff_i") = raneff_i);
}

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
double log_lik_lmm_cpp(arma::vec y, arma::mat Z, arma::mat X,
                      arma::vec group_indicator, arma::vec beta, arma::mat Omega,
                      double sigma2, int J)
{

  // Containers
  double log_lik_lmm = 0.0;

  for ( int j = 1; j < (J+1); j++ ) {
    arma::uvec rows_j = find(group_indicator==j);
    arma::mat Z_j=Z.rows(rows_j);
    arma::mat X_j=X.rows(rows_j);
    arma::vec y_j=y(rows_j);
    int n_j = X_j.n_rows;
    arma::mat diag_sigma(n_j,n_j);
    diag_sigma.eye();
    arma::mat G_j=Z_j * Omega * Z_j.t()+diag_sigma/pow(sigma2,-1);
    arma::vec mu_j=X_j*beta;
    log_lik_lmm+=dmvnrm_arma(y_j.t(),mu_j.t(),G_j,true);
  }

  return log_lik_lmm;
}
