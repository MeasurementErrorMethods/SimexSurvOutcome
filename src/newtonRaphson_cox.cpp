
#include <RcppArmadillo.h>
using namespace Rcpp;

#include <iostream>
using namespace std;

#include "newtonRaphson_cox.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat newtonRaphson(arma::mat pred_mat, arma::mat curr_beta,
                        arma::mat failure_mat) {

  arma::mat eta = pred_mat * curr_beta;
  arma::vec haz = exp(eta);

  arma::vec rsk = as<NumericVector>(wrap(haz));

  std::reverse(rsk.begin(), rsk.end());

  arma::vec rsk_cumsum = cumsum(rsk);

  std::reverse(rsk_cumsum.begin(), rsk_cumsum.end());

  arma::colvec haz_colvec = arma::conv_to<arma::colvec>::from(haz);
  arma::rowvec rsk_cumsum_recip = arma::conv_to<arma::rowvec>::from(pow(rsk_cumsum, -1));

  arma::mat P = haz * rsk_cumsum_recip;

  arma::mat upper_tri_P = trimatl(P);

  arma::mat W = -upper_tri_P * diagmat(failure_mat) * upper_tri_P.t();
  int W_rows = W.n_rows;

  arma::mat ones_mat(upper_tri_P.n_rows, upper_tri_P.n_cols);
  ones_mat.ones();

  arma::mat intermed = upper_tri_P * diagmat(failure_mat) * (ones_mat - upper_tri_P).t();
  arma::mat intermed_diag = intermed.diag(0);

  for (int i = 0; i < W_rows; i++) {
    W(i,i) = intermed_diag(i);
  }

  arma::mat new_beta = arma::inv(pred_mat.t() * W * pred_mat) * pred_mat.t() *
    (failure_mat - (upper_tri_P * failure_mat)) + curr_beta;

  arma::umat diff = all(abs(new_beta - curr_beta) < 0.0001, 0);

  if (diff(0,0) == 1) {
    return new_beta;
  }

  return newtonRaphson(pred_mat, new_beta, failure_mat);

}

