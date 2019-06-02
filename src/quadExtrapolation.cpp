
#include <RcppArmadillo.h>
using namespace Rcpp;

#include <iostream>
using namespace std;

#include "fastLM.h"

#include "quadExtrapolation.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat quadExtrapolation(arma::colvec lambda_mat, arma::colvec coef1_mat,
                            arma::colvec coef2_mat, arma::colvec coef3_mat) {

  // fit quadratic function to new naive estimates
  int num_simex = lambda_mat.size();
  arma::vec int_vec(num_simex);
  int_vec.ones();

  // create matrix of lambda predictors
  arma::mat lambda_pred_mat(num_simex, 3);
  lambda_pred_mat.col(0) = int_vec;
  lambda_pred_mat.col(1) = lambda_mat;
  lambda_pred_mat.col(2) = pow(lambda_mat, 2);

  // run linear regression
  List quad_fit1 = fastLM(lambda_pred_mat, coef1_mat);
  arma::vec fit1 = quad_fit1(0);

  List quad_fit2 = fastLM(lambda_pred_mat, coef2_mat);
  arma::vec fit2 = quad_fit2(0);

  List quad_fit3 = fastLM(lambda_pred_mat, coef3_mat);
  arma::vec fit3 = quad_fit3(0);

  // create matrix of coefficients
  arma::mat coef_mat(3, 3);
  coef_mat.col(0) = fit1;
  coef_mat.col(1) = fit2;
  coef_mat.col(2) = fit3;

  // extrapolate to lambda_simex = -1
  arma::vec pred_vec(3);
  pred_vec(0) = 1.0;
  pred_vec(1) = -1.0;
  pred_vec(2) = 1.0;

  // get predictions at lambda_simex = -1
  arma::mat simex_betas = coef_mat.t() * pred_vec;

  return simex_betas;

}
