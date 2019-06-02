// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include <iostream>
using namespace std;

#include "newtonRaphson_cox.h"

#include "quadExtrapolation.h"

// [[Rcpp::export]]
arma::mat simexOutcome(arma::mat pred_mat, arma::vec lambdavec,
                       arma::vec fail_times, arma::mat fail_ind,
                       arma::mat init_b, int B,
                       double simex_sd, int n) {

  int num_lambda = lambdavec.size();

  // generate SIMEX estimates
  arma::mat simex_mat(num_lambda * B, 4);

  int j = 0;

  for (int i = 0; i < num_lambda; i++) {

    double lambda = lambdavec.at(i);

    for (int k = 0; k < B; k++) {

      arma::vec simex_err = sqrt(lambda) * rnorm(n, 0, simex_sd);
      arma::vec time_simex = fail_times % arma::exp(simex_err);

      // order failure indicators by time_simex
      arma::vec fail_ind_ordered = fail_ind(arma::sort_index(time_simex));

      // order design mat by time_simex
      arma::mat pred_mat_ordered = pred_mat.rows(arma::sort_index(time_simex));

      // newton-raphson to get cox model estimates
      arma::mat intermed_cox_results = newtonRaphson(pred_mat_ordered,
                                                     init_b,
                                                     fail_ind_ordered);

      // fill vectors with results
      simex_mat(j, 0) = lambda;
      simex_mat(j, 1) = intermed_cox_results(0,0);
      simex_mat(j, 2) = intermed_cox_results(1,0);
      simex_mat(j, 3) = intermed_cox_results(2,0);

      j += 1;
    }
  }

  // run quadratic extrapolation and get predictions at lambda = -1
  arma::mat corrected_betas = quadExtrapolation(simex_mat.col(0),
                                                simex_mat.col(1),
                                                simex_mat.col(2),
                                                simex_mat.col(3));

  return corrected_betas;

}

