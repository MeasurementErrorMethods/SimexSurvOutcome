
#' Calculates SIMEX estimates for resampled data
#'
#' Runs the bootstrap to calculate standard
#' errors for the SIMEX procedure
#'
#' @param pred_mat Matrix of predictors
#'
#' @param lambavec Vector of lambda params
#' for SIMEX algorithm
#'
#' @param fail_times Vector of error-prone
#' failure times
#' 
#' @param init_b Matrix of starting values
#' for optimization
#'
#' @param fail_ind Matrix of failure indicators
#'
#' @param B Parameter for SIMEX algorithm
#'
#' @param simex_sd Standard deviation of
#' multiplicative error in event-time
#' 
#' @param n Number of samples
#'
#' @return Matrix of corrected betas from 
#' resampled data
#'
#' @useDynLib SimexSurvOutcome
#' @importFrom Rcpp sourceCpp
#' @exportPattern "^[[:alpha:]]+"
#'
#' @rdname bootstrap_simex
#' @export
bootstrap_simex <- function(pred_mat, lambdavec, fail_times,
                            fail_ind, init_b, B, simex_sd, n) {
  
  # resample n samples with replacement
  boot_samps <- sample(1:n, n, replace = TRUE)
  
  # get resampled samples
  boot_pred_mat <- pred_mat[boot_samps,]
  boot_fail_times <- fail_times[boot_samps]
  boot_fail_ind <- fail_ind[boot_samp,]
  
  # get SIMEX estimates
  boot_betas_simex <- simexOutcome(pred_mat, lambdavec, fail_times,
                                   fail_ind, init_b, B, simex_sd, n)
  
  return(boot_betas_simex)
  
}