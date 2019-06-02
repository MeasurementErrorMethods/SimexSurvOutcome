
#' Calculates SIMEX estimates for error in a time-to-event
#'
#' The SIMEX method adds additional error to the error-prone
#' time-to-event that has either known standard deviation or
#' is calculated from a validation subset. Then, a quadratic
#' function is fit to the error-prone estimates with
#' additional error and we extrapolate the function to
#' where theoretically, there is no error
#'
#' @param pred_mat Matrix of predictors
#'
#' @param lambavec Vector of lambda params
#' for SIMEX algorithm
#'
#' @param fail_times Vector of error-prone
#' failure times
#'
#' @param fail_ind Matrix of failure indicators
#'
#' @param B Parameter for SIMEX algorithm
#'
#' @param simex_sd Standard deviation of
#' multiplicative error in event-time
#'
#' @return Matrix of corrected betas
#'
#' @useDynLib SimexSurvOutcome
#' @importFrom Rcpp sourceCpp
#' @exportPattern "^[[:alpha:]]+"
#'
#' @rdname simex_surv_outcome
#' @export
simex_surv_outcome <- function(pred_mat, lambdavec, fail_times,
                               fail_ind, B = 50, simex_sd) {

  # prepare data for analysis
  pred_mat <- as.matrix(pred_mat)
  lambdavec <- c(lambdavec)
  fail_times <- c(fail_times)
  fail_ind <- as.matrix(fail_ind)

  # create necessary variables
  init_b <- matrix(c(0, 0, 0), nrow = 3, ncol = 1)
  n <- dim(pred_mat)[1]
  simex_sd <- as.numeric(simex_sd)

  # get SIMEX estimates
  betas_simex <- simexOutcome(pred_mat, lambdavec, fail_times,
                              fail_ind, init_b, B, simex_sd, n)

  return(betas_simex)

}



