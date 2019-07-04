
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
#' @param num_boot Number of bootstrap
#' replicates desired
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
                               fail_ind, B = 50, simex_sd,
                               num_boot = 250) {
  
  if (B%%1 != 0) {
    stop("B must be an integer")
  }
  
  if (num_boot%%1 != 0) {
    stop("num_boot must be an integer")
  }
  
  if (simex_sd <= 0) {
    stop("simex_sd must be greater than 0")
  }
  
  if (any(lambdavec < 0)) {
    stop("All lambda values must be at least 0")
  }

  # prepare data for analysis
  pred_mat <- as.matrix(pred_mat)
  lambdavec <- c(lambdavec)
  fail_times <- c(fail_times)
  fail_ind <- as.matrix(fail_ind)
  
  if (all(sapply(list(dim(pred_mat)[1],
                      length(fail_times)) 
                 function(x) x == dim(fail_ind)[1])) == FALSE) {
    stop("X_mat, time_star, and delta must all contain the same
         number of observations")
  }

  # create necessary variables
  init_b <- matrix(c(0, 0, 0), nrow = 3, ncol = 1)
  n <- dim(pred_mat)[1]
  simex_sd <- as.numeric(simex_sd)
  
  message("Running SIMEX")

  # get SIMEX estimates
  betas_simex <- simexOutcome(pred_mat, lambdavec, fail_times,
                              fail_ind, init_b, B, simex_sd, n)
  
  message("Calculating standard errors for SIMEX")
  
  simex_boot <- replicate(num_boot, bootstrap_simex(pred_mat,
                                                    lambdavec,
                                                    fail_times,
                                                    fail_ind,
                                                    init_b,
                                                    B, 
                                                    simex_sd,
                                                    n))
  
  message("Done!")
  
  simex_mat <- matrix(c(betas_simex[1,],
                        sd(simex_boot[1,], na.rm = TRUE),
                        betas_simex[2,],
                        sd(simex_boot[2,], na.rm = TRUE),
                        betas_simex[3,],
                        sd(simex_boot[3,], na.rm = TRUE)),
                      nrow = 3, ncol = 2, byrow = TRUE)
  
  rownames(simex_mat) <- colnames(pred_mat)
  colnames(simex_mat) <- c("Estimate", "Standard Error")
  
  return(as.data.frame(simex_mat))

}



