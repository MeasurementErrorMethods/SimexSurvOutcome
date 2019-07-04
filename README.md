
# SimexSurvOutcome

This package implements the Simulation and
Extrapolation (SIMEX) method to correct for bias in 
time-to-event data in the presence of multiplicative
measurement errors in the censored event time. The 
assumed model is the Cox proportional hazards model.

The method assumes that you have an estimate of
the variance of the measurement error on the additive
(ie. log) scale. If not, this package can be used to 
perform sensitivity analyses of the Cox model
to varying magnitudes of the measurement error 
variance. 

We demonstrate how to use the method below. For
further details, see our [paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.7554).

## Installation

To install and load this package in R from GitHub, run the following commands:
  
```{r}
install.packages("devtools")
library(devtools)
devtools::install_github("ericoh17/SimexSurvOutcome")
```

## Load package
```R
library(SimexSurvOutcome)
```

## Getting started

First, we will read in a simulated dataset. 

```{r}
data(example_dat)
```

Then we extract separate dataframes containing the
covariates, error-prone censored event time, and
the failure indicator. 

```{r}
cov <- c("X1", "X2", "X3")
X_mat <- example_dat[, cov]

time_star <- example_dat[, "time_star"]

delta <- example_dat[, "delta"]
```

Due to this being a simulated dataset,
we know the true variance of the measurement
error on the log scale. For real data, this 
parameter could be varied. 
```{r}
var_error <- 0.5
```

## Settings

Set the number of SIMEX replicates, B, and
vector containing additional error to add
to the error-prone censored event time.
See our [paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.7554)
for further discussion.

```{r}
B <- 50
lambda <- seq(0.1, 0.9, 0.1)
```

Set the number of bootstrap replicates to
perform to calculate standard errors:

```R
num_boot <- 250
```


## Running SIMEX

We call the function `simex_surv_outcome`
with the above arguments.

```{r}
res <- simex_surv_outcome(X_mat,
                          lambda,
                          time_star,
                          delta,
                          B,
                          var_error,
                          num_boot)
```

`res` is a dataframe containing the estimates 
and standard errors for each of the covariates in `X_mat`. 
Confidence intervals can then be calculated using 
your favorite method. 

