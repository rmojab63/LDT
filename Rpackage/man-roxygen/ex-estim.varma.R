# Example 1 (simulation, ARMA):
num_eq <- 1L
num_ar <- 2L
num_ma <- 1L
num_exo <- 1L
sample <- sim.varma(num_eq, arList = num_ar, maList = num_ma, exoCoef = num_exo, nObs = 100)
# estimate:
fit <- estim.varma(sample$y, sample$x, params = c(num_ar, 0, num_ma, 0, 0, 0))
# split coefficient matrix:
get.varma.params(fit$estimations$coefs, numAR = num_ar, numMA = num_ma, numExo = num_exo)

# Example 2 (simulation, VARMA):
num_eq <- 3L
num_ar <- 2L
num_ma <- 1L
num_ma <- 1L
num_exo <- 2L
sample <- sim.varma(num_eq, arList = num_ar, maList = num_ma, exoCoef = num_exo, nObs = 100)
# estimate:
fit <- estim.varma(sample$y, sample$x, params = c(num_ar, 0, num_ma, 0, 0, 0))
# split coefficient matrix:
get.varma.params(fit$estimations$coefs, numAR = num_ar, numMA = num_ma, numExo = num_exo)
