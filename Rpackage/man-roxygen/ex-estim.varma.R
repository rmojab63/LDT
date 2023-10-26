# Example 1 (simulation, ARMA):
num_eq <- 1L
num_ar <- 2L
num_ma <- 1L
num_exo <- 1L
sample <- sim.varma(num_eq, arList = num_ar, maList = num_ma, exoCoef = num_exo, nObs = 110)
# estimate:
fit <- estim.varma(data = get.data(cbind(sample$y, sample$x)[1:100,],
                                   endogenous = num_eq,
                                   newData = sample$x[101:110,, drop=FALSE]),
                   params = c(num_ar, 0, num_ma, 0, 0, 0),
                   maxHorizon = 10,
                   simFixSize = 5,
                   simHorizons = c(1:10))
print(fit)

pred <- predict(fit, actualCount = 10)
plot(pred, simMetric = "mape")

# split coefficient matrix:
get.varma.params(fit$estimations$coefs, numAR = num_ar, numMA = num_ma, numExo = num_exo)

# Example 2 (simulation, VARMA):
num_eq <- 3L
num_ar <- 2L
num_ma <- 1L
num_ma <- 1L
num_exo <- 2L
sample <- sim.varma(num_eq, arList = num_ar, maList = num_ma, exoCoef = num_exo, nObs = 110)
# estimate:
fit <- estim.varma(data = get.data(cbind(sample$y, sample$x)[1:100,],
                                   endogenous = num_eq,
                                   newData = sample$x[101:110,]),
                   params = c(num_ar, 0, num_ma, 0, 0, 0),
                   maxHorizon = 10,
                   simFixSize = 5,
                   simHorizons = c(1:10))

pred <- predict(fit, actualCount = 10)
plot(pred, simMetric = "mape")

# split coefficient matrix:
get.varma.params(fit$estimations$coefs, numAR = num_ar, numMA = num_ma, numExo = num_exo)
