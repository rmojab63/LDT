num_y <- 2L
num_x <- 3L
n_obs = 100
data <- sim.sur(sigma = num_y, coef = num_x, nObs = n_obs)

# Use systemfit to estimate:
exp_names <- paste0(colnames(data$x), collapse = " + ")
fmla <- lapply(1:num_y, function(i) as.formula(paste0("Y", i, " ~ -1 +", exp_names)))
fit <- systemfit::systemfit(fmla, data = data.frame(data$y, data$x), method = "SUR")
summary(fit)

# Use estim.sur function in this package:
fit2 <- estim.sur(data$y, data$x, addIntercept = FALSE)
fit2$estimations$coefs # coefficient matrix
