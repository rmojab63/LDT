# Example 1 (simulation, small model):
sample <- sim.sur(sigma = 2L, coef = 3L, nObs = 100)
#    Use systemfit to estimate:
exp_names <- paste0(colnames(sample$x), collapse = " + ")
fmla <- lapply(1:ncol(sample$y), function(i) as.formula(paste0("Y", i, " ~ -1 +", exp_names)))
fit <- systemfit::systemfit(fmla, data = data.frame(sample$y, sample$x), method = "SUR")
#    Use estim.sur function in this package:
fit2 <- estim.sur(sample$y, sample$x, addIntercept = FALSE)
coefs <- fit2$estimations$coefs # coefficient matrix
table <- get.coefs.table(fit2)


# Example 2 (simulation, large model with significancy search):
num_obs <- 100
num_x_ir <- 20
sample <- sim.sur(sigma = 2L, coef = 3L, nObs = num_obs)
x_ir <- matrix(rnorm(num_obs * num_x_ir), ncol = num_x_ir) # irrelevant sample
data_x <- data.frame(sample$x, x_ir)
colnames(data_x) <- c(colnames(sample$x), paste0("z", 1:num_x_ir))

fit3 <- estim.sur(sample$y, data_x, addIntercept = FALSE,
                  searchSigMaxIter = 100, searchSigMaxProb = 0.05)
coefs <- fit2$estimations$coefs # coefficient matrix, probably with lots of zero restrictions

# Example 3 (simulation, large model with PCA):
# I use the previous example's data
fit3 <- estim.sur(sample$y, data_x, addIntercept = FALSE,
                  pcaOptionsX = get.options.pca(2,4))
coefs <- fit3$estimations$coefs # has intercept and the first exogenous variable and 4 PCs
table <- get.coefs.table(fit3, expList = 2) # doesn't display PC coefficients
