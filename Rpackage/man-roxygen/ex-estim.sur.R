# Example 1 (simulation, small model):
set.seed(123)
sample <- sim.sur(sigma = 2L, coef = 3L, nObs = 100)
print(sample$coef)
print(sample$sigma)

data <- data.frame(sample$y, sample$x)

#    Use systemfit to estimate:
exp_names <- paste0(colnames(sample$x), collapse = " + ")
fmla <- lapply(1:ncol(sample$y), function(i) as.formula(paste0("Y", i, " ~ -1 +", exp_names)))
fit <- systemfit::systemfit(fmla, data = data, method = "SUR")
print(fit)

#    Use 'ldt::estim.sur' function
fit <- estim.sur(data = get.data(cbind(sample$y, sample$x),
                                  endogenous = ncol(sample$y),
                                  addIntercept = FALSE))

# or, by using formula list:
fit <- estim.sur(data = get.data(data = data,
                                 equations = fmla,
                                 addIntercept = FALSE))

print(fit)
print(fit$estimations$sigma)
plot_data <- plot(fit, equation = 1)


# Example 2 (simulation, large model with significancy search):
num_obs <- 100
sample <- sim.sur(sigma = 2L, coef = 3L, nObs = num_obs)
print(sample$coef)

#   create irrelevant data:
num_x_ir <- 20
x_ir <- matrix(rnorm(num_obs * num_x_ir), ncol = num_x_ir)
data_x <- data.frame(sample$x, x_ir)
colnames(data_x) <- c(colnames(sample$x), paste0("z", 1:num_x_ir))

fit <- estim.sur(data = get.data(cbind(sample$y, data_x),
                                 endogenous = ncol(sample$y),
                                 addIntercept = FALSE),
                 searchSigMaxIter = 100,
                 searchSigMaxProb = 0.05)

print(fit$estimations$coefs)
# coefficient matrix, with lots of zero restrictions

# Example 3 (simulation, large model with PCA):
#   by using data of the previous example
fit <- estim.sur(data = get.data(cbind(sample$y, data_x),
                                 endogenous = ncol(sample$y),
                                 addIntercept = FALSE),
                 pcaOptionsX = get.options.pca(2,4))
print(fit$estimations$coefs)
#  coefficients are: intercept and the first exogenous variable and 4 PCs

