# Example 1 (simulation, small model):
set.seed(123)
sample <- sim.bin(3L, 100)
print(sample$coef)

data <- data.frame(sample$y, sample$x)

#   Estimate using glm
fit <- glm(Y ~ X1 + X2, data = data, family = binomial())
print(fit)

#   Estimate using 'ldt::estim.bin'
fit <- estim.bin(data = get.data(data = data,
                                 equations = list(Y ~ X1 + X2)),
                  linkFunc = "logit")
print(fit)
plot_data <- plot(fit, type = 1)
#   See 'plot.ldt.estim()' function documentation


# Example 2 (simulation, large model with PCA analysis):
sample <- sim.bin(30L, 100, probit = TRUE)
data <- data.frame(sample$y, sample$x)
colnames(data) <- c(colnames(sample$y),colnames(sample$x))
pca_options <- get.options.pca(ignoreFirst = 1, exactCount = 3)
fit <- estim.bin(data = get.data(cbind(sample$y, sample$x),
                                  endogenous = ncol(sample$y),
                                  addIntercept = FALSE),
                  linkFunc = "probit",
                  pcaOptionsX = pca_options)
print(fit)
plot_data <- plot(fit, type = 2)
