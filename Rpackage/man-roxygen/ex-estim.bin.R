# Example 1 (simulation, small model):
sample <- sim.bin(3L, 100)
data <- data.frame(sample$y, sample$x)
#   Estimate using glm
fit <- glm(Y ~ X1 + X2, data = data, family = binomial())
#   Estimate using estim.bin in this package
fit1 <- estim.bin(sample$y, sample$x[,2:ncol(sample$x), drop=FALSE], linkFunc = "logit",
                 newX = sample$x[,2:3, drop=FALSE])

# Example 2 (simulation, large model with PCA analysis):
sample <- sim.bin(30L, 100, probit = TRUE)
data <- data.frame(sample$y, sample$x)
colnames(data) <- c(colnames(sample$y),colnames(sample$x))
pca_options <- get.options.pca(ignoreFirst = 1, exactCount = 3)
fit1 <- estim.bin(sample$y, sample$x[,2:ncol(sample$x), drop=FALSE], linkFunc = "probit",
                 pcaOptionsX = pca_options)

