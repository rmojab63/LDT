# Generate data from a logit model with 3 variables
sample <- sim.bin(3L, 100)

data <- data.frame(sample$y, sample$x)
colnames(data) <- c(colnames(sample$y),colnames(sample$x))

# Estimate using glm
fit <- glm(Y ~ X1 + X2, data = data, family = binomial())
# you can compare fit$coefficients and sample$coef

# Estimate using estim.bin in this package
fit1 <- estim.bin(sample$y, sample$x[,2:3, drop=FALSE], probType = "logit",
                 newX = sample$x[,2:3, drop=FALSE])

# Generate data from a probit model with specified coefficients
sample1 <- sim.bin(c(-1, 0.5, 2), 100, probit = TRUE)
data <- data.frame(sample1$y, sample1$x)
colnames(data) <- c(colnames(sample1$y),colnames(sample1$x))

# Estimate the probit model using glm
fit <- glm(Y ~ X1 + X2, data = data, family = binomial( link = "probit"))

# Estimate using estim.bin in this package
fit1 <- estim.bin(sample$y, sample$x[,2:3, drop=FALSE], probType = "probit")

