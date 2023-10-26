# We simulate some data for this example:
# sample data:
n = 50 # number of observations
num_x_r <- 3L # number of relevant explanatory variables
num_x_ir <- 20 # (relatively large) number of irrelevant explanatory variables
set.seed(340)
sample <- sim.bin(num_x_r, n)
x_ir <- lapply(1:num_x_ir, function(x) rnorm(n))

# prepare data:
data <- data.frame(sample$y, sample$x, x_ir)
colnames(data) <- c("Y", colnames(sample$x), paste0("z", 1:num_x_ir))

# Use glm function to estimate and analyse:
fit <- glm(Y ~ . - Y, data = data, family = binomial())
summary(fit)

# You can also use this package estimation function:
data0 <- get.data(data,
                equations = list(Y ~ . - Y),
                addIntercept = FALSE)
fit <- estim.bin(data = data0)
# format and print coefficients:
print(fit)

# Alternatively, You can define a binary choice model set:
x_sizes = c(1:3) # assuming we know the number of relevant explanatory variables is less than 3
metric_options <- get.search.metrics(typesIn = c("sic")) # We use SIC for searching
search_res <- search.bin(data = data0,
                         combinations = get.combinations(sizes = x_sizes),
                         metrics = metric_options)
print(search_res)

# Use summary function to estimate the best model:
search_sum <- summary(search_res, y = sample$y, x = data[,3:ncol(data)])

# format and print coefficients:
s_fit <- summary(search_res)
print(s_fit$results[[1]]$value)

# Try a step-wise search for creating a larger model set:
search_res <- search.bin(data = data0,
                         combinations = get.combinations(
                           sizes = list(c(1, 2, 3), c(4)),
                           stepsNumVariables = c(NA, 7)),
                         metrics = metric_options)
# combinations argument is different

print(search_res)
# Use summary like before.

