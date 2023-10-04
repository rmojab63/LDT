# We simulate some data for this example:

# sample data:
n = 50 # number of observations
num_x_r <- 3L # number of relevant explanatory variables
num_x_ir <-
  20 # (relatively large) number of irrelevant explanatory variables
set.seed(340)
sample <- sim.bin(num_x_r, n)
x_ir <- lapply(1:num_x_ir, function(x) rnorm(n))

# prepare data:
data <- data.frame(sample$y, sample$x, x_ir)
colnames(data) <- c("Y", colnames(sample$x), paste0("z", 1:num_x_ir))

# Use glm function to estimate and analyse:
model1 <- glm(Y ~ . - Y, data = data, family = binomial())
summary(model1)

# You can also use this package estimation function:
model2 <- estim.bin(sample$y, data[, 3:ncol(data), drop = FALSE])
# format and print coefficients:
coefs2 <- data.frame(model2$estimations[1:4])
colnames(coefs2) <- names(model2$estimations)[1:4]
print(coefs2)

# Alternatively, You can define a search process:
x_sizes = c(1:4) # assuming we know the number of relevant explanatory variables is less than 4
metric_options <-
  get.search.metrics(typesIn = c("sic")) # We use SIC for searching
search_res <- search.bin(sample$y, data[, 3:ncol(data)],
                        xSizes = x_sizes, metrics = metric_options)
print(search_res$sic$target1$model$bests$best1$exogenous) # print best model's explanatory indexes

# Use summary function to estimate the best model:
search_sum <-
  summary(search_res, y = sample$y, x = data[, 3:ncol(data)])
# format and print coefficients:
model3 <- search_sum$sic$target1$model$bests$best1
coefs3 <- data.frame(model3$estimations[1:4])
colnames(coefs3) <- names(model3$estimations)[1:4]
print(coefs3)

# Try a step-wise search (you can estimate larger models, faster):
x_sizes_steps = list(c(1, 2, 3), c(4), c(5))
counts_steps = c(NA, 10, 9)
search_items <- get.search.items(bestK = 10)
search_step_res <- search.bin.stepwise(
  sample$y,
  data[, 3:ncol(data)],
  xSizeSteps = x_sizes_steps,
  countSteps = counts_steps,
  metrics = metric_options,
  items = search_items
)
print(search_step_res$sic$target1$model$bests$best1$exogenous)
# Use summary like before.
