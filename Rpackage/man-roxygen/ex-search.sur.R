num_y <- 2L # number of equations
num_x_r <- 3L # number of relevant explanatory variables
num_x_ir <-
  20 # (relatively large) number of irrelevant explanatory variables
num_obs = 100  # number of observations

# create random data
sample <- sim.sur(sigma = num_y, coef = num_x_r, nObs = num_obs)
x_ir <- matrix(rnorm(num_obs * num_x_ir), ncol = num_x_ir) # irrelevant data

# prepare data for estimation
data <- data.frame(sample$y, sample$x, x_ir)
colnames(data) <- c(colnames(sample$y), colnames(sample$x), paste0("z", 1:num_x_ir))

# Use systemfit to estimate and analyse:
exp_names <- paste0(colnames(data)[(num_y + 1):(length(colnames((data))))], collapse = " + ")
fmla <- lapply(1:num_y, function(i) as.formula(paste0("Y", i, " ~ -1 + ", exp_names)))
fit <- systemfit::systemfit(fmla, data = data, method = "SUR")
summary(fit)

# You can also use this package estimation function:
model2 <-
  estim.sur(data[, 1:num_y], data[, (num_y + 1):(length(data))], addIntercept = FALSE)
# format and print coefficients:
for (j in c(1:num_y)) {
  coefs2 <-
    data.frame(lapply(c(1:4), function(c)
      model2$estimations[[c]][, j]))
  colnames(coefs2) <-
    lapply(c(1:4), function(c)
      names(model2$estimations[c]))
  print(paste0("------------ Equation: ", j))
  print(coefs2)
}

# Alternatively, You can define a search process:
x_sizes = c(1:4) # assuming we know the number of relevant explanatory variables is less than 4
num_targets = 2
metric_options <-
  get.search.metrics(typesIn = c("sic")) # We use SIC for searching
search_res <-
  search.sur(
    sample$y,
    data[, 2:(length(data))],
    numTargets = num_targets,
    xSizes = x_sizes,
    metrics = metric_options
  )
# best model's explanatory indexes for the first and second variables:
print(search_res$sic$target1$model$bests$best1$exoIndices)
print(search_res$sic$target2$model$bests$best1$exoIndices)

# Use summary function to estimate the best models:
search_sum <-
  summary(search_res, y = data[, 1:num_y], x = data[, (num_y + 1):(length(data))])
# format and print coefficients:
for (j in c(1:num_targets)) {
  model3 <- search_sum$sic[[j]]$model$bests$best1
  coefs2 <-
    data.frame(lapply(c(1:4), function(c)
      model3$estimations[[c]][, j]))
  colnames(coefs2) <-
    lapply(c(1:4), function(c)
      names(model3$estimations[c]))
  print(paste0("------------ Equation: ", j))
  print(coefs2)
}

# Try a step-wise search (you can estimate larger models, faster):
x_sizes_steps = list(c(1, 2, 3), c(4), c(5))
counts_steps = c(NA, 10, 9)
search_items <- get.search.items(bestK = 10)
search_step_res <-
  search.sur.stepwise(
    y = data[, 1:num_y],
    x = data[, 2:(length(data))],
    xSizeSteps = x_sizes_steps,
    countSteps = counts_steps,
    metrics = metric_options,
    items = search_items
  )
print(search_step_res$sic$target1$model$bests$best1$exoIndices)
# Use summary like before.
