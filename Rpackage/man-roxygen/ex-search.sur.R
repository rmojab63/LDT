num_y <- 2L # number of equations
num_x_r <- 3L # number of relevant explanatory variables
num_x_ir <-
  10 # (relatively large) number of irrelevant explanatory variables
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
fit <- estim.sur(data = get.data(data, endogenous = num_y, addIntercept = FALSE))
print(fit)

# Alternatively, You can define an SUR model set:
x_sizes = c(1:3) # assuming we know the number of relevant explanatory variables is less than 3
num_targets = 2
metric_options <- get.search.metrics(typesIn = c("sic")) # We use SIC for searching
search_res <- search.sur(data = get.data(data, endogenous = num_y, addIntercept = FALSE),
                         combinations = get.combinations(numTargets = num_targets,
                                                         sizes = x_sizes,
                                                         innerGroups = list(c(1), c(2))),
                         metrics = metric_options)
print(search_res)

# Use summary function to estimate the best models:
search_sum <- summary(search_res)

# Print the best model:
print(search_sum$results[[2]]$value)
#   see 'estim.sur' function

# Using a step-wise search to build a larger model set:
x_sizes_steps = list(c(1, 2, 3), c(4))
counts_steps = c(NA, 7)
search_step_res <- search.sur(data = get.data(data, endogenous = num_y, addIntercept = FALSE),
                              combinations = get.combinations(numTargets = num_targets,
                                                              sizes = x_sizes_steps,
                                                              stepsNumVariables = counts_steps,
                                                              innerGroups = list(c(1,2))),
                              metrics = metric_options)
# combinations argument is different

print(search_step_res)

