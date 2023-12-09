# We simulate some data for this example:
set.seed(340)
n = 100
num_eq <- 3L
num_ar <- 2L
num_ma <- 1L
num_ma <- 1L
num_exo <- 2L
sample <- sim.varma(num_eq, arList = num_ar, maList = num_ma, exoCoef = num_exo, nObs = n)

# (relatively large) number of irrelevant explanatory variables:
num_y_ir <- 10
y_ir <- lapply(1:num_y_ir, function(x) rnorm(n))

# prepare data:
data <- data.frame(sample$y, y_ir, sample$x)
colnames(data) <- c(colnames(sample$y), paste0("w", 1:num_y_ir), colnames(sample$x))


# Create a VARMA model set:
y_sizes = 3 # assuming we know the number of relevant endogenous variables
metric_options <- get.search.metrics(typesIn = c("aic")) # We use SIC for searching
search_res <- search.varma(data = get.data(data, endogenous = num_eq + num_y_ir),
                           combinations = get.combinations(sizes = y_sizes,
                                                           numTargets = 3),
                           metrics = metric_options,
                           maxHorizon = 0)
print(search_res)


