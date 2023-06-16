
#' Search for Best SUR Models
#'
#' Use this function to create a Seemingly Unrelated Regression model set and search for the best models (and other information) based on in-sample and out-of-sample evaluation measures.
#'
#' @param y A matrix of endogenous data with variables in the columns.
#' @param x A matrix of exogenous data with variables in the columns.
#' @param numTargets An integer for the number of targets variables.
#' If for example 2, the first two variable in the first columns of \code{y} will be target.
#' Information is saved just for the target variables.
#' It must be positive and cannot be larger than the number of endogenous variables.
#' @param xSizes An integer vector specifying the number of exogenous variables in the regressions.
#' E.g., \code{c(1,2)} means the model set contains all regressions with 1 and 2 exogenous variables.
#' If \code{NULL}, \code{c(1)} is used.
#' @param xPartitions A list of integer vectors that partition the indexes of the exogenous variables.
#' No regression is estimated with two variables in the same partition.
#' If \code{NULL}, each variable is placed in its own group, and the size of the model set is maximized.
#' @param numFixXPartitions Number of partitions at the beginning of \code{xPartitions} to be included in all regressions.
#' @param yGroups A list of integer vectors that determine different combinations of the indexes of the endogenous variables to be used as endogenous variables in the SUR regressions.
#' @param searchSigMaxIter An integer for the maximum number of iterations in searching for significant coefficients. Use 0 to disable the search.
#' @param searchSigMaxProb A number for the maximum value of type I error to be used in searching for significant coefficients. If p-value is less than this, it is interpreted as significant.
#' @param measureOptions A list of options for measuring performance.
#' Use [get.options.measure] function to get them.
#' @param modelCheckItems A list of options for excluding a subset of the model set.
#' See and use [get.items.modelcheck] function to get them.
#' @param searchItems A list of options for specifying the purpose of the search.
#' See and use [get.items.search] function to get them.
#' @param searchOptions A list of extra options for performing the search.
#' See and use [get.options.search] function to get them.
#'
#' @return A nested list with the following members:
#' \tabular{ll}{
#' \code{counts} \tab Information about the expected number of models, number of estimated models, failed estimations, and some details about the failures.\cr
#' \code{...} \tab Results reported separately for each measure, then for each target variable, then for each requested type of output. This part of the output is highly nested and items are reported based on the arguments of the search.\cr
#' \code{info} \tab General information about the search process, some arguments, elapsed time, etc.
#' }
#'
#' @export
#'
#' @examples
#' # We simulate some data for this example:
#'
#' n = 50 # number of observations
#' num_x_relevant <- 3 # number of relevant explanatory variables
#' num_x_irrelevant <- 20 # (relatively large) number of irrelevant explanatory variables
#' num_y <- 3 # number of endogenous data
#'
#' # create relevant and irrelevant explanatory variables:
#' x_relevant <- matrix(rnorm(n * num_x_relevant), ncol = num_x_relevant)
#' x_irrelevant <- matrix(rnorm(n * num_x_irrelevant), ncol = num_x_irrelevant)
#'
#' # calculate the dependent variable
#' beta <- matrix(rnorm((num_x_relevant + 1) * num_y), ncol = num_y) # coefficients (including intercepts)
#' Sigma <- crossprod(matrix(rnorm(num_y^2), ncol = num_y)) # errors covariance
#' errors <- MASS::mvrnorm(n, mu = rep(0, num_y), Sigma = Sigma)
#' y <- cbind(rep(1, n), x_relevant) %*% beta + errors
#'
#' # prepare data for estimation
#' data <- data.frame(y, x_relevant, x_irrelevant)
#' colnames(data) <- c(paste0("y", 1:num_y), paste0("x", 1:num_relevant), paste0("z", 1:num_irrelevant))
#'
#' # Use systemfit to estimate and analyse:
#' exp_names <- paste0(colnames(data)[(num_y+1):(length(colnames((data))))], collapse = " + ")
#' fmla <- lapply(1:num_y, function(i) as.formula(paste0("y", i, " ~ 1 + ", exp_names)))
#' fit <- systemfit::systemfit(fmla, data = data, method = "SUR")
#' summary(fit)
#'
#' # You can also use this package estimation function:
#' model2 <- estim.sur(data[,1:num_y], data[,(num_y+1):(length(data))])
#' # format and print coefficients:
#' for (j in c(1:num_y)){
#'   coefs2 <- data.frame(lapply(c(1:4), function(c)model2$estimations[[c]][,j]))
#'   colnames(coefs2) <- lapply(c(1:4), function(c)names(model2$estimations[c]))
#'   print(paste0("------------ Equation: ", j))
#'   print(coefs2)
#' }
#'
#' # Alternatively, You can define a search process:
#' x_sizes = c(1:4) # assuming we know the number of relevant explanatory variables is less than 4
#' num_targets = 2
#' measure_options <- get.options.measure(typesIn = c("sic")) # We use SIC for searching
#' search_res <- search.sur(y, data[,2:(length(data))], numTargets = num_targets,
#'                          xSizes = x_sizes, measureOptions = measure_options)
#' print(search_res$sic$target1$model$bests$best1$exoIndices) # best model's explanatory indexes for the first variable
#' print(search_res$sic$target2$model$bests$best1$exoIndices) # best model's explanatory indexes for the second variable
#'
#' # Use summary function to estimate the best models:
#' search_sum <- summary(search_res, y = data[,1:num_y], x = data[,(num_y+1):(length(data))])
#' # format and print coefficients:
#' for (j in c(1:num_targets)){
#'   model3 <- search_sum$sic[[j]]$model$bests$best1
#'   coefs2 <- data.frame(lapply(c(1:4), function(c)model3$estimations[[c]][,j]))
#'   colnames(coefs2) <- lapply(c(1:4), function(c)names(model3$estimations[c]))
#'   print(paste0("------------ Equation: ", j))
#'   print(coefs2)
#' }
#'
#' # Try a step-wise search (you can estimate larger models, faster):
#' x_sizes_steps = list(c(1, 2, 3), c(4), c(5))
#' counts_steps = c(NA, 10, 9)
#' search_items <- get.items.search(bestK = 10)
#' search_step_res <- search.sur.stepwise(y = data[,1:num_y], x = data[,2:(length(data))],
#'                                        xSizesSteps = x_sizes_steps, countSteps = counts_steps,
#'                                        measureOptions = measure_options,
#'                                        searchItems = search_items)
#' print(search_step_res$sic$target1$model$bests$best1$exoIndices)
#' # Use summary like before.
#'
search.sur <- function(y, x, numTargets = 1, xSizes = NULL,
                      xPartitions = NULL, numFixXPartitions = 0,
                      yGroups = NULL, searchSigMaxIter = 0,
                      searchSigMaxProb = 0.1, measureOptions = get.options.measure(),
                      modelCheckItems = get.items.modelcheck(), searchItems = get.items.search(),
                      searchOptions = get.options.search()){

  y = as.matrix(y)
  x = as.matrix(x)
  numTargets = as.integer(numTargets)
  xSizes = if (is.null(xSizes)) NULL else as.integer(xSizes)
  numFixXPartitions = as.integer(numFixXPartitions)
  searchSigMaxIter = as.integer(searchSigMaxIter)
  searchSigMaxProb = as.numeric(searchSigMaxProb)

  if (is.null(xPartitions) == FALSE){
    xPartitions = as.list(xPartitions)
    for (i in c(1:length(xPartitions)))
      xPartitions[[i]] = as.integer(xPartitions[[i]])
  }

  if (is.null(yGroups) == FALSE){
    yGroups = as.list(yGroups)
    for (i in c(1:length(yGroups)))
      yGroups[[i]] = as.integer(yGroups[[i]])
  }

  if (is.null(measureOptions))
    measureOptions = get.options.measure()
  else
    measureOptions <- CheckMeasureOptions(measureOptions)

  if (is.null(modelCheckItems))
    modelCheckItems = get.items.modelcheck()
  else
    modelCheckItems <- CheckModelCheckItems(modelCheckItems)

  if (is.null(searchItems))
    searchItems = get.items.search()
  else
    searchItems <- CheckSearchItems(searchItems)

  if (is.null(searchOptions))
    searchOptions = get.options.search()
  else
    searchOptions <- CheckSearchOptions(searchOptions)

  res <- .SearchSur(y, x, numTargets, xSizes, xPartitions, numFixXPartitions,
                    yGroups, searchSigMaxIter, searchSigMaxProb, measureOptions,
                    modelCheckItems, searchItems, searchOptions)
  res
}


#' Estimate a SUR Model
#'
#' Use this function to estimate a Seemingly Unrelated Regression model.
#'
#' @param y A matrix of endogenous data with variables in the columns.
#' @param x A matrix of exogenous data with variables in the columns.
#' @param addIntercept If \code{TRUE}, intercept is added automatically to \code{x}.
#' @param searchSigMaxIter An integer for the maximum number of iterations in searching for significant coefficients. Use 0 to disable the search.
#' @param searchSigMaxProb A number for the maximum value of type I error to be used in searching for significant coefficients. If p-value is less than this, it is interpreted as significant.
#' @param restriction A \code{km x q} matrix of restrictions where \code{m=ncols(y)}, \code{k=ncols(x)} and \code{q} is the number of unrestricted coefficients.
#' @param newX A matrix with new exogenous data to be used in the projections. Its number of columns must be equal to \code{x}. It can be \code{NULL}.
#' @param pcaOptionsY A list of options to use principal components of the \code{y}, instead of the actual values. Set \code{NULL} to disable. Use [get.options.pca()] for initialization.
#' @param pcaOptionsX A list of options to use principal components of the \code{x}, instead of the actual values. Set \code{NULL} to disable. Use [get.options.pca()] for initialization.
#' @param simFixSize An integer that determines the number of pseudo out-of-sample simulations. Use zero to disable the simulation.
#' @param simTrainFixSize An integer representing the number of data points in the training sample in the pseudo out-of-sample simulation. If zero, \code{trainRatio} will be used.
#' @param simTrainRatio A number representing the size of the training sample relative to the available size, in the pseudo out-of-sample simulation. It is effective if \code{trainFixSize} is zero.
#' @param simSeed A seed for the random number generator. Use zero for a random value.
#' @param simMaxConditionNumber A number for the maximum value for the condition number in the simulation.
#' @param printMsg Set to \code{TRUE} to enable printing some details.
#'
#' @return A nested list with the following items:
#' \tabular{ll}{
#' \code{counts} \tab Information about different aspects of the estimation such as the number of observation, number of exogenous variables, etc.\cr
#' \code{estimations} \tab Estimated coefficients, standard errors, z-statistics, p-values, etc.\cr
#' \code{measures} \tab Value of different goodness of fit and out-of-sample performance measures. \cr
#' \code{projections} \tab Information on the projected values, if \code{newX} is provided.\cr
#' \code{info} \tab Some other general information.
#' }
#'
#' @details
#' The main purpose of exporting this method is to show the inner calculations of the search process in [search.sur] function. See the details of this function for more information.
#'
#' @export
#' @examples
#' See the example in the 'search.sur' function.
#'
#' @seealso [search.sur], [search.sur.stepwise]
estim.sur <- function(y, x, addIntercept = TRUE,
                     searchSigMaxIter = 0, searchSigMaxProb = 0.1,
                     restriction = NULL, newX = NULL,
                     pcaOptionsY = NULL, pcaOptionsX = NULL,
                     simFixSize = 0, simTrainFixSize = 0,
                     simTrainRatio = 0.75, simSeed = 0,
                     simMaxConditionNumber = Inf, printMsg = FALSE){
  y = as.matrix(y)
  x = as.matrix(x)
  addIntercept = as.logical(addIntercept)
  searchSigMaxIter = as.integer(searchSigMaxIter)
  searchSigMaxProb = as.numeric(searchSigMaxProb)
  restriction = if (is.null(restriction)) NULL else as.matrix(restriction)
  newX = if (is.null(newX)) NULL else as.matrix(newX)
  simFixSize = as.integer(simFixSize)
  simTrainRatio = as.numeric(simTrainRatio)
  simTrainFixSize = as.integer(simTrainFixSize)
  simSeed = as.integer(simSeed)
  simMaxConditionNumber = as.numeric(simMaxConditionNumber)
  printMsg = as.logical(printMsg)

  if (is.null(pcaOptionsY) == FALSE)
    pcaOptionsY = CheckPcaOptions(as.list(pcaOptionsY))
  if (is.null(pcaOptionsX) == FALSE)
    pcaOptionsX = CheckPcaOptions(as.list(pcaOptionsX))

  res <- .EstimSur(y, x, addIntercept,
                   searchSigMaxIter, searchSigMaxProb,
                   restriction, newX, pcaOptionsY, pcaOptionsX,
                   simFixSize, simTrainRatio,
                   simTrainFixSize, simSeed,
                   simMaxConditionNumber, printMsg)
  res
}

# get estimation from search result
GetEstim_sur <- function(searchRes, endoIndices,
                         exoIndices, y, x, printMsg, ...) {
  y <- y[, endoIndices, drop = FALSE]
  x <- if (is.null(exoIndices) || is.null(x)) {
    NULL
  } else {
    x[, c(exoIndices), drop = FALSE]
  }

  M <- estim.sur(
    y = y,
    x = x,
    addIntercept = FALSE,
    searchSigMaxIter = searchRes$info$searchSigMaxIter,
    searchSigMaxProb = searchRes$info$searchSigMaxProb,
    restriction = NULL,
    newX = NULL,
    pcaOptionsY = NULL,
    pcaOptionsX = NULL,
    simFixSize = searchRes$info$measureOptions$simFixSize,
    simTrainRatio = searchRes$info$measureOptions$trainRatio,
    simTrainFixSize = searchRes$info$measureOptions$trainFixSize,
    simSeed = abs(searchRes$info$measureOptions$seed),
    simMaxConditionNumber = searchRes$info$modelCheckItems$maxConditionNumber,
    printMsg = printMsg
  )

  return(M)
}


#' Step-wise Search for Best SUR Models
#'
#' For a large model set, use this function to find the best seemingly unrelated regression models.
#' It selects a subset of variables from smaller models and moves to the bigger ones.
#'
#' @param y A matrix of endogenous data with variables in the columns.
#' @param x A matrix of exogenous data with variables in the columns.
#' @param xSizesSteps A list of model dimensions to be estimated in each step.
#' Its size determines the number of steps.
#' @param countSteps A vector to determine the number of variables to be used in each step.
#' \code{NA} means all variables. Variables are selected based on best estimations.
#' All variables in the best models (all measures and targets) are selected until the corresponding suggested number is reached.
#' Select an appropriate value for \code{bestK} in the options.
#' @param savePre A directory for saving and loading the progress.
#' Each step's result is saved in a file (name=\code{paste0(savePre,i)} where \code{i} is the index of the step.
#' @param ... other arguments to pass to [search.sur] function such as the \code{w} argument.
#' Note that \code{xSizes} is ineffective here.
#'
#' @return Similar to [search.sur] function.
#' @export
#'
#' @examples
#' See the example in the 'search.sur' function.
#'
#' @seealso [search.sur], [estim.sur]
search.sur.stepwise <- function(y, x, xSizesSteps = list(c(1, 2), c(3, 4), c(5), c(6:10)),
                        countSteps = c(NA, 40, 30, 20),
                        savePre = NULL, ...) {
  Search_s("sur", x, xSizesSteps, countSteps, savePre, y = y, ...)
}
