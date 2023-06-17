
#' Search for Best VARMA Models
#'
#' Use this function to create a Vector Autoregressive Moving Average model set and search for the best models (and other information) based on in-sample and out-of-sample evaluation measures.
#'
#' @param y A matrix of endogenous data with variables in the columns.
#' @param x A matrix of exogenous data with variables in the columns. It can be \code{NULL}.
#' @param numTargets An integer representing the number of target variables.
#' For example, if 2, the first two variables in the first columns of \code{y} will be targets.
#' Information is saved just for the target variables.
#' It must be positive and cannot be larger than the number of endogenous variables.
#' @param ySizes An integer vector specifying the number of endogenous variables in the regressions.
#' For example, \code{c(1,2)} means the model set contains all regressions with 1 and 2 endogenous variables.
#' If \code{NULL}, \code{c(1)} is used.
#' @param yPartitions A list of integer vectors that partition the indexes of the endogenous variables.
#' No regression is estimated with two variables in the same partition.
#' If \code{NULL}, each variable is placed in its own partition, and the size of the model set is maximized.
#' @param xGroups A list of integer vectors that determine different combinations of the indexes of the exogenous variables to be used as exogenous variables in the SUR regressions.
#' @param maxParams An integer vector that determines the maximum values for the parameters of the VARMA model: \code{(p,d,q,P,D,Q)}. If \code{NULL}, \code{c(2,0,0,0,0,0)} is used.
#' @param seasonsCount An integer value representing the number of observations per unit of time.
#' @param maxHorizon An integer value representing the maximum value for the prediction horizon if \code{type1} is \code{TRUE} in the \code{modelCheckItems} argument. Also, it is used as the maximum prediction horizon in checking predictions.
#' @param newX A matrix of new exogenous data for out-of-sample prediction. It must have the same number of columns as the \code{x} argument.
#' @param simUsePreviousEstim If \code{TRUE}, parameters are initialized only in the first step of the simulation. The initial values of the n-th simulation (with one more observation) are the estimations from the previous step.
#' @param olsStdMultiplier A number used as a multiplier for the standard deviation of OLS, used for restricting maximum likelihood estimation.
#' @param lbfgsOptions A list containing L-BFGS optimization options.
#' Use [get.options.lbfgs] function for initialization.
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
#' Note that the output does not contain any estimation results,
#' but minimum required data to estimate the models (Use \code{summary()} function to get the estimation).
#'
#' @export
#'
#' @examples
#' # We simulate some data for this example:
#'
#' n = 100 # number of observations
#' num_y_relevant <- 2 # number of relevant dependent variables
#' num_y_irrelevant <- 5 # (relatively large) number of irrelevant dependent variables
#' lags_ar <- 1 # number of AR lags
#' lags_ma <- 2 # number of MA lags
#' exoCoef <- 0 # (number of) exogenous data
#' intercept = TRUE
#'
#' # generate dynamic data
#' sim_data <- sim.varma(num_y_relevant, lags_ar, lags_ma, exoCoef, n, n/2, intercept)
#'
#' # prepare endogenous variables by adding irrelevant data
#' y <- cbind(sim_data$y, matrix(rnorm(n * num_y_irrelevant), ncol = num_y_irrelevant))
#' colnames(y) <- c(paste0("Y",c(1:num_y_relevant)), paste0("W",c(1:num_y_irrelevant)))
#'
#' # You might use MTS package for estimation and analysis:
#' fixed_mas <- lapply(seq_len(lags_ma),
#'                     function(x) matrix(!diag(num_y_relevant),
#'                     num_y_relevant, num_y_relevant)) #for diagonal MAs
#' # the following might fail because it might not support restrictions on MA or
#' # it might not support exogenous data
#' varma_fit <- try(MTS::VARMA(y, p = lags_ar, q = lags_ma,
#'                             fixed = list(ar = NULL, ma = fixed_mas), include.mean = intercept))
#'
#'
#' # You can also use this package estimation function:
#' model2 <- estim.varma(y = y, x = sim_data$x, params = c(lags_ar,0,lags_ma,0,0,0))
#'
#'
#' # Alternatively, You can define a search process:
#' y_sizes = c(1:2) # assuming we know the number of relevant dependent variables is less than 2
#' max_params <- c(lags_ar,0,lags_ma,0,0,0) # assuming we know something about the dynamics
#' num_targets = 2
#' measure_options <- get.options.measure(typesIn = c("sic")) # We use SIC for searching
#' search_res <- search.varma(y = y, x = sim_data$x, maxParams = max_params,
#'                            numTargets = num_targets,
#'                            ySizes = y_sizes, measureOptions = measure_options)
#' # best model's dependent indexes for the first two variable:
#' print(search_res$sic$target1$model$bests$best1$depIndices)
#' print(search_res$sic$target2$model$bests$best1$depIndices)
#'
#' # Use summary function to estimate the best models:
#' search_sum <- summary(search_res, y = y, x = sim_data$x)
#'
#'
#' # Try a step-wise search (you can estimate larger models, faster):
#' x_sizes_steps = list(c(1), c(2))
#' counts_steps = c(NA, 10)
#' search_items <- get.items.search(bestK = 20)
#' search_step_res <- search.varma.stepwise(y = y, x = sim_data$x,
#'                      maxParams = max_params, numTargets = num_targets,
#'                      ySizeSteps = x_sizes_steps, countSteps = counts_steps,
#'                      measureOptions = measure_options,
#'                      searchItems = search_items)
#' print(search_step_res$sic$target1$model$bests$best1$depIndices)
#'
#' @seealso [estim.varma], [search.varma.stepwise]
search.varma <- function(y, x = NULL, numTargets = 1,
                        ySizes = NULL, yPartitions = NULL,
                        xGroups = NULL, maxParams = c(1,0,0,0,0,0),
                        seasonsCount = 0, maxHorizon = 0,
                        newX = NULL, simUsePreviousEstim = TRUE,
                        olsStdMultiplier = 2.0, lbfgsOptions = get.options.lbfgs(),
                        measureOptions = get.options.measure(),
                        modelCheckItems = get.items.modelcheck(), searchItems = get.items.search(),
                        searchOptions = get.options.search()){

  y = as.matrix(y)
  x = if (is.null(x)) NULL else as.matrix(x)
  numTargets = as.integer(numTargets)
  ySizes = if (is.null(ySizes)) c(1L) else as.integer(ySizes)
  maxParams = if (is.null(maxParams)) c(2L,0L,0L,0L,0L,0L) else as.integer(maxParams)
  seasonsCount = as.integer(seasonsCount)
  maxHorizon = as.integer(maxHorizon)
  newX = if (is.null(newX)) NULL else as.matrix(newX)
  simUsePreviousEstim = as.logical(simUsePreviousEstim)
  olsStdMultiplier = as.numeric(olsStdMultiplier)

  if (is.null(yPartitions) == FALSE){
    yPartitions = as.list(yPartitions)
    for (i in c(1:length(yPartitions)))
      yPartitions[[i]] = as.integer(yPartitions[[i]])
  }

  if (is.null(xGroups) == FALSE){
    xGroups = as.list(xGroups)
    for (i in c(1:length(xGroups)))
      xGroups[[i]] = as.integer(xGroups[[i]])
  }

  if (is.null(lbfgsOptions))
    lbfgsOptions = get.options.lbfgs()
  else
    lbfgsOptions = CheckLbfgsOptions(lbfgsOptions)

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

  res <- .SearchVarma(y, x, numTargets, ySizes, yPartitions,
                      xGroups, maxParams, seasonsCount, maxHorizon,
                      newX, simUsePreviousEstim,
                      olsStdMultiplier, lbfgsOptions,
                      measureOptions,modelCheckItems,searchItems,searchOptions)
  res
}

#' Estimate a VARMA Model
#'
#' Use this function to estimate a Vector Autoregressive Moving Average model.
#'
#' @param y A matrix of endogenous data with variables in the columns.
#' @param x A matrix of exogenous data with variables in the columns. It can be \code{NULL}.
#' @param params An integer vector that determines the values for the parameters of the VARMA model: \code{(p,d,q,P,D,Q)}. If \code{NULL}, \code{c(1,0,0,0,0,0)} is used.
#' @param seasonsCount An integer value representing the number of observations per unit of time.
#' @param addIntercept If \code{TRUE}, an intercept is automatically added to x.
#' @param lbfgsOptions A list containing L-BFGS optimization options.
#' Use [get.options.lbfgs] function for initialization.
#' @param olsStdMultiplier A number used as a multiplier for the standard deviation of OLS, used for restricting maximum likelihood estimation.
#' @param pcaOptionsY A list of options to use principal components of \code{y}, instead of the actual values. Set to \code{NULL} to disable. Use [get.options.pca()] for initialization.
#' @param pcaOptionsX A list of options to use principal components of \code{x}, instead of the actual values. Set to \code{NULL} to disable. Use [get.options.pca()] for initialization.
#' @param maxHorizon An integer representing the maximum prediction horizon. Set to zero to disable prediction.
#' @param newX A matrix containing new exogenous variables to be used in predictions. Its columns must be the same as \code{x}.
#' @param simFixSize An integer that determines the number of pseudo out-of-sample simulations. Use zero to disable simulation.
#' @param simHorizons An integer vector representing the prediction horizons to be used in pseudo out-of-sample simulations. See also [get.options.measure()].
#' @param simUsePreviousEstim If \code{TRUE}, parameters are initialized only in the first step of the simulation. The initial values of the n-th simulation (with one more observation) are the estimations from the previous step.
#' @param simMaxConditionNumber A number representing the maximum value for the condition number in simulation.
#' @param printMsg Set to \code{TRUE} to enable printing some details.
#'
#' @return A nested list with the following items:
#' \tabular{ll}{
#' \code{counts} \tab Information about different aspects of the estimation such as the number of observation, number of exogenous variables, etc.\cr
#' \code{estimations} \tab Estimated coefficients, standard errors, z-statistics, p-values, etc.\cr
#' \code{measures} \tab Value of different goodness of fit and out-of-sample performance measures. \cr
#' \code{prediction} \tab Information on the predicted values.\cr
#' \code{simulation} \tab Information on the simulations. \cr
#' \code{info} \tab Some other general information.
#' }
#'
#' @details
#' The main purpose of exporting this method is to show the inner calculations of the search process in [search.varma] function. See the details of this function for more information.
#'
#' @export
#' @examples
#' See the example in the 'search.varma' function.
#' @seealso [search.varma], [search.varma.stepwise]
estim.varma <- function(y, x = NULL, params = NULL,
                       seasonsCount = 0, addIntercept = TRUE,
                       lbfgsOptions = get.options.lbfgs(), olsStdMultiplier = 2,
                       pcaOptionsY = NULL, pcaOptionsX = NULL,
                       maxHorizon = 0, newX = NULL, simFixSize = 0,
                       simHorizons = NULL, simUsePreviousEstim = TRUE,
                       simMaxConditionNumber = Inf, printMsg = FALSE){

  y = as.matrix(y)
  x = if (is.null(x)) NULL else as.matrix(x)
  params = if (is.null(params)) c(1L,0L,0L,0L,0L,0L) else as.integer(params)
  seasonsCount = as.integer(seasonsCount)
  addIntercept = as.logical(addIntercept)
  olsStdMultiplier = as.numeric(olsStdMultiplier)
  maxHorizon = as.integer(maxHorizon)
  newX = if (is.null(newX)) NULL else as.matrix(newX)
  simFixSize = as.integer(simFixSize)
  simHorizons = if (is.null(simHorizons)) NULL else as.integer(simHorizons)
  simUsePreviousEstim = as.logical(simUsePreviousEstim)
  simMaxConditionNumber = as.numeric(simMaxConditionNumber)
  printMsg = as.logical(printMsg)

  if (is.null(lbfgsOptions))
    lbfgsOptions = get.options.lbfgs()
  else
    lbfgsOptions = CheckLbfgsOptions(lbfgsOptions)

  if (is.null(pcaOptionsY) == FALSE)
    pcaOptionsY = CheckPcaOptions(as.list(pcaOptionsY))
  if (is.null(pcaOptionsX) == FALSE)
    pcaOptionsX = CheckPcaOptions(as.list(pcaOptionsX))

  res <- .EstimVarma(y, x, params, seasonsCount, addIntercept,
                     lbfgsOptions, olsStdMultiplier,
                     pcaOptionsY, pcaOptionsX,
                     maxHorizon, newX, simFixSize, simHorizons,
                     simUsePreviousEstim, simMaxConditionNumber, printMsg)
  res
}


# get estimation from search result
GetEstim_varma <- function(searchRes, endoIndices,
                           exoIndices, y, x, printMsg,
                           params, newX = NULL, ...) {
  M <- estim.varma(y[, endoIndices, drop = FALSE],
                  x = if (is.null(exoIndices) || is.null(x)) {
                    NULL
                  } else {
                    x[, exoIndices, drop = FALSE]
                  },
                  params = params,
                  seasonsCount = searchRes$info$seasonsCount,
                  addIntercept = FALSE,
                  lbfgsOptions = searchRes$info$lbfgsOptions,
                  olsStdMultiplier = searchRes$info$olsStdMultiplier,
                  pcaOptionsY = NULL,
                  pcaOptionsX = NULL,
                  maxHorizon = searchRes$info$maxHorizon,
                  newX = if (is.null(exoIndices) || is.null(newX)) {
                    NULL
                  } else {
                    as.matrix(newX[, exoIndices])
                  },
                  simFixSize = searchRes$info$measureOptions$simFixSize,
                  simHorizons = searchRes$info$measureOptions$simHorizons,
                  simUsePreviousEstim = searchRes$info$simUsePreviousEstim,
                  simMaxConditionNumber = searchRes$info$modelCheckItems$maxConditionNumber,
                  printMsg = printMsg
  )

  return(M)
}


#' Step-wise Search for Best VARMA Models
#'
#' For a large model set, use this function to find the best Vector Autoregressive Moving Average models.
#' It selects a subset of variables from smaller models and moves to the bigger ones.
#'
#' @param y A matrix of endogenous data with variables in the columns.
#' @param ySizeSteps A list of model dimensions to be estimated in each step.
#' Its size determines the number of steps.
#' @param countSteps An integer vector to determine the number of variables to be used in each step.
#' \code{NA} means all variables. Variables are selected based on best estimations.
#' All variables in the best models (all measures and targets) are selected until the corresponding suggested number is reached.
#' Select an appropriate value for \code{bestK} in the options.
#' @param savePre A directory for saving and loading the progress.
#' Each step's result is saved in a file (name=\code{paste0(savePre,i)} where \code{i} is the index of the step.
#' @param ... other arguments to pass to [search.varma] function such as the \code{x} argument.
#' Note that \code{ySizes} is ineffective here.
#'
#' @return Similar to [search.varma] function.
#' @export
#'
#' @examples
#' See the example in the 'search.varma' function.
#'
#' @seealso [search.varma], [estim.varma]
search.varma.stepwise <- function(y, ySizeSteps = list(c(1, 2), c(3, 4), c(5), c(6:10)),
                          countSteps = c(NA, 40, 30, 20),
                          savePre = NULL, ...) {
  Search_s("varma", y, ySizeSteps, countSteps, savePre, ...)
}



#' Generate data for a multivariate VARMA model
#'
#' This function generates a multivariate time series using a VARMA process.
#'
#' @param sigma A positive definite matrix representing the covariance matrix of the white noise series or an integer representing the dimension of a random covariance matrix to generate.
#' @param arList A list of matrices representing the AR coefficients of the VARMA model or an integer representing the number of random AR coefficients to generate.
#' @param maList A list of matrices representing the MA coefficients of the VARMA model or an integer representing the number of random MA coefficients to generate. For identification purposes, it generates diagonal matrices.
#' @param exoCoef A matrix representing the coefficients of the exogenous variables or an integer representing the number of random exogenous coefficients to generate.
#' @param nObs An integer representing the number of observations to generate.
#' @param nBurn An integer representing the number of burn-in observations to remove from the generated time series.
#' @param intercept A numeric vector representing the intercept of the VARMA model or a logical value indicating whether to generate a random intercept.
#'
#' @return A list containing the generated data, exogenous data, errors, and all input arguments.
#'
#' @export
#' @importFrom stats rnorm 
#'
#' @examples
#' ar1 <- matrix(c(0.7,0.2,-0.4,0.3),2,2)
#' ar2 <- matrix(c(-0.4,0.1,0.2,-0.3),2,2)
#' ma1 <- matrix(c(0.5,-0.1,0.3,0.4),2,2)
#' Sigma <- matrix(c(1,0.3,0.3,1),2,2)
#' B <- matrix(c(0.5,-0.3),2)
#'
#' result <- sim.varma(Sigma, list(ar1, ar2), list(ma1), exoCoef = B ,
#'                     nObs =100, nBurn =10 ,intercept = c(1,-1))
#'
#' # Plot the y series
#' # matplot(result$y,type = "l")
#'
sim.varma <- function(sigma, arList, maList, exoCoef, nObs, nBurn, intercept) {

  # Set the dimension of the VARMA model
  p <- ifelse(is.numeric(sigma),sigma,ncol(sigma))

  # Generate random AR coefficients if arList is an integer
  if (is.numeric(arList)) {
    arList <- lapply(seq_len(arList), function(x) matrix(rnorm(p^2)*0.1, p, p))
  }

  # Generate random MA coefficients if maList is an integer
  if (is.numeric(maList)) {
    maList <- lapply(seq_len(maList), function(x) diag(rnorm(1)*0.1, p, p))
    # in order to get identifiable model, I make MAs diagonal
    # see https://doi.org/10.1080/07350015.2021.1904960
  }

  # Generate a random covariance matrix if sigma is an integer
  if (is.numeric(sigma)) {
    A <- matrix(rnorm(p^2), p, p)
    sigma <- crossprod(A)
  }

  # Generate random exogenous coefficients if exoCoef is an integer
  if (is.numeric(exoCoef)) {
    if (exoCoef == 0)
      exoCoef = NULL
    else
      exoCoef <- matrix(rnorm(p * exoCoef), p, exoCoef)
  }

  # Generate a random intercept if intercept is TRUE
  if (is.logical(intercept) && intercept) {
    intercept <- rnorm(p)
  }

  # Check the dimensions of the input arguments
  if (length(arList) > 0 && any(sapply(arList, nrow) != p | sapply(arList, ncol) != p)) {
    stop("Invalid dimensions of AR coefficients")
  }
  if (length(maList) > 0 && any(sapply(maList, nrow) != p | sapply(maList, ncol) != p)) {
    stop("Invalid dimensions of MA coefficients")
  }
  if (!is.null(exoCoef) && (nrow(exoCoef) != p || ncol(exoCoef) < 1)) {
    stop("Invalid dimensions of exogenous coefficients")
  }
  if (!is.null(intercept) && length(intercept) != p) {
    stop("Invalid dimensions of intercept")
  }

  # Compute the Cholesky decomposition of the covariance matrix
  L <- chol(sigma)

  # Generate the white noise series
  e <- matrix(rnorm((nObs + nBurn) * p), nObs + nBurn, p) %*% t(L)

  # Generate the exogenous data
  q <- ifelse(is.null(exoCoef),0,ncol(exoCoef))
  x <- matrix(rnorm((nObs + nBurn) * q), nObs + nBurn, q)

  # Initialize the y series
  y <- matrix(0, nObs + nBurn, p)

  # Set the maximum lag order
  p_ar <- length(arList)
  p_ma <- length(maList)
  p_max <- max(p_ar, p_ma)

  # Generate the y series using the VARMA process
  for (t in (p_max+1):(nObs+nBurn)) {
    y[t,] <- e[t,] + ifelse(is.null(exoCoef),0,x[t,] %*% t(exoCoef))
    for (j in seq_along(maList)) {
      y[t,] <- y[t,] + maList[[j]] %*% e[t-j,]
    }
    for (j in seq_along(arList)) {
      y[t,] <- y[t,] + arList[[j]] %*% y[t-j,]
    }
    if (!is.null(intercept)) {
      y[t,] <- y[t,] + intercept
    }
  }

  # Remove the burn-in observations
  y <- y[(nBurn+1):(nObs+nBurn),]
  colnames(y) <- paste0("Y",c(1:ncol(y)))

  if (length(x) == 0)
    x <- NULL
  else
  {
    x <- x[(nBurn+1):(nObs+nBurn),]
    colnames(x) <- paste0("X",c(1:ncol(x)))
  }

  e <- e[(nBurn+1):(nObs+nBurn),]
  colnames(e) <- paste0("E",c(1:ncol(e)))

  return(list(y = y,
              x = x,
              error = e,
              sigma = sigma,
              arList = arList,
              maList = maList,
              exoCoef = exoCoef,
              intercept = intercept,
              nObs = nObs,
              nBurn = nBurn))
}
