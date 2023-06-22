
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
#' \item{counts}{Information about the expected number of models, number of estimated models, failed estimations, and some details about the failures.}
#' \item{...}{Results reported separately for each measure, then for each target variable, then for each requested type of output. This part of the output is highly nested and items are reported based on the arguments of the search.}
#' \item{info}{General information about the search process, some arguments, elapsed time, etc.}
#'
#' Note that the output does not contain any estimation results,
#' but minimum required data to estimate the models (Use \code{summary()} function to get the estimation).
#'
#' @export
#' @example man-roxygen/ex-search.varma.R
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
#' \item{counts}{Information about different aspects of the estimation such as the number of observation, number of exogenous variables, etc.}
#' \item{estimations}{Estimated coefficients, standard errors, z-statistics, p-values, etc.}
#' \item{measures}{Value of different goodness of fit and out-of-sample performance measures. }
#' \item{prediction}{Information on the predicted values.}
#' \item{simulation}{Information on the simulations. }
#' \item{info}{Some other general information.}
#'
#' @details
#' The main purpose of exporting this method is to show the inner calculations of the search process in [search.varma] function. See the details of this function for more information.
#'
#' @export
#' @example man-roxygen/ex-estim.varma.R
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
#' # See the example in the 'search.varma' function.
#'
#' @seealso [search.varma], [estim.varma]
search.varma.stepwise <- function(y, ySizeSteps = list(c(1, 2), c(3, 4), c(5), c(6:10)),
                                  countSteps = c(NA, 40, 30, 20),
                                  savePre = NULL, ...) {
  Search_s("varma", y, ySizeSteps, countSteps, savePre, ...)
}



#' Generate Random Sample from a VARMA Model
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
#' @return A list with the following items:
#'   \item{y}{The simulated endogenous data.}
#'   \item{x}{The simulated exogenous data.}
#'   \item{e}{The simulated white noise series.}
#'   \item{sigma}{The covariance matrix of the white noise series.}
#'   \item{arList}{The list of autoregressive coefficients.}
#'   \item{maList}{The list of moving average coefficients.}
#'   \item{exoCoef}{The matrix of exogenous coefficients.}
#'   \item{intercept}{The intercept vector.}
#'   \item{nObs}{The number of observations generated.}
#'   \item{nBurn}{The number of burn-in observations removed.}
#'
#' A list containing the generated data, exogenous data, errors, and all input arguments.
#'
#' @export
#' @importFrom stats rnorm
#' @example man-roxygen/ex-sim.varma.R
#'
sim.varma <- function(sigma = 2L, arList = 1L, maList = 0L, exoCoef = 0L, nObs = 100, nBurn = 10, intercept = TRUE) {

  nObs = as.integer(nObs)
  if (nObs <= 0)
    stop("nObs must be a positive integer.")

  nBurn = as.integer(nBurn)
  if (nBurn < 0)
    stop("nBurn cannot be a negative integer.")

  # Set the dimension of the VARMA model
  p <- ifelse(is.integer(sigma) && length(sigma) == 1,sigma,ncol(sigma))

  # Generate random AR coefficients if arList is an integer
  if (is.integer(arList) && length(arList) == 1) {
    arList <- lapply(seq_len(arList), function(x) matrix(rnorm(p^2)*0.1, p, p))
  }

  # Generate random MA coefficients if maList is an integer
  if (is.integer(maList) && length(maList) == 1) {
    maList <- lapply(seq_len(maList), function(x) diag(rnorm(1)*0.1, p, p))
    # in order to get identifiable model, I make MAs diagonal
    # see https://doi.org/10.1080/07350015.2021.1904960
  }

  # Generate a random covariance matrix if sigma is an integer
  if (is.integer(sigma) && length(sigma) == 1) {
    A <- matrix(rnorm(p^2), p, p)
    sigma <- crossprod(A)
  }

  # Generate random exogenous coefficients if exoCoef is an integer
  if (is.integer(exoCoef) && length(exoCoef) == 1) {
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
  y <- y[(nBurn+1):(nObs+nBurn),, drop = FALSE]
  colnames(y) <- paste0("Y",c(1:ncol(y)))

  if (length(x) == 0)
    x <- NULL
  else
  {
    x <- x[(nBurn+1):(nObs+nBurn),, drop = FALSE]
    colnames(x) <- paste0("X",c(1:ncol(x)))
  }

  e <- e[(nBurn+1):(nObs+nBurn),, drop = FALSE]
  colnames(e) <- paste0("E",c(1:ncol(e)))



  # Set row and column names for arList items
  if (!is.null(arList)) {
    for (i in seq_along(arList)) {
      colnames(arList[[i]]) <- paste0("Y", seq_len(ncol(arList[[i]])), "(-", i, ")")
      rownames(arList[[i]]) <- paste0("Eq", seq_len(nrow(arList[[i]])))
    }
  }

  if (!is.null(maList)) {
    # Set row and column names for maList items
    for (i in seq_along(maList)) {
      colnames(maList[[i]]) <- paste0("E", seq_len(ncol(maList[[i]])), "(-", i, ")")
      rownames(maList[[i]]) <- paste0("Eq", seq_len(nrow(maList[[i]])))
    }
  }

  # Set row and column names for sigma
  if (!is.null(sigma)) {
    colnames(sigma) <- paste0("Eq", seq_len(ncol(sigma)))
    rownames(sigma) <- paste0("Eq", seq_len(nrow(sigma)))
  }

  # Set row and column names for exoCoef
  if (!is.null(exoCoef)) {
    colnames(exoCoef) <- paste0("X", seq_len(ncol(exoCoef)))
    rownames(exoCoef) <- paste0("Eq", seq_len(nrow(exoCoef)))
  }

  # Set name for intercept
  if (!is.null(intercept)) {
    names(intercept) <- "Intercept"
  }


  return(list(y = y,
              x = x,
              e = e,
              sigma = sigma,
              arList = arList,
              maList = maList,
              exoCoef = exoCoef,
              intercept = intercept,
              nObs = nObs,
              nBurn = nBurn))
}



varma.to.latex.mat <- function(arList, int, exoCoef, maList, numFormat = "%.2f") {

  numEq <- nrow(arList[[1]])
  numAR <- length(arList)
  numMA <- length(maList)
  numExo <- ncol(exoCoef)

  # Initialize the LaTeX string
  latex_str <- ""

  # Add the dependent variable vector
  latex_str <- paste0(latex_str, "\\begin{bmatrix}", paste(paste0("Y_", seq_len(numEq)), collapse = "\\\\"), "\\end{bmatrix}")

  # Add the equal sign
  latex_str <- paste0(latex_str, " = ")

  # Add the intercept vector
  if (!is.null(int)) {
    latex_str <- paste0(latex_str, "\\begin{bmatrix}", paste(sprintf(numFormat, int), collapse = "\\\\"), "\\end{bmatrix}")

    # Add a plus sign if there are more terms
    if (numAR > 0 || numExo > 0 || numMA > 0) {
      latex_str <- paste0(latex_str, " + ")
    }
  }

  # Add the AR matrices
  for (lag in seq_len(numAR)) {
    ar_mat <- matrix(sprintf(numFormat, arList[[lag]]), nrow = nrow(arList[[lag]]), ncol = ncol(arList[[lag]]))
    latex_str <- paste0(
      latex_str,
      "\\begin{bmatrix}",
      paste(apply(ar_mat, 1, paste, collapse = " & "), collapse = "\\\\"),
      "\\end{bmatrix}"
    )

    # Add the lagged dependent variable vector
    latex_str <- paste0(
      latex_str,
      "\\begin{bmatrix}",
      paste(paste0("Y_", seq_len(numEq), "(-", lag, ")"), collapse = "\\\\"),
      "\\end{bmatrix}"
    )

    # Add a plus sign if there are more terms
    if (lag < numAR || numExo > 0 || numMA > 0) {
      latex_str <- paste0(latex_str, " + ")
    }
  }

  # Add the exogenous matrix
  if (numExo > 0) {
    exo_mat <- matrix(sprintf(numFormat, exoCoef), nrow = nrow(exoCoef), ncol = ncol(exoCoef))
    latex_str <- paste0(
      latex_str,
      "\\begin{bmatrix}",
      paste(apply(exo_mat, 1, paste, collapse = " & "), collapse = "\\\\"),
      "\\end{bmatrix}"
    )

    # Add the exogenous variable vector
    latex_str <- paste0(
      latex_str,
      "\\begin{bmatrix}",
      paste(paste0("X_", seq_len(numExo)), collapse = "\\\\"),
      "\\end{bmatrix}"
    )

    # Add a plus sign if there are more terms
    if (numMA > 0) {
      latex_str <- paste0(latex_str, " + ")
    }
  }

  # Add the MA matrices
  for (lag in seq_len(numMA)) {
    ma_mat <- matrix(sprintf(numFormat, maList[[lag]]), nrow = nrow(maList[[lag]]), ncol = ncol(maList[[lag]]))
    latex_str <- paste0(
      latex_str,
      "\\begin{bmatrix}",
      paste(apply(ma_mat, 1, paste, collapse = " & "), collapse = "\\\\"),
      "\\end{bmatrix}"
    )

    # Add the lagged error vector
    latex_str <- paste0(
      latex_str,
      "\\begin{bmatrix}",
      paste(paste0("E_", seq_len(numEq), "(-", lag, ")"), collapse = "\\\\"),
      "\\end{bmatrix}"
    )

    # Add a plus sign if there are more terms
    if (lag < numMA) {
      latex_str <- paste0(latex_str, " + ")
    }
  }

  # Add the error vector
  latex_str <- paste0(latex_str, " + \\begin{bmatrix}", paste(paste0("E_", seq_len(numEq)), collapse = "\\\\"), "\\end{bmatrix}")

  return(latex_str)
}

varma.to.latex.eqs <- function(arList, int, exoCoef, maList, numFormat = "%.2f") {

  numEq <- nrow(arList[[1]])
  numAR <- length(arList)
  numMA <- length(maList)
  numExo <- ncol(exoCoef)

  # Initialize the LaTeX string
  latex_str <- ""

  # Add the equations
  for (eq in seq_len(numEq)) {
    # Add the dependent variable
    latex_str <- paste0(latex_str, "Y_", eq)

    # Add the intercept
    if (!is.null(int)) {
      latex_str <- paste0(latex_str, " = ", round(int[eq], digits = 2))
    }

    # Add the AR terms
    for (lag in seq_len(numAR)) {
      coefs <- round(arList[[lag]][eq, ], digits = 2)
      signs <- ifelse(coefs >= 0, " + ", " - ")
      terms <- paste0(signs, abs(coefs), " Y_", seq_len(numEq), "(-", lag, ")")
      latex_str <- paste0(latex_str, paste(terms, collapse = ""))
    }

    # Add the exogenous terms
    for (exo in seq_len(numExo)) {
      coef <- round(exoCoef[eq, exo], digits = 2)
      sign <- ifelse(coef >= 0, " + ", " - ")
      term <- paste0(sign, abs(coef), " X_", exo)
      latex_str <- paste0(latex_str, term)
    }

    # Add the MA terms
    for (lag in seq_len(numMA)) {
      coefs <- round(maList[[lag]][eq, ], digits = 2)
      signs <- ifelse(coefs >= 0, " + ", " - ")
      terms <- paste0(signs, abs(coefs), " E_", seq_len(numEq), "(-", lag, ")")
      latex_str <- paste0(latex_str, paste(terms, collapse = ""))
    }

    # Add the error term
    latex_str <- paste0(latex_str, " + E_", eq)

    # Add a line break if this is not the last equation
    if (eq < numEq) {
      latex_str <- paste0(latex_str, "\\\\")
    }
    latex_str <- paste0("\\begin{aligned} ", latex_str, " \\end{aligned}")
  }

  return(latex_str)
}


#' Split VARMA parameter into its Components
#'
#' Use this function to extract AR, MA, intercept, and exogenous coefficients from the VARMA estimation.
#'
#' @param coef A matrix of coefficients with dimensions \code{numEq} x \code{numAR * numEq + numMA * numEq + numExo + ifelse(intercept, 1, 0)}.
#' @param numAR A non-negative integer scalar specifying the number of AR lags.
#' @param numMA A non-negative integer scalar specifying the number of MA lags.
#' @param numExo A non-negative integer scalar specifying the number of exogenous variables.
#' @param intercept A logical scalar indicating whether an intercept is included in the model.
#' @param numFormat A character string that determines how to format the numbers, to be used as the argument of the \code{sprintf} function.
#' If \code{NULL}, conversion to latex or html representations are disabled.
#'
#' @return A list with the following items:
#' \itemize{
#'   \item \code{arList}: A list of length \code{numAR} containing the AR coefficients for each lag.
#'   \item \code{intercept}: A numeric vector of length \code{numEq} containing the intercept, or \code{NULL} if \code{intercept = FALSE}.
#'   \item \code{exoCoef}: A matrix of dimensions \code{numEq} x \code{numExo} containing the exogenous coefficients, or \code{NULL} if \code{numExo = 0}.
#'   \item \code{maList}: A list of length \code{numMA} containing the MA coefficients for each lag.
#'   \item \code{eqsLatex}: A character string containing the Latex representation of the equations of the system.
#'   \item \code{eqsLatexSys}: A character string containing the Latex representation of the system in matrix form.
#' }
#' @export
#' @examples
#' # see 'search.varma' or 'estim.varma' functions.
#'
#' @seealso [estim.varma]
get.varma.params <- function(coef, numAR = 1, numMA = 0,
                             numExo = 0, intercept = TRUE,
                             numFormat = "%.2f") {
  coef <- as.matrix(coef)
  numEq <- nrow(coef)
  numAR <- as.integer(numAR)
  numExo <- as.integer(numExo)
  intercept <- as.logical(intercept)

  if (numEq <= 0)
    stop("'numEq' must be a positive integer scalar")
  if (numAR < 0)
    stop("'numAR' must be a non-negative integer scalar")
  if (numMA < 0)
    stop("'numMA' must be a non-negative integer scalar")
  if (numExo < 0)
    stop("'numExo' must be a non-negative integer scalar")

  expected_ncol <- numAR * numEq + numMA * numEq + numExo + ifelse(intercept, 1, 0)
  if (ncol(coef) != expected_ncol) {
    stop("The number of columns of 'coef' is not consistent with the other arguments")
  }

  # Set row names
  if (is.null(rownames(coef))) {
    rownames(coef) <- paste0("Eq", seq_len(numEq))
  }

  # Set column names
  if (is.null(colnames(coef))) {
    colnames(coef) <- c(
      paste0("Y", rep(seq_len(numEq), times = numAR), "(-", rep(seq_len(numAR), each = numEq), ")"),
      ifelse(intercept, "Intercept", NULL),
      paste0("X", seq_len(numExo)),
      paste0("E", rep(seq_len(numEq), times = numMA), "(-", rep(seq_len(numMA), each = numEq), ")")
    )
  }

  start <- 1
  if (numAR > 0) {
    ar <- lapply(1:numAR, function(i) coef[, (start + (i - 1) * numEq):(start + i * numEq - 1), drop = FALSE])
    start <- start + numAR * numEq
  } else {
    ar <- NULL
  }
  if (intercept) {
    int <- coef[, start, drop = FALSE]
    start <- start + 1
  } else {
    int <- NULL
  }
  if (numExo > 0) {
    exog <- coef[, start:(start + numExo - 1), drop = FALSE]
    start <- start + numExo
  } else {
    exog <- NULL
  }
  if (numMA > 0) {
    ma <- lapply(1:numMA, function(i) coef[, (start + (i - 1) * numEq):(start + i * numEq - 1), drop = FALSE])
  } else {
    ma <- NULL
  }


  return(list(arList = ar, integer = int,
              exoCoef = exog, maList = ma,
              eqsLatex = ifelse(is.null(numFormat), NULL, varma.to.latex.eqs(ar, int, exog, ma, as.character(numFormat))),
              eqsLatexSys = ifelse(is.null(numFormat), NULL, varma.to.latex.mat(ar, int, exog, ma, as.character(numFormat)))))
}
