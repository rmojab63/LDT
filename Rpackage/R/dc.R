

#' Search for Best Discrete-Choice Models
#'
#' Use this function to create a discrete-choice model set and search for the best models (and other information) based on in-sample and out-of-sample evaluation measures.
#'
#' @param y A matrix of endogenous data with variable in the column.
#' @param x A matrix of exogenous data with variables in the columns.
#' @param w Weights of the observations in \code{y}.
#' \code{NULL} means equal weights for all observations.
#' @param xSizes An integer vector specifying the number of exogenous variables in the regressions.
#' E.g., \code{c(1,2)} means the model set contains all regressions with 1 and 2 exogenous variables.
#' If \code{NULL}, \code{c(1)} is used.
#' @param xPartitions A list of integer vectors that partition the indexes of the exogenous variables.
#' No regression is estimated with two variables in the same partition.
#' If \code{NULL}, each variable is placed in its own group, and the size of the model set is maximized.
#' @param costMatrices A list of numeric matrices where each one determines how to score the calculated probabilities.
#' Given the number of choices \code{n}, a frequency cost matrix is an \code{m x n+1} matrix.
#' The first column determines the thresholds.
#' Elements in the \code{j}-th column determine the costs corresponding to the \code{j-1}-th choice in \code{y}.
#' It can be \code{NULL} if it is not selected in \code{measureOptions}.
#' @param searchLogit If \code{TRUE}, logit regressions are added to the model set.
#' @param searchProbit If \code{TRUE}, probit regressions are added to the model set.
#' @param optimOptions A list for Newton optimization options.
#' Use [get.options.newton] function to get the options.
#' @param aucOptions A list for AUC calculation options.
#' Use [get.options.roc] function to get the options.
#' @param measureOptions A list of options for measuring performance.
#' Use [get.options.measure] function to get them.
#' @param modelCheckItems A list of options for excluding a subset of the model set.
#' See and use [get.items.modelcheck] function to get them.
#' @param searchItems A list of options for specifying the purpose of the search.
#' See and use [get.items.search] function to get them.
#' @param searchOptions A list of extra options for performing the search.
#' See and use [get.options.search] function to get them.
#'
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
#' @importFrom stats glm
#'
#' @example man-roxygen/ex-search.dc.R
#'
#' @seealso [estim.dc], [search.dc.stepwise]
search.dc <- function(y, x, w = NULL, xSizes = NULL,
                      xPartitions = NULL, costMatrices = NULL,
                      searchLogit = TRUE, searchProbit = FALSE,
                      optimOptions = get.options.newton(), aucOptions = get.options.roc(),
                      measureOptions = get.options.measure(),
                      modelCheckItems = get.items.modelcheck(),
                      searchItems = get.items.search(),
                      searchOptions = get.options.search()){

  y = as.matrix(y)
  x = as.matrix(x)
  w = if (is.null(w)) NULL else as.matrix(w)
  xSizes = if (is.null(xSizes)) NULL else as.integer(xSizes)
  searchLogit = as.logical(searchLogit)
  searchProbit = as.logical(searchProbit)


  if (is.null(xPartitions) == FALSE){
    xPartitions = as.list(xPartitions)
    for (i in c(1:length(xPartitions)))
      xPartitions[[i]] = as.integer(xPartitions[[i]])
  }

  if (is.null(costMatrices) == FALSE){
    costMatrices = as.list(costMatrices)
    for (i in c(1:length(costMatrices)))
      costMatrices[[i]] = as.matrix(costMatrices[[i]])
  }

  if (is.null(optimOptions))
    optimOptions = get.options.newton()
  else
    optimOptions = CheckNewtonOptions(optimOptions)

  if (is.null(aucOptions))
    aucOptions = get.options.roc()
  else
    aucOptions = CheckRocOptions(aucOptions)

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

  res <- .SearchDc(y, x, w, xSizes, xPartitions, costMatrices,
                   searchLogit, searchProbit,
                   optimOptions, aucOptions, measureOptions ,
                   modelCheckItems, searchItems,
                   searchOptions)
  res
}


#' Estimate a Discrete Choice Model
#'
#' Use this function to estimate a discrete choice model.
#'
#' @param y A matrix of endogenous data with variable in the column.
#' @param x A matrix of exogenous data with variables in the columns.
#' @param w Weights of the observations in \code{y}.
#' \code{NULL} means equal weights for all observations.
#' @param distType A character string that shows the distribution assumption. It can be \code{logit} or \code{probit}.
#' @param newX A numeric matrix where for each row in it, probabilities are projected and reported. It can be \code{NULL}.
#' @param pcaOptionsX A list of options to use principal components of the \code{x}, instead of the actual values. Set \code{NULL} to disable. Use [get.options.pca()] for initialization.
#' @param costMatrices A list of numeric matrices where each one determines how to score the calculated probabilities. See and use [search.dc] for more information and initialization.
#' @param aucOptions A list of options for AUC calculation. See and use \code{[get.options.roc()]} for more information and initialization.
#' @param simFixSize An integer that determines the number of pseudo out-of-sample simulations. Use zero to disable the simulation.
#' @param simTrainFixSize An integer representing the number of data points in the training sample in the pseudo out-of-sample simulation. If zero, \code{trainRatio} will be used.
#' @param simTrainRatio A number representing the size of the training sample relative to the available size, in the pseudo out-of-sample simulation. It is effective if \code{trainFixSize} is zero.
#' @param simSeed A seed for the random number generator. Use zero for a random value.
#' @param weightedEval If \code{TRUE}, weights will be used in evaluations.
#' @param printMsg Set to \code{TRUE} to enable printing some details.
#'
#' @return A nested list with the following items:
#' \item{counts}{Information about different aspects of the estimation such as the number of observation, number of exogenous variables, etc.}
#' \item{estimations}{Estimated coefficients, standard errors, z-statistics, p-values, etc.}
#' \item{measures}{Value of different goodness of fit and out-of-sample performance measures.}
#' \item{projections}{Information on the projected values, if \code{newX} is provided.}
#' \item{info}{Some other general information.}
#'
#' @details
#' The main purpose of exporting this method is to show the inner calculations of the search process in [search.dc] function. See the details of this function for more information.
#'
#' @export
#' @example man-roxygen/ex-estim.dc.R
#' @seealso [search.dc], [search.dc.stepwise]
estim.dc <- function(y, x, w = NULL,
                     distType = c("logit", "probit"), newX = NULL,
                     pcaOptionsX = NULL, costMatrices = NULL,
                     aucOptions = get.options.roc(), simFixSize = 200,
                     simTrainFixSize = 0, simTrainRatio = 0.5, simSeed = 0, weightedEval = FALSE, printMsg = FALSE)
{

  y = as.matrix(y)
  x = as.matrix(x)
  w = if (is.null(w)) NULL else as.matrix(w)
  newX = if (is.null(newX)) NULL else as.matrix(newX)
  simFixSize = as.integer(simFixSize)
  simTrainRatio = as.numeric(simTrainRatio)
  simTrainFixSize = as.integer(simTrainFixSize)
  simSeed = as.integer(simSeed)
  weightedEval = as.logical(weightedEval)
  printMsg = as.logical(printMsg)

  distType <- match.arg(as.character(distType), c("logit", "probit"))

  if (is.null(pcaOptionsX) == FALSE)
    pcaOptionsX = CheckPcaOptions(as.list(pcaOptionsX))

  if (is.null(costMatrices) == FALSE){
    costMatrices = as.list(costMatrices)
    for (i in c(1:length(costMatrices)))
      costMatrices[[i]] = as.matrix(costMatrices[[i]])
  }

  if (is.null(aucOptions))
    aucOptions = get.options.roc()
  else
    aucOptions = CheckRocOptions(aucOptions)

  res <- .EstimDc(y, x, w, distType, newX,
                  pcaOptionsX, costMatrices,
                  aucOptions, simFixSize, simTrainRatio,
                  simTrainFixSize, simSeed,
                  weightedEval, printMsg)
  res
}

# get estimation from search result
GetEstim_dc <- function(searchRes, endoIndices, exoIndices, y, x, printMsg, w, distType, ...) {
  M <- estim.dc(y,
                x = if (is.null(exoIndices) || is.null(x)) NULL else x[, exoIndices, drop = FALSE],
                w = w,
                distType = distType,
                newX = if (is.null(exoIndices) || is.null(searchRes$info$newX)) {
                  NULL
                } else {
                  as.matrix(searchRes$info$newX[, exoIndices])
                },
                pcaOptionsX = NULL,
                costMatrices = searchRes$info$costMatrices,
                simFixSize = searchRes$info$measureOptions$simFixSize,
                simTrainRatio = searchRes$info$measureOptions$trainRatio,
                simTrainFixSize = searchRes$info$measureOptions$trainFixSize,
                simSeed = abs(searchRes$info$measureOptions$seed),
                printMsg = printMsg
  )

  return(M)
}

#' Step-wise Search for Best Discrete-Choice Models
#'
#' For a large model set, use this function to find the best discrete choice models.
#' It selects a subset of variables from smaller models and moves to the bigger ones.
#'
#' @param y A matrix of endogenous data with variable in the column.
#' @param x A matrix of exogenous data with variables in the columns.
#' @param xSizeSteps A list of model dimensions to be estimated in each step. Its size determines the number of steps.
#' @param countSteps A integer vector to determine the number of variables to be used in each step.
#' \code{NA} means all variables. Variables are selected based on best estimations.
#' All variables in the best models (all measures and targets) are selected until the corresponding suggested number is reached.
#' Select an appropriate value for \code{bestK} in the options.
#' @param savePre A directory for saving and loading the progress.
#' Each step's result is saved in a file (name=\code{paste0(savePre,i)} where \code{i} is the index of the step.
#' @param ... other arguments to pass to [search.dc()] function such as the \code{w} argument.
#' Note that \code{xSizes} is ineffective here.
#'
#' @return Similar to [search.dc] function.
#' @export
#'
#' @examples
#' # See the example in the 'search.dc' function.
#'
#' @seealso [search.dc], [estim.dc]
search.dc.stepwise <- function(y, x, xSizeSteps = list(c(1, 2), c(3, 4), c(5), c(6:10)),
                               countSteps = c(NA, 40, 30, 20),
                               savePre = NULL, ...) {
  Search_s("dc", x, xSizeSteps, countSteps, savePre, y = y, ...)
}


dc.to.latex.eqs <- function(coef, probit, numFormat = "%.2f") {

  if (is.null(xNames))
    xNames <- paste0("X_",c(1:length(coef)))

  terms <- character(length(coef))
  terms[1] <- sprintf0(numFormat, coef[1])
  for (i in seq_along(coef)[-1]) {
    # avoid +- signs
    if (coef[i] >= 0) {
      sign <- " + "
    } else {
      sign <- " - "
    }
    terms[i] <- paste0(sign, sprintf(numFormat, abs(coef[i])), " ", xNames[i])
  }

  formula_str <- paste(terms, collapse = "")

  if (probit) {
    formula_str <- paste0("P(Y = 1 | ", paste(xNames[-1], collapse = ", "), ") = \\Phi(", formula_str, ")")
  } else {
    formula_str <- paste0("P(Y = 1 | ", paste(xNames[-1], collapse = ", "), ") = \\frac{1}{1 + e^{-(", formula_str, ")}}")
  }

  return(formula_str)
}

#' Generate Random Sample from a DC Model
#'
#' This function generates a random sample from an discrete choice regression model.
#'
#' @param coef Either a single integer specifying the number of variables in the model,
#'   or a numeric vector of coefficients for the regression.
#' @param nObs The number of observations to generate.
#' @param probit Logical value indicating whether to generate data from a probit model
#'   (if \code{TRUE}) or a logit model (if \code{FALSE}).
#' @param pPos The percentage of positive observations (\code{y=1}) in the dependent variable y.
#'   Must be between 0 and 1.
#' @param sampleFactor The factor used to control the size of the initial sample.
#'   A larger value generates a larger initial sample, which can increase the accuracy
#'   of the generated sample but also takes more time and memory.
#'
#' @return A list with the following items:
#'   \item{y}{The dependent variable.}
#'   \item{x}{The independent variables.}
#'   \item{coef}{The coefficients of the regression.}
#'   \item{probit}{Logical value indicating whether data was generated from a probit model.}
#'   \item{pPos}{The percentage of negative observations in y.}
#'
#' @export
#' @importFrom stats pnorm rnorm rbinom
#' @example man-roxygen/ex-sim.dc.R
#'
#' @seealso [estim.dc], [search.dc]
sim.dc <- function(coef = 2L, nObs = 100, probit = FALSE, pPos = 0.5, sampleFactor = 2) {

  if (nObs <= 0)
    stop("nObs must be a positive integer")

  probit <- as.logical(probit)

  pPos <- as.numeric(pPos)
  if (pPos < 0 || pPos > 1)
    stop("pPos must be between 0 and 1")

  sampleFactor = as.numeric(sampleFactor)
  if (sampleFactor < 1)
    stop("sampleFactor must be larger than 1")

  # Set the number of variables
  if (length(coef) == 1 && is.integer(coef)) {
    nVar <- coef
    coef <- rnorm(nVar)
  } else {
    nVar <- length(coef)
  }

  # Generate the independent variables
  if (nVar == 0) {
    x <- matrix(1, ncol = 1, nrow = nObs * sampleFactor)
    colnames(x) <- c("Intercept")
  } else {
    x <- cbind(1, matrix(rnorm(nObs * sampleFactor * (nVar - 1)), ncol = nVar - 1))
    colnames(x) <- c("Intercept", paste0("X", seq_len(nVar - 1)))
  }

  # Calculate the probability of choosing option 1
  if (probit) {
    p1 <- pnorm(x %*% coef)
  } else {
    p1 <- 1 / (1 + exp(-x %*% coef))
  }

  # Generate the dependent variable
  y <- matrix(as.numeric(rbinom(n = nObs * sampleFactor, size = 1, prob = p1)), ncol = 1)
  colnames(y) <- "Y"

  # Remove observations until we have nPos positive observations
  nPos <- round(nObs * pPos)
  posToRemove <- sum(y) - nPos
  if (posToRemove > 0) {
    posIndexes <- sample(which(y == 1), posToRemove)
    y <- y[-posIndexes, , drop = FALSE]
    x <- x[-posIndexes, , drop = FALSE]
  } else if (posToRemove < 0) {
    stop("Failed to generate enough positive observations. Try increasing the sampleFactor argument.")
  }

  # Remove excess negative observations
  negToRemove <- sum(!y) - (nObs - nPos)
  if (negToRemove > 0) {
    negIndexes <- sample(which(y == 0), negToRemove)
    y <- y[-negIndexes, , drop = FALSE]
    x <- x[-negIndexes, , drop = FALSE]
  } else if (negToRemove < 0) {
    stop("Failed to generate enough negative observations. Try increasing the sampleFactor argument.")
  }

  # Return the results as a list
  list(y = y, x = x, coef = coef, probit = probit, pPos = pPos)
}
