

#' Discrete Choice Search
#'
#' @param y (numeric matrix) endogenous data with variable in the column.
#' @param x (numeric matrix) exogenous data with variables in the columns.
#' @param w (numeric vector) weights of the observations in \code{y}. null means equal weights.
#' @param xSizes (nullable vector) Number of exogenous variables in the regressions. E.g., c(1,2) means the model set contains all the regressions with 1 and 2 exogenous variables. If null, c(1) is used.
#' @param xPartitions (nullable list of vector) a partition over the indexes of the exogenous variables. No regression is estimated with two variables in the same group. If null, each variable is placed in its own group and the size of the model set is maximized.
#' @param costMatrices (list of numeric matrix) each frequency cost matrix determines how to score the calculated probabilities. Given the number of choices 'n', a frequency cost matrix is a 'm x n+1' matrix. The first column determines the thresholds. Cells in the j-th column determines the costs corresponding to the (j-1)-th choice in \code{y}. It can be null if it is not selected in \code{measureOptions}.
#' @param searchLogit (bool) if \code{TRUE}, logit regressions are added to the model set.
#' @param searchProbit (bool) if \code{TRUE}, probit regressions are added to the model set.
#' @param optimOptions (list) Newton optimization options. see \code{[get.newton.options()]}.
#' @param aucOptions (list) AUC calculation options. see \code{[get.roc.options()]}.
#' @param measureOptions (list) see \code{[get.measure.options()]}.
#' @param modelCheckItems (list) see \code{[get.modelcheck.items()]}.
#' @param searchItems (list) see \code{[get.search.items()]}.
#' @param searchOptions (list) see \code{[get.search.options()]}.
#'
#' @return An object of class \code{ldtsearch}. It does not contain any estimation results,
#' but minimum required data to estimate the models (Use \code{[summary()]} for this goal).
#' An object of class \code{ldtsearch} has a nested structure.
#'
#' @export
search.dc <- function(y, x, w = NULL, xSizes = NULL,
                     xPartitions = NULL, costMatrices = NULL,
                     searchLogit = TRUE, searchProbit = FALSE,
                     optimOptions = get.newton.options(), aucOptions = get.roc.options(),
                     measureOptions = get.measure.options(),
                     modelCheckItems = get.modelcheck.items(),
                     searchItems = get.search.items(),
                     searchOptions = get.search.options()){

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
    optimOptions = get.newton.options()
  else
    optimOptions = CheckNewtonOptions(optimOptions)

  if (is.null(aucOptions))
    aucOptions = get.roc.options()
  else
    aucOptions = CheckRocOptions(aucOptions)

  if (is.null(measureOptions))
    measureOptions = get.measure.options()
  else
    measureOptions <- CheckMeasureOptions(measureOptions)

  if (is.null(modelCheckItems))
    modelCheckItems = get.modelcheck.items()
  else
    modelCheckItems <- CheckModelCheckItems(modelCheckItems)

  if (is.null(searchItems))
    searchItems = get.search.items()
  else
    searchItems <- CheckSearchItems(searchItems)

  if (is.null(searchOptions))
    searchOptions = get.search.options()
  else
    searchOptions <- CheckSearchOptions(searchOptions)

  res <- .SearchDc(y, x, w, xSizes, xPartitions, costMatrices,
                   searchLogit, searchProbit,
                   optimOptions, aucOptions, measureOptions ,
                   modelCheckItems, searchItems,
                   searchOptions)
  res
}


#' Estimates an Discrete Choice Model
#'
#' @param y (numeric matrix) Data with dependent variable in the column. Given the number of choices 'n', it must contain 0,1,...,n-1 and 'sum(y==i)>0' for i=0,...,n-1.
#' @param x (numeric matrix) Exogenous data with variables in the columns.
#' @param w (numeric vector) Weights of the observations in \code{y}. Null means equal weights.
#' @param distType (string) Distribution assumption. It can be \code{logit} or \code{probit}.
#' @param newX (numeric matrix) If not null, probabilities are projected for each row of this matrix.
#' @param pcaOptionsX (list) A list of options in order to use principal components of the \code{x}, instead of the actual values. set null to disable. Use [get.pca.options()] for initialization.
#' @param costMatrices (list of matrices) Each cost table determines how you score the calculated probabilities.
#' @param aucOptions (nullable list) AUC calculation options. see \code{[get.roc.options()]}.
#' @param simFixSize (int) Number of pseudo out-of-sample simulations. Use zero to disable the simulation. (see [get.measure.options()]).
#' @param simTrainRatio (double) Size of the training sample as a ratio of the number of the observations. It is effective only if \code{simTrainFixSize} is zero.
#' @param simTrainFixSize (int) A fixed size for the training sample. If zero, \code{simTrainRatio} is used.
#' @param simSeed (int) A seed for the pseudo out-of-sample simulation.
#' @param weightedEval (bool) If TRUE, weights will be used in evaluations.
#' @param printMsg (bool) Set FALSE to disable printing the details.
#'
#' @return An object of class \code{ldtestimdc}.
#'
#' @export
estim.dc <- function(y, x, w = NULL,
                    distType = c("logit", "probit"), newX = NULL,
                    pcaOptionsX = NULL, costMatrices = NULL,
                    aucOptions = get.roc.options(), simFixSize = 200, simTrainRatio = 0.5,
                    simTrainFixSize = 0, simSeed = 0, weightedEval = FALSE, printMsg = FALSE)
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
    aucOptions = get.roc.options()
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

#' Step-wise Discrete Choice Search
#'
#' A helper class to deal with large model sets. It selects a subset of variables
#' from smaller models and moves to the bigger ones.
#'
#' @param x exogenous data
#' @param xSizes a list of model dimension to be estimated in each step.
#' @param counts a list of suggested number of variables to be used in each step.
#' \code{NA} means all variables. Variables are selected based on best estimations
#' (select an appropriate value for \code{searchItems$bestK}). All variables in the
#' best models (all measures and targets) are selected until corresponding
#' suggested number is reached.
#' @param savePre if not \code{NULL}, it saves and tries to load the progress of search
#' step in a file (name=\code{paste0(savePre,i)} where \code{i} is the index of the step).
#' @param ... other arguments to pass to [search.dc()] function such as endogenous data.
#' Note that \code{xSizes} is treated differently.
#'
#' @return A combined \code{LdtSearch} object
#' @export
search.dc_s <- function(x, xSizes = list(c(1, 2), c(3, 4), c(5), c(6:10)),
                       counts = c(NA, 40, 30, 20),
                       savePre = NULL, ...) {
  Search_s("dc", x, xSizes, counts, savePre, ...)
}
