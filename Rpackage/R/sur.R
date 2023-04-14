
#' SUR Search
#'
#' @param y (numeric matrix) endogenous data with variables in the columns.
#' @param x (numeric matrix) exogenous data with variables in the columns.
#' @param numTargets (int) determines the number of variable in the first columns of \code{y} for which the information is saved. It must be positive and cannot be larger than the number of endogenous variables.
#' @param xSizes (nullable integer vector) Number of exogenous variables in the regressions. E.g., c(1,2) means the model set contains all the regressions with 1 and 2 exogenous variables. If null, c(1) is used.
#' @param xPartitions (nullable list of integer vector) a partition over the indexes of the exogenous variables. No regression is estimated with two variables in the same group. If \code{NULL}, each variable is placed in its own group and the size of the model set is maximized.
#' @param numFixXPartitions (int) number of partitions at the beginning of \code{xPartitions} to be included in all regressions.
#' @param yGroups (nullable list of integer vector) different combinations of the indexes of the endogenous variables to be used as endogenous variables in the SUR regressions.
#' @param searchSigMaxIter (int) maximum number of iterations in searching for significant coefficients. Use 0 to disable the search.
#' @param searchSigMaxProb (double) maximum value of type I error to be used in searching for significant coefficients. If p-value is less than this, it is interpreted as significant.
#' @param measureOptions (list) see \code{[GetMeasureOptions()]}.
#' @param modelCheckItems (list) see \code{[GetModelCheckItems()]}.
#' @param searchItems (list) see \code{[GetSearchItems()]}.
#' @param searchOptions (list) see \code{[GetSearchOptions()]}.
#'
#' @return An object of class \code{ldtsearch}. It does not contain any estimation results,
#' but minimum required data to estimate the models (Use \code{[summary()]} for this goal).
#' An object of class \code{ldtsearch} has a nested structure.
#'
#' @export
SurSearch <- function(y, x, numTargets = 1, xSizes = NULL,
                      xPartitions = NULL, numFixXPartitions = 0,
                      yGroups = NULL, searchSigMaxIter = 0,
                      searchSigMaxProb = 0.1, measureOptions = GetMeasureOptions(),
                      modelCheckItems = GetModelCheckItems(), searchItems = GetSearchItems(),
                      searchOptions = GetSearchOptions()){

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
    measureOptions = GetMeasureOptions()
  else
    measureOptions <- CheckMeasureOptions(measureOptions)

  if (is.null(modelCheckItems))
    modelCheckItems = GetModelCheckItems()
  else
    modelCheckItems <- CheckModelCheckItems(modelCheckItems)

  if (is.null(searchItems))
    searchItems = GetSearchItems()
  else
    searchItems <- CheckSearchItems(searchItems)

  if (is.null(searchOptions))
    searchOptions = GetSearchOptions()
  else
    searchOptions <- CheckSearchOptions(searchOptions)

  res <- .SurSearch(y, x, numTargets, xSizes, xPartitions, numFixXPartitions,
                    yGroups, searchSigMaxIter, searchSigMaxProb, measureOptions,
                    modelCheckItems, searchItems, searchOptions)
  res
}


#' Estimates an SUR Model
#'
#' @param y (numeric matrix) Endogenous data with variables in the columns.
#' @param x (numeric matrix) Exogenous data with variables in the columns.
#' @param addIntercept (bool) If TRUE, intercept is added automatically to x.
#' @param searchSigMaxIter (int) Maximum number of iterations in searching for significant coefficients. Use 0 to disable the search.
#' @param searchSigMaxProb (double) Maximum value of type I error to be used in searching for significant coefficients. If p-value is less than this, it is interpreted as significant and removed in the next iteration (if any exists).
#' @param restriction (nullable numeric matrix) A km x q matrix in which m=ncols(y), k=ncols(x) and q is the number of unrestricted coefficients.
#' @param newX (nullable numeric matrix) Data of new exogenous variables to be used in the predictions. Its columns must be the same as \code{x}. If null, projection is disabled.
#' @param pcaOptionsY (nullable list) A list of options in order to use principal components of the \code{y}, instead of the actual values. Set null to disable. Use \code{[GetPcaOptions()]} for initialization.
#' @param pcaOptionsX (nullable list) Similar to \code{pcaOptionsY} but for \code{x}. see \code{pcaOptionsY}.
#' @param simFixSize (int) Number of pseudo out-of-sample simulations. Use zero to disable the simulation. See also \code{GetMeasureOptions()]}.
#' @param simTrainRatio (double) Size of the training sample as a ratio of the number of the observations. It is effective only if \code{simTrainFixSize} is zero.
#' @param simTrainFixSize (int) A fixed size for the training sample. If zero, \code{simTrainRatio} is used.
#' @param simSeed (int) A seed for the pseudo out-of-sample simulation.
#' @param simMaxConditionNumber (double) Maximum value for the condition number in the simulation.
#' @param printMsg (bool) Set TRUE to enable printing details.
#'
#' @return An object of class \code{ldtestimsur}.
#'
#' @export
SurEstim <- function(y, x, addIntercept = TRUE,
                     searchSigMaxIter = 0, searchSigMaxProb = 0.1,
                     restriction = NULL, newX = NULL,
                     pcaOptionsY = NULL, pcaOptionsX = NULL,
                     simFixSize = 0, simTrainRatio = 0.75,
                     simTrainFixSize = 0, simSeed = 0,
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

  res <- .SurEstim(y, x, addIntercept,
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

  M <- SurEstim(
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


#' Step-wise SUR Search
#'
#' A helper class to deal with large model sets.
#' It selects a subset of variables from smaller models and
#' moves to the bigger ones.
#'
#' @param x exogenous data
#' @param xSizes a list of model dimension to be estimated in each step.
#' @param counts a list of suggested number of variables to be used
#' in each step. \code{NA} means all variables. Variables are selected
#' based on best estimations (select an appropriate value for
#' \code{searchItems$bestK}). All variables in the best models (all measures
#' and targets) are selected until corresponding suggested number is reached.
#' @param savePre if not \code{NULL}, it saves and tries to load the progress
#' of search step in a file (name=\code{paste0(savePre,i)} where \code{i}
#' is the index of the step).
#' @param ... other arguments to pass to [SurSearch()] function such
#' as endogenous data. Note that \code{xSizes} is treated differently.
#'
#' @return A combined \code{LdtSearch} object
#' @export
SurSearch_s <- function(x, xSizes = list(c(1, 2), c(3, 4), c(5), c(6:10)),
                        counts = c(NA, 40, 30, 20),
                        savePre = NULL, ...) {
  Search_s("sur", x, xSizes, counts, savePre, ...)
}
