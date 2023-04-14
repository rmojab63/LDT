
#' VARMA Search
#'
#' @param y (numeric vector) Endogenous data with variables in the columns.
#' @param x (nullable numeric matrix) Exogenous data with variables in the columns. It can be null.
#' @param numTargets (int) Number of variables in the first columns of \code{y}, regarded as targets. It must be positive and cannot be larger than the number of endogenous variables.
#' @param ySizes (nullable integer vector) Determines the number of endogenous variables (or equations) in the regressions.
#' @param yPartitions (nullable list of int vector) A partition over the indexes of the endogenous variables. No regression is estimated with two variables in the same group. If \code{NULL}, each variable is placed in its own group.
#' @param xGroups (nullable list of int vector) different combinations of the indexes of the exogenous variables to be used as exogenous variables in the SUR regressions.
#' @param maxParams (integer vector, length=6) Maximum values for the parameters of the VARMA model (p,d,q,P,D,Q). If null, c(1,1,1,0,0,0) is used.
#' @param seasonsCount (integer) number of observations per unit of time
#' @param maxHorizon (integer) maximum value for the prediction horizon if \code{type1} is \code{TRUE} in \code{checkItems}. Also, it is used as the maximum prediction horizon in checking the predictions.
#' @param newX (matrix) New exogenous data for out-of-sample prediction. It must have the same number of columns as \code{x}.
#' @param simUsePreviousEstim (logical) if \code{TRUE}, parameters are initialized in just the first step of the simulation. The initial values of the n-th simulation (with one more observation) is the estimations in the previous step.
#' @param olsStdMultiplier (numeric) a multiplier for the standard deviation of OLS, used for restricting the maximum likelihood estimation.
#' @param lmbfgsOptions (list) Optimization options. see \code{[GetLmbfgsOptions()]}. Use null for default values.
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
VarmaSearch <- function(y, x = NULL, numTargets = 1,
                        ySizes = NULL, yPartitions = NULL,
                        xGroups = NULL, maxParams = c(1,0,0,0,0,0),
                        seasonsCount = 0, maxHorizon = 0,
                        newX = NULL, simUsePreviousEstim = TRUE,
                        olsStdMultiplier = 2.0, lmbfgsOptions = GetLmbfgsOptions(),
                        measureOptions = GetMeasureOptions(),
                        modelCheckItems = GetModelCheckItems(), searchItems = GetSearchItems(),
                        searchOptions = GetSearchOptions()){

  y = as.matrix(y)
  x = if (is.null(x)) NULL else as.matrix(x)
  numTargets = as.integer(numTargets)
  ySizes = if (is.null(ySizes)) NULL else as.integer(ySizes)
  maxParams = as.integer(maxParams)
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

  if (is.null(lmbfgsOptions))
    lmbfgsOptions = GetLmbfgsOptions()
  else
    lmbfgsOptions = CheckLmbfgsOptions(lmbfgsOptions)

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

  res <- .VarmaSearch(y, x, numTargets, ySizes, yPartitions,
                      xGroups, maxParams, seasonsCount, maxHorizon,
                      newX, simUsePreviousEstim,
                      olsStdMultiplier, lmbfgsOptions,
                      measureOptions,modelCheckItems,searchItems,searchOptions)
  res
}

#' Estimates an VARMA Model
#'
#' @param y (matrix) endogenous data with variables in the columns.
#' @param x (matrix) exogenous data with variables in the columns.
#' @param params (integer vector, length=6) parameters of the VARMA model (p,d,q,P,D,Q).
#' @param seasonsCount (integer) number of observations per unit of time
#' @param addIntercept (logical) if \code{TRUE}, intercept is added automatically to x.
#' @param lmbfgsOptions (list) optimization options. See \code{[GetLmbfgsOptions()]}.
#' @param olsStdMultiplier (numeric) a multiplier for the standard deviation of OLS, used for restricting the maximum likelihood estimation.
#' @param pcaOptionsY (list) a list of options in order to use principal components of the \code{y}, instead of the actual values. set \code{NULL} to disable. Use \code{[GetPcaOptions()]} for initialization.
#' @param pcaOptionsX (list) similar to \code{pcaOptionsY} but for \code{x}. see \code{pcaOptionsY}.
#' @param maxHorizon (integer) maximum prediction horizon. Set zero to disable.
#' @param newX (matrix) data of new exogenous variables to be used in the predictions. Its columns must be the same as \code{x}.
#' @param simFixSize (integer) number of pseudo out-of-sample simulations. Use zero to disable the simulation. see also \code{[GetMeasureOptions()]}.
#' @param simHorizons (integer vector) prediction horizons to be used in pseudo out-of-sample simulations. see also \code{[GetMeasureOptions()]}.
#' @param simUsePreviousEstim (logical) if \code{TRUE}, parameters are initialized in just the first step of the simulation. The initial values of the n-th simulation (with one more observation) is the estimations in the previous step.
#' @param simMaxConditionNumber (numeric) maximum value for the condition number in the pseudo out-of-sample simulations.
#' @param printMsg (logical) set \code{FALSE} to disable printing the details.
#'
#' @return An object of class \code{ldtestimvarma}.
#'
#' @export
VarmaEstim <- function(y, x = NULL, params = NULL,
                       seasonsCount = 0, addIntercept = TRUE,
                       lmbfgsOptions = GetLmbfgsOptions(), olsStdMultiplier = 2,
                       pcaOptionsY = NULL, pcaOptionsX = NULL,
                       maxHorizon = 0, newX = NULL, simFixSize = 0,
                       simHorizons = NULL, simUsePreviousEstim = TRUE,
                       simMaxConditionNumber = Inf, printMsg = FALSE){

  y = as.matrix(y)
  x = if (is.null(x)) NULL else as.matrix(x)
  params = as.integer(params)
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

  if (is.null(lmbfgsOptions))
    lmbfgsOptions = GetLmbfgsOptions()
  else
    lmbfgsOptions = CheckLmbfgsOptions(lmbfgsOptions)

  if (is.null(pcaOptionsY) == FALSE)
    pcaOptionsY = CheckPcaOptions(as.list(pcaOptionsY))
  if (is.null(pcaOptionsX) == FALSE)
    pcaOptionsX = CheckPcaOptions(as.list(pcaOptionsX))

  res <- .VarmaEstim(y, x, params, seasonsCount, addIntercept,
                     lmbfgsOptions, olsStdMultiplier,
                     pcaOptionsY, pcaOptionsX,
                     maxHorizon, newX, simFixSize, simHorizons,
                     simUsePreviousEstim, simMaxConditionNumber, printMsg)
  res
}


# get estimation from search result
GetEstim_varma <- function(searchRes, endoIndices,
                           exoIndices, y, x, printMsg,
                           params, newX = NULL, ...) {
  M <- VarmaEstim(y[, endoIndices, drop = FALSE],
                  x = if (is.null(exoIndices) || is.null(x)) {
                    NULL
                  } else {
                    x[, exoIndices, drop = FALSE]
                  },
                  params = params,
                  seasonsCount = searchRes$info$seasonsCount,
                  addIntercept = FALSE,
                  lmbfgsOptions = searchRes$info$lmbfgsOptions,
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


#' Step-wise VARMA Search
#'
#' A helper class to deal with large model sets.
#' It selects a subset of variables from smaller
#' models and moves to the bigger ones.
#'
#' @param y endogenous data
#' @param ySizes a list of model dimension to be estimated in each step.
#' @param counts a list of suggested number of variables to be used in
#' each step. \code{NA} means all variables. Variables are selected based
#' on best estimations (select an appropriate value for
#' \code{searchItems$bestK}). All variables in the best models (all
#' measures and targets) are selected until corresponding suggested
#' number is reached.
#' @param savePre if not \code{NULL}, it saves and tries to load the
#' progress of search step in a file (name=\code{paste0(savePre,i)}
#' where \code{i} is the index of the step).
#' @param ... other arguments to pass to [VarmaSearch()] function such
#' as endogenous data. Note that \code{ySizes} is treated differently.
#'
#' @return A combined \code{LdtSearch} object
#' @export
VarmaSearch_s <- function(y, ySizes = list(c(1, 2), c(3, 4), c(5), c(6:10)),
                          counts = c(NA, 40, 30, 20),
                          savePre = NULL, ...) {
  Search_s("varma", y, ySizes, counts, savePre, ...)
}
