
#' Get Options for ROC and AUC Calculations
#'
#' Use this function to get the required options for [search.dc], [estim.dc], or [s.roc] functions.
#'
#' @param lowerThreshold A number representing the lower bound for calculating partial AUC.
#' @param upperThreshold A number representing the upper bound for calculating partial AUC.
#' @param epsilon A small number used to ignore small floating point differences when comparing scores.
#' @param pessimistic If \code{TRUE}, sequences of equally scored instances are treated differently and a pessimistic measure is calculated (see Fawcett (2006) An introduction to ROC analysis, fig. 6).
#' @param costs The cost of each observation. If \code{NULL}, the cost of all observations will be 1.
#' @param costMatrix A \code{2x2} cost matrix in which: (1,1) is the cost of TN,
#' (2,2) is the cost of TP, (1,2) is the cost of FP and (2,1) is the cost of FN. The first
#' column is multiplied by the corresponding value in the costs vector (see
#' Fawcett (2006), ROC graphs with instance-varying costs).
#'
#' @details
#' See details of [s.roc] function.
#'
#' @return A list with the given options.
#' @examples
#' See [s.roc] function.
#'
#'
#' @export
#' @seealso [search.dc], [estim.dc], [s.roc]
get.options.roc <- function(lowerThreshold = 0, upperThreshold = 1, epsilon = 1e-12,
                            pessimistic = FALSE, costs = NULL, costMatrix = NULL){

  O = list(lowerThreshold = lowerThreshold, upperThreshold = upperThreshold,
           epsilon=epsilon, pessimistic = pessimistic, costs=costs,
           costMatrix = costMatrix)
  O = CheckRocOptions(O)
  O
}

CheckRocOptions <- function(O){
  O$lowerThreshold = as.numeric(O$lowerThreshold);
  O$upperThreshold = as.numeric(O$upperThreshold);
  O$epsilon = as.numeric(O$epsilon);
  O$pessimistic = as.logical(O$pessimistic);

  if (O$lowerThreshold < 0)
    stop("Invalid ROC option. 'lowerThreshold' cannot be negative.")
  if (O$upperThreshold > 1)
    stop("Invalid ROC option. 'upperThreshold' cannot be larger than 1.")
  if (O$epsilon < 0 || O$epsilon > 1)
    stop("Invalid ROC option. 'epsilon' must be positive and less than 1.")

  if (is.null(O$costs) == FALSE) {
    O$costs = as.numeric(O$costs)
    O$costMatrix = as.matrix(O$costMatrix)

    if (is.null(O["costMatrix"]))
      stop("Invalid ROC option. 'costMatrix' cannot be NULL when there is a cost vector.")
    if (nrow(O$costMatrix) != 2 || ncol(O$costMatrix) != 2)
      stop("Invalid ROC option. 'costMatrix' must be a 2x2 matrix.")
  }
  O
}

#' Options for Nelder-Mead Optimization
#'
#'
#'
#' @param maxIterations (int) Maximum number of iterations.
#' @param epsilon (double) A small value to test convergence.
#' @param alpha (double) the reflection coefficient.
#' @param beta (double) the contraction coefficient.
#' @param gamma (double) the expansion coefficient.
#' @param scale (double) A scale in initializing the simplex.
#'
#' @return A list with the given options.
#'
#'  TODO: export
#'
get.options.neldermead <- function(maxIterations = 100, epsilon = 1e-8,
                                   alpha = 1, beta = 0.5, gamma = 2,
                                   scale = 1){
  O = list(maxIterations = maxIterations,
           epsilon = epsilon, alpha = alpha,
           beta = beta, gamma = gamma,
           scale = scale)
  O = CheckNelderMeadOptions(O)
  O
}


CheckNelderMeadOptions <- function(O){
  O$maxIterations <- as.integer(O$maxIterations)
  O$epsilon <- as.numeric(O$epsilon)
  O$alpha <- as.numeric(O$alpha)
  O$beta <- as.numeric(O$beta)
  O$gamma <- as.numeric(O$gamma)
  O$scale <- as.numeric(O$scale)


  if (O$maxIterations <= 0)
    stop("Invalid Nelder-Mead option: 'maxIterations' must be positive.")

  if (O$epsilon < 0)
    stop("Invalid Nelder-Mead option: 'epsilon' cannot be negative.")

  if (O$alpha <= 0)
    stop("Invalid Nelder-Mead option: 'alpha' must be positive.")

  if (O$beta <= 0)
    stop("Invalid Nelder-Mead option: 'beta' must be positive.")

  if (O$gamma <= 0)
    stop("Invalid Nelder-Mead option: 'gamma' must be positive.")

  if (O$scale <= 0)
    stop("Invalid Nelder-Mead option: 'scale' must be positive.")

  O
}



#' Get Options for PCA
#'
#' Use this function to get PCA options in [estim.dc], [estim.sur], [estim.varma], or [s.pca] functions.
#'
#' @param ignoreFirst A number representing the number of variables to exclude at the beginning of data matrices (such as intercept) from PCA.
#' @param exactCount A number that determines the number of components to be used. If zero, the number of components is determined by the \code{cutoffRate}.
#' @param cutoffRate A number between 0 and 1 that determines the cutoff rate for the cumulative variance ratio in order to determine the number of PCA components. It is not used if \code{exactCount} is positive.
#' @param max A number representing the maximum number of components when \code{cutoffRate} is used.
#'
#' @details
#' See details of [s.pca] function.
#'
#' @return A list with the given options.
#' @examples
#' See [s.pca] function.
#'
#' @export
#' @seealso [estim.dc], [estim.sur], [estim.varma], [s.pca]
get.options.pca <- function(ignoreFirst = 1, exactCount = 0, cutoffRate = 0.8, max = 1000){
  O = list(
    ignoreFirst = ignoreFirst, exactCount = exactCount,
    cutoffRate = cutoffRate, max = max)
  O = CheckPcaOptions(O)
  O
}

CheckPcaOptions <- function(O){
  O$ignoreFirst <- as.integer(O$ignoreFirst)
  O$exactCount <- as.integer(O$exactCount)
  O$cutoffRate <- as.numeric(O$cutoffRate)
  O$max <- as.integer(O$max)

  if (O$ignoreFirst < 0)
    stop("Invalid Pca option. 'ignoreFirst' cannot be negative.")

  if (O$exactCount < 0)
    stop("Invalid Pca option. 'exactCount' cannot be negative.")

  if (O$cutoffRate <= 0 || O$cutoffRate > 1)
    stop("Invalid Pca option. 'cutoffRate' must be positive and less than 1.")

  if (O$max <= 0)
    stop("Invalid Pca option. 'max' must be positive.")

  O
}


#' Get Options for L-BFGS Optimization
#'
#' Use this function to get optimization options in [estim.varma] or [search.varma] functions.
#'
#' @param maxIterations A positive integer representing the maximum number of iterations.
#' @param factor A number that determines the condition for stopping the iterations. Use, for example, 1e12 for low accuracy, 1e7 (default) for moderate accuracy, and 1e1 for extremely high accuracy. The default is 1e7.
#' @param projectedGradientTol A number used to stop the iteration using the projected gradient. The default is zero.
#' @param maxCorrections The maximum number of variable metric corrections allowed in the limited memory matrix. The default is 5.
#'
#' @return A list with the given options.
#'
#' @export
get.options.lmbfgs <- function(maxIterations = 100, factor = 1e7,
                               projectedGradientTol = 0, maxCorrections = 5){
  O = list(maxIterations = maxIterations,
           factor = factor,
           projectedGradientTol = projectedGradientTol,
           maxCorrections = maxCorrections)
  O = CheckLmbfgsOptions(O)
  O
}


CheckLmbfgsOptions <- function(O){
  O$maxIterations <- as.integer(O$maxIterations)
  O$factor <- as.numeric(O$factor)
  O$projectedGradientTol <- as.numeric(O$projectedGradientTol)
  O$maxCorrections <- as.integer(O$maxCorrections)


  if (O$maxIterations <= 0)
    stop("Invalid L-BFGS option. 'maxIterations' must be positive.")

  if (O$factor <= 0)
    stop("Invalid L-BFGS option. 'factor' must be positive.")

  if (O$projectedGradientTol < 0)
    stop("Invalid L-BFGS option. 'projectedGradientTol' cannot be negative.")

  if (O$maxCorrections <= 0)
    stop("Invalid L-BFGS option. 'maxCorrections' must be positive.")

  O
}


#' Get Options for Newton Optimization
#'
#' Use this function to get optimization options in [estim.dc] or [search.dc] functions.
#'
#' @param maxIterations An integer representing the maximum number of iterations.
#' @param functionTol A small value used to test the convergence of the objective function.
#' @param gradientTol A small value used to test the convergence of the gradient.
#' @param useLineSearch If \code{TRUE}, line search is used.
#'
#' @return A list with the given options.
#'
#' @export
get.options.newton <- function(maxIterations = 100, functionTol = 1e-4,
                               gradientTol = 0, useLineSearch = TRUE){
  O = list(maxIterations = maxIterations,
           functionTol = functionTol,
           gradientTol = gradientTol,
           useLineSearch = useLineSearch)
  O = CheckNewtonOptions(O)
  O
}


CheckNewtonOptions <- function(O){
  O$maxIterations <- as.integer(O$maxIterations)
  O$functionTol <- as.numeric(O$functionTol)
  O$gradientTol <- as.numeric(O$gradientTol)
  O$useLineSearch <- as.logical(O$useLineSearch)

  if (O$maxIterations <= 0)
    stop("Invalid Newton option. 'maxIterations' must be positive.")

  if (O$functionTol < 0)
    stop("Invalid Newton option. 'functionTol' cannot be negative.")

  if (O$gradientTol < 0)
    stop("Invalid Newton option. 'gradientTol' cannot be negative.")

  O
}


#' Specify the Purpose of the Model Search Process
#'
#' Use this function to list the required items and information that should be saved and retrieved from the model set search process in \code{search.?} functions.
#'
#' @param model If \code{TRUE}, some information about the models is saved.
#' @param type1 If \code{TRUE} and implemented, extra information is saved. This can be the coefficients in the SUR search or predictions in the VARMA search.
#' @param type2 If \code{TRUE} and implemented, extra information is saved. This is similar to \code{type1}. **It is reserved for future updates.**
#' @param bestK The number of best items to be saved in \code{model}, \code{type1}, or \code{type2} information.
#' @param all If \code{TRUE}, all models' information is saved.
#' @param inclusion If \code{TRUE}, inclusion weights are saved.
#' @param cdfs Weighted average of the CDFs at each given point is calculated (for \code{type1} and \code{type2} cases).
#' @param extremeMultiplier A number that determines the multiplier in the extreme bound analysis (for \code{type1} and \code{type2} cases). Use zero to disable it.
#' @param mixture4 If \code{TRUE}, the first four moments of the average distributions are calculated in \code{type1} and \code{type2} cases.
#'
#' @return A list with the given options.
#'
#' @export
get.items.search <- function(model = TRUE, type1 = FALSE, type2 = FALSE,
                             bestK = 1, all = FALSE, inclusion = FALSE,
                             cdfs = numeric(0), extremeMultiplier = 0,
                             mixture4 = FALSE){
  O = list(model = model, type1 = type1,
           type2 = type2, bestK = bestK,
           all = all, inclusion = inclusion,
           cdfs = cdfs,
           extremeMultiplier = extremeMultiplier,
           mixture4 = mixture4)
  O = CheckSearchItems(O)
  O
}


CheckSearchItems <- function(O){
  O$model <- as.logical(O$model)
  O$type1 <- as.logical(O$type1)
  O$bestK <- as.integer(O$bestK)
  O$all <- as.logical(O$all)
  O$inclusion <- as.logical(O$inclusion)
  O$cdfs <- as.numeric(O$cdfs)
  O$extremeMultiplier <- as.numeric(O$extremeMultiplier)
  O$mixture4 <- as.logical(O$mixture4)

  if (O$bestK < 0)
    stop("Invalid Search item. 'bestK' cannot be negative.")

  if (O$extremeMultiplier < 0)
    stop("Invalid Search item. 'extremeMultiplier' cannot be negative.")

  O
}


#' Get Extra Options for Model Search Process
#'
#' Use this function to determine how the model search is performed.
#'
#' @param parallel If \code{TRUE}, a parallel search algorithm is used. This generally changes the speed and memory usage.
#' @param reportInterval An integer representing the time interval (in seconds) for reporting progress (if any significant change has occurred). Set to zero to disable reporting.
#' @param printMsg Set to \code{TRUE} to enable printing details.
#'
#' @return A list with the given options.

#'
#' @export
get.options.search <- function(parallel = FALSE, reportInterval = 2,
                               printMsg = FALSE){
  O = list(parallel = parallel,
           reportInterval = reportInterval,
           printMsg = printMsg)
  O = CheckSearchOptions(O)
  O
}

CheckSearchOptions <- function(O){
  if (O$parallel && .SupportsParallel() == FALSE) {
    O$parallel = FALSE
    warning("Warning: 'parallel' option is not available.")
  }
  O$reportInterval = as.integer(O$reportInterval)

  if (O$reportInterval < 0)
    stop("Invalid Search option. 'reportInterval' cannot be negative.")

  O
}


#' Set Options to Exclude a Model Subset
#'
#' Use this function to determine which models should be skipped in the search process.
#'
#' @param estimation If \code{TRUE}, the model is estimated with all data and is ignored if this estimation fails. If \code{FALSE}, you might get a 'best model' that cannot be estimated.
#' @param maxConditionNumber A number used to ignore an estimation that has a high condition number (if implemented in the search).
#' @param minObsCount An integer used to ignore an estimation where the number of observations (after dealing with \code{NA}) is low. Use 0 to disable this check.
#' @param minDof An integer used to ignore an estimation with low degrees of freedom (equation-wise). Use 0 to disable this check.
#' @param minOutSim An integer used to ignore estimations with a low number of out-of-sample simulations (if implemented in the search).
#' @param minR2 A number used to ignore estimations with a low value for 'R2' (if implemented in the search).
#' @param maxAic A number used to ignore estimations with a high 'AIC' (if implemented in the search).
#' @param maxSic A number used to ignore estimations with a high 'SIC' (if implemented in the search).
#' @param prediction If \code{TRUE}, model data is predicted given all data and is ignored if this process fails. If \code{FALSE}, you might get a 'best model' that cannot be used for prediction.
#' @param predictionBoundMultiplier A positive number used to create a bound and check predictions.
#' The bound is created by multiplying this value by the average growth rate of the data.
#' A model is ignored if its prediction lies outside of this bound. Use zero to disable this check.
#'
#' @return A list with the given options.
#'
#' @export
get.items.modelcheck <- function( estimation = TRUE, maxConditionNumber = Inf,
                                  minObsCount = 0, minDof = 0, minOutSim = 0,
                                  minR2 = -Inf, maxAic = Inf, maxSic = Inf,
                                  prediction = FALSE, predictionBoundMultiplier = 4){

  O = list(
    estimation = estimation,
    maxConditionNumber = maxConditionNumber,
    minObsCount = minObsCount, minDof = minDof,
    minOutSim = minOutSim, maxSic = maxSic,
    minR2 = minR2, maxAic = maxAic,
    maxSic = maxSic, prediction = prediction,
    predictionBoundMultiplier = predictionBoundMultiplier)
  O = CheckModelCheckItems(O)
  O
}


CheckModelCheckItems <- function(O){

  O$estimation = as.logical(O$estimation)
  O$maxConditionNumber = as.numeric(O$maxConditionNumber)
  O$minObsCount = as.integer(O$minObsCount)
  O$minDof = as.integer(O$minDof)
  O$minOutSim = as.integer(O$minOutSim)
  O$maxSic = as.numeric(O$maxSic)
  O$minR2 = as.numeric(O$minR2)
  O$maxAic = as.numeric(O$maxAic)
  O$maxSic = as.numeric(O$maxSic)
  O$prediction = as.logical(O$prediction)
  O$predictionBoundMultiplier = as.numeric(O$predictionBoundMultiplier)

  if (O$minObsCount < 0)
    stop("Invalid model-Check option. 'minObsCount' cannot be negative.")

  if (O$minDof < 0)
    stop("Invalid model-Check option. 'minDof' cannot be negative.")

  if (O$minOutSim < 0)
    stop("Invalid model-Check option. 'minOutSim' cannot be negative.")

  if (O$maxConditionNumber < 0)
    stop("Invalid model-Check option. 'maxConditionNumber' cannot be negative.")

  if (O$predictionBoundMultiplier < 0)
    stop("Invalid model-Check option. 'predictionBoundMultiplier' cannot be negative.")

  O
}

#' Get Options for Measuring Performance
#'
#' Use this function to get measuring options in \code{search.?} functions.
#'
#' @param typesIn A list of evaluation measures when the model is estimated using all available data. It can be \code{aic}, \code{sic}, \code{frequencyCostIn}, or \code{aucIn}. \code{NULL} means no measure.
#' @param typesOut A list of evaluation measures in a pseudo out-of-sample simulation. It can be \code{sign}, \code{direction}, \code{rmse}, \code{scaledRmse}, \code{mae}, \code{scaledMae}, \code{crps}, \code{frequencyCostOut}, or \code{aucOut}. Null means no measure.
#' @param simFixSize An integer that determines the number of pseudo out-of-sample simulations. Use zero to disable the simulation.
#' @param trainFixSize An integer representing the number of data points in the training sample in the pseudo out-of-sample simulation. If zero, \code{trainRatio} will be used.
#' @param trainRatio A number representing the size of the training sample relative to the available size, in the pseudo out-of-sample simulation. It is effective if \code{trainFixSize} is zero.
#' @param seed A seed for the random number generator. Use zero for a random value. It can be negative to get reproducible results between the \code{search.?} function and the \code{estim.?} function.
#' @param horizons An array of integers representing the prediction horizons to be used in pseudo out-of-sample simulations, if the model supports time-series prediction. If \code{NULL}, \code{c(1)} is used.
#' @param weightedEval If \code{TRUE}, weights are used in evaluating discrete-choice models.
#'
#' @return A list with the given options.
#'
#' @export
get.options.measure <- function(typesIn = character(0), typesOut = character(0),
                                simFixSize = 10, trainRatio = 0.75,
                                trainFixSize = 0, seed = 0,
                                horizons = c(1L), weightedEval = FALSE){
  O = list(
    typesIn = typesIn, typesOut = typesOut,
    simFixSize = simFixSize, trainRatio = trainRatio,
    trainFixSize = trainFixSize, seed = seed,
    horizons = horizons, weightedEval = weightedEval)
  O = CheckMeasureOptions(O)
  O
}


CheckMeasureOptions <- function(O){

  O$typesIn = as.character(O$typesIn)
  O$typesOut = as.character(O$typesOut)
  O$weightedEval = as.logical(O$weightedEval)


  if (length(O$typesIn) == 0 && length(O$typesOut) == 0)
    stop("Invalid Measure option. Both 'typesIn' and 'typesOut' are empty.")

  if (length(O$typesOut) > 0) {
    O$horizons = as.integer(O$horizons)
    O$simFixSize = as.integer(O$simFixSize)
    O$trainRatio = as.numeric(O$trainRatio)
    O$trainFixSize = as.integer(O$trainFixSize)
    O$seed = as.numeric(O$seed)

    for (h in O$horizons) {
      if (h <= 0)
        stop("Invalid Measure option. zero or negative value in 'horizons'.")
    }

    if (O$simFixSize < 0)
      stop("Invalid Measure option. 'simFixSize' cannot be negative.")

    if (O$trainRatio < 0 || O$trainRatio > 1)
      stop("Invalid Measure option. 'trainRatio' cannot be negative or cannot be larger than 1.")

    if (O$trainFixSize < 0)
      stop("Invalid Measure option. 'trainFixSize' cannot be negative.")

    # if (O$seed < 0)  It can be negative for similar distribution
    # of the seeds in the searchers

    if (O$trainRatio == 0 && O$trainFixSize == 0)
      stop("Invalid Measure option. Both 'trainRatio' and 'trainFixSize' are zero.")
  }
  O
}
