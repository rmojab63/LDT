
#' Options for ROC and AUC
#'
#' @param lowerThreshold (double) Lower bound for calculating partial AUC.
#' @param upperThreshold (double) Upper bound for calculating partial AUC.
#' @param epsilon (double) A value to ignore small floating podifferences in comparing scores.
#' @param pessimistic (bool) If TRUE, sequences of equally scored instances are treated differently and a pessimistic measure is calculated (see Fawcett (2006) An introduction to roc analysis, fig. 6).
#' @param costs (numeric vector) cost of each observations. If null, cost of all observations will be 1.
#' @param costMatrix (numeric matrix) a 2x2 cost matrix in which: (1,1) is cost of TN,
#' (2,2) is cost of TP, (1,2) is cost of FP and (2,1) is cost of FN. First
#' column is multiplied by the corresponding value in costs vector (see
#' Fawcett (2006, ROC graphs with instance-varying costs).
#'
#' @return A list with the given options.
#'
#' @export
GetRocOptions <- function(lowerThreshold = 0, upperThreshold = 1, epsilon = 1e-12,
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
#' @param maxIterations (int) Maximum number of iterations.
#' @param epsilon (double) A small value to test convergence.
#' @param alpha (double) the reflection coefficient.
#' @param beta (double) the contraction coefficient.
#' @param gamma (double) the expansion coefficient.
#' @param scale (double) A scale in initializing the simplex.
#'
#' @return A list with the given options.
#'
#' @export
GetNelderMeadOptions <- function(maxIterations = 100, epsilon = 1e-8,
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



#' Options for PCA
#'
#' @param ignoreFirst (int) Excludes variables at the beginning of data matrices (such as intercept) from PCA.
#' @param exactCount (int) Determines the number of components to be used. If zero, number of components are determined by the \code{cutoffRate}.
#' @param cutoffRate (between 0 and 1) Determines the cutoff rate for cumulative variance ratio in order to determine the number of PCA components. It is not used if \code{exactCount} is positive.
#' @param max (int) Maximum number of components when \code{cutoffRate} is used.
#'
#' @return A list with the given options.
#'
#' @export
GetPcaOptions <- function(ignoreFirst = 1, exactCount = 0, cutoffRate = 0.8, max = 1000){
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


#' Options for LMBFGS Optimization
#'
#' @param maxIterations (int) A positive integer for maximum number of iterations.
#' @param factor (double) A condition for stopping the iterations. The iteration will stop when (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} < \code{factor}*epsmch where epsmch is the machine precision, which is automatically generated by the code. Use e.g., 1e12 for low accuracy, 1e7 (default) for moderate accuracy and 1e1 for extremely high accuracy. default is 1e7
#' @param projectedGradientTol (double) The iteration will stop when \code{max{|proj g_i | i = 1, ..., n} < projectedGradientTol} where \code{pg_i} is the ith component of the projected gradient. default is zero.
#' @param maxCorrections (int) Maximum number of variable metric corrections allowed in the limited memory Matrix. default is 5.
#'
#' @return A list with the given options.
#'
#' @export
GetLmbfgsOptions <- function(maxIterations = 100, factor = 1e7,
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
    stop("Invalid LMBFGS option. 'maxIterations' must be positive.")

  if (O$factor <= 0)
    stop("Invalid LMBFGS option. 'factor' must be positive.")

  if (O$projectedGradientTol < 0)
    stop("Invalid LMBFGS option. 'projectedGradientTol' cannot be negative.")

  if (O$maxCorrections <= 0)
    stop("Invalid LMBFGS option. 'maxCorrections' must be positive.")

  O
}


#' Options for Newton Optimization
#'
#' @param maxIterations (int) Maximum number of iterations.
#' @param functionTol (double) A small value to test convergence of the objective function.
#' @param gradientTol (double) A small value to test convergence of the gradient.
#' @param useLineSearch (bool) If TRUE, it uses line search.
#'
#' @return A list with the given options.
#'
#' @export
GetNewtonOptions <- function(maxIterations = 100, functionTol = 1e-4,
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


#' Options for 'Search Items'
#'
#' @description Creates a with predefined items which determines the information to be saved and retrieved.
#'
#' @param model (bool) If TRUE, information about the models is saved.
#' @param type1 (bool) If TRUE and implemented, extra information is saved. This can be the coefficients in the SUR search or predictions in VARMA search.
#' @param type2 (bool) If TRUE and implemented, extra information is saved. This is similar to \code{type1}. **It is reserved for future updates.**
#' @param bestK (int) Number of best items to be saved in \code{model}, \code{type1}, or \code{type2} information.
#' @param all (bool) If TRUE, all available information is saved.
#' @param inclusion (bool) If TRUE, inclusion weights are saved in \code{model}.
#' @param cdfs (nullable numeric vector) Weighted average of the CDFs at each given pois calculated (for \code{type1} and \code{type2}).
#' @param extremeMultiplier (double) Determined the multiplier in the extreme bound analysis (for \code{type1} and \code{type2}).
#' @param mixture4 (bool) If TRUE, the first 4 moments of the average distributions are calculated in \code{type1} and \code{type2}.
#'
#' @return A list with the given options.
#'
#' @export
GetSearchItems <- function(model = TRUE, type1 = FALSE, type2 = FALSE,
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


#' Options for 'Search Options'
#'
#' @description Creates a with predefined Search O.
#'
#' @param parallel (bool) If TRUE, it uses a parallel search. It generally changes the speed and memory usage.
#' @param reportInterval (int) Time interval (in seconds) for reporting the progress (if the change is significant). Set zero to disable.
#' @param printMsg (bool) Set FALSE to disable printing the details.
#'
#' @return A list with the given options.
#'
#' @export
GetSearchOptions <- function(parallel = FALSE, reportInterval = 2,
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


#' Options for 'Model Check Items'
#'
#' @param estimation (bool) If TRUE, model is estimated with all data. If FALSE, you might get a 'best model' that cannot be estimated.
#' @param maxConditionNumber (double) Maximum value for the condition number (if implemented in the search).
#' @param minObsCount (int) Minimum value for the number of observations. Use 0 to disable.
#' @param minDof (int) Minimum value for the degrees of freedom (equation-wise). Use 0 to disable.
#' @param minOutSim (int) Minimum value for the number of valid out-of-sample simulations (if implemented in the search).
#' @param minR2 (double) Minimum value for R2 (if implemented in the search).
#' @param maxAic (double) Maximum value for AIC (if implemented in the search).
#' @param maxSic (double) Maximum value for SIC (if implemented in the search).
#' @param prediction (bool) If TRUE, model data is predicted given all data. If FALSE, you might get a 'best model' that cannot be used in prediction.
#' @param predictionBoundMultiplier (double) If positive, a bound is created by multiplying this value to the average growth rate. A model is ignored, if its prediction lies outside of this bound.
#'
#' @return A list with the given options.
#'
#' @export
GetModelCheckItems <- function( estimation = TRUE, maxConditionNumber = Inf,
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

#' Options for 'Measuring Performance'
#'
#' @param typesIn (nullable string vector) Evaluations when model is estimated using all available data. It can be \code{aic}, \code{sic}, \code{frequencyCostIn}, \code{aucIn}. Null means no measure.
#' @param typesOut (nullable string vector) Evaluations in an pseudo out-of-sample simulation. It can be \code{sign}, \code{direction}, \code{rmse}, \code{scaledRmse}, \code{mae}, \code{scaledMae}, \code{crps}, \code{frequencyCostOut}, \code{aucOut}. Null means no measure.
#' @param simFixSize (int) Number of pseudo out-of-sample simulations. Use zero to disable the simulation.
#' @param trainFixSize (int) Number of data-points in the training sample in the pseudo out-of-sample simulation. If zero, \code{trainRatio} will be used.
#' @param trainRatio (double) Number of data-points, as a ratio of the available size, in the training sample in the pseudo out-of-sample simulation.
#' @param seed (int) A seed for random number generator. Use zero for a random value.
#' @param horizons (nullable integer vector) prediction horizons to be used in pseudo out-of-sample simulations, if model supports time-series prediction. If null, c(1) is used.
#' @param weightedEval (bool) If TRUE, weights are used in evaluationg discrete-choice models
#' @return A list with the given options.
#'
#' @export
GetMeasureOptions <- function(typesIn = character(0), typesOut = character(0),
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
