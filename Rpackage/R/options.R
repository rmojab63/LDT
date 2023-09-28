
#' Get Options for ROC and AUC Calculations
#'
#' Use this function to get the required options for [search.bin], [estim.bin], or [s.roc] functions.
#'
#' @param lowerThreshold A number representing the lower bound for calculating partial AUC.
#' @param upperThreshold A number representing the upper bound for calculating partial AUC.
#' @param epsilon A small number used to ignore small floating point differences when comparing scores.
#' @param pessimistic If \code{TRUE}, sequences of equally scored instances are treated differently and a pessimistic metric is calculated (see Fawcett (2006) An introduction to ROC analysis, fig. 6).
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
#' # See 's.roc' function.
#'
#'
#' @export
#' @seealso [search.bin], [estim.bin], [s.roc]
get.options.roc <- function(lowerThreshold = 0, upperThreshold = 1, epsilon = 1e-12,
                            pessimistic = FALSE, costs = NULL, costMatrix = NULL){
  stopifnot(is.numeric(lowerThreshold) && length(lowerThreshold) == 1 && lowerThreshold >= 0)
  stopifnot(is.numeric(upperThreshold) && length(upperThreshold) == 1 && lowerThreshold <= 1)
  stopifnot(is.numeric(epsilon) && length(epsilon) == 1 && epsilon >= 0 && epsilon <= 1)

  if (!is.null(costs)) {
    stopifnot(is.numeric(costs))
    stopifnot(!is.null(costMatrix) && is.matrix(costMatrix) &&
                nrow(costMatrix) == 2 && ncol(costMatrix) == 2)
  }

  res = list(lowerThreshold = lowerThreshold,
             upperThreshold = upperThreshold,
             epsilon=epsilon,
             pessimistic = pessimistic,
             costs=costs,
             costMatrix = costMatrix)
  class(res) <- c("ldt.options.roc", "list")
  res
}



#' Options for Nelder-Mead Optimization
#'
#' Use this function to get the required options when Nelder-Mead optimization is needed such as [s.gld.from.moments] function.
#'
#' @param maxIterations (int) Maximum number of iterations.
#' @param tolerance A small number to determine the convergence.
#' The algorithm terminates when the difference between the best and worst points in the simplex is less than this value.
#' @param reflection A number for reflection coefficient.
#' It controls how far the worst point is reflected through the centroid of the remaining points.
#' @param expansion A number that determines the expansion coefficient.
#' It controls how far the reflected point is expanded along the line connecting it to the centroid.
#' @param contraction A number that determines the contraction coefficient.
#' It controls how far the worst point is contracted towards the centroid.
#' @param shrink A number that determines the shrink coefficient.
#' It controls how much the simplex is shrunk towards the best point when all other moves are rejected.
#'
#' @return A list with the given options.
#' @export
get.options.neldermead <- function(maxIterations = 100, tolerance = 1e-6,
                                   reflection = 1, expansion = 2, contraction = 0.5,
                                   shrink = 1){
  stopifnot(is.numeric(maxIterations) && length(maxIterations) == 1 && maxIterations > 0)
  stopifnot(is.numeric(tolerance) && length(tolerance) == 1 && tolerance >= 0)
  stopifnot(is.numeric(reflection) && length(reflection) == 1 && reflection > 0)
  stopifnot(is.numeric(expansion) && length(expansion) == 1 && expansion > 0)
  stopifnot(is.numeric(contraction) && length(contraction) == 1 && contraction > 0)
  stopifnot(is.numeric(shrink) && length(shrink) == 1 && shrink > 0)

  res = list(maxIterations = maxIterations,
           tolerance = tolerance, reflection = reflection,
           expansion = expansion, contraction = contraction,
           shrink = shrink)
  class(res) <- c("ldt.options.neldermead", "list")
  res
  res
}


#' Get Options for PCA
#'
#' Use this function to get PCA options in [estim.bin], [estim.sur], [estim.varma], or [s.pca] functions.
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
#' # See 's.pca' function.
#'
#' @export
#' @seealso [estim.bin], [estim.sur], [estim.varma], [s.pca]
get.options.pca <- function(ignoreFirst = 1, exactCount = 0, cutoffRate = 0.8, max = 1000){

  stopifnot(is.numeric(ignoreFirst) && length(ignoreFirst) == 1 && ignoreFirst >= 0)
  stopifnot(is.numeric(exactCount) && length(exactCount) == 1 && exactCount >= 0)
  stopifnot(is.numeric(cutoffRate) && length(cutoffRate) == 1 && cutoffRate > 0 && cutoffRate < 1)
  stopifnot(is.numeric(max) && length(max) == 1 && max > 0)

  res = list(
    ignoreFirst = ignoreFirst,
    exactCount = exactCount,
    cutoffRate = cutoffRate,
    max = max)
  class(res) <- c("ldt.options.pca", "list")
  res
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
get.options.lbfgs <- function(maxIterations = 100, factor = 1e7,
                               projectedGradientTol = 0, maxCorrections = 5){
  stopifnot(is.numeric(maxIterations) && length(maxIterations) == 1 && maxIterations > 0)
  stopifnot(is.numeric(factor) && length(factor) == 1 && factor > 0)
  stopifnot(is.numeric(projectedGradientTol) && length(projectedGradientTol) == 1 && projectedGradientTol >= 0)
  stopifnot(is.numeric(maxCorrections) && length(maxCorrections) == 1 && maxCorrections > 0)

  res = list(maxIterations = maxIterations,
           factor = factor,
           projectedGradientTol = projectedGradientTol,
           maxCorrections = maxCorrections)
  class(res) <- c("ldt.options.lbfgs", "list")
  res
}


#' Get Options for Newton Optimization
#'
#' Use this function to get optimization options in [estim.bin] or [search.bin] functions.
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

  stopifnot(is.numeric(maxIterations) && length(maxIterations) == 1 && maxIterations > 0)
  stopifnot(is.numeric(functionTol) && length(functionTol) == 1 && functionTol >= 0)
  stopifnot(is.numeric(gradientTol) && length(gradientTol) == 1 && gradientTol >= 0)
  stopifnot(is.logical(useLineSearch) && length(useLineSearch) == 1)


  res = list(maxIterations = maxIterations,
           functionTol = functionTol,
           gradientTol = gradientTol,
           useLineSearch = useLineSearch)
  class(res) <- c("ldt.options.newton", "list")
  res
}

