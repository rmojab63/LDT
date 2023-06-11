

#' Converts a Measure to Weight
#'
#' @param value (double) the measure
#' @param measureName (string) measure name
#'
#' @return the weight
#' @export
#'
#' @examples
#' weight <- GetWeightFromMeasure(-3.4, "sic")
GetWeightFromMeasure <- function(value, measureName)
{
  res <- .GetWeightFromMeasure(value, measureName)
  res
}


#' Converts a Measure to Weight
#'
#' @param value (double) the measure
#' @param measureName (string) measure name
#'
#' @return the measure
#' @export
#'
#' @examples
#' weight <- GetWeightFromMeasure(-3.4, "sic")
#' measure <- GetMeasureFromWeight(weight, "sic")
GetMeasureFromWeight <- function(value, measureName)
{
  res <- .GetMeasureFromWeight(value, measureName)
  res
}

#' ROC curve for a binary case
#'
#' It does not draw the ROC, but calculatates the required points. It also
#' Calculates the AUC with different options
#'
#' @param y (numeric vector, \code{Nx1}) Actual values
#' @param scores (numeric vector, \code{Nx1}) Calculated probabilities for the negative observations
#' @param weights (numeric vector, \code{Nx1}) Weights of the observations. Use \code{NULL} for equal weights.
#' @param options (list) More options. See [GetRocOptions()] function for details.
#' @param printMsg (bool) Set TRUE to report some details.
#'
#' @return A list with the following items:
#' \item{N}{(integer) Number of observations}
#' \item{AUC}{(numeric) Value of AUC}
#' \item{Points}{(numeric matrix) Points for ploting ROC}
#'
#' @export
#'
#' @examples
#' y <- c(1, 0, 1, 0, 1, 1, 0, 0, 1, 0)
#' scores <- c(0.1, 0.2, 0.3, 0.5, 0.5, 0.5, 0.7, 0.8, 0.9, 1)
#' res1 = GetRoc(y,scores, printMsg = FALSE)
#' costs <- c(1,2,1,4,1,5,1,1,0.5,1)
#' costMatrix = matrix(c(0.02,-1,-3,3),2,2)
#' opt <- GetRocOptions(costs = costs, costMatrix = costMatrix)
#' res2 = GetRoc(y,scores,NULL,options = opt, printMsg = FALSE)
#' #plot(res1$Points)
#' #lines(res2$Points)
#'
GetRoc <- function(y, scores, weights = NULL, options = GetRocOptions(),
            printMsg = FALSE)
{
  res <- .GetRoc(y, scores, weights, options, printMsg)
  res
}


#' Gets the GLD-FKML Parameters from the moments
#'
#' @description Calculates the parameters of the generalized lambda distribution (FKML), given the first four moments of the distribution.
#'
#' @details
#' The type of the distribution is determined by one or two restrictions:
#' - **type 0:** general
#' - **type 1:** symmetric 'type 0'
#' - **type 2:** uni-modal continuous tail: L3<1 & L4<1
#' - **type 3:** symmetric 'type 2' L3==L4
#' - **type 4:** uni-modal continuous tail finite slope  L3<=0.5 &  L4<=5
#' - **type 5:** symmetric 'type 4' L3==L4
#' - **type 6:** uni-modal truncated density curves: L3>=2 & L4>=2 (includes uniform distribution)
#' - **type 7:** symmetric 'type 6' L3==L4
#' - **type 8:** S shaped L3>2 & 1<L4<2 or 1<L3<2 & L4>2
#' - **type 9:** U shaped 1<L3<=2 and 1<L4<=2
#' - **type 10:** symmetric 'type 9' L4==L4
#' - **type 11:** monotone L3>1 & L4<=1
#'
#' @param mean (double) mean of the distribution.
#' @param variance (double) variance of the distribution.
#' @param skewness (double) skewness of the distribution.
#' @param excessKurtosis (double) excess kurtosis of the distribution.
#' @param type (int) The type of the distribution.
#' @param start (numeric vector, length=2) starting value for L3 and L4. Use null for c(0,0).
#' @param nelderMeadOptions (list) The optimization parameters. Use null for default.
#' @param printMsg (bool) If \code{TRUE}, details are printed.
#'
#' @return a vector with the parameters of the GLD distribution.
#' @export
#'
#' @examples
#' res = GetGldFromMoments(0,1,0,0,0,c(0,0))
GetGldFromMoments <- function(mean = 0, variance = 1,
                                skewness = 0, excessKurtosis = 0,
                                type = 0, start = NULL,
                                nelderMeadOptions = GetNelderMeadOptions(),
                                printMsg = FALSE)
{
  res <- .GetGldFromMoments(mean , variance, skewness, excessKurtosis,
                  type, start, nelderMeadOptions, printMsg)
  res
}

#' Gets GLD Quantile
#'
#' @param data (numeric vector) data
#' @param L1 (double) First parameter
#' @param L2 (double) Second parameter
#' @param L3 (double) Third parameter
#' @param L4 (double) Fourth parameter
#'
#' @return (numeric vector) result
#' @export
GldQuantile <- function(data, L1, L2, L3, L4)
{
  res <- .GldQuantile(data, L1, L2, L3, L4)
  res
}

#' Gets GLD Density Quantile
#'
#' @param data (numeric vector) data
#' @param L1 (double) First parameter
#' @param L2 (double) Second parameter
#' @param L3 (double) Third parameter
#' @param L4 (double) Fourth parameter
#'
#' @return (numeric vector) result
#' @export
GldDensityQuantile <- function(data, L1, L2, L3, L4)
{
  res <- .GldDensityQuantile(data, L1, L2, L3, L4)
  res
}

#' Combines Two Distributions Defined by their First 4 Moments
#'
#' @param mix1 (list) First distribution which is defined by a list with mean, variance, skewness, kurtosis, sumWeights, count
#' @param mix2 (list) Second distribution (similar to \code{mix1}).
#'
#' @return (list) A list similar to \code{mix1}
#' @export
#'
#' @examples
#' #see its \code{test_that} function
GetCombination4Moments <- function(mix1, mix2)
{
  res <- .GetCombination4Moments(mix1, mix2)
  res
}

#' Principle Component Analysis
#'
#' @param x (numeric matrix) data with variables in columns.
#' @param center (bool) if \code{TRUE}, it demeans the variables.
#' @param scale (bool) if \code{TRUE}, it scales the variables to unit variance.
#' @param newX (numeric matrix) data to be used in projection. Its structure must be similar to the \code{x}.
#'
#' @return (list) results
#' \item{removed0Var}{(integer vector) Zero-based indices of removed columns with zero variances.}
#' \item{directions}{(numeric matrix) Directions}
#' \item{stds}{(integer vector) Standard deviation of the principle components}
#' \item{stds2Ratio}{(integer vector) stds^2/sum(stds^2)}
#' \item{projections}{(numeric matrix) Projections if \code{newX} is given.}
#'
#' @export
#'
GetPca <- function(x, center = TRUE, scale = TRUE, newX = NULL)
{
  res <- .GetPca(x, center, scale, newX)
  res
}



#' Gets Distances Between Variables
#'
#' @param data (numeric matrix) Data with variables in the columns.
#' @param distance (string) Determines how distances are calculated. It can be \code{correlation}, \code{absCorrelation}, \code{euclidean}, \code{manhattan}, \code{maximum}.
#' @param correlation (string) If \code{distance} is correlation, it determines the type of the correlation. It can be \code{pearson}, \code{spearman}.
#' @param checkNan (bool) If false, \code{NAN}s are not omitted.
#'
#' @return A symmetric matrix (lower triangle as a vector).
#'
#' @export
GetDistance <- function(data,
                        distance = "correlation",
                        correlation = "pearson",
                        checkNan = TRUE) {
  res <- .GetDistance(data, distance, correlation, checkNan)
  res
}

#' Hierarchical Clustering
#'
#'
#' @param distances (numeric vector) Determines the distances. This must be the lower triangle of a (symmetric) distance matrix (without the diagonal).
#' @param numVariables (int) Determines the number of variables. This should hold: '2 * length(\code{distances}) = \code{numVariables}(\code{numVariables} - 1)'.
#' @param linkage (string) Determines how Distances are calculated in a left-right node merge. It can be \code{single}, \code{complete}, \code{uAverage}, \code{wAverage}, \code{ward}.
#'
#' @return A list:
#' \item{merge}{(integer matrix)}
#' \item{height}{(numeric vector)}
#' \item{order}{(integer vector)}
#'
#' @export
ClusterH <- function(distances, numVariables,
                     linkage = "single"){
  res <- .ClusterH(distances, numVariables, linkage)
  res
}


#' Groups Variables with Hierarchical Clustering
#'
#' @details The results might be different from R's 'cutree' function. I don't know how 'cutree' works, but here I iterate over the nodes and whenever a split occurs, I add a group until the required number of groups is reached.
#'
#' @param data (numeric matrix) Data with variables in the columns.
#' @param nGroups (int) Number of groups
#' @param threshold (double) A threshold for omitting variables. If distance between two variables in a group is less than this value, the second one will be omitted. Note that a change in the order of the columns might change the results.
#' @param distance (string)  Determines how distances are calculated. It can be \code{correlation}, \code{absCorrelation}, \code{euclidean}, \code{manhattan}, \code{maximum}.
#' @param linkage (string) Determines how Distances are calculated in a left-right node merge. It can be \code{single}, \code{complete}, \code{uAverage}, \code{wAverage}, \code{ward}.
#' @param correlation (string) If \code{distance} is correlation, it determines the type of the correlation. It can be \code{pearson}, \code{spearman}.
#'
#' @return A list:
#' \item{groups}{(of integer vectors) indexes of variables in each group.}
#' \item{removed}{(integer vector) indexes of removed variables.}
#'
#' @export
ClusterHGroup <- function(data, nGroups = 2, threshold = 0,
                          distance = "correlation",
                          linkage = "single",
                          correlation = "pearson")

{
  res <- .ClusterHGroup(data, nGroups, threshold,
                        distance, linkage, correlation)
  res
}
