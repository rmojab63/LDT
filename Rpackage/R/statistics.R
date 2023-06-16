

#' Converts a Measure to Weight
#'
#' This function converts a measure to its weight equivalent.
#'
#' @param value Numeric value of the measure.
#' @param measureName Character string specifying the name of the measure.
#' See [get.options.measure] function for the list of available options.
#'
#' @details
#' Given a collection of models for the data, a measure is not
#' generally a measure of the relative quality of a model. This function
#' converts the value of a measure to such a number.
#'
#' The main purpose of exporting this statistics helper method is to show the inner calculations of the package.
#'
#' @return A numeric value representing the converted measure.
#' @export
#'
#' @examples
#' weight <- s.weight.from.measure(-3.4, "sic")
#' measure <- s.measure.from.weight(weight, "sic")
#'
#' @seealso [s.measure.from.weight]
s.weight.from.measure <- function(value, measureName)
{
  value = as.numeric(value)
  measureName = as.character(measureName)
  res <- .GetWeightFromMeasure(value, measureName)
  res
}


#' Converts a Weight to Measure
#'
#' This function converts a weight to its measure equivalent.
#'
#' @param value Numeric value of the weight.
#' @param measureName Character string specifying the name of the measure.
#' See \code{\link{get.options.measure}} function for the list of available options.
#'
#' @details
#' The main purpose of exporting this statistics helper method is to show the inner calculations of the package.
#'
#' @return A numeric value representing the converted weight.
#' @export
#'
#' @examples
#' weight <- s.weight.from.measure(-3.4, "sic")
#' measure <- s.measure.from.weight(weight, "sic")
#'
#' @seealso [s.weight.from.measure]
s.measure.from.weight <- function(value, measureName)
{
  value = as.numeric(value)
  measureName = as.character(measureName)

  res <- .GetMeasureFromWeight(value, measureName)
  res
}

#' Gets ROC Curve Data for Binary Classification
#'
#' This function calculates the required points for plotting the ROC curve and the AUC.
#'
#' @param y A numeric vector (\code{Nx1}) representing the actual values.
#' @param scores A numeric vector (\code{Nx1}) representing the calculated probabilities for the negative observations.
#' @param weights A numeric vector (\code{Nx1}) representing the weights of the observations.
#' Use \code{NULL} for equal weights.
#' @param options A list from \code{\link{get.options.roc}} function for more options.
#' @param printMsg Set to \code{TRUE} to enable printing some details.
#'
#' @details
#' This is generally a statistics helper method in this package and it shows the inner calculations.
#'
#'
#' @return A list with the following items:
#' \tabular{ll}{
#' \code{n} \tab Number of observations. \cr
#' \code{auc} \tab Value of AUC. \cr
#' \code{points} \tab Points for plotting ROC.
#' }
#'
#' @export
#'
#' @examples
#' y <- c(1, 0, 1, 0, 1, 1, 0, 0, 1, 0)
#' scores <- c(0.1, 0.2, 0.3, 0.5, 0.5, 0.5, 0.7, 0.8, 0.9, 1)
#' res1 <- s.roc(y,scores, printMsg = FALSE)
#' costs <- c(1,2,1,4,1,5,1,1,0.5,1)
#' costMatrix = matrix(c(0.02,-1,-3,3),2,2)
#' opt <- get.options.roc(costs = costs, costMatrix = costMatrix)
#' res2 <- s.roc(y,scores,NULL,options = opt, printMsg = FALSE)
s.roc <- function(y, scores, weights = NULL,
                  options = get.options.roc(), printMsg = FALSE)
{
  y <- as.numeric(y)
  scores <- as.numeric(scores)
  if (length(y) != length(scores))
    stop("Inconsistent length between actual data and the scores.")
  if (is.null(weights) == FALSE){
    weights <- as.numeric(weights)
    if (length(y) != length(weights))
      stop("Inconsistent length between actual data and the weights.")
  }
  if (is.null(options))
    options <- get.options.roc()
  options <-as.list(options)
  printMsg <- as.logical(printMsg)

  res <- .GetRoc(y, scores, weights, options, printMsg)
  res
}


#' Gets the GLD Parameters from the moments
#'
#' Calculates the parameters of the generalized lambda distribution (FKML), given the first four moments of the distribution.
#'
#' @details
#' The type of the distribution is determined by one or two restrictions:
#' - **type 0:** general
#' - **type 1:** symmetric 'type 0'
#' - **type 2:** uni-modal continuous tail: p3<1 & p4<1
#' - **type 3:** symmetric 'type 2' p3==p4
#' - **type 4:** uni-modal continuous tail finite slope  p3<=0.5 &  p4<=5
#' - **type 5:** symmetric 'type 4' p3==p4
#' - **type 6:** uni-modal truncated density curves: p3>=2 & p4>=2 (includes uniform distribution)
#' - **type 7:** symmetric 'type 6' p3==p4
#' - **type 8:** S shaped p3>2 & 1<p4<2 or 1<p3<2 & p4>2
#' - **type 9:** U shaped 1<p3<=2 and 1<p4<=2
#' - **type 10:** symmetric 'type 9' p4==p4
#' - **type 11:** monotone p3>1 & p4<=1
#'
#' @param mean (double) mean of the distribution.
#' @param variance (double) variance of the distribution.
#' @param skewness (double) skewness of the distribution.
#' @param excessKurtosis (double) excess kurtosis of the distribution.
#' @param type (int) The type of the distribution.
#' @param start (numeric vector, length=2) starting value for p3 and p4. Use null for c(0,0).
#' @param nelderMeadOptions (list) The optimization parameters. Use null for default.
#' @param printMsg (bool) If \code{TRUE}, details are printed.
#'
#'
#' @return a vector with the parameters of the GLD distribution.
#' export (TODO: the Nelder-Mead algorithm in the c++ code need revision)
#'
#' @examples
#' res = s.gld.from.moments(0,1,0,0,0,c(0,0))
s.gld.from.moments <- function(mean = 0, variance = 1,
                               skewness = 0, excessKurtosis = 0,
                               type = 0, start = NULL,
                               nelderMeadOptions = get.options.neldermead(),
                               printMsg = FALSE)
{
  res <- .GetGldFromMoments(mean , variance, skewness, excessKurtosis,
                            type, start, nelderMeadOptions, printMsg)
  res
}

#' GLD Quantile Function
#'
#' This function calculates the quantiles of a Generalized Lambda Distribution (FKML).
#'
#' @param probs A numeric vector of probabilities.
#' @param p1 Numeric value representing the first parameter of the distribution (location of the distribution).
#' @param p2 Numeric value representing the second parameter of the distribution (scale of the distribution).
#' @param p3 Numeric value representing the third parameter of the distribution (skewness of the distribution).
#' @param p4 Numeric value representing the fourth parameter of the distribution (kurtosis of the distribution).
#'
#' @details
#' It is a helper statistics method in this package and is generally used to plot density function of a GLD distribution.
#' See the example of \code{\link{s.gld.density.quantile}} function for more details.
#'
#' @return A numeric vector representing the quantiles for each probability in \code{probs}.
#' @export
#' @examples
#' res = s.gld.quantile(c(0.1,0.5,0.95), 0,1,0,0) # standard normal distribution
#'
#' @seealso [s.gld.density.quantile]
s.gld.quantile <- function(probs, p1, p2, p3, p4)
{
  probs <- as.numeric(probs)
  p1 <- as.numeric(p1)
  p2 <- as.numeric(p2)
  p3 <- as.numeric(p3)
  p4 <- as.numeric(p4)

  res <- .GldQuantile(probs, p1, p2, p3, p4)
  res
}

#' GLD Density-Quantile Function
#'
#' This function calculates the densities of a Generalized Lambda Distribution (FKLM) given a vector of probabilities.
#'
#' @param probs A numeric vector representing the probabilities.
#' @param p1 Numeric value representing the first parameter (location) of the distribution.
#' @param p2 Numeric value representing the second parameter (scale) of the distribution.
#' @param p3 Numeric value representing the third parameter (skewness) of the distribution.
#' @param p4 Numeric value representing the fourth parameter (kurtosis) of the distribution.
#'
#' @details
#' It is a helper statistics method in this package and is generally used to plot density function of a GLD distribution.
#'
#' @return A numeric vector representing the densities for each probability in \code{probs}.
#' @export
#' @examples
#' # In this example we use this function and plot the density function for
#' # standard normal distribution:
#' probs <- seq(0.1,0.9,0.1)
#' x = s.gld.quantile(probs, 0,1,0,0)
#' y = s.gld.density.quantile(probs, 0,1,0,0)
#' # plot(x,y)
#' # lines(x,y)
#' @seealso [s.gld.quantile]
s.gld.density.quantile <- function(probs, p1, p2, p3, p4)
{
  probs <- as.numeric(probs)
  p1 <- as.numeric(p1)
  p2 <- as.numeric(p2)
  p3 <- as.numeric(p3)
  p4 <- as.numeric(p4)

  res <- .GldDensityQuantile(probs, p1, p2, p3, p4)
  res
}

#' Combines Two Distributions Given their First 4 Moments
#'
#' This function combines two distributions and calculates the first 4 moments of the combined distribution.
#'
#' @param mix1 A list representing the first distribution with elements \code{mean}, \code{variance}, \code{skewness}, \code{kurtosis}, \code{sumWeights}, and \code{count}.
#' @param mix2 A list representing the second distribution (similar to \code{mix1}).
#'
#' @details
#' Let \code{A} and \code{B} be two distributions. By combining these two distributions we mean
#' if we get a sample of size \code{n} from both distributions and the combined distribution,
#' the first four moments of the combined samples equals the first four moments of the
#' sample from the combined distribution, as \code{n} goes to infinity.
#'
#'
#' @return A list similar to \code{mix1}.
#' @export
#'
#' @examples
#' n <- 1000 # sample size (increase it for more accurate result)
#' sample1 <- rchisq(n,3)
#' sample2 <- rchisq(n,5)
#'
#' d1 <- list(mean = mean(sample1),
#'            variance = var(sample1),
#'            skewness = moments::skewness(sample1),
#'            kurtosis = moments::kurtosis(sample1),
#'            count=length(sample1),
#'            sumWeights = length(sample1))
#' d2 <- list(mean = mean(sample2),
#'            variance = var(sample2),
#'            skewness = moments::skewness(sample2),
#'            kurtosis = moments::kurtosis(sample2),
#'            count=length(sample2),
#'            sumWeights = length(sample2))
#' c <- s.combine.by.moments4(d1,d2)
#'
#' # we can compare the results:
#' combined <- c(x1,x2)
#' mean_c = mean(combined)
#' variance_c = var(combined)
#' skewness_c = moments::skewness(combined)
#' kurtosis_c = moments::kurtosis(combined)
#'
s.combine.by.moments4 <- function(mix1, mix2)
{
  mix1 <- as.list(mix1)
  mix2 <- as.list(mix2)

  res <- .CombineByMoments4(mix1, mix2)
  res
}

#' Principal Component Analysis
#'
#' This function performs PCA on the columns of a matrix.
#'
#' @param x A numeric matrix with variables in the columns.
#' @param center Logical value indicating whether to demean the columns of \code{x}.
#' @param scale Logical value indicating whether to scale the columns of \code{x} to unit variance.
#' @param newX A numeric matrix to be used in projection.
#' Its structure must be similar to \code{x}.
#'
#' @details
#' The main purpose of exporting this statistics helper method is to show the inner calculations of the package.
#'
#' @return A list with the following items:
#' \tabular{ll}{
#' \code{removed0Var} \tab An integer vector showing the zero-based indices of removed columns with zero variances.\cr
#' \code{directions} \tab Directions matrix.\cr
#' \code{stds} \tab An integer vector showing the standard deviation of the principal components.\cr
#' \code{stds2Ratio} \tab Shows \code{stds^2/sum(stds^2)}.\cr
#' \code{projections} \tab Projections matrix if \code{newX} is provided.
#' }
#' @export
#'
#' @examples
#' data <- data.frame( x = rnorm(100), y = rnorm(100), z = rep(0, 100))
#' res <- s.pca(data)
#'
#' # Note that one of the columns has zero variance.
#' res_invalid <- try(stats::prcomp(data, center = TRUE,
#'                    scale. = TRUE)) # We should remove 'z' first
#'
s.pca <- function(x, center = TRUE,
                  scale = TRUE, newX = NULL)
{
  x <- as.matrix(x)
  center <- as.logical(center)
  scale <- as.logical(scale)
  if (is.null(newX) == FALSE){
    newX <- as.matrix(newX)
    if (ncol(x) != ncol(newX))
      stop("Inconsistent number of columns between 'x' and 'newX'.")
  }

  res <- .GetPca(x, center, scale, newX)
  res
}


#' Gets the Distances Between Variables
#'
#' This function calculates the distances between the columns of a numeric matrix.
#'
#' @param data A numeric matrix with variables in the columns.
#' @param distance Character string specifying the type of distance.
#' It can be \code{correlation}, \code{absCorrelation}, \code{euclidean}, \code{manhattan}, or \code{maximum}.
#' @param correlation Character string specifying the type of correlation if \code{distance} is correlation.
#' It can be \code{pearson} or \code{spearman}.
#' @param checkNan Logical value indicating whether to check for \code{NA}s (and omit them if any exist).
#'
#' @details
#' The main purpose of exporting this statistics helper method is to show the inner calculations of the package.
#'
#' @return A symmetric matrix (lower triangle as a vector).
#'
#' @export
#' @examples
#' n <- 10
#' data <- data.frame(x = rnorm(n), y = rnorm(n), z = rnorm(n))
#' distances <- s.distance(data)
#'
s.distance <- function(data, distance = "correlation",
                       correlation = "pearson", checkNan = TRUE) {

  data <- as.matrix(data)
  distance <- as.character(distance)
  correlation <- as.character(correlation)
  checkNan <- as.logical(checkNan)

  res <- .GetDistance(data, distance, correlation, checkNan)
  res
}

#' Hierarchical Clustering
#'
#' This function performs hierarchical clustering on a group of variables, given their distances from each other.
#'
#' @param distances Lower triangle of a symmetric distance matrix (without the diagonal).
#' This can be the output of \code{\link{s.distance}} function.
#' @param linkage Character string specifying the method for calculating the distance in a left-right node merge.
#' It can be \code{single}, \code{complete}, \code{uAverage}, \code{wAverage}, or \code{ward}.
#'
#' @details
#' The main purpose of exporting this statistics helper method is to show the inner calculations of the package.
#'
#'
#' @return A list with the following items:
#' \tabular{ll}{
#' \code{merge} \tab An integer matrix representing the merge matrix. \cr
#' \code{height} \tab A numeric vector representing the heights. \cr
#' \code{order} \tab An integer vector representing the orders.
#' }
#'
#' @export
#' @examples
#' n <- 10
#' data <- data.frame(x = rnorm(n), y = rnorm(n), z = rnorm(n))
#' distances <- s.distance(data)
#' clusters <- s.cluster.h(distances)
#'
s.cluster.h <- function(distances, linkage = "single"){

  distances <- as.numeric(distances)
  linkage <- as.character(linkage)

  res <- .ClusterH(distances, linkage)
  res
}


#' Groups Variables with Hierarchical Clustering
#'
#' This function groups the columns of a numeric matrix based on the hierarchical clustering algorithm.
#'
#' @param data A numeric matrix with variables in the columns.
#' @param nGroups Integer value specifying the number of required groups.
#' @param threshold Numeric value specifying a threshold for omitting variables.
#' If the distance between two variables in a group is less than this value, the second one will be omitted.
#' Note that a change in the order of the columns might change the results.
#' @param distance Character string specifying how distances are calculated.
#' It can be \code{correlation}, \code{absCorrelation}, \code{euclidean}, \code{manhattan}, or \code{maximum}.
#' See \code{\link{s.distance}} function.
#' @param linkage Character string specifying how distances are calculated in a left-right node merge.
#' It can be \code{single}, \code{complete}, \code{uAverage}, \code{wAverage}, or \code{ward}.
#' See \code{\link{s.cluster.h}} function.
#' @param correlation Character string specifying the type of correlation if \code{distance} is correlation.
#' It can be \code{pearson} or \code{spearman}. See \code{\link{s.distance}} function.
#'
#' @details
#' The results might be different from R's 'cutree' function.
#' (I don't know how 'cutree' works) Here this function iterates over the nodes and
#' whenever a split occurs, it adds a group until the required number of groups is reached.
#'
#' @return A list with the following items:
#' \tabular{ll}{
#' \code{groups} \tab A list of integer vectors representing the indexes of variables in each group. \cr
#' \code{removed} \tab An integer vector representing the indexes of removed variables.
#' }
#'
#' @export
s.cluster.h.group <- function(data, nGroups = 2, threshold = 0,
                             distance = "correlation",
                             linkage = "single",
                             correlation = "pearson")

{
  data <- as.matrix(data)
  nGroups <- as.integer(nGroups)
  threshold <- as.numeric(threshold)
  distance <- as.character(distance)
  linkage <- as.character(linkage)
  correlation <- as.character(correlation)

  res <- .ClusterHGroup(data, nGroups, threshold,
                        distance, linkage, correlation)
  res
}
