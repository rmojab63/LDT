

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




#' Calculate Long-run Growth
#'
#'
#' @param data (integer vector) data
#' @param trimStart (integer) if the number of leading NAs is larger
#' than this number, it returns NA. Otherwise, it finds the first
#' number and continue the calculations.
#' @param trimEnd (integer) if the number of trailing NAs is larger
#' than this number, it returns NA. Otherwise, it finds the first number
#' and continue the calculations.
#' @param cont (logical) if \code{TRUE} it will use the continuous formula.
#' @param skipZero (logical) if \code{TRUE} leading and trailing
#' zeros are skipped.
#' @param isPercentage (logical) if the unit of measurement in \code{data}
#' is percentage (e.g., growth rate) use \code{TRUE}. Long-run growth rate
#' is calculated by arithmetic mean for continuous case, and geometric mean
#' otherwise. If missing data exists, it returns NA.
#' @param ... additional arguments
#'
#' @return the growth rate (percentage)
#' @export
#'
#' @examples
#' y <- c(NA, 0, c(60, 70, 80, 90), 0, NA, NA)
#' g <- LongrunGrowth(y, 2, 3, skipZero = TRUE, isPercentage = TRUE, cont = TRUE)
#'
LongrunGrowth <- function(data, trimStart = 0, trimEnd = 0,
                          cont = FALSE, skipZero = TRUE, isPercentage = FALSE, ...) {
  data <- as.numeric(data)
  trimStart <- as.integer(trimStart)
  trimEnd <- as.integer(trimEnd)

  N <- length(data)
  I <- 1
  J <- N
  start <- data[[I]]
  end <- data[[J]]

  is.na.zero <- function(d, skipz) {
    return(is.na(d) || (skipz && d == 0))
  }

  if (is.na.zero(start, skipZero) && trimStart > 0) {
    for (i in c(1:(trimStart + 1))) {
      I <- i
      start <- data[[I]]
      if (is.na.zero(start, skipZero) == FALSE) {
        break
      }
    }
  }
  if (is.na.zero(start, skipZero)) {
    return(NA)
  }

  if (is.na.zero(end, skipZero) && trimEnd > 0) {
    for (i in c(1:trimEnd)) {
      J <- N - i
      end <- data[[J]]
      if (is.na.zero(end, skipZero) == FALSE) {
        break
      }
    }
  }
  if (is.na.zero(end, skipZero)) {
    return(NA)
  }
  if (J < I) {
    return(NA)
  }

  if (isPercentage) {
    data <- data[I:J]
    if (cont) {
      return(mean(data))
    } else {
      d <- 1
      for (i in data) {
        d <- d * (i / 100 + 1)
      }
      return((d^(1 / (J - I + 1)) - 1) * 100)
    }
  } else {
    if (cont) {
      return(log(end / start) / (J - I) * 100)
    } else {
      return(((end / start)^(1 / (J - I)) - 1) * 100)
    }
  }
}
