

#' Convert Data to Daily Frequency
#'
#' Use this function to convert a time-series data (currently implemented: \code{Date-List}, \code{Daily-In-Week}) to a time-series data with daily frequency.
#'
#' @param variable A variable.
#' @param aggregateFun Function to aggregate data within each interval (see details).
#'
#' @details
#' In some cases, conversion sorts the dates and fills any gaps between them with \code{NA}.
#' However, in other cases, conversion requires aggregation. For example,
#' when aggregating hourly data over a period of \code{k} hours to generate daily data,
#' we expect \code{k} numbers in each interval.
#' The aggregate function can be set to calculate the mean,
#' variance, median, etc., or any function that takes
#' the vector of \code{k} values and returns a number.
#'
#' @return A variable with daily frequency, with data sorted from the original variable and missing dates filled with \code{NA}.
#' @export
#'
#' @examples
#' startFreq <- f.list.date(c("20220904","20220901"), "20220901")
#' v <- variable(c(4,1), startFreq)
#' w <- convert.to.daily(v)
convert.to.daily <- function(variable, aggregateFun = NULL){

res <- .ConvertTo_Daily(variable, aggregateFun)
res

}


#' Convert Data to Multi-Day Frequency
#'
#' Use this function to convert a time-series data (currently implemented: daily) to a time-series data with multi-day frequency.
#'
#' @param variable A variable.
#' @param k Number of days in multi-day frequency, must be positive.
#' @param aggregateFun Function to aggregate data within each interval.
#' @param fromEnd If the number of observations is not divisible by \code{k}, this argument matters. If \code{TRUE}, the last observation is the combination of \code{k} observations. Otherwise, the last observation may be created from fewer observations.
#'
#' @details
#' See the details section of the \code{\link{convert.to.daily}} function.
#'
#'
#' @return A variable with multi-day frequency.
#' @export
#'
#' @examples
#' startFreq <- f.daily(c(2022, 9, 1))
#' v <- variable(c(1,2,3,4,5,6,7,8), startFreq)
#' w <- convert.to.multidaily(v, 3, function(x)mean(x, na.rm=TRUE))
convert.to.multidaily <- function(variable, k, aggregateFun, fromEnd = TRUE){

  k = as.integer(k)
  if (k<=0)
    stop("Invalid 'k'. It must be positive.")

  fromEnd = as.logical(fromEnd)

  if (is.function(aggregateFun) == FALSE && is.character(aggregateFun) == FALSE)
    stop("Invalid 'aggregateFun'. It must be an R function or a string that represents a valid aggregator.");

  res = .ConvertTo_MultiDaily(variable,k,aggregateFun,fromEnd)
  res
}


#' Convert Data to Weekly Frequency
#'
#' Use this function to convert time-series data (currently implemented: daily) to time-series data with weekly frequency.
#'
#' @param variable A variable.
#' @param weekStart Determines the start day of the week, can be \code{sun}, \code{mon}, \code{tue}, \code{wed}, \code{thu}, \code{fri}, or \code{sat}.
#' @param aggregateFun Function to aggregate data within each interval.
#'
#' @details
#' See the details section of the \code{\link{convert.to.daily}} function.
#'
#' @return A variable with weekly frequency.
#' @export
#'
#' @examples
#' startFreq <- f.daily(c(2022, 9, 1))
#' v <- variable(c(1,2,3,4,5,6,7,8), startFreq)
#' w <- convert.to.weekly(v, "mon", function(x)mean(x, na.rm=TRUE))
#'
convert.to.weekly <- function(variable, weekStart, aggregateFun){
  weekStart = as.character(weekStart)
  if (is.function(aggregateFun) == FALSE && is.character(aggregateFun) == FALSE)
    stop("Invalid 'aggregateFun'. It must be an R function or a string that represents a valid aggregator.");

  res = .ConvertTo_Weekly(variable,weekStart,aggregateFun)
  res
}


#' Convert Data to Year-Based Frequency
#'
#' Use this function to convert time-series data (currently implemented: daily) to time-series data with year-based frequency such as monthly, quarterly, yearly, etc.
#'
#' @param variable A variable.
#' @param x Determines the number of partitions in each year, for example, use 12 for monthly data.
#' @param aggregateFun Function to aggregate data within each interval.
#'
#'
#' @details
#' See the details section of the \code{\link{convert.to.daily}} function.
#'
#' @return A variable with year-based frequency.
#' @export
#'
#' @examples
#' startFreq <- f.daily(c(2023,1,1))
#' v <- variable(c(1:(365*2)), startFreq)
#' w <- convert.to.XxYear(v,12,function(x)mean(x))
#'
convert.to.XxYear <- function(variable, x, aggregateFun){
  x = as.integer(x)
  if (x < 0)
    stop("Invalid 'x'. It should be a positive integer.")
  if (is.function(aggregateFun) == FALSE && is.character(aggregateFun) == FALSE)
    stop("Invalid 'aggregateFun'. It must be an R function or a string that represents a valid aggregator.");

  res = .ConvertTo_XxYear(variable,aggregateFun,x)
  res
}
