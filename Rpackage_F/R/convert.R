

#' Converts \code{Date-List} Data to Daily Data
#'
#' Use this to convert a time-series data with \code{Date-List} frequency to a time-series data with daily frequency.
#'
#' @param variable A variable with \code{Date-List} frequency.
#'
#' @details
#' It generally sorts the dates and fills any gaps between them with \code{NA}.
#'
#'
#' @return A variable with the daily frequency.
#' @export
#'
#' @examples
#' startFreq <- f.list.date(c("20220904","20220901"), "20220901")
#' v <- Variable(c(4,1), "var", startFreq, list())
#' w <- c.datelist.daily(v)
c.datelist.daily <- function(variable){

res <- .c_DateList_to_Daily(variable)
res

}


#' Converts Daily Data to Multi-Day Data
#'
#' Use this to convert a time-series data with daily frequency to a time-series data with multi-day frequency.
#'
#' @param variable A variable with daily frequency.
#' @param k Number of days in multi-day frequency. It should be positive.
#' @param aggregateFun A function that aggregates the data within each interval.
#' @param fromEnd If the number of observations is not divisible by k, this argument matters. If TRUE, the last observation is the combination of k observations. Otherwise, the last observation might be created from a lower number of observations.
#'
#' @details
#' When aggregating daily frequency over a period of \code{k} days,
#' we expect \code{k} numbers in each interval.
#' The aggregate function can be set to calculate the mean,
#' variance, etc. of the values within that interval.
#' Alternatively, you can select the first or last value.
#' Note that you may need to handle \code{NA} values.
#'
#'
#' @return A variable with the multi-day frequency.
#' @export
#'
#' @examples
#' startFreq <- f.daily(2022, 9, 1)
#' v <- Variable(c(1,2,3,4,5,6,7,8), "var", startFreq, list())
#' w <- c.daily.multidaily(v, 3, function(x)mean(x, na.rm=TRUE))
c.daily.multidaily <- function(variable, k, aggregateFun, fromEnd = TRUE){

  k = as.integer(k)
  if (k<=0)
    stop("Invalid 'k'. It must be positive.")
  aggregateFun = as.function(aggregateFun)
  fromEnd = as.logical(fromEnd)

  res = .c_Daily_to_MultiDaily(variable,k,aggregateFun,fromEnd)
  res

}
