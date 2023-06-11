

#' @useDynLib tdata
#' @importFrom Rcpp sourceCpp
NULL


#' Data for Vignette
#'
#' This is oil price data from 2010 retrieved by using the following code:
#' \code{oil_price <- Quandl::Quandl("OPEC/ORB", start_date="2010-01-01")}
#' It is saved due to the fact that CRAN checks may fail if the vignette relies on an external API call.
#'
#' @docType data
#' @format A \code{data.frame} with 2 columns: \code{Date} and \code{Value}
"oil_price"
