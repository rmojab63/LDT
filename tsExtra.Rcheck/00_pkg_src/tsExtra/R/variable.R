
#' Creates a Variable
#'
#' @param data Data of the variable
#' @param name Name of the variable
#' @param startFrequency Frequency of the first data-point. It is an \code{ldtf} object. See \code{F_?} functions.
#' @param fields Named of any other fields
#'
#' @return An object of class \code{ldtv}.
#' @export
#' @examples
#' v1 = ldt::Variable(c(1,2,3,2,3,4,5),"V1",F_Monthly(2022,12),
#'      list(c("key1","value1"), c("key2", "value2")))
Variable <- function(data, name, startFrequency,
         fields) {
  res <- .Variable(data, name, startFrequency, fields)
  res
}

#' Converts a Variable to String
#'
#' @param w The variable
#'
#' @return String representation of the variable in compact form
#' @export
VariableToString <- function(w) {
  res <- .VariableToString(w)
  res
}



#' Binds a of Variables
#'
#' @param varA of variables ((i.e., \code{ldtv} objects)) with similar frequency class
#' @param interpolate If \code{TRUE} missing observations are interpolated.
#' @param adjustLeadLags If \code{TRUE} leads and lags are adjusted with respect to the first variable.
#' @param numEndo (integer) If \code{adjustLeadLags} is \code{TRUE}, this must be
#' the number of endogenous variables. The rest is exogenous.
#' @param horizon (integer) If \code{adjustLeadLags} is \code{TRUE} and there is exogenous variables,
#' this determines the required length of out-of-sample data. It creates lag of exogenous variables
#' or omits \code{NaN}s to make data available.
#'
#' @return (list) results
#' \item{data}{(numeric matrix) Final data after the requested fixes. It is a matrix with variables in the columns and frequencies as the row names.}
#' \item{info}{(integer matrix) Information about the columns of the final data. E.g., Range of data, missing data, lags/leads, etc.}
#' \item{startFrequency}{Start frequency of data}
#' \item{startClass}{frequency class of the data}
#'
#' @export
#' @examples
#' v1 = ldt::Variable(c(1,2,3,2,3,4,5),"V1",F_Monthly(2022,12), list())
#' v2 = ldt::Variable(c(10,20,30,20,30,40,50),"V2",F_Monthly(2022,8), list())
#' L = ldt::BindVariables(list(v1,v2))
BindVariables <- function(varList, interpolate = FALSE,
              adjustLeadLags = FALSE, numExo = 0,
              horizon = 0)
{
  res <- .BindVariables(varList, interpolate, adjustLeadLags, numExo, horizon)
  res
}
