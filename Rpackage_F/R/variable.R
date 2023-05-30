
#' Creates a Variable
#'
#' Use this to create data array with a frequency. It can have a name and other named fields.
#'
#' @param data Data of the variable.
#' @param name Name of the variable.
#' @param startFrequency Frequency of the first element.
#' @param fields A list that contains named fields.
#'
#' @return An object of class \code{ldtv}. This is also a list with the following memebers:
#' \itemize{
#' \item **data:** Determines the \code{data}.
#' \item **name:** Determines the \code{name}.
#' \item **startFrequency:** Determines the \code{startFrequency}.
#' \item **fields:** Determines the \code{fields}.
#' }
#' @export
#'
#' @examples
#' data <- c(1,2,3,2,3,4,5)
#' start_f <- f.monthly(2022,12)
#' fields <- list(c("key1","value1"), c("key2", "value2"))
#' v1 = Variable(data,"V1",start_f, fields)
Variable <- function(data, name, startFrequency,
         fields) {
  data = as.numeric(data)
  name = as.character(name)

  res <- .Variable(data, name, startFrequency, fields)
  res
}

#' Converts a Variable to String
#'
#' Use this to convert the variable in a compact form.
#'
#' @param w An object of class \code{ldtv}.
#' @param ... Additional arguments.
#'
#' @details
#' The returned character will have just one line, items are separated by \code{tab} or \code{semi-colon}.
#'
#' @return A character that represents the variable.
#' @export
#'
#' @examples
#' # define the variable:
#' data <- c(1,2,3,2,3,4,5)
#' start_f <- f.monthly(2022,12)
#' fields <- list(c("key1","value1"), c("key2", "value2"))
#' v1 = Variable(data,"V1",start_f, fields)
#'
#' #string representation:
#' v1_str <- as.character(v1)
#'
as.character.ldtv <- function(w, ...) {
  res <- .VariableToString(w)
  res
}

#' Prints a Variable
#'
#' Use this to print an \code{ldtv} object.
#'
#' @param x An \code{ldtv} object
#' @param ... additional arguments
#'
#' @return \code{NULL}
#' @export
print.ldtv <- function(x, ...) {
  if (is.null(x)) {
    stop("argument is null.")
  }
  if (any(class(x) == "ldtv") == FALSE) {
    stop("Invalid class")
  }
  if (any(class(x$startFrequency) == "ldtf") == FALSE) {
    stop("Invalid frequency class")
  }

  n <- length(x$fields)
  if (n > 0) {
    fields <- paste0(
      unlist(sapply(
        c(1:n),
        function(i) {
          paste0(
            x$fields[[i]][1], " = ",
            x$fields[[i]][2]
          )
        }
      )),
      sep = "\n", strrep(" ", 8), collapse = ""
    )
  } else {
    fields <- NULL
  }
  s <- .ToString_F0(x$startFrequency)
  cat("Variable:\n",
      "    Name = ", x$name, "\n",
      "    Length = ", length(x$data), "\n",
      "    Frequency Class = ", s$classType, ": ", s$class, "\n",
      "    Start Frequency = ", s$value, "\n",
      "    Fields:\n", strrep(" ", 8), if (is.null(fields)) "" else fields,
      sep = ""
  )
  return(NULL)
}


#' Converts Variable to \code{data.frame}
#'
#' Use this to convert an \code{ldtv} object to a \code{data.frame}. You can use the result for plotting.
#'
#' @param x An \code{ldtv} object
#' @param ... additional arguments
#'
#' @return A \code{data.frame} in which row names are set from the frequency of the variable.
#' @export
#' @examples
#' # define the variable:
#' data <- c(1,2,3,2,3,4,5)
#' start_f <- f.monthly(2022,12)
#' fields <- list(c("key1","value1"), c("key2", "value2"))
#' v1 = Variable(data,"V1",start_f, fields)
#'
#' # convert it to data.frame
#' df1 <- to.data.frame(v1)
#'
to.data.frame <- function(x, ...) {
  if (is.null(x)) {
    stop("argument is null.")
  }
  if (any(class(x) == "ldtv") == FALSE) {
    stop("Invalid class.")
  }
  if (is.null(x$data)) {
    stop("Variable's data array is null.")
  }
  if (is.null(x$startFrequency)) {
    stop("Variable's frequency is null.")
  }

  df <- as.data.frame(x$data)
  colnames(df) <- if (is.null(x$name)) "V" else x$name

  freqs <- get.seq0(x$startFrequency, length(x$data))
  rownames(df) <- freqs

  df
}


#' Binds Variables
#'
#' Use this to bind variables with the same class of frequency together.
#'
#' @param var A list of variables ((i.e., \code{ldtv} objects)) with similar frequency class.
#' @param interpolate If \code{TRUE}, missing observations are interpolated.
#' @param adjustLeadLags If \code{TRUE}, leads and lags are adjusted with respect to the first variable.
#' @param numEndo (integer) If \code{adjustLeadLags} is \code{TRUE}, this must be
#' the number of endogenous variables. The rest is exogenous.
#' @param horizon (integer) If \code{adjustLeadLags} is \code{TRUE} and there is exogenous variables,
#' this determines the required length of out-of-sample data. It creates lag of exogenous variables
#' or omits \code{NaN}s to make data available.
#'
#' @return A list with the following members:
#' \itemize{
#' \item **data:** (numeric matrix) Final data after the requested fixes. It is a matrix with variables in the columns and frequencies as the row names.
#' \item **info:** (integer matrix) Information about the columns of the final data. E.g., Range of data, missing data, lags/leads, etc.
#' \item **startFrequency:** Start frequency of data.
#' \item **startClass:** frequency class of the data.
#' }
#'
#' @export
#' @examples
#' v1 = Variable(c(1,2,3,2,3,4,5),"V1",f.monthly(2022,12), list())
#' v2 = Variable(c(10,20,30,20,30,40,50),"V2",f.monthly(2022,8), list())
#' L = BindVariables(list(v1,v2))
BindVariables <- function(varList, interpolate = FALSE,
              adjustLeadLags = FALSE, numExo = 0,
              horizon = 0)
{
  adjustLeadLags = as.logical(adjustLeadLags)
  interpolate = as.logical(interpolate)
  numExo = as.integer(numExo)
  horizon = as.integer(horizon)
  if (numExo < 0)
    stop("Invalid 'numExo'. It cannot be negative.")
  if (horizon < 0)
    stop("Invalid 'horizon'. It cannot be negative.")


  res <- .BindVariables(varList, interpolate, adjustLeadLags, numExo, horizon)
  res
}