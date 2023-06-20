
#' Create a Variable
#'
#' Use this function to create a variable, which is a data array with frequencies. It can have a name and other named fields.
#'
#' @param data The data of the variable.
#' @param startFrequency The frequency of the first element.
#' @param name The name of the variable.
#' @param fields A list that contains named fields.
#'
#' @return An object of class \code{ldtv}, which is also a list with the following members:
#' \itemize{
#' \item **data**: Determines the \code{data}.
#' \item **name**: Determines the \code{name}.
#' \item **startFrequency**: Determines the \code{startFrequency}.
#' \item **fields**: Determines the \code{fields}.
#' }
#' @export
#'
#' @examples
#' data <- c(1,2,3,2,3,4,5)
#' start_f <- f.monthly(2022,12)
#' fields <- list(c("key1","value1"), c("key2", "value2"))
#' v1 = variable(data, start_f, "V1", fields)
variable <- function(data, startFrequency = NULL, name = NULL, fields = NULL) {
  data = as.numeric(data)
  if (is.null(name) == FALSE)
    name = as.character(name)
  if (is.null(startFrequency))
    startFrequency <- f.cross.section(1)

  res <- .Variable(data, name, startFrequency, fields)
  res
}

#' Convert a Variable to Character String
#'
#' Use this function to convert a variable to a compact form.
#'
#' @param x An object of class \code{ldtv}.
#' @param ... Additional arguments.
#'
#' @details
#' The returned character will have just one line, with items separated by \code{tab} or \code{semi-colon}.
#'
#' @return A character that represents the variable.
#' @export
#'
#' @examples
#' # define the variable:
#' data <- c(1,2,3,2,3,4,5)
#' start_f <- f.monthly(2022,12)
#' fields <- list(c("key1","value1"), c("key2", "value2"))
#' v1 = variable(data,start_f, "V1", fields)
#'
#' #string representation:
#' v1_str <- as.character(v1)
#'
as.character.ldtv <- function(x, ...) {
  res <- .VariableToString(x)
  res
}

#' Print a Variable
#'
#' Use this to print a variable.
#'
#' @param x A variable which is an object of class \code{ldtv}.
#' @param ... Additional arguments
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
      "    Name = ", if (is.null(x$name)) "" else x$name, "\n",
      "    Length = ", length(x$data), "\n",
      "    Frequency Class = ", s$classType, ": ", s$class, "\n",
      "    Start Frequency = ", s$value, "\n",
      "    Fields:", if (is.null(fields)) " NULL" else paste0("\n", strrep(" ", 8), fields),
      sep = ""
  )
  return(NULL)
}


#' Convert Variable to Data Frame
#'
#' Use this function to convert a variable to a data frame. You can use the result for plotting.
#'
#' @param x An \code{ldtv} object.
#' @param ... Additional arguments.
#'
#' @return A data frame in which row names are set from the frequency of the variable.
#' @export
#' @examples
#' # Define the variable:
#' data <- c(1,2,3,2,3,4,5)
#' start_f <- f.monthly(2022,12)
#' fields <- list(c("key1","value1"), c("key2", "value2"))
#' v1 = variable(data,start_f,"V1", fields)
#'
#' # convert it to data.frame
#' df1 <- as.data.frame(v1)
#'
as.data.frame.ldtv <- function(x, ...) {
  if (is.null(x)) {
    stop("argument is null.")
  }
  if (any(class(x) == "ldtv") == FALSE) {
    stop("Invalid class.")
  }
  if (is.null(x$data)) {
    stop("variable's data array is null.")
  }
  if (is.null(x$startFrequency)) {
    stop("variable's frequency is null.")
  }

  df <- as.data.frame(x$data)
  colnames(df) <- if (is.null(x$name)) "V" else x$name

  freqs <- get.seq0(x$startFrequency, length(x$data))
  rownames(df) <- freqs

  df
}


#' Bind Variables and Create a Data.frame
#'
#' Use this function to bind variables with the same class of frequency together.
#'
#' @param varList A list of variables (i.e., \code{ldtv} objects) with similar frequency class.
#' @param interpolate If \code{TRUE}, missing observations are interpolated.
#' @param adjustLeadLags If \code{TRUE}, leads and lags are adjusted with respect to the first variable.
#' @param numExo An integer representing the number of exogenous variables.
#' @param horizon An integer representing the required length of out-of-sample data if \code{adjustLeadLags} is \code{TRUE} and there are exogenous variables.
#' It creates lags of exogenous variables or omits \code{NaN}s to make data available.
#'
#' @return A list with the following members:
#' \item{data}{A numeric matrix representing the final data after the requested fixes. It is a matrix with variables in the columns and frequencies as the row names.}
#' \item{info}{An integer matrix containing information about the columns of the final data, such as range of data, missing data, lags/leads, etc.}
#' \item{startFrequency}{The start frequency of the data.}
#' \item{startClass}{The frequency class of the data.}
#'
#' @export
#' @examples
#' v1 = variable(c(1,2,3,2,3,4,5),f.monthly(2022,12),"V1")
#' v2 = variable(c(10,20,30,20,30,40,50),f.monthly(2022,8),"V2")
#' L = bind.variables(list(v1,v2))
bind.variables <- function(varList, interpolate = FALSE,
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
  varList = as.list(varList)

  res <- .BindVariables(varList, interpolate, adjustLeadLags, numExo, horizon)
  res
}
