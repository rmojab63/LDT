

#' Scenarios for Removing \code{NA}s
#'
#' Use this to remove \code{NA} from a matrix. You can optimize the size of information by this method.
#'
#' @param data A matrix that contains \code{NA}.
#' @param countFun A function to determine how strategies are
#' sorted. Default function counts the number of observations.
#' You might want to give columns a higher level of importance
#' for example by using \code{nRows*nCols^1.5}.
#' @param rowIndices Indices of the sorted rows to search.
#' Use it to create jumps for large number of rows (E.g.,
#' if the first sorted strategies suggest small number of
#' columns and you are looking for other strategies).
#' Use \code{NULL} to disable it.
#' @param colIndices Similar to \code{rowIndices} for columns.
#' @param printMsg If \code{TRUE}, it prints the progress.
#'
#' @details
#' When a matrix has \code{NA}, one can omit columns
#' with \code{NA}, rows with \code{NA}, or a combination of
#' these two. The total number of observations is a function
#' of the order. This function tries all combinations and returns the results.
#'
#'
#' @return A list of lists, each with the following elements:
#' \tabular{ll}{
#'    \code{nRows} \tab number of rows in the matrix.\cr
#'   \code{nCols} \tab number of cols in the matrix.\cr
#'    \code{colFirst} \tab whether to remove columns or rows first.\cr
#'    \code{colRemove} \tab indexes of the columns to be removed.\cr
#'    \code{rowRemove} \tab indexes of the rows to be removed.
#'   }
#'
#' @export
#'
#' @examples
#' data <- matrix(c(NA, 2, 3, 4, NA, 5, NA, 6, 7, NA, 9, 10, 11, 12, 13, 14, 15, NA, 16, 17), 4, 5)
#' res <- remove.na.strategies(data)
remove.na.strategies <- function(data, countFun = function(nRows, nCols) nRows * nCols,
                                 rowIndices = NULL, colIndices = NULL, printMsg = FALSE) {
  data <- as.matrix(data)
  Nrow = nrow(data)
  Ncol = ncol(data)
  data_na <- is.na(data)
  c_cols <- sort(colSums(data_na), index.return = TRUE, decreasing = TRUE)


  rowIndices = as.integer(rowIndices)
  colIndices = as.integer(colIndices)

  if (length(colIndices) == 0)
    colIndices = c(1:length(c_cols$ix))
  cmaxCount = length(colIndices)
  if (printMsg)
    cat("Searching columns: \n")
  result <- list()
  j <- 0
  for (i in colIndices) { # remove i columns with most NAs
    j = j + 1
    colRemove <- c(c_cols$ix[1:i])
    if (length(colRemove) > 0) {
      mat <- data[, -colRemove, drop = FALSE]
    }
    rowRemove <- which(rowSums(is.na(mat)) > 0)
    if (length(rowRemove) == Nrow){
      if (printMsg)
        cat("    breaking the loop. All rows are removed after removing the columns\n")
      break # all rows will be removed
    }
    result[[length(result) + 1]] <- list(
      nRows = Nrow - length(rowRemove),
      nCols = Ncol - length(colRemove), colFirst = TRUE,
      colRemove = sort(colRemove), rowRemove = as.integer(sort(rowRemove))
    )
    if (printMsg)
      cat(paste0("    ", j, "/", cmaxCount,
                 ", index=", i,
                 ", nrow=", Nrow - length(rowRemove),
                 ", ncol=",Ncol - length(colRemove), "\n"))
  }

  c_rows <- sort(rowSums(data_na), index.return = TRUE, decreasing = TRUE)

  if (length(rowIndices) == 0)
    rowIndices = c(1:length(c_rows$ix))
  rmaxCount = length(rowIndices)

  if (printMsg)
    cat("Searching rows: \n")
  j <- 0
  for (i in rowIndices) { # remove i row with most NAs
    j <- j + 1
    rowRemove <- c(c_rows$ix[1:i])
    if (length(rowRemove) > 0) {
      mat <- data[-rowRemove, , drop = FALSE]
    }
    colRemove <- which(colSums(is.na(mat)) > 0)
    if (length(colRemove) == Ncol){
      if (printMsg)
        cat("    breaking the loop. All columns are removed after removing the rows\n")
      break
    }
    result[[length(result) + 1]] <- list(
      nRows = Nrow - length(rowRemove),
      nCols = Ncol - length(colRemove), colFirst = FALSE,
      colRemove = as.integer(sort(colRemove)), rowRemove = sort(rowRemove)
    )

    if (printMsg)
      cat(paste0("    ", j, "/", rmaxCount,
                 ", index=", i,
                 ", nrow=", Nrow - length(rowRemove),
                 ", ncol=",Ncol - length(colRemove), "\n"))
  }

  inx <- sort(sapply(result, function(r) countFun(r$nRows, r$nCols)),
              index.return = TRUE, decreasing = TRUE
  )

  result <- result[inx$ix]
  result
}





#' Calculate Long-run Growth
#'
#' Use this to calculate long-run growth of a time-series data.
#'
#' @param data (numeric vector) Determines the data of the series.
#' @param continuous (logical) if \code{TRUE}, it will use the continuous formula.
#' @param isPercentage (logical) If the unit of measurement in \code{data} is percentage (e.g., growth rate), use \code{TRUE}. The long-run growth rate is calculated by the arithmetic mean for the continuous case and the geometric mean otherwise. If missing data exists, it returns \code{NA}.
#' @param trimStart If the number of leading \code{NA}s is larger than this number, the function returns NA. Otherwise, it finds the first non-NA value and continues the calculations.
#' @param trimEnd Similar to \code{trimStart}, but for the end of the series.
#' @param skipZero If \code{TRUE}, leading and trailing zeros are skipped, similar to \code{NA}.
#'
#' @details
#' A variable can have continuous growth (\eqn{y(t)=y(0) (1+g_1)(1+g_2)\ldots (1+g_t)})
#' or discrete growth (\eqn{y(t)=y(0)e^{g_1}e^{g_2}\ldots e^{g_t}}) over \eqn{t} periods.
#' \eqn{y(0)} is the first value and \eqn{y(n)} is the last value.
#' By long-run growth rate, we mean a number such as \eqn{g} such that if
#' we start from \eqn{y(0)} and the variable growth is \eqn{g} every period,
#' we reach \eqn{y(t)} after t periods. This number summarizes all \eqn{g_i}s,
#' however, it is not generally the average of these rates.
#'
#' @return Long-run growth rate (percentage)
#' @export
#'
#' @examples
#' y <- c(60, 70, 80, 95)
#' g <- get.longrun.growth(y, isPercentage = TRUE, continuous = FALSE)
#' # Note that 'g' is different from 'mean(y)'.
#'
get.longrun.growth <- function(data, continuous = FALSE, isPercentage = FALSE,
                               trimStart = 0, trimEnd = 0, skipZero = TRUE) {
  data <- as.numeric(data)
  continuous = as.logical(continuous)
  isPercentage = as.logical(isPercentage)

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
    if (continuous) {
      return(mean(data))
    } else {
      d <- 1
      for (i in data) {
        d <- d * (i / 100 + 1)
      }
      return((d^(1 / (J - I + 1)) - 1) * 100)
    }
  } else {
    if (continuous) {
      return(log(end / start) / (J - I) * 100)
    } else {
      return(((end / start)^(1 / (J - I)) - 1) * 100)
    }
  }
}

#' Get Descriptive Statistics
#'
#' Use this function to get descriptive statistics of a numeric vector.
#'
#' @param data A numeric vector that contains the data.
#' @param type Its a character array that determines the type of the descriptive statistics.
#' @param skipNAN If \code{TRUE}, it checks for and skips any \code{NA}s. If you are sure that data does not have \code{NA} set it to be \code{False}.
#'
#' @return An array with different sizes based on the inputs
#'
get.descriptive <- function(data, type, skipNAN){

  # TODO: test

  type = as.character(type)
  data = as.numeric(data)
  skipNAN = as.logical(skipNAN)

  res = .Get_Descriptive(data, type, skipNAN)
  res
}
