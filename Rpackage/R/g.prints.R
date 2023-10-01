
# Define the print method for 'options' class
rec.print.list <- function(x, indent = "") {
  for(i in seq_along(x)) {
    if (is.null(x[[i]])){
      cat(indent, names(x)[i], ": (NULL)\n", sep = "")
    }
    else if (is.list(x[[i]])) {
      cat(indent, names(x)[i], ":\n", sep = "")
      rec.print.list(x[[i]], paste0(indent, "    "))  # Recursive call with increased indent
    }
    else if (is.matrix(x[[i]])){
      mat = x[[i]]
      rows <- apply(mat, 1, function(row) {
        if (ncol(mat) > 4) {
          row_string <- paste(row[1:4], collapse = ", ")
          row_string <- paste0(row_string, ", ...; ")
        } else {
          row_string <- paste(row, collapse = ", ")
          row_string <- paste0(row_string, "; ")
        }
        return(row_string)
      })
      mat_str = paste(rows, collapse = "")

      if (nrow(mat) > 4) {
        mat_str = paste(mat_str, " ...")
      }

      cat(indent, names(x)[i], " (", nrow(x[[i]]) ,"x", ncol(x[[i]]), "): ", mat_str, "\n", sep = "")
    }
    else if (is.vector(x[[i]]) && length(x[[i]]) > 1) {
      if (length(x[[i]]) > 10)
        cat(indent, names(x)[i], " (", length(x[[i]]) ,"x1): ", paste(x[[i]][1:10], collapse = ", "), "...\n", sep = "")
      else
        cat(indent, names(x)[i], " (", length(x[[i]]) ,"x1): ", paste(x[[i]], collapse = ", "), "\n", sep = "")
    }
    else if (is.vector(x[[i]]) && length(x[[i]])  == 0) {
      cat(indent, names(x)[i], " (0x1): \n", sep = "")
    }
    else {
      cat(indent, names(x)[i], ": ", x[[i]], "\n", sep = "")
    }
  }
}


#' Prints an \code{ldt.list} object
#'
#' @param x An object of class \code{ldt.list}
#' @param ... Additional arguments
#'
#' @return This function has no output.
#' @export
print.ldt.list <- function(x, ...) {
  if (is.null(x))
    stop("argument is null.")
  if (!is.list(x))
    stop("Invalid class. A list is expected.")
  if (length(x) == 0)
    cat("An empty list")
  else{
    cat("A list with the following elements:\n\n")
    rec.print.list(x, " ")
  }
  cat("\n")
}


#' Prints an \code{ldt.search} object
#'
#' An \code{ldt.search} object is an output from one of the \code{search.?} functions (see \code{search.sur}, \code{search.varma}, or \code{search.bin}).
#'
#' @param x An object of class \code{ldt.search}
#' @param ... Additional arguments
#'
#' @return This function has no output.
#' @export
print.ldt.search <- function(x, ...) {
  if (is.null(x))
    stop("argument is null.")
  if (!is(x, "ldt.search"))
    stop("Invalid class. An 'ldt.search' object is expected.")

  is_summary = is(x, "ldt.search.summary")

  cat("LDT search result:\n")
  cat(" Method in the search process:", attr(x, "method"), "\n")
  cat(" Expected number of models: ", prettyNum(x$counts$expectedCount, big.mark = ","),
      ", searched: ", prettyNum(x$counts$searchedCount, big.mark = ","), " ",
      ", failed: ", prettyNum(x$counts$failedCount, big.mark = ","), " (",
      round(x$counts$failedCount / x$counts$searchedCount * 100,1), "%)\n",
      sep = ""
  )
  et <- difftime(x$info$endTime, x$info$startTime, units = "mins")
  cat(" Elapsed time:", et, "minutes \n")
  cat(" Length of results:", length(x$results), "\n")
  cat("--------\n")
  if (x$counts$failedCount > 0){
    cat(" Failures:\n")
    i <- 0
    for (fd in x$counts$failedDetails){
      i <- i + 1
      cat("   ", i,". ",fd$message, ": ", fd$count
          ," (", round((fd$count/x$counts$searchedCount)*100,1) ,"%)\n", sep = "")
    }
    cat("--------\n")
  }

  evalNames <- unique(sapply(x$results, function(y)y$evalName))
  targNames <- unique(sapply(x$results, function(y)y$targetName))

  for (tName in targNames){
    cat(" Target (", tName, "):\n" ,sep = "")
    for (eName in evalNames){

      values <- x$results[which(sapply(x$results, function(y)
        {y$targetName == tName && y$evalName == eName}))]

      cat("   Evaluation (", eName , "):\n" ,sep = "")
      indent = "     "
      if (x$info$items$model){
        if (x$info$items$bestK > 0){

          best <- values[which(sapply(values, function(y){y$typeName == "best model" && y$info == 0}))][[1]]
          cat(indent, "Best model:\n")
          if (is_summary){
            cat("\n+++++++++++++++ MODEL SUMMARY +++++++++++++++++\n")
            print(best$value)
            cat("\n+++++++++++++++      END      +++++++++++++++++\n")
          }
          else{
            rec.print.list(list(
              endogenous = colnames(x$info$data$data)[best$value$depIndices+1],
              exogenous = if (is.null(best$value$exoIndices)) {NULL} else {colnames(x$info$data$data)[best$value$exoIndices+1]},
              metric = best$value$metric
            ), indent = c(indent, "  "))
          }

        }

        if (x$info$items$inclusion){

          inclusion <- values[which(sapply(values, function(y){y$typeName == "best model" && y$info == 0}))][[1]]
          cat(indent, "Inclusion weights:\n")
          rec.print.list(list(
            maximum = "TODO"
          ), indent = c(indent, "  "))
        }

      }

      if (x$info$items$type1){

        if (x$info$items$bestK){
          #TODO: print for best coefficients
        }

        if (x$info$items$extremeMultiplier){
          #TODO:
        }

        if (x$info$items$cdfs){
          #TODO:
        }

        if (x$info$items$mixture4){
          #TODO:
        }

      }

    }


    cat("--------\n")
  }

  if (x$info$items$bestK > 1)
    cat(" ** results for ", x$info$items$bestK, " best model(s) are saved\n", sep="")
  if (x$info$items$all)
    cat(" ** results for all estimated models are saved\n")

}


