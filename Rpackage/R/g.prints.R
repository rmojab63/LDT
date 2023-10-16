
paste00 <- function(x, collapse = ", ") {
  paste(sprintf(paste0("%.", options()$digits, "g"), x), collapse = collapse)
}

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
          row_string <- paste00(row[1:4], collapse = ", ")
          row_string <- paste0(row_string, ", ...; ")
        } else {
          row_string <- paste00(row, collapse = ", ")
          row_string <- paste0(row_string, "; ")
        }
        return(row_string)
      })
      mat_str = paste(rows, collapse = "")

      if (nrow(mat) > 4) {
        mat_str = paste(mat_str, " ...")
      }

      cat(indent, names(x)[i], ": (", nrow(x[[i]]) ,"x", ncol(x[[i]]), ") ", mat_str, "\n", sep = "")
    }
    else if (is.character(x[[i]])) {
      if (length(x[[i]]) == 1)
        cat(indent, names(x)[i], ": ", x[[i]], "\n", sep = "")
      else if (length(x[[i]]) > 10)
        cat(indent, names(x)[i], ": (", length(x[[i]]) ,"x1) ", paste0(x[[i]][1:10], collapse = ", "), "...\n", sep = "")
      else
        cat(indent, names(x)[i], ": (", length(x[[i]]) ,"x1) ", paste0(x[[i]], collapse = ", "), "\n", sep = "")
    }
    else if (is.vector(x[[i]]) && length(x[[i]]) > 1) {
      if (length(x[[i]]) > 10)
        cat(indent, names(x)[i], ": (", length(x[[i]]) ,"x1) ", paste00(x[[i]][1:10], collapse = ", "), "...\n", sep = "")
      else
        cat(indent, names(x)[i], ": (", length(x[[i]]) ,"x1) ", paste00(x[[i]], collapse = ", "), "\n", sep = "")
    }
    else if (is.vector(x[[i]]) && length(x[[i]])  == 0) {
      cat(indent, names(x)[i], ": (0x1) \n", sep = "")
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
#' Prints the main results in an \code{ldt.search} object.
#' This includes information about the first best models and significant coefficients.
#'
#'
#' @param x An object of class \code{ldt.search}
#' @param ... Additional arguments
#'
#' @details
#' An \code{ldt.search} object is an output from one of the \code{search.?} functions (see \code{search.sur}, \code{search.varma}, or \code{search.bin}).
#'
#'
#' @return This function has no output.
#' @export
print.ldt.search <- function(x, ...) {
  if (is.null(x))
    stop("argument is null.")
  if (!is(x, "ldt.search"))
    stop("Invalid class. An 'ldt.search' object is expected.")

  is_summary = is(x, "ldt.search.summary")
  method <- tolower(attr(x, "method"))
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
  if (x$counts$failedCount == x$counts$searchedCount){
    cat("All Failed!\n")
  }
  else{

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
                endogenous = best$value$endogenous,
                exogenous = best$value$exogenous,
                metric = best$value$metric
              ), indent = c(indent, "  "))
            }

          }

          if (x$info$items$inclusion){
            inclusion <- values[which(sapply(values, function(y){y$typeName == "inclusion"}))][[1]]
            cat(indent, "Inclusion weights average:\n")
            mi <- which.max(inclusion$value[,1])
            rec.print.list(list(
              "maximum value" = max(inclusion$value[,1]),
              "name" = rownames(inclusion$value)[mi],
              "count" = inclusion$value[,2][mi]
            ), indent = c(indent, "  "))
          }

        }

        if (x$info$items$type1){

          if (x$info$items$bestK){
            best_coefs <- values[which(sapply(values, function(y){y$info == 0 && startsWith(y$typeName, "best item for")}))]
            cat(indent, "Best significant coeffiecients [mean-1.95*std, mean, mean+1.95*std]:\n")
            lst <- lapply(best_coefs, function(b)c(b$value$mean - 1.95 * sqrt(b$value$var),
                                                   b$value$mean,
                                                   b$value$mean + 1.95 * sqrt(b$value$var)))
            names(lst) <- sapply(best_coefs, function(b)strsplit(sub(".*'([^']*)'.*", "\\1", b$typeName), split = "'")[[1]])
            lst <- lst[which(sapply(lst,function(b)b[1]>0 || b[3]<0))]
            rec.print.list(lst, indent = c(indent, "  "))
          }

          if (x$info$items$extremeMultiplier > 0){
            extremeB <- values[which(sapply(values, function(y){startsWith(y$typeName, "extreme")}))]
            cat(indent, "Extreme bounds (significant):\n")
            sig <- extremeB$value[which(extremeB$value[,1]>0 | extremeB$value[,2]<0),,drop=FALSE]
            if (!is.null(sig)  && nrow(sig) > 0){
              lst <- lapply(1:nrow(sig), function(e)c(sig[e,1],sig[e,2]))
              names(lst) <- rownames(sig)
              rec.print.list(lst, indent = c(indent, "  "))
            }
            else
              cat(c(indent, "  "), "(none)\n")
          }

          if (length(x$info$items$cdfs) > 0){
            cdfs <- values[which(sapply(values, function(y){ y$info == 0 && startsWith(y$typeName, "cdf")}))]
            cat(indent, "CDF (significant):\n")
            if (length(cdfs) == 0)
              cat(c(indent, "  "), "CDF at 0 is not available.")
            else{
              cdfs <- cdfs[[1]]
              sig <- cdfs$value[which(cdfs$value[,1]<0.1 | cdfs$value[,1]>0.9),,drop=FALSE]
              if (!is.null(sig) && nrow(sig) > 0){
                lst <- lapply(1:nrow(sig), function(e)c(sig[e,1]))
                names(lst) <- rownames(sig)
                rec.print.list(lst, indent = c(indent, "  "))
              }
              else
                cat(c(indent, "  "), "(none)\n")
            }
          }

          if (x$info$items$mixture4){
            mixture <- values[which(sapply(values, function(y){startsWith(y$typeName, "mixture")}))][[1]]
            cat(indent, "Mixture significant [mean-1.95*std, mean, mean+1.95*std]:\n")
            lst <- lapply(1:nrow(mixture$value), function(i)c(mixture$value[i,1] - 1.95 * sqrt(mixture$value[i,2]),
                                                              mixture$value[i,1],
                                                              mixture$value[i,1] + 1.95 * sqrt(mixture$value[i,2])))
            names(lst) <- rownames(mixture$value)
            lst <- lst[which(sapply(lst,function(b)b[1]>0 || b[3]<0))]
            rec.print.list(lst, indent = c(indent, "  "))
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
}

#' Prints an \code{ldt.estim} object
#'
#' Prints the main results in an \code{ldt.estim} object.
#'
#' @param x An object of class \code{ldt.estim}
#' @param ... Additional arguments
#'
#' @details
#' An \code{ldt.search} object is an output from one of the \code{search.?} functions (see \code{search.sur}, \code{search.varma}, or \code{search.bin}).
#'
#'
#' @return This function has no output.
#' @export
print.ldt.estim <- function(x, ...) {

  if (is.null(x))
    stop("argument is null.")
  if (!is(x, "ldt.estim"))
    stop("Invalid class. An 'ldt.estim' object is expected.")

  method <- tolower(attr(x, "method"))
  cat("LDT '", attr(x, "method") ,"' estimation result\n", sep = "")
  cat("Model: ", do.call(paste0("estim.", method, ".model.string"), list(obj = x)), "\n")

  dist <- "t"
  if (method == "binary")
    dist <- "z"

  coefs <- x$estimations$coefs
  stds <- x$estimations$std
  tzstats <- x$estimations[[paste0(dist, "stats")]]
  pValues <- x$estimations$pValues

  max_chars <- max(nchar(rownames(coefs)))
  max_widths <- sapply(list(coefs, stds, tzstats, pValues), function(mat) max(nchar(format(mat))))

  dep_vars <- colnames(coefs)

  j <- 0
  for (dep_var in dep_vars) {
    j <- j + 1
    cat(paste0("\nDependent variable: ", dep_var, "\n"))
    cat("----------------------------------------------------\n")
    if (!is.null(x$estimations$resid)){
      cat("Residuals:\n")
      residuals_summary <- summary(x$estimations$resid[,dep_var])
      print(residuals_summary)
      cat("\n")
    }

    df <- data.frame(coefs[,dep_var], stds[,dep_var], tzstats[,dep_var], pValues[,dep_var])
    rownames(df) <- rownames(coefs)

    ms <- max(ifelse(df[,4] < .001, 3,
                     ifelse(df[,4] < .01, 2,
                            ifelse(df[,4] < .1, 1, 0))), na.rm = TRUE)
    stars <- logical(nrow(df))
    for (i in seq_len(nrow(df))){
      if (!is.null(x$estimations$isRestricted) && x$estimations$isRestricted[i,dep_var] == 1)
        stars[[i]] <- paste0("r", strrep(" ", ms - 1))
      else if (df[i,4] < .001)
        stars[[i]] <- "***"
      else if(df[i,4] < .01)
        stars[[i]] <- paste0("**", strrep(" ", ms - 2))
      else if(df[i,4] < .05)
        stars[[i]] <- paste0("*", strrep(" ", ms - 1))
      else if(df[i,4] < .1)
        stars[[i]] <- paste0(".", strrep(" ", ms - 1))
      else
        stars[[i]] <- strrep(" ", ms)
    }

    colnames(df) <- c("Estimate", "Std. Error", paste(dist, "stats"), paste0("Pr(>|t|)", strrep(" ", ms+1)))

    formatted_df <- format(df, digits = 4, nsmall = 4, justify = "left")
    formatted_df[,4] =  paste0(format(df[,4], digits = 4, nsmall = 4), " " , stars)
    formatted_df[df$estimations$isRestricted[,dep_var] == 1, 3:4] <- gsub("NaN", "  -", formatted_df[df[, 2] == 0, 3:4])
    print(formatted_df)

    cat("---\n")
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n r restricted\n\n")


    numRes <- ifelse (is.null(x$estimations$isRestricted), 0, sum(x$estimations$isRestricted[,dep_var]))
    cat("Number of restrictions in the equation: ", numRes,  "\n")
    cat("Number of estimated coefficients in the equation: ", nrow(x$estimations$coefs) - numRes,  "\n")
    if ("r2" %in% rownames(x$metrics))
      cat("R-squared: ", x$metrics["r2",dep_var], "\n")
    #cat("F-statistic: ", x$metrics["f",dep_var], " on ", x$counts$fProbD1, "and", x$counts$fProbD2, " DF,  p-value: ", x$metrics["fProb",dep_var]) #TODO: check its validity, esp. the degrees of freedom in multivariate case
    cat("Akaike Information Criterion: ", x$metrics["aic",dep_var], "\n")
    cat("Schwartz Information Criterion: ", x$metrics["sic",dep_var], "\n")
  }
  cat("\n")

  cat("----------------------------------------------------\n")

  cat("Number of observations: ", nrow(x$estimations$Y), "\n")
  omitted_dp <- nrow(x$info$data$data) - nrow(x$estimations$Y)
  if (omitted_dp != 0)
    cat("Number of omitted data-points: ", omitted_dp, "\n")
  cat("Number of equations in the system: ", x$info$data$numEndo, "\n")

  if (!is.null(x$projection)){
    cat(" ** Projection results are available.\n")
  }
  if (!is.null(x$simulation)){
    cat(" ** Simulation results are available (number of out-of-sample simulations: ", x$simulation$validCounts,")\n")
  }
  if (!is.null(x$info$pcaOptionsY)){
    cat(" ** Principle components are used for endogenous data.\n")
  }
  if (!is.null(x$info$pcaOptionsX)){
    cat(" ** Principle components are used for exogenous data.\n")
  }
  if (!is.null(x$info$searchSigMaxIter) && x$info$searchSigMaxIter > 0){
    cat(" ** Zero restrictions are automatically selected (alpha=", x$info$searchSigMaxProb,").\n")
  }
}


#' Prints an \code{ldt.varma.prediction} object
#'
#' An \code{ldt.varma.prediction} object is the output of [predict.ldt.estim.varma()] function.
#'
#' @param x An object of class \code{ldt.varma.prediction}
#' @param ... Additional arguments
#'
#' @return This function has no output.
#' @export
print.ldt.varma.prediction <- function(x, ...) {
  if (is.null(x))
    stop("argument is null.")
  if (!is(x, "ldt.varma.prediction"))
    stop("Invalid class. An 'ldt.varma.prediction' object is expected.")
  if (is.null(x$means))
    stop("Invalid data. Predictions are not available.")

  cat("LDT VARMA prediction:\n")
  rec.print.list(list(
    `Maximum horizon` = nrow(x$means) - x$actualCount,
    `Number of variables` = ncol(x$means),
    `Has variance` = !is.null(x$vars),
    `Horizon 1 at` = rownames(x$means)[x$actualCount+1],
    `Prediction at horizon 1` = x$means[x$actualCount+1,],
    `Variance at horizon 1` = if (is.null(x$vars)) {NULL} else {x$vars[x$actualCount+1,]}
  ), "  ")
}

#' Prints an \code{ldt.estim.projection} object
#'
#' An \code{ldt.estim.projection} object is the output of [predict.ldt.estim()] function.
#'
#' @param x An object of class \code{ldt.estim.projection}
#' @param ... Additional arguments
#'
#' @return This function has no output.
#' @export
print.ldt.estim.projection <- function(x, ...) {
  if (is.null(x))
    stop("argument is null.")
  if (!is(x, "ldt.estim.projection"))
    stop("Invalid class. An 'ldt.estim.projection' object is expected.")
  if (is.null(x$means) || nrow(x$means) == 0)
    stop("Invalid data. Predictions are not available.")

  cat("LDT ", x$method," model projection:\n")
  rec.print.list(list(
    `Length` = nrow(x$means),
    `Number of variables` = ncol(x$means),
    `Has variance` = !is.null(x$vars),
    `First exogenous observation` = x$newX[1,],
    `First prediction` = x$means[1,],
    `First variance` = if (is.null(x$vars)) {NULL} else {x$vars[1,]}
  ), "  ")
}

