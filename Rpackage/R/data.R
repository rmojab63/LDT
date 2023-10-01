
#' Transform and Prepare Data for Analysis
#'
#' This function prepares a data matrix for analysis by applying a Box-Cox transformation to the endogenous variables, adding an intercept column, and optionally adding new rows with exogenous data.
#'
#' @param data A data.frame or a numeric matrix as the main data source.
#' @param endogenous A single number for the number of endogenous variables at the first columns or a list of names indicating the endogenous variables. The rest will be exogenous variables.
#' @param equations A formula or a list of formula objects representing the equations to be used instead of \code{endogenous}. If given, the final data is a matrix where the response variables are in the first columns, and the predictor variables are in the subsequent columns.
#' @param lambdas A numeric vector, a single number, NA, or NULL indicating the lambda parameter(s) for the Box-Cox transformation.
#' Use \code{NULL} for no transformation, \code{NA} for estimating the lambda parameter for each variable,
#' a single number for equal lambda parameter for all variables, and a numeric vector for distinct lambda parameters for the corresponding variables.
#' @param newData A data.frame or a numeric matrix representing new data for exogenous variables. Its structure should be similar to \code{data} without endogenous or response variables.
#' @param addIntercept A logical value indicating whether to add an intercept column to the final matrix.
#' @param weights A numeric vector or a column matrix representing weights of observations. It is not implemented in all applications.
#' @param ... additional parameters for \code{MASS::boxcox} function.
#'
#' @return A list of items suitable for being used in the \code{ldt::search.?} functions. It has the following members:
#' \item{data}{The final data matrix. Endogenous variables are in the first columns, then comes the weights (if given), then the intercept (if added), and finally the exogenous variables.}
#' \item{numEndo}{The number of endogenous variables in the data.}
#' \item{numExo}{The number of exogenous variables in the data (including 'intercept' if it is added).}
#' \item{obsCount}{The number of observations in the original data.}
#' \item{newX}{The matrix of new observations for the exogenous variables.}
#' \item{newObsCount}{The number of observations in the new data.}
#' \item{lambdas}{The lambda parameters used in the Box-Cox transformation.}
#' \item{hasIntercept}{Indicates whether an intercept column in added in the final matrix.}
#' \item{hasWeight}{Indicates whether there is a weight column in the final matrix.}
#'
#' @examples
#' # Example 1:
#' data <- matrix(1:24, ncol = 6)
#' result <- get.data(data, endogenous = 1)
#' print(result$data)
#'
#' # Example 2:
#' data <- matrix(1:24, ncol = 6,
#'                dimnames = list(NULL,c("V1", "V2", "V3", "V4", "V5", "V6")))
#' result <- get.data(data, endogenous = c("V6", "V1"))
#' print(result$data)
#'
#' # Example 3:
#' data <- data.frame(matrix(1:24, ncol = 6))
#' colnames(data) <- c("X1", "X2", "Y2", "X3", "Y1", "X4")
#' equations <- list(
#'    Y1 ~ X2 + X1,
#'    Y2 ~ X4 + X3)
#' result <- get.data(data, equations = equations)
#' print(result$data)
#'
#' @export
get.data <- function(data, endogenous = 1, equations = NULL,
                     weights = NULL, lambdas = NULL, newData = NULL,
                     addIntercept = TRUE, ...) {

  stopifnot(is.matrix(data) || is.data.frame(data))
  if (!is.null(newData))
    stopifnot(is.matrix(newData) || is.data.frame(newData))

  if (is.null(equations)) { # use endogenous
    data <- as.matrix(data)

    stopifnot(is.number(endogenous) || is.character(endogenous))
    numEndo <- ifelse(is.numeric(endogenous), as.integer(endogenous), length(endogenous))
    if (numEndo == 0)
      stop("Number of endogenous variables cannot be zero.")
    if (numEndo > ncol(data))
      stop("Invalid number of endogenous variables. It is larger than the number of columns in data.")
    numExo <- ncol(data) - numEndo

    # make sure it is a numeric matrix:
    data <- matrix(as.numeric(data), ncol = ncol(data), dimnames = list(NULL, colnames(data)))

    # column names might be missing:
    if (any(nchar(colnames(data)) == 0))
      stop("data has columns with empty names.")
    if (is.null(colnames(data))) {
      if (is.character(endogenous))
        stop("data has no column name while endogenous is a list of names.")

      colnames(data) <- c(paste0("Y", 1:numEndo), paste0("X", 1:numExo))
      if (numExo > 0){
        if (!is.null(newData)){
          if (ncol(newData) != numExo)
            stop("newData does not have the same number of columns as the exogenous part of data.")
          colnames(newData) <- paste0("X", 1:numExo) # override column names (?!)
        } # else, newData is null
      } # there is no exogenous variable
    } # else, data has column names
    else{
      if (is.null(newData) == FALSE && is.null(colnames(newData)))
        stop("data has column names, while column names in newData is missing")
    }
    if (is.character(endogenous)){ # reorder columns
      idx <- match(endogenous, colnames(data))
      data <- data[, c(idx, setdiff(seq_len(ncol(data)), idx))]
    }

    if (addIntercept){
      data <- cbind(data[,1:numEndo,drop=FALSE],rep(1,nrow(data)), data[,(numEndo+1):ncol(data),drop=FALSE])
      colnames(data)[numEndo + 1] <- "(Intercept)"
      if (!is.null(newData)){
        newData <- cbind(rep(1,nrow(newData)), newData)
        colnames(newData)[1] <- "(Intercept)"
      }
      numExo = numExo + 1
    }

  }
  else{

    data0 <- eqList2Matrix(equations, data, addIntercept)
    data <- data0$result
    numEndo <- data0$numResponse
    if (numEndo == 0)
      stop("Number of endogenous variables cannot be zero.")
    numExo <- ncol(data) - numEndo # it counts the intercept

    # deal with newData
    if (!is.null(newData)){
      newData <- as.data.frame(newData)
      # add endogenous data if missing
      responseNames <- colnames(data)[1:numEndo]
      for (name in responseNames) {
        if (!name %in% names(newData))
          newData[[name]] <- 1.0 # we will ignore this part
      }
      newData <- tryCatch(eqList2Matrix(equations, newData, addIntercept)$result,
                          error = function(err){
                            stop("Compared with the exogenous variables, newData either lacks some necessary columns or does not have the correct column names. Additional information: ", err)
                          })
    }
  }

  # make sure data and newData are numeric matrices
  data <- matrix(as.numeric(data), ncol = ncol(data), dimnames = list(NULL, colnames(data)))
  if (!is.null(newData))
    newData <- matrix(as.numeric(newData), ncol = ncol(newData), dimnames = list(NULL, colnames(newData)))
  obsCount <- nrow(data)

  if (!is.null(lambdas)) {
    bcRes <- boxCoxTransform(data[, 1:numEndo,drop=FALSE], lambdas, ...)
    data[, 1:numEndo] <- bcRes$data
    lambdas <- bcRes$lambda
  }

  if(!is.null(weights)){
    stopifnot(is.numeric(weights) || (is.matrix(weights) && ncol(weights) == 1))
    stopifnot(length(weights) == nrow(data))
    data <- cbind(data[,1:numEndo,drop=FALSE], weights, data[,(numEndo+1):ncol(data),drop=FALSE])
    colnames(data)[numEndo+1] <- "(weights)"
  }

  if (!is.null(newData)) {
    exoNames <- colnames(data[, (numEndo+1+ifelse(is.null(weights),0,1)):ncol(data),drop=FALSE])
    if (is.null(colnames(newData)) ||
        !all(exoNames %in% colnames(newData)))
      stop("Error:Compared with the exogenous variables, newData either lacks some necessary columns or does not have the correct column names.")

    newData <- newData[,exoNames, drop=FALSE] # reorder and filter columns
  }

  res <- list(data = data,
              numEndo = numEndo,
              numExo = numExo,
              obsCount = obsCount,
              newX = newData, # it contains just exogenous data, if not NULL
              newObsCount = ifelse(is.null(newData),0, nrow(newData)),
              lambdas = lambdas,
              hasIntercept = addIntercept,
              hasWeight = !is.null(weights))
  class(res) <- c("ldt.search.data")
  return(res)
}

#' Append \cpde{newX} to \code{data$data} matrix.
#'
#' Use it for VARMA estimation
#'
#' @param data The output of [get.data] function.
#'
#' @return The input \code{data} with updated data matrix
get.data.append.newX <- function(data){
  new_rows <- cbind(matrix(NA,
                           ncol = data$numEndo + ifelse(data$hasWeight,1,0),
                           nrow = nrow(data$newX)), data$newX)
  colnames(new_rows) <- colnames(data$data)
  data$data <- rbind(data$data, new_rows)
  data
}

#' Remove Rows with Missing Observations from Data
#'
#' @param data Output of [get.data] function
#' @param warn If true,  warning message about the indices of the removed rows is shown
#'
#' @return The input \code{data} but with updated \code{data$data} and \code{data$obsCount}
get.data.keep.complete <- function(data, warn = TRUE){
  incomplete_rows <- which(apply(data$data, 1, function(x) any(is.na(x) | is.nan(x))))

  if (length(incomplete_rows) > 0) {
    if (warn)
      warning(paste("Removed rows with missing observations: ", paste(incomplete_rows, collapse = ", ")))
    data$data <- data$data[-incomplete_rows, ]
    data$obsCount <- nrow(data$data)
  }
  data
}

#' Convert a List of Equations to a Matrix
#'
#' This function takes a list of equations and a data frame, and returns a matrix where the response variables are in the first columns, and the predictor variables are in the subsequent columns.
#'
#' The function checks for duplicate response variables across equations and throws an error if any are found. It also ensures that predictor variables are not duplicated in the final matrix.
#'
#' @param equations A formula or a list of formula objects representing the equations.
#' @param data A matrix or a data frame containing the variables in the equations.
#' @param addIntercept Logical. If \code{TRUE}, an intercept column is added after the response variables. Default is \code{FALSE}.
#'
#' @return
#' \item{result}{A matrix with response variables in the first columns, predictor variables in subsequent columns, and optionally an intercept column. The matrix does not contain duplicate columns.}
#' \item{numResponse}{Number of response variables in the first columns.}
#'
#' @examples
#' data <- data.frame(income = c(50000, 60000, 80000, 85000, 65000),
#'                    age = c(25, 32, 47, 51, 36),
#'                    education = c(16, 18, 20, 20, 16),
#'                    savings = c(20000, 25000, 30000, 35000, 40000))
#' equations <- list(as.formula("income ~ age + education"),
#'                   as.formula("savings ~ age + education"))
#' matrix_data <- eqList2Matrix(equations, data, addIntercept = TRUE)
#' print(matrix_data)
#'
eqList2Matrix <- function(equations, data, addIntercept = FALSE) {

  stopifnot(inherits(equations, "formula") || is.list(equations))
  if (inherits(equations, "formula"))
    equations <- list(equations)
  for(eq. in equations)
    stopifnot(inherits(eq., "formula"))
  stopifnot(is.data.frame(data) || is.matrix(data))
  stopifnot(is.logical(addIntercept))

  data <- as.data.frame(data)

  res <- data.frame() # storage
  response_vars <- c() # names of response variables

  for (eq in equations) {
    rv <- all.vars(eq)[1] # response variable
    if (rv %in% response_vars)
      stop("Duplicate response variable '", rv, "' found in equations.")
    response_vars <- c(response_vars, rv)

    model_mat <- model.matrix(eq, data)
    if ("(Intercept)" %in% colnames(model_mat)) # remove intercept
      model_mat <- model_mat[, colnames(model_mat) != "(Intercept)", drop = FALSE]

    # Add the response variable as a column of the matrix
    rvdata <- model.response(model.frame(eq, data))
    model_mat <- cbind(rvdata, model_mat)
    colnames(model_mat)[1] <- rv

    # Add the predictor variables to the matrix, excluding those already present
    newvars <- setdiff(colnames(model_mat), colnames(res))
    if (nrow(res) == 0) {
      res <- model_mat
    } else {
      res <- cbind(res, model_mat[, newvars, drop = FALSE])
    }
  }

  res <- res[, c(response_vars, setdiff(colnames(res), response_vars))]

  if (addIntercept) {
    intercept <- rep(1, nrow(res))
    res <- cbind(res[,1:length(response_vars)], intercept, res[,-(1:length(response_vars)), drop = FALSE])
    colnames(res)[length(response_vars)+1] <- "(Intercept)"
  }

  return(list(result = res, numResponse = length(response_vars)))
}

#' Box-Cox Transformation of Numeric Matrix
#'
#' This function applies the Box-Cox transformation to the columns of a numeric matrix.
#'
#' @param data A numeric matrix to be transformed.
#' @param lambda A numeric vector, a single number, NA, or NULL indicating the lambda parameter(s) for the Box-Cox transformation.
#' Use \code{NULL} for no transformation, \code{NA} for estimating the lambda parameter for each variable,
#' a single number for equal lambda parameter for all variables, and a numeric vector for distinct lambda parameters for the corresponding variables.
#' @param ... additional parameters for \code{MASS::boxcox} function.
#'
#' @return
#' \code{data}{transformed data}
#' \code{lambda}{final lambda vector used in the calculations.}
#'
#' @examples
#' data <- matrix(rnorm(40), ncol = 2)
#' result <- boxCoxTransform(data, c(0.5, 0.5))
#'
#' @importFrom MASS boxcox
boxCoxTransform <- function(data, lambda, ...) {

  if (!is.matrix(data) || !is.numeric(data))
    stop("Data should be a numeric matrix")

  if (!is.null(lambda) && !is.numeric(lambda) && length(lambda) != 1 && length(lambda) != ncol(data))
    stop("lambda should be NULL, NA, a single number, or a vector of numbers with length equal to 'ncol(data)'.")

  if (is.null(lambda))
    return(list(data = data, lambda = NULL))

  lambda_used <- numeric(ncol(data))
  res <- data

  for (i in 1:ncol(data)) {
    if ((length(lambda) == 1 && is.na(lambda)) ||
        (length(lambda) == ncol(data) && is.na(lambda[i]))) { # estimate
      model <- boxcox(data[,i]~1, plotit = FALSE, ...)
      li <- model$x[which.max(model$y)]
    }
    else if (length(lambda) == 1)
      li <- lambda
    else if (length(lambda) == ncol(data))
      li <- lambda[i]
    else
      stop("Invalid lambda argument")

    lambda_used[i] <- li

    if (li == 0)
      res[,i] <- log(data[,i])
    else
      res[,i] <- (data[,i]^li - 1) / li
  }

  return(list(data = res, lambda = lambda_used))
}



#' Define Combinations for Search Process
#'
#' Assuming a data matrix with variables in the columns and specific column names is generated by [get.data()] function,
#' this function defines a structure for a two-level nested loop used in a search process.
#' The outer loop is defined by a vector of sizes and all the combinations of the variables are generated automatically.
#' The inner loop is defined by a list of predefined combinations of the variables.
#' Each variable can belong to either endogenous or exogenous variables based on their usage.
#'
#' @param sizes A numeric vector that determines the sizes of outer loop combinations.
#' For example, if the outer loop belongs to the endogenous variables, \code{c(1, 2)} means all models with 1 and 2 equations.
#' If the outer loop belongs to exogenous variables, \code{c(1,2)} means all regressions with 1 and 2 exogenous variables.
#' @param partitions A list of numeric vectors or character vectors that partitions the outer loop variables. No model is estimated with two variables from the same partition.
#' @param numFixPartitions A single number that determines the number of partitions at the beginning of \code{partitions} to be included in all models.
#' @param innerGroups A list of numeric vectors or character vectors that determines different combinations of the variables for the inner loop.
#' For example, if the inner loop belongs to exogenous data, \code{list(c(1), c(1, 2))} means estimating all models with just the first exogenous variable and all models with the first and second exogenous variables.
#' @param numTargets An integer for the number of target variables at the first columns of the data matrix. Results of a search process are specific to these variables. A model is not estimated if it does not contain a target variable.
#'
#' @details
#' The indexation in \code{ldt} is a two-level nested loop over different combinations of endogenous and exogenous variables, similar to the following code:
#' \preformatted{
#' for (endo in list(c(1), c(1, 2)))
#'   for (exo in list(c(1), c(1, 2)))
#'     Estimate a model by using \code{endo} and \code{exo} indexation
#' }
#' However, it is not memory efficient to predefine both loops and \code{ldt} uses a running algorithm to define the outer loop.
#' It asks for the desired size of endogenous or exogenous variables in the model (i.e., \code{sizes})
#' and creates the outer groups by using all possible combinations of the variables.
#' Of course, the \code{partitions} and \code{numFixPartitions} can be used to restrict this set.
#' However, for the inner loop, the desired combination of the variables (endogenous or exogenous) must be provided.
#'
#' Note that given \code{m} as the number of variables, you can generate all possible combinations by the following code:
#' \preformatted{
#' m <- 4
#' combinations <- unlist(lapply(1:m, function(i) {
#'  t(combn(1:m, i, simplify = FALSE))
#' }), recursive = FALSE)
#' }
#' and use it as the \code{innerGroups} argument. Of course, this might result in a large model set.
#'
#' Also note that in \code{ldt}, if the data matrix does not have column names,
#' default names for the endogenous variables are \code{Y1, Y2, ...}
#' and default names for the exogenous variables are \code{X1, X2, ...}. See [get.data()] function.
#'
#' @return A list of items suitable for being used in the \code{ldt::search.?} functions.
#' @export
#'
get.combinations <- function(sizes = c(1, 2),
                             partitions = NULL,
                             numFixPartitions = 0,
                             innerGroups = list(c(1)),
                             numTargets = 1){
  stopifnot(is.numeric(sizes))
  stopifnot(is.numeric(numTargets) && length(numTargets) == 1 && numTargets >= 0)

  stopifnot(is.null(partitions) || is.list(partitions))
  if (!is.null(partitions))
    for (partition in partitions)
      stopifnot(is.numeric(partition) || is.character(partition))

  stopifnot(is.numeric(numFixPartitions) && length(numFixPartitions) == 1)

  stopifnot(is.list(innerGroups))
  for (innerGroup in innerGroups)
    stopifnot(is.numeric(innerGroup) || is.character(innerGroup))

  res <- list(sizes = sizes,
              partitions = partitions,
              numFixPartitions = numFixPartitions,
              innerGroups = innerGroups,
              numTargets = numTargets)
  class(res) <- "ldt.search.combinations"
  return(res)
}


#' Get Numeric Indices in a Combination
#'
#' This function takes the output of the [get.combinations] function and a numeric matrix with given column names.
#' It converts all character vectors in \code{innerGroups} or \code{partitions} to numeric vectors based on the index of the columns.
#'
#' @param combinations A list returned by the [get.combinations] function.
#' @param data A list returned by [get.data] function.
#' @param outerOverEndogenous \code{TRUE} if outer loop is defined over the endogenous variables. \code{FALSE} for exogenous.
#'
#' @return A list similar to the input \code{combinations}, but with all character vectors in \code{innerGroups} or \code{partitions} converted to numeric vectors based on the index of the columns in the \code{data} matrix.
#' It sums the exogenous indexes with the number of endogenous variables and returns zero-based indexation for C code.
#'
get.indexation <- function(combinations, data, outerOverEndogenous) {

  stopifnot(!is.null(data))
  if (!(is(data, "ldt.search.data")))
    stop("Invalid class. Use 'get.data()' function to generate 'data'.")
  stopifnot(!is.null(combinations))
  if (!(is(combinations, "ldt.search.combinations")))
    stop("Invalid class. Use 'get.combinations()' function to generate 'combinations'.")
  stopifnot(is.matrix(data$data) && is.numeric(data$data) && !is.null(colnames(data$data)))

  if (outerOverEndogenous == FALSE && data$numExo == 0)
    stop("Number of exogenous variables cannot be zero.") #TODO: check for intercept
  # we check for non-zero number of endogenous variables in 'get.data'

  if (combinations$numTargets > data$numEndo)
    stop("Invalid number of target variables (=",combinations$numTargets,
         "). It cannot be larger than the number of endogenous variables (=", data$numEndo, ").")

  colnames_to_index <- function(names, colnames) {
    indices <- match(names, colnames)
    if (any(is.na(indices))) {
      stop("Some column names are not found in the data matrix.")
    }
    return(as.integer(indices))
  }

  if (!is.null(combinations$partitions)) {
    combinations$partitions <-
      lapply(combinations$partitions,
             function(x){
               j <- if (is.character(x)) {colnames_to_index(x, colnames(data$data))}
               else {as.integer(x)}
               if (any(j) < 1)
                 stop("Zero or negative index in partitions.")
               if (outerOverEndogenous && any(j > data$numEndo))
                 stop("Invalid index in partitions. Make sure they are less than the number of endogenous variables.")
               if (outerOverEndogenous == FALSE && any(j > data$numExo))
                 stop("Invalid index in partitions. Make sure they are less than the number of exogenous variables.")
               j
             })

  }
  else{ # put each variable in its own partition
    combinations$partitions <-
      lapply(1:ifelse(outerOverEndogenous, data$numEndo, data$numExo),
             function(x)c(x))
  }

  combinations$innerGroups <-
    lapply(combinations$innerGroups,
           function(x){
             j <- if (is.character(x))
               {colnames_to_index(x, colnames(data$data))}
             else
               {as.integer(x)}
             if (any(j < 1))
               stop("Zero or negative index in innerGroups.")
             if (outerOverEndogenous && any(j > data$numExo))
               stop("Invalid index in innerGroups. Make sure they are less than the number of exogenous variables.")
             if (outerOverEndogenous  == FALSE && any(j > data$numEndo))
               stop("Invalid index in innerGroups. Make sure they are less than the number of endogenous variables.")
             j
           })

  # zero-based indexation and moving exogenous indices:
  for (i in c(1:length(combinations$partitions))){
    combinations$partitions[[i]] <- combinations$partitions[[i]] - 1
    if (outerOverEndogenous == FALSE)
      combinations$partitions[[i]] <- combinations$partitions[[i]] + data$numEndo
  }
  for (i in c(1:length(combinations$innerGroups)))
    combinations$innerGroups[[i]] <- combinations$innerGroups[[i]] - 1

  return(combinations)
}
