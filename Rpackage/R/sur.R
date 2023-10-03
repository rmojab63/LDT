
#' Search for Best SUR Models
#'
#' Use this function to create a Seemingly Unrelated Regression model set and search for the best models (and other information) based on in-sample and out-of-sample evaluation metrics.
#'
#' @param data A list that determines data and other required information for the search process.
#' Use [get.data()] function to generate it from a \code{matrix} or a \code{data.frame}.
#' This function is designed to help you specify endogenous and exogenous variables by a list of equations.
#' It also helps to add intercept or define Box-Cox transformations for the endogenous variables.
#' @param combinations A list that determines the combinations of endogenous and exogenous variables in the search process.
#' Use [get.combinations()] function to define it.
#' This is a two-level nested loop and in this function, the outer loop is defined over the exogenous variables.
#' @param metrics A list of options for measuring performance.
#' Use [get.search.metrics] function to get them.
#' @param modelChecks A list of options for excluding a subset of the model set.
#' See and use [get.search.modelchecks] function to get them.
#' @param items A list of options for specifying the purpose of the search.
#' See and use [get.search.items] function to get them.
#' @param options A list of extra options for performing the search.
#' See and use [get.search.options] function to get them.
#' @param searchSigMaxIter Maximum number of iterations in searching for significant coefficients. Use 0 to disable the search.
#' @param searchSigMaxProb Maximum value of type I error to be used in searching for significant coefficients. If p-value is less than this, it is interpreted as significant.
#'
#' @return A nested list with the following members:
#' \item{counts}{Information about the expected number of models, number of estimated models, failed estimations, and some details about the failures.}
#' \item{...}{Results reported separately for each metric, then for each target variable, then for each requested type of output. This part of the output is highly nested and items are reported based on the arguments of the search.}
#' \item{info}{General information about the search process, some arguments, elapsed time, etc.}
#'
#' Note that the output does not contain any estimation results,
#' but minimum required data to estimate the models (Use \code{summary()} function to get the estimation).
#'
#' @export
#'
#' @example man-roxygen/ex-search.sur.R
#'
#' @seealso [estim.sur], [search.sur.stepwise]
search.sur <- function(data = get.data(),
                       combinations = get.combinations(),
                       metrics = get.search.metrics(),
                       modelChecks = get.search.modelchecks(),
                       items = get.search.items(),
                       options = get.search.options(),
                       searchSigMaxIter = 0,
                       searchSigMaxProb = 0.1){
  if (data$hasWeight)
    stop("SUR search does not support weighted observations.")

  combinations <- get.indexation(combinations, data, FALSE) # it also check for inconsistencies, etc.

  stopifnot(is.zero.or.positive.number(searchSigMaxIter))
  stopifnot(is.positive.number(searchSigMaxProb) && searchSigMaxProb < 1)
  if (searchSigMaxIter > 0 && searchSigMaxProb == 0)
    stop("searchSigMaxProb cannot be zero when search is enabled.")

  if (is.null(metrics))
    metrics = get.search.metrics()
  if (is.null(modelChecks))
    modelChecks = get.search.modelchecks()
  if (is.null(items))
    items = get.search.items()
  if (is.null(options))
    options = get.search.options()


  startTime <- Sys.time()
  res <- .SearchSur(data, combinations, metrics, modelChecks, items, options,
                    as.integer(searchSigMaxIter), searchSigMaxProb)
  endTime <- Sys.time()

  res$info$data <- data
  res$info$combinations <- combinations
  res$info$metrics <- metrics
  res$info$options <- options
  res$info$modelChecks <- modelChecks
  res$info$items <- items
  res$info$startTime <- startTime
  res$info$endTime <- endTime

  res$info$searchSigMaxIter <- searchSigMaxIter
  res$info$searchSigMaxProb <- searchSigMaxProb

  class(res) <- c("ldt.search.sur", "ldt.search", "list")
  attr(res, "method") <- "SUR"
  res
}


#' Estimate a SUR Model
#'
#' Use this function to estimate a Seemingly Unrelated Regression model.
#'
#' @param data A list that determines data and other required information for the search process.
#' Use [get.data()] function to generate it from a \code{matrix} or a \code{data.frame}.
#' This function is designed to help you specify endogenous and exogenous variables by a list of equations.
#' It also helps to add intercept or define Box-Cox transformations for the endogenous variables.
#' @param searchSigMaxIter An integer for the maximum number of iterations in searching for significant coefficients. Use 0 to disable the search.
#' @param searchSigMaxProb A number for the maximum value of type I error to be used in searching for significant coefficients. If p-value is less than this, it is interpreted as significant.
#' @param restriction A \code{km x q} matrix of restrictions where \code{m} is the number of endogenous data, \code{k} is the number of exogenous data, and \code{q} is the number of unrestricted coefficients.
#' @param pcaOptionsY A list of options to use principal components of the endogenous data, instead of the actual values. Set \code{NULL} to disable. Use [get.options.pca()] for initialization.
#' @param pcaOptionsX A list of options to use principal components of the exogenous data, instead of the actual values. Set \code{NULL} to disable. Use [get.options.pca()] for initialization.
#' @param simFixSize An integer that determines the number of out-of-sample simulations. Use zero to disable the simulation.
#' @param simTrainFixSize An integer representing the number of data points in the training sample in the out-of-sample simulation. If zero, \code{trainRatio} will be used.
#' @param simTrainRatio A number representing the size of the training sample relative to the available size, in the out-of-sample simulation. It is effective if \code{trainFixSize} is zero.
#' @param simSeed A seed for the random number generator. Use zero for a random value.
#' @param simMaxConditionNumber A number for the maximum value for the condition number in the simulation.
#'
#' @return A nested list with the following items:
#' \item{counts}{Information about different aspects of the estimation such as the number of observation, number of exogenous variables, etc.}
#' \item{estimations}{Estimated coefficients, standard errors, z-statistics, p-values, etc.}
#' \item{metrics}{Value of different goodness of fit and out-of-sample performance metrics. }
#' \item{projections}{Information on the projected values, if \code{newX} is provided.}
#' \item{info}{Some other general information.}
#'
#' @details
#' As described in section 10.2 in \insertCite{greene2020econometric;textual}{ldt}, this type of statistical model consists of multiple regression equations, where each equation may have a different set of independent variables and the disturbances between the equations are assumed to be correlated. The general form with \eqn{m} equations can be written as \eqn{y_i=z_i'\gamma_i+v_i} and \eqn{E(v_i v_j)=\sigma_{ij}^2} for \eqn{i=1,\ldots m}. Assuming that a sample of \eqn{N} independent observations is available, we can stack the observations and use the following system for estimation:
#' \deqn{
#' Y = X B + V, \quad \operatorname{vec}B = R\gamma,
#' }
#'
#' where the columns of \eqn{Y:N \times m} contain the dependent variables for each equation and the columns of \eqn{X: N\times k} contain the explanatory variables, with \eqn{k} being the number of unique explanatory variables in all equations. Note that $X$ combines the \eqn{z_i} variables, and the restrictions imposed by \eqn{R:mk\times q} and \eqn{\gamma:q\times 1} determine a set of zero constraints on \eqn{B: k \times m}, resulting in a system of equations with different sets of independent variables.
#'
#' Imposing restrictions on the model using the \eqn{R} matrix is not user-friendly, but it is suitable for use in this package, as users are not expected to specify such restrictions, but only to provide a list of potential regressors. Note that in this package, most procedures, including significance search, are supposed to be automated.
#'
#' The unrestricted estimators (i.e., \eqn{\hat{B}=(X'X)^{-1}X'Y}, and \eqn{\hat{\Sigma}=(\hat{V}'\hat{V})/N} where \eqn{\hat{V}=Y-X\hat{B}}) are used to initialize the feasible GLS estimators:
#'
#' \deqn{
#'     \tilde{B} = RW^{-1}R'[\hat{V}-1 \otimes x']\operatorname{vec}Y, \quad \tilde{\Sigma}=(\tilde{V}'\tilde{V})/N,
#' }
#' where \eqn{W = R'[\hat{V}^{-1} \otimes X'X]R} and \eqn{\tilde{V}=Y-X\tilde{B}}. The properties of these estimators are discussed in proposition 5.3 in \insertCite{lutkepohl2005new;textual}{ldt}. See also section 10.2 in \insertCite{greene2020econometric;textual}{ldt}. The maximum likelihood value is calculated by \eqn{-\frac{N}{2}(m(\ln 2\pi+1)+\ln|\tilde{\Sigma}|)}. The condition number is calculated by multiplying 1-norm of $W$ and its inverse (e.g., see page 94 in \insertCite{trefethen1997numerical;textual}{ldt}). Furthermore, given an out-of-sample observation such as \eqn{x:k\times 1}, the prediction is \eqn{y^f = \tilde{B}'x}, and its variance is estimated by the following formula:
#'
#' \deqn{
#'      \operatorname{var}y^f = \tilde{V} + (x' \otimes I_m)R W^{-1}R'(x \otimes I_m).
#' }
#'
#'
#' Note that the focus in \code{ldt} is model uncertainty and for more sophisticated implementations of the FGLS estimator, you may consider using other packages such as \code{systemfit}.
#'
#' Finally, note that the main purpose of exporting this method is to show the inner calculations of the search process in [search.sur] function.
#'
#' @references
#'   \insertAllCited{}
#' @importFrom Rdpack reprompt
#'
#' @export
#' @example man-roxygen/ex-estim.sur.R
#'
#' @seealso [search.sur], [search.sur.stepwise]
estim.sur <- function(data, searchSigMaxIter = 0,
                      searchSigMaxProb = 0.1,
                      restriction = NULL,
                      pcaOptionsY = NULL,
                      pcaOptionsX = NULL,
                      simFixSize = 0,
                      simTrainFixSize = 0,
                      simTrainRatio = 0.75,
                      simSeed = 0,
                      simMaxConditionNumber = Inf){
  data <- get.data.keep.complete(data)
  stopifnot(is.zero.or.positive.number(searchSigMaxIter))
  stopifnot(is.positive.number(searchSigMaxProb) && searchSigMaxProb < 1)
  if (searchSigMaxIter > 0 && searchSigMaxProb == 0)
    stop("searchSigMaxProb cannot be zero when search is enabled.")
  if (searchSigMaxIter > 0 && !is.null(restriction))
    stop("restriction must be null when  significant search is enabled.")

  if (!is.null(restriction))
    stopifnot(is.matrix(restriction))

  stopifnot(is.zero.or.positive.number(simFixSize))
  stopifnot(is.zero.or.positive.number(simTrainRatio) && simTrainRatio <= 1)
  stopifnot(is.zero.or.positive.number(simTrainFixSize))
  stopifnot(is.zero.or.positive.number(simSeed))
  if (simSeed == 0)
    simSeed = runif(1,10,10e4) # set it here such that it is reported in the info section
  stopifnot(is.number(simMaxConditionNumber))

  if (!is.null(pcaOptionsX))
    stopifnot(is.list(pcaOptionsX))
  if (!is.null(pcaOptionsY))
    stopifnot(is.list(pcaOptionsY))

  res <- .EstimSur(data, as.integer(searchSigMaxIter), searchSigMaxProb,
                   restriction, pcaOptionsY, pcaOptionsX,
                   as.integer(simFixSize), simTrainRatio,
                   as.integer(simTrainFixSize), as.integer(simSeed),
                   simMaxConditionNumber)


  res$info$data <- data
  res$info$searchSigMaxIter <- searchSigMaxIter
  res$info$searchSigMaxProb <- searchSigMaxProb
  res$info$pcaOptionsY <- pcaOptionsY
  res$info$pcaOptionsX <- pcaOptionsX
  res$info$simFixSize <- simFixSize
  res$info$simTrainFixSize <- simTrainFixSize
  res$info$simTrainRatio <- simTrainRatio
  res$info$simSeed <- simSeed
  res$info$simMaxConditionNumber <- simMaxConditionNumber

  class(res) <- c("ldt.estim.sur", "ldt.estim", "list")
  attr(res, "method") <- "SUR"
  res
}


estim.sur.from.search <- function(searchResult, endoIndices, exoIndices, ...){
  search_data <- searchResult$info$data
  exoIndices <- exoIndices + 1 # adjust for zero-based
  endoIndices <- endoIndices + 1
  data <- get.data(search_data$data[,c(endoIndices, exoIndices), drop = FALSE],
                   endogenous = length(endoIndices),
                   weights = NULL,
                   lambdas = search_data$lambdas,
                   addIntercept = FALSE,...)


  estim.sur(
    data = data,
    searchSigMaxIter = searchResult$info$searchSigMaxIter,
    searchSigMaxProb = searchResult$info$searchSigMaxProb,
    simFixSize = searchResult$info$metrics$simFixSize,
    simTrainRatio = searchResult$info$metrics$trainRatio,
    simTrainFixSize = searchResult$info$metrics$trainFixSize,
    simSeed = abs(searchResult$info$metrics$seed),
    simMaxConditionNumber = searchResult$info$modelChecks$maxConditionNumber
  )
}


# get estimation from search result
GetEstim_sur <- function(searchRes, endoIndices,
                         exoIndices, y, x, printMsg, ...) {
  y <- y[, endoIndices, drop = FALSE]
  x <- if (is.null(exoIndices) || is.null(x)) {
    NULL
  } else {
    x[, c(exoIndices), drop = FALSE]
  }

  M <- estim.sur(
    y = y,
    x = x,
    addIntercept = FALSE,
    searchSigMaxIter = searchRes$info$searchSigMaxIter,
    searchSigMaxProb = searchRes$info$searchSigMaxProb,
    restriction = NULL,
    newX = NULL,
    pcaOptionsY = NULL,
    pcaOptionsX = NULL,
    simFixSize = searchRes$info$metrics$simFixSize,
    simTrainRatio = searchRes$info$metrics$trainRatio,
    simTrainFixSize = searchRes$info$metrics$trainFixSize,
    simSeed = abs(searchRes$info$metrics$seed),
    simMaxConditionNumber = searchRes$info$modelChecks$maxConditionNumber,
    simTransform = searchRes$info$metrics$transform,
    printMsg = printMsg
  )

  return(M)
}


#' Step-wise Search for Best SUR Models
#'
#' For a large model set, use this function to find the best seemingly unrelated regression models.
#' It selects a subset of variables from smaller models and moves to the bigger ones.
#'
#' @param y A matrix of endogenous data with variables in the columns.
#' @param x A matrix of exogenous data with variables in the columns.
#' @param xSizeSteps A list of model dimensions to be estimated in each step.
#' Its size determines the number of steps.
#' @param countSteps An integer vector to determine the number of variables to be used in each step.
#' \code{NA} means all variables. Variables are selected based on best estimations.
#' All variables in the best models (all metrics and targets) are selected until the corresponding suggested number is reached.
#' Select an appropriate value for \code{bestK} in the options.
#' @param savePre A directory for saving and loading the progress.
#' Each step's result is saved in a file (name=\code{paste0(savePre,i)} where \code{i} is the index of the step.
#' @param ... other arguments to pass to [search.sur] function such as the \code{numTargets} argument.
#' Note that \code{xSizes} is ineffective here.
#'
#' @return Similar to [search.sur] function.
#' @export
#'
#' @examples
#' # See the example in the 'search.sur' function.
#'
#' @seealso [search.sur], [estim.sur]
search.sur.stepwise <- function(y, x, xSizeSteps = list(c(1), c(2)),
                                countSteps = c(NA, 10),
                                savePre = NULL, ...) {

  if (is.null(savePre) == FALSE && (is.character(savePre) == FALSE || length(savePre) > 1))
    stop("'savePre' must be a character string.")

  Search_s("sur", data = x, sizes = xSizeSteps, counts = countSteps,
          savePre = savePre, y = y, ...)
}


sur.to.latex.eqs <- function(sigma, coef, intercept, numFormat = "%.2f") {
  num_y <- ncol(coef)
  num_x <- nrow(coef) - ifelse(intercept, 1, 0)
  eqs <- character(num_y)
  for (i in seq_len(num_y)) {
    b <- coef[, i]
    if (intercept) {
      eqs[i] <- paste0("Y_", i, " = ", sprintf0(numFormat, b[1]))
      if (num_x > 0) {
        eqs[i] <- paste0(eqs[i], " + ", paste0(sprintf0(numFormat, b[-1]), " X_", seq_len(num_x), collapse = " + "))
      }
    } else {
      eqs[i] <- paste0("Y_", i, " = ", paste0(sprintf0(numFormat, b), " X_", seq_len(num_x), collapse = " + "))
    }
    eqs[i] <- paste0(eqs[i], " + E_", i, ", \\sigma_",i,"^2 = ", sprintf0(numFormat, sigma[[i,i]]))
  }
  eqs_latex <- paste(eqs, collapse = " \\\\ ")

  return(eqs_latex)
}

sur.to.latex.mat <- function(sigma, coef, intercept = TRUE, numFormat = "%.2f",
                             num_x_break = 3, y_label = "Y", x_label= "X", e_label = "E") {
  num_y <- ncol(coef)
  num_x <- nrow(coef)

  y_vec <- latex.variable.vector(num_y, y_label)

  coef_t <- t(coef)
  x_mat <- latex.matrix(mat = coef_t, numFormat = numFormat)

  x_vec <- latex.variable.vector(num_x, x_label, intercept, num_x_break)

  e_vec <- latex.variable.vector(num_y, e_label)

  s_mat <- latex.matrix(mat = sigma, numFormat = numFormat)

  eq_latex <- paste0(y_vec, " = ", x_mat," ", x_vec," + ", e_vec, ", \\Sigma = ", s_mat)



  return(eq_latex)
}


#' Generate Random Sample from an SUR Model
#'
#' This function generates a random sample from an Seemingly Unrelated Regression model.
#'
#' @param sigma covariance matrix of the errors.
#' If it is an integer value, it specifies the number of equations in the SUR model and covariance matrix is generated randomly.
#' @param coef Coefficients of the model.
#' If it is an integer value, it specifies the number of independent variables in each equation of the SUR model and coefficient matrix is generated randomly.
#' @param nObs Number of observations to generate.
#' @param intercept If \code{TRUE}, an intercept is included in the model as the first exogenous variable.
#' @param numFormat A character string that determines how to format the numbers, to be used as the argument of the \code{sprintf} function.
#' If \code{NULL}, conversion to latex or html representations are disabled.
#'
#' @return A list with the following items:
#'   \item{y}{matrix, the generated dependent variable.}
#'   \item{x}{matrix, the generated independent variable.}
#'   \item{e}{matrix, the generated errors.}
#'   \item{sigma}{matrix, the covariance matrix of the disturbances.}
#'   \item{coef}{matrix, the coefficients used in the model.}
#'   \item{intercept}{logical, whether an intercept was included in the model.}
#'   \item{eqsLatex}{character string, Latex representation of the equations of the system.}
#'   \item{eqsLatexSys}{character string, Latex representation of the system in matrix form.}
#'
#' @export
#' @importFrom stats rnorm
#' @example man-roxygen/ex-sim.sur.R
#' @seealso [sim.varma],[estim.sur],[search.sur]
sim.sur <- function(sigma = 1L, coef = 1L,
                    nObs = 100, intercept = TRUE,
                    numFormat = "%.2f") {

  nObs = as.integer(nObs)
  if (nObs <= 0)
    stop("nObs must be a positive integer.")

  intercept = as.logical(intercept)

  if (is.null(coef)) {
    stop("coef must be provided")
  }

  if (is.integer(coef) && length(coef) == 1) {
    num_x <- coef - ifelse(intercept, 1, 0)
    if (is.integer(sigma) && length(sigma) == 1) {
      num_y <- sigma
    } else {
      num_y <- nrow(sigma)
    }
    coef <- matrix(rnorm((num_x + ifelse(intercept, 1, 0)) * num_y), ncol = num_y)
  } else {
    num_y <- ncol(coef)
    num_x <- nrow(coef) - ifelse(intercept, 1, 0)
  }

  if (is.null(sigma)) {
    stop("sigma must be provided")
  }

  if (is.integer(sigma) && length(sigma) == 1) {
    sigma <- crossprod(matrix(rnorm(num_y^2), ncol = num_y))
  } else {
    if (!is.matrix(sigma) || nrow(sigma) != ncol(sigma)) {
      stop("sigma must be a square matrix")
    }
    if (nrow(sigma) != num_y) {
      stop("The number of rows and columns in sigma must be equal to num_y")
    }
  }

  if (num_x == 0){
    if (intercept)
      x <- matrix(rep(1,nObs), ncol = 1)
    else
      x <- matrix()
  }
  else{
    x <- matrix(rnorm(nObs * num_x), ncol = num_x)
    if (intercept){
      x <- cbind(rep(1,nObs),x)
    }
  }



  errors <- rand.mnormal(nObs, mu = rep(0, num_y), sigma = sigma)
  e <- errors$sample
  if (length(x) != 0)
    y <- x %*% coef + e
  else
    y <- e

  colnames(y) <- paste0("Y",c(1:ncol(y)))
  if (num_x == 0){
    if (intercept)
      colnames(x) <- c("Intercept")
  }
  else {
    if (intercept)
      colnames(x) <- c("Intercept", paste0("X",c(1:num_x)))
    else
      colnames(x) <- paste0("X",c(1:num_x))
  }
  colnames(e) <- paste0("E",c(1:ncol(e)))

  colnames(sigma) <- colnames(y)
  rownames(sigma) <- colnames(y)
  rownames(coef) <- colnames(x)
  colnames(coef) <- colnames(y)

  result <- list(y = y, x = x, e = e,
                 sigma = sigma, coef = coef,
                 intercept = intercept,
                 eqsLatex = ifelse(is.null(numFormat), NULL, sur.to.latex.eqs(sigma, coef, intercept, as.character(numFormat))),
                 eqsLatexSys = ifelse(is.null(numFormat), NULL, sur.to.latex.mat(sigma, coef, intercept, as.character(numFormat))))

  return(result)
}
