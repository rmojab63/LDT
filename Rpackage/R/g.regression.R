
checkEquation <- function(object, equations, justOne){
  if (justOne && is.null(equations))
    stop("'equation' cannot be null.")
  if (justOne && length(equations) != 1)
    stop("Length of 'equation' must be one.")
  if (justOne == FALSE && is.null(equations))
    equations <- c(1L:object$info$data$numEndo)
  stopifnot(is.numeric(equations) || is.character(equations))
  if (is.numeric(equations))
    equations <- as.integer(equations)
  else
    equations <- match(equations, colnames(object$info$data$data))
  if (any(equations < 1 | equations > object$info$data$numEndo)){
    if (justOne)
      stop("'equation' is not found.")
    else
      stop("At least one element in 'equations' is not valid or not found.")
  }
  equations
}

#' Extract Dependent Variable(s) Data
#'
#' This function extracts data of a dependent variable(s) from an estimated model.
#'
#' @param object An object of class \code{ldt.estim}.
#' @param equations A number, a numeric array or a string array specifying the indices or names of the endogenous variables in the equations. \code{NULL} means all equations.
#' @param ... Additional arguments.
#'
#' @return A data frame containing the endogenous data.
#'
#' @export
endogenous <- function(object, equations = NULL, ...){
  equations <- checkEquation(object, equations, FALSE)
  object$info$data$data[, equations, drop=FALSE]
}

#' Extract Exogenous Variable(s) Data
#'
#' This function extracts data of an exogenous variable(s) in an equation from an estimated model.
#' It takes zero restrictions imposed into account.
#'
#' @param object An object of class \code{ldt.estim}.
#' @param equation A number or a string specifying the equation with exogenous data.
#' @param ... Additional arguments.
#'
#' @return A matrix containing the exogenous data.
#'
#' @export
exogenous <- function(object, equation = 1, ...){
  equation <- checkEquation(object, equation, TRUE)
  exogenous <- c((object$info$data$numEndo+1+ifelse(object$info$data$hasWeight,1,0)):ncol(object$info$data$data))
  if (!is.null(object$estimations$isRestricted)){
    not_restricted <- as.logical(1 - object$estimations$isRestricted[,equation])
    exogenous <- exogenous[not_restricted]
  }
  object$info$data$data[, exogenous, drop=FALSE]
}


#' Extract Coefficients Matrix
#'
#' This function extracts coefficient matrix from an \code{ldt.estim} object.
#'
#' @param object An object of class \code{ldt.estim}
#' @param equations A number, a numeric array or a string array specifying the indices or names of the endogenous variables in the equations.
#' \code{NULL} means all equations.
#' @param removeZeroRest If \code{TRUE} and model is restricted, zero restrictions are removed.
#' @param ... Additional arguments.
#'
#' @return If zero restrictions are not removed, it is a matrix containing the coefficients of the system.
#' Each column of the matrix belongs to an equation. Explanatory variables are in the rows.
#' Otherwise, coefficients of different equations are reported in a list.
#' @export
coef.ldt.estim <- function(object, equations = NULL, removeZeroRest = FALSE, ...){
  equations <- checkEquation(object, equations, FALSE)
  c <- object$estimations$coefs[, equations]
  if (removeZeroRest == FALSE ||
      is.null(object$estimations$isRestricted) ||
      sum(object$estimations$isRestricted) == 0)
    return(c)

  result <- list()
  for (j in 1:ncol(c)){
    notRest <- as.logical(1 - object$estimations$isRestricted[,j])
    result[[colnames(c)[j]]] <- c[,j,drop=FALSE][notRest]
  }
  result
}


hatvalues0 <- function(object, equation){
  X <- exogenous(object, equation = equation)
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  diag(H)
}

cooks.distance0 <- function(object, equation){
  h <- hatvalues0(object, equation)
  r <- resid(object, equation)
  mse <- sum(r^2)/length(r)
  c <- coef(object, equations = equation, removeZeroRest = TRUE)
  p <- length(c[1])
  if (object$info$data$hasWeight){
    stop("Cook's distance is not implemented for weighted observations.") # couldn't find a reference
  }
  else
    r^2 * h / (p*mse*(1-h^2))
}

#' Extract Residuals Data
#'
#' This function returns residuals from or calculates the standardized residuals for an \code{ldt.estim} object.
#'
#' @param object An object of class \code{ldt.estim}.
#' @param equation A number, a numeric array or a string array specifying the equations with residual data. If \code{NULL}, residuals in all equations are returned.
#' @param type Determines the type of the calculations. If more than one type is specified, the first item is used. See details.
#' @param ...
#'
#' @details
#' The standardized residuals have identical variance.
#' In order to calculate the standardized residuals, each residual is divided by \eqn{s\sqrt{w_i(1-h_{ii})}} where \eqn{s} is the standard error of residuals and \eqn{h_{ii}} is the leverage of \eqn{i}-th observation. \eqn{w_i} is the weight of the \eqn{i}-th observation if data is weighted, and 1 otherwise.
#' Note that while the residuals are estimated in a system, the \eqn{h_{ii}} is calculated in a univariate context as the \eqn{i}-th diagonal of \eqn{X(X'X)^{-1}X'} matrix, where \eqn{X} is the exogenous variables in the corresponding equation.
#'
#' @return A matrix containing the residuals data.
#' @export
residuals.ldt.estim <- function(object, equations = NULL, standardized = FALSE, ...){
  stopifnot(is.logical(standardized))
  equations <- checkEquation(object, equations, FALSE)
  resid <- object$estimations$resid[, equations, drop=FALSE]
  if (standardized) {
    for (i in equations){
      hii <- hatvalues0(object, equation = i)
      s <- sd(resid)
      if (object$info$data$hasWeight){
        weights <- object$info$data[,object$info$data$numEndo + 1]
        resid[,i] <- resid[,i] / (s * weights * sqrt(1 - hii))
      }
      else{
        resid[,i] <- resid[,i] / (s * sqrt(1 - hii))
      }
    }
  }
  resid
}

#' Extract Fitted Data
#'
#' This function calculates and returns fitted values for an \code{ldt.estim} object.
#'
#' @param object An object of class \code{ldt.estim}.
#' @param equation A number, a numeric array or a string array specifying the equations with residual data. If \code{NULL}, residuals in all equations are returned.
#' @param ... Additional arguments.
#'
#' @return A matrix containing the exogenous data.
#'
#' @export
fitted.ldt.estim <- function(object, equations = NULL, ...){
  equations <- checkEquation(object, equations, FALSE)
  resid <- object$estimations$resid[, equations, drop=FALSE]
  y <- object$info$data$data[, equations, drop=FALSE]
  y - resid
}

#' Extract Maximum Log-Likelihood
#'
#' This function extracts maximum log-likelihood from an \code{ldt.estim} object.
#'
#' @param object An object of class \code{ldt.estim}
#' @param ... Additional arguments.
#'
#' @return The value of the maximum log-likelihood for the whole system.
#' @export
logLik.ldt.estim <- function(object, ...){
  #equation <- checkEquation(object, equation, TRUE) It is not equation-specific
  object$metrics["logL", 1]
}

#' Akaike Information Criterion
#'
#' This function extracts Akaike information criterion from an \code{ldt.estim} object.
#'
#' @param object An object of class \code{ldt.estim}
#' @param ... Additional arguments.
#'
#' @return The value of AIC for the whole system.
#' @export
AIC.ldt.estim <- function(object, ...){
  #equation <- checkEquation(object, equation, TRUE) It is not equation-specific
  object$metrics["aic", 1]
}

#' Bayesian Information Criterion
#'
#' This function extracts Baysian information criterion from an \code{ldt.estim} object.
#'
#' @param object An object of class \code{ldt.estim}
#' @param ... Additional arguments.
#'
#' @return The value of BIC for the whole system.
#' @export
BIC.ldt.estim <- function(object, ...){
  #equation <- checkEquation(object, equation, TRUE) It is not equation-specific
  object$metrics["sic", 1]
}



#' Extract Prediction Results
#'
#' This function extracts predicted mean and its variance from an \code{ldt.estim} object.
#' new data must be provided while estimating the model.
#'
#' @param object An object of class \code{ldt.estim}
#' @param ... Additional arguments.
#'
#' @return A list containing the predicted (projected) mean and variancesf.
#' @export
predict.ldt.estim <- function(object, ...){
  object$projection
}


