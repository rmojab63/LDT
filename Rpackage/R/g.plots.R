
#' Plot Diagnostics for \code{ldt.estim} Object
#'
#' This function creates diagnostic plots for estimated regression models of \code{ldt.estim} class.
#'
#' @param x An object of type \code{ldt.estim}.
#' @param equation A number or a name of endogenous variable specifying an equation in the estimated system.
#' @param type One of these numbers: 1, 2, 3, or 5. See \code{which} argument in [plot.lm] documentation.
#' @param textArgs A list of additional arguments to customize the "abline" function.
#' @param ablineArgs A list of additional arguments to customize the "text" function used for labeling influential observations.
#' @param ... additional arguments to be passed to "plot" (or "qqnorm" function for \code{type=2}, or "barplot" for \code{type=4}).
#'
#' @details
#' This function is designed to be similar to [plot.lm] function.
#' However, note that an \code{ldt.estim} object might be a system estimation.
#'
#' Some plots use standardized residuals. Note that they are not calculated in a system estimation context. See [residuals.ldt.estim] documentation for a description.
#' Cook's distance is also calculated equation-wise. Its formula is:
#' \deqn{
#' d = \frac{r_i^2}{k*var(r)}\frac{h_{ii}}{(1-h_{ii})^2}
#' }
#' where \eqn{r_i} and \eqn{h_{ii}} are residual and leverage in \eqn{i}-th observation, respectively. \eqn{var(r)} is variance of residuals and \eqn{k} is the number of estimated coefficients in the equation.
#' Note that Cook's distance is not implemented for weighted observations.
#'
#'
#' @return This function creates diagnostic plots for regression models.
#' It also returns a list with \code{x} and \code{y} data used in plot functions.
#' @export
plot.ldt.estim <- function(x,
                           equation = 1,
                           type = c(1,2,3,4,5,6),
                           ablineArgs = list(col = "lightblue"),
                           textArgs = list(pos = 3, cex = 0.7, col = "red"),
                           ...) {
  equation <- checkEquation(x, equation, TRUE)
  stopifnot(is.numeric(type))
  if (length(type) > 1)
    type = type[1]
  stopifnot(type %in% c(1, 2, 3, 4, 5, 6))

  fitted <- fitted(x, equation = equation)
  if (type == 1)
    residuals <- resid(x, equation = equation, standardized = FALSE)
  else
    residuals <- resid(x, equation = equation, standardized = TRUE)

  args <- list(...)

  x_data <- NULL
  y_data <- NULL
  if (type == 1) { # Residuals vs Fitted
    x_data <- fitted
    y_data <- residuals
    args <- modifyList(list(main = "Residuals vs Fitted", xlab = "Fitted Values", ylab = "Residuals"), args)

    do.call(plot, c(list(x=x_data, y = y_data), args))
    do.call(abline, c(list(h = 0), ablineArgs))

    high_data <- order(abs(y_data), decreasing = TRUE)[1:(length(y_data)/10)]
    if (length(high_data) > 0)
      do.call(text, c(list(x = x_data[high_data], y = y_data[high_data], labels = high_data), textArgs))
  }
  else if (type == 2) { # Residual Q-Q
    y_data <- residuals
    args <- modifyList(list(main = "Normal Q-Q", ylab = "Standardized Residuals"), args)

    points <- do.call(qqnorm, c(list(y_data), args))
    qqline(y_data)

    high_data <- order(abs(y_data), decreasing = TRUE)[1:(length(y_data)/10)]
    if (length(high_data) > 0)
      do.call(text, c(list(x = points$x[high_data], y = y_data[high_data], labels = high_data), textArgs))
  }
  else if (type == 3) { # Scale-Location
    x_data <- fitted
    y_data <- sqrt(abs(residuals))
    args <- modifyList(list(main = "Scale-Location", xlab = "Fitted Values", ylab = expression(sqrt(abs("Standardized Residuals"))),
                            ylim = c(0,max(y_data)*1.1)), args)

    do.call(plot, c(list(x=x_data, y = y_data), args))
    do.call(abline, c(list(h = 0), ablineArgs))

    high_data <- order(abs(y_data), decreasing = TRUE)[1:(length(y_data)/10)]
    if (length(high_data) > 0)
      do.call(text, c(list(x = x_data[high_data], y = y_data[high_data], labels = high_data), textArgs))
  }
  else if (type == 4) {
    x_data <- as.numeric(cooks.distance0(x, equation = equation))
    args <- modifyList(list(main = "Cook's distance", xlab = "Observation Number", ylab = "Cook's distance"), args)
    midpoints <- do.call(barplot, c(list(x_data), args))

    high_data <- x_data > 7/length(x_data)
    if (length(high_data) > 0)
      do.call(text, c(list(x = midpoints, y = x_data, labels = ifelse(high_data, 1:length(x_data), "")), textArgs))
  }
  else if (type == 5) {
    x_data <- hatvalues0(x, equation = equation)
    y_data <- residuals
    args <- modifyList(list(main = "Standardized Residuals vs Leverage", xlab = "Leverage", ylab = "Standardized Residuals"), args)

    do.call(plot, c(list(x=x_data, y = y_data), args))
    do.call(abline, c(list(h = 0), ablineArgs))

    high_data <- order(abs(y_data), decreasing = TRUE)[1:(length(y_data)/10)]
    if (length(high_data) > 0)
      do.call(text, c(list(x = x_data[high_data], y = y_data[high_data], labels = high_data), textArgs))
  }
  else if (type == 6) {
    x_data <- hatvalues0(x, equation = equation)
    y_data <- cooks.distance0(x, equation = equation)
    args <- modifyList(list(main = "Cook's dist vs Leverage/(1-Leverage)", xlab = "Leverage/(1-Leverage)", ylab = "Cook's distance"), args)

    do.call(plot, c(list(x=x_data, y = y_data), args))
    do.call(abline, c(list(h = 0), ablineArgs))

    high_y <- 7 / length(y_data)
    high_x <- 2 * mean(x_data)
    high_data <- which(y_data > high_y & x_data > high_x)
    if (length(high_data) > 0)
      do.call(text, c(list(x = x_data[high_data], y = y_data[high_data], labels = ifelse(high_data, 1:length(x_data), "")), textArgs))
  }
  else {
    print("Invalid 'type' argument.")
  }


  return(list(x = x_data, y = y_data))
}

