
#' Plot Estimation Results
#'
#' This function plots estimated coefficients, which may be point estimates, bounds, intervals, or distributions.
#'
#' @param points A list where each element is a point estimate to be drawn
#' as a shape. Each element must also be a list with the following elements:
#' 1.\code{value}: estimated coefficient,
#' 2.\code{y}: vertical position (default=0),
#' 3.\code{pch}: \code{pch} of the point (default="1"), and other usual attributes.
#' @param bounds A list where each element is bound estimation information (e.g.
#' output of an extreme bound analysis) to be drawn as a rectangle. Each element is defined by the following items:
#' 1.\code{xmin}: where bound starts on x-axis,
#' 2.\code{xmax}: where bound ends on x-axis,
#' 3.\code{ymin}: where bound starts on y-axis (default=-0.1),
#' 4.\code{ymax}: where bound ends on y-axis (default=+0.1),
#' 5.\code{alpha}: alpha for \code{rgb} color, and other usual attributes.
#' @param intervals A list where each element is interval estimation information. Each element is a
#' list similar to \code{bounds} but with a value to be drawn at the middle of the interval. In other words,
#' each element is defined by 1.\code{value}: center of the interval, 2.\code{xmin}: see \code{bound}, 3.\code{xmax}: see \code{bound}, etc.
#' @param distributions A list where each element is a distribution estimate to be drawn
#' by its density function. Each element is defined by a
#' 1.\code{type} which is the name of the distribution. Other elements of the list depend on the type of distribution.
#' For \code{type=normal}, elements must be 2.\code{mean}, 3.\code{var}, 4.\code{sdMultiplier}. For
#' \code{type=GLD}, they can be 2.\code{p1},..., 5.\code{p4}, 6.\code{quantiles}. For \code{type==cdfs}
#' they can be 2.\code{xs}, 3.\code{cdfs}, 4.\code{smoothFun}.
#' @param newPlot If \code{TRUE}, initializes a new plot.
#' @param xlim Array of size two for x-axis limits. If \code{NULL}, auto-generated.
#' @param ylim Array of size two for y-axis limits. If \code{NULL}, auto-generated.
#' @param boundFun Function to control \code{xlim} and \code{ylim} in the plot, with arguments being the computed bounds.
#' @param legendsTitle List of titles for legends.
#' @param legendSize Numeric value for size of legend (width or height) in lines of text (passed to \code{oma}).
#'
#' @param ... Additional properties for plot or legend, such as \code{xlab} and \code{ylab}.
#'
#' @return This function does not return any value.
#' @export
#' @importFrom graphics legend par
#' @importFrom grDevices rgb
#' @importFrom stats dnorm
#' @examples
#' points <- list()
#' points$one <- list(value = 1, label = "Point 1")
#' points$two <- list(value = 2, label = "Point 2", col = "red", pch = 22, cex = 4)
#' coefs.plot(points = points)
#'
#' bounds <- list()
#' bounds$one <- list(xmin = -1, xmax = 0.5, label = "Bound 1")
#' bounds$two <- list(
#'   xmin = 0, xmax = 1, ymin = 0.2, ymax = 0.3,
#'   label = "Bound 2", alpha = 0.2, col = rgb(0, 0, 1.0, alpha = 0.3)
#' )
#' coefs.plot(points = points, bounds = bounds)
#'
#' intervals <- list()
#' intervals$one <- list(value = 2, xmin = 0, xmax = 3, label = "Interval 1")
#' intervals$two <- list(
#'   value = 1.5, xmin = 1, xmax = 2, y = 4,
#'   label = "Interval 2", col = "blue", lwd = 3, pch = 11, cex = c(1.2, 3, 1.2)
#' )
#' coefs.plot(points = points, bounds = bounds, intervals = intervals)
#'
#' distributions <- list()
#' distributions$one <- list(type = "normal", mean = 0, var = 1, label = "Distribution 1")
#' distributions$two <- list(
#'   type = "gld", p1 = 0, p2 = 1.5, p3 = 1.2,
#'   p4 = 1.2, label = "Distribution 2", col = "blue", lwd = 3
#' )
#' distributions$three <- list(
#'   type = "cdfs", xs = seq(-2, 2, 0.1),
#'   cdfs = pnorm(seq(-2, 2, 0.1)), label = "Distribution 3",
#'   col = rgb(1, 0, 0, alpha = 0.5), lwd = 8
#' )
#' coefs.plot(
#'   points = points, bounds = bounds, intervals = intervals,
#'   distributions = distributions, legendsTitle = NULL, legendSize = 7
#' )
#'
coefs.plot <- function(points = NULL, bounds = NULL, intervals = NULL, distributions = NULL,
                       newPlot = TRUE, xlim = NULL, ylim = NULL,
                       boundFun = function(b, type) ifelse(type == "xmin" || type == "ymin",0.9 * b,1.1 * b),
                       legendsTitle = c("Point", "Bound", "Interval", "Density"),
                       legendSize = 5, ...) {

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  dots <- list(...)

  def_bound_ymin <- -0.1
  def_bound_ymax <- 0.1
  numPoints <- 100


  points <- as.list(points)
  bounds <- as.list(bounds)
  intervals <- as.list(intervals)
  distributions <- as.list(distributions)

  if (is.null(legendsTitle)) {
    legendsTitle <- c(NULL, NULL, NULL, NULL)
  }

  dists <- list()
  if (length(distributions) > 0) {
    for (g in distributions) {
      type <- g$type
      if (is.null(type)) {
        stop("missing: type of the distribution, e.g., normal, gld, cdfs.")
      }

      if (type == "normal") {
        sd <- sqrt(g$var)
        sdint <- !is.null(g$sdmultiplier, 2) * sd
        x <- seq(g$mean - sdint, g$mean + sdint, length.out = numPoints)
        dists[[length(dists) + 1]] <- list(x = x, y = dnorm(x, g$mean, sd))
      } else if (type == "gld") {
        quantiles <- !is.null(g$quantiles, seq(0.01, 0.99, length.out = numPoints))
        dists[[length(dists) + 1]] <- list(
          x = s.gld.quantile(quantiles, g$p1, g$p2, g$p3, g$p4),
          y = s.gld.density.quantile(quantiles, g$p1, g$p2, g$p3, g$p4)
        )
      } else if (type == "cdfs") {
        y <- diff(g$cdfs) / diff(g$xs)
        if (is.null(g$smoothFun) == FALSE) {
          y <- g$smoothFun(y)
        }

        dists[[length(dists) + 1]] <- list(
          x = g$xs,
          y = c(y, NA)
        )
      } else {
        stop("not supported distribution.")
      }
    }
  }



  if (is.null(xlim)) {
    xlim <- numeric(2) * NA
  }
  if (is.null(ylim)) {
    ylim <- numeric(2) * NA
  }
  if (is.na(xlim[[1]])) { # TODO: for normal, gld, first get the data
    xlim[[1]] <- boundFun(min(as.numeric(c(
      sapply(points, function(p) p$value),
      sapply(bounds, function(p) p$xmin),
      sapply(intervals, function(p) p$xmin),
      sapply(dists, function(p) min(p$x, na.rm = TRUE))
    )), na.rm = TRUE), "xmin")
  }
  if (is.na(xlim[[2]])) {
    xlim[[2]] <- boundFun(max(as.numeric(c(
      sapply(points, function(p) p$value),
      sapply(bounds, function(p) p$xmax),
      sapply(intervals, function(p) p$xmax),
      sapply(dists, function(p) max(p$x, na.rm = TRUE))
    )), na.rm = TRUE), "xmax")
  }
  if (is.na(ylim[[1]])) {
    ylim[[1]] <- boundFun(min(as.numeric(c(
      sapply(points, function(p) !is.null(p$y, 0)),
      sapply(bounds, function(p) !is.null(p$ymin, def_bound_ymin)),
      sapply(intervals, function(p) !is.null(p$y, (def_bound_ymin + def_bound_ymax) / 2)),
      sapply(dists, function(p) min(p$y, na.rm = TRUE))
    )), na.rm = TRUE), "ymin")
  }

  if (is.na(ylim[[2]])) {
    ylim[[2]] <- boundFun(max(as.numeric(c(
      sapply(points, function(p) !is.null(p$y, 0)),
      sapply(bounds, function(p) !is.null(p$ymax, def_bound_ymax)),
      sapply(intervals, function(p) !is.null(p$y, (def_bound_ymin + def_bound_ymax) / 2)),
      sapply(dists, function(p) max(p$y, na.rm = TRUE))
    )), na.rm = TRUE), "ymax")
  }

  legendPos <- "right" # TODO: other options

  xjust <- 0
  yjust <- 1
  leg_x_f <- NULL
  leg_y_f <- NULL
  if (newPlot) {
    L <- list(rect = list(top = ylim[[2]], left = par("usr")[2], h = 0))
    if (legendPos == "bottom") {
      par(oma = c(legendSize, 0, 0, 0))
    } else if (legendPos == "left") {
      par(oma = c(0, legendSize, 0, 0))

      leg_x_f <- function(p, rect) p[1] - 3
      leg_y_f <- function(p, rect) {
        rect$top - rect$h
      }
    } else if (legendPos == "top") {
      par(oma = c(0, 0, legendSize, 0))
    } else if (legendPos == "right") {
      par(oma = c(0, 0, 0, legendSize))
      leg_x_f <- function(p, rect) p[2]
      leg_y_f <- function(p, rect) rect$top - rect$h
    }


    plot(NULL, xlim = xlim, ylim = ylim, ylab = dots$ylab, xlab = dots$xlab)
  }


  # draw points
  if (length(points) > 0) {
    lgn_lgn <- character(length(points))
    lgn_col <- character(length(points))
    lgn_pch <- integer(length(points))
    lgn_cex <- integer(length(points))

    i <- 0
    for (g in points) {
      i <- i + 1
      pch <- !is.null(g$pch, 4)
      col <- !is.null(g$col, "black")
      cex <- !is.null(g$cex, 3)
      graphics::points(
        x = g$value, y = !is.null(g$y, 0), type = "p",
        pch = pch, col = col, cex = cex
      )
      lgn_lgn[[i]] <- g$label
      lgn_col[[i]] <- col
      lgn_pch[[i]] <- pch
      lgn_cex[[i]] <- cex
    }
    L <- legend(leg_x_f(par("usr"), L$rect), leg_y_f(par("usr"), L$rect),
                xpd = NA, xjust = xjust, yjust = yjust,
                legend = lgn_lgn, bty = "n", title = legendsTitle[[1]],
                col = lgn_col, pch = lgn_pch, x.intersp = 1.4,
                pt.cex = 1.5, fill = NA, border = NA
    ) #
  }

  if (length(bounds) > 0) {
    lgn_lgn <- character(length(bounds))
    lgn_fill <- character(length(bounds))

    i <- 0
    for (g in bounds) {
      i <- i + 1

      col <- !is.null(g$col, rgb(0, 0, 0.0, alpha = 0.2))
      density <- !is.null(g$density, NULL)
      border <- !is.null(g$border, NA)
      graphics::rect(
        xleft = g$xmin, xright = g$xmax,
        ybottom = !is.null(g$ymin, def_bound_ymin),
        ytop = !is.null(g$ymax, def_bound_ymax), density = density,
        col = col, border = border
      )
      lgn_lgn[[i]] <- g$label
      lgn_fill[[i]] <- col
    }

    L <- legend(leg_x_f(par("usr"), L$rect), leg_y_f(par("usr"), L$rect),
                xpd = NA, xjust = xjust, yjust = yjust,
                legend = lgn_lgn, bty = "n", title = legendsTitle[[2]],
                col = lgn_col,
                fill = lgn_fill, border = lgn_fill,
                pch = NA, x.intersp = 1.4
    )
  }

  if (length(intervals) > 0) {
    lgn_lgn <- character(length(intervals))
    lgn_col <- character(length(intervals))
    lgn_lwd <- integer(length(intervals))
    lgn_lty <- integer(length(intervals))
    lgn_pch <- integer(length(intervals))

    i <- 0
    for (g in intervals) {
      i <- i + 1

      col <- !is.null(g$col, "black")
      lty <- !is.null(g$lty, 1)
      lwd <- !is.null(g$lwy, 1)
      pch <- !is.null(g$pch, 8)
      y <- !is.null(g$y, 0)
      cex <- !is.null(g$cex, c(1, 1, 1))

      graphics::points(
        x = c(g$xmin, g$value, g$xmax), y = c(y, y, y), type = "b",
        col = col, lwd = lwd, lty = lty, pch = c(15, pch, 15), cex = cex
      )

      lgn_lgn[[i]] <- g$label
      lgn_col[[i]] <- col
      lgn_lty[[i]] <- lty
      lgn_lwd[[i]] <- lwd
      lgn_pch[[i]] <- pch
    }

    L <- legend(leg_x_f(par("usr"), L$rect), leg_y_f(par("usr"), L$rect),
                xpd = NA, xjust = xjust, yjust = yjust,
                legend = lgn_lgn, bty = "n", title = legendsTitle[[3]],
                col = lgn_col, pch = lgn_pch, lty = lgn_lty, lwd = lgn_lwd
    )
  }

  if (length(distributions) > 0) {
    lgn_lgn <- character(length(distributions))
    lgn_col <- character(length(distributions))
    lgn_lwd <- integer(length(distributions))
    lgn_lty <- integer(length(distributions))

    i <- 0
    for (g in distributions) {
      i <- i + 1
      dist <- dists[[i]]
      col <- !is.null(g$col, "black")
      lty <- !is.null(g$lty, 1)
      lwd <- !is.null(g$lwd, 1)

      graphics::lines(
        x = dist$x, y = dist$y, type = "l",
        col = col, lwd = lwd, lty = lty
      )

      lgn_lgn[[i]] <- g$label
      lgn_col[[i]] <- col
      lgn_lty[[i]] <- lty
      lgn_lwd[[i]] <- lwd
    }

    L <- legend(leg_x_f(par("usr"), L$rect), leg_y_f(par("usr"), L$rect),
                xpd = NA, xjust = xjust, yjust = yjust,
                legend = lgn_lgn, bty = "n", title = legendsTitle[[4]],
                col = lgn_col, lty = lgn_lty, lwd = lgn_lwd,
                pch = NA
    ) # set pch for better placement
  }
}




#' Create a fan plot from a matrix of distribution parameters
#'
#' @param data A matrix with n columns where each row represents a prediction as a probability distribution.
#' @param dist The type of distribution to use: "normal" or "log-normal".
#' @param quantiles A list of quantiles to plot.
#' @param gradient Logical value indicating whether to add a gradient plot from the median to the minimum and maximum \code{quantiles}.
#' @param gradPointShape The shape of the point to use when gradient is TRUE.
#' @param limExpand An array of size 2 where the elements can be \code{NA} or a multiplier to  expand the range of the x and y axes.
#' @param limExpand A multiplier to to expand the range of the y-axis.
#' @param new.plot Logical value indicating whether to create a new plot or use the current plot.
#' @param boundColor The color to use for the bounds around the median.
#' @param midColor The color to use for the median or points.
#' @param ... Other arguments to pass to \code{plot} function if \code{new.plot} is \code{TRUE}.
#'
#' @examples
#' data <- matrix(c(0, 1, 1, 4), ncol = 2)
#' rownames(data) <- c("A", "B")
#' fan.plot(data)
#' fan.plot(data, dist = "log-normal")
#' fan.plot(data, gradient = TRUE)
#'  fan.plot(data, gradient = TRUE)
#' fan.plot(data, gradient = FALSE, new.plot = FALSE,
#'         boundColor = adjustcolor("red", alpha.f = 0.2),
#'         midColor = NA)
#'
#' @export
#' @importFrom graphics points lines polygon axis
#' @importFrom grDevices colorRampPalette
#' @importFrom stats qlnorm qnorm
fan.plot <- function(data, dist = "normal",
                     quantiles = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95),
                     gradient = FALSE, gradPointShape = 19,
                     limExpand = c(0.1,0.1), new.plot = TRUE,
                     boundColor = "blue", midColor = "black", ...) {
  if (dist == "normal") {
    if (ncol(data) != 2) {
      stop("For normal distribution, the matrix must have two columns: mean and variance.")
    }
    means <- data[, 1]
    sds <- sqrt(data[, 2])
    medians <- means
    quants <- t(apply(data, 1, function(x) qnorm(quantiles, x[1], sqrt(x[2]))))

  } else if (dist == "log-normal") {
    if (ncol(data) != 2) {
      stop("For log-normal distribution, the matrix must have two columns: meanlog and sdlog.")
    }
    meanlogs <- data[, 1]
    sdlogs <- data[, 2]
    medians <- exp(meanlogs + sdlogs^2 / 2)
    quants <- t(apply(data, 1, function(x) qlnorm(quantiles, x[1], x[2])))

  } else {
    stop(paste("Distribution", dist, "is not supported."))
  }

  quants <- cbind(quants[, quantiles <= 0.5], medians, quants[, quantiles > 0.5])

  ymin <- min(quants)
  ymax <- max(quants)

  if (length(limExpand) != 2)
    stop("Invalid 'limExpand'. An array of size 2 is expected")
  if ( is.na(limExpand[[2]]) == FALSE && limExpand[[2]] > 0) {
    yrange <- ymax - ymin
    ymin <- ymin - limExpand[[2]] * yrange
    ymax <- ymax + limExpand[[2]] * yrange
  }

  xmin = 1
  xmax = length(medians)
  if ( is.na(limExpand[[1]]) == FALSE && limExpand[[1]] > 0) {
    xrange <- xmax - xmin
    xmin <- xmin - limExpand[[1]] * xrange
    xmax <- xmax + limExpand[[1]] * xrange
  }

  cols <- colorRampPalette(c("white", boundColor, "white"), alpha = TRUE)(ncol(quants) + 1)
  cols <- cols[2:(length(cols)-1)] # remove whites

  if (new.plot) {
    plot(medians, col = midColor,
         ylim = c(ymin, ymax),
         xlim = c(xmin, xmax),
         ...)
    axis(1, at = seq_along(medians), labels = rownames(data))
  }

  if (!gradient) {
    for (i in seq_len(ncol(quants) - 1)) {
      polygon(c(seq_along(medians), rev(seq_along(medians))),
              c(quants[, i], rev(quants[, i + 1])),
              col = cols[i], border = NA)
    }

    lines(medians, col = midColor)

  } else {
    gradient_cols <- colorRampPalette(c(boundColor, "white"), alpha = TRUE)(100)

    for (i in seq_along(medians)) {
      x <- c(i - 0.4, i + 0.4)
      y1 <- seq(medians[i], quants[i, ncol(quants)], length.out = 100)
      y2 <- seq(medians[i], quants[i, 1], length.out = 100)

      for (j in seq_along(y1)) {
        polygon(c(x[1], x[2], x[2], x[1]), c(y1[j], y1[j], y1[j + 1], y1[j + 1]), col = gradient_cols[j], border = NA)
        polygon(c(x[1], x[2], x[2], x[1]), c(y2[j], y2[j], y2[j + 1], y2[j + 1]), col = gradient_cols[j], border = NA)
      }
    }

    points(seq_along(medians), medians, pch = gradPointShape, col = midColor)

  }
}


