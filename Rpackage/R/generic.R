
#' Prints the Output of a Search Process
#'
#' @param x An output from one of the \code{search.?} functions (see \code{search.sur}, \code{search.varma}, or \code{search.dc}).
#' @param ... Additional arguments
#'
#' @return This function has no output.
#' @export
print.ldtsearch <- function(x, ...) {
  if (is.null(x)) {
    stop("argument is null.")
  }
  if (any(class(x) == "ldtsearch") == FALSE) {
    stop("Invalid class")
  }

  cat("method:", attr(x, "method"), "\n")
  cat("expected: ", prettyNum(x$counts$expectedCount, big.mark = ","), ", searched: ",
    prettyNum(x$counts$searchedCount, big.mark = ","), " (",
    round(x$counts$searchedCount / x$counts$expectedCount * 100,1), "%)",
    ", failed: ", prettyNum(x$counts$failedCount, big.mark = ","), " (",
    round(x$counts$failedCount / x$counts$searchedCount * 100,1), "%)\n",
    sep = ""
  )
  cat("--------\n")
  if (x$counts$failedCount > 0){
    cat("Failures:\n")
    i <- 0
    for (fd in x$counts$failedDetails){
      i <- i + 1
      cat(i,". ",fd$message, ": ", fd$count
          ," (", round((fd$count/x$counts$failedCount)*100,1) ,"%)\n", sep = "")
    }
    cat("--------\n")
  }

  nms <- names(x)
  i <- 0
  for (m in x) {
    i <- i + 1
    if (i == 1 || i > length(x) - 1 || is.null(m)) {
      next()
    }
    cat(i - 1, ". ", nms[[i]], ":\n", sep = "")
    j <- 0
    for (t in m) {
      j <- j + 1
      if (is.null(t)) {
        next()
      }
      cat(" ", t$name, " ", sep = "")
      if (is.null(t$model) || is.null(t$model$bests) || length(t$model$bests) > 0) {
        cat("(best=", round(s.measure.from.weight(t$model$bests[[1]]$weight, nms[[i]]), 3), ")\n", sep = "")
      } else {
        cat("(best model is missing)\n", sep = "")
      }
    }
  }
  return(NULL)
}


#' Gets Estimation from Search Result
#'
#' This function retrieves the estimation result from the output of a search process.
#' It uses the search options to re-estimate the same model estimated in the search process,
#' given a list of required indices.
#'
#' @param searchRes Output from one of the \code{search.?} functions (see [search.sur], [search.varma], or [search.dc]).
#' @param endoIndices Endogenous indices corresponding to the columns of the dependent variable.
#' @param exoIndices Exogenous indices corresponding to the columns of the exogenous variable.
#' @param y Data for dependent variables.
#' @param x Data for exogenous variables.
#' @param printMsg Argument passed to estimation methods
#' @param ... Additional arguments
#'
#' @return Estimation result similar to output of [estim.sur], [estim.varma], or [estim.dc].
#' @export
#' @seealso [estim.sur], [estim.varma], [estim.dc]
h.get.estim <- function(searchRes, endoIndices, exoIndices, y, x, printMsg, ...) {
  if (is.null(endoIndices) == FALSE && is.null(y) == FALSE) {
    y <- as.matrix(y)
  }
  if (is.null(exoIndices) == FALSE && is.null(x) == FALSE) {
    x <- as.matrix(x)
  }

  method <- attr(searchRes, "method")
  # tryCatch({
  if (method == "sur") {
    GetEstim_sur(searchRes, endoIndices, exoIndices, y, x, printMsg = printMsg, ...)
  } else if (method == "dc") {
    GetEstim_dc(searchRes, endoIndices, exoIndices, y, x, printMsg = printMsg, ...)
  } else if (method == "varma") {
    GetEstim_varma(searchRes, endoIndices, exoIndices, y, x, printMsg = printMsg, ...)
  } else {
    stop("not implemented method")
  }
  # },error = function(e) e)
}

getMeasureFrom <- function(model, measureName, tarIndex, method) {
  rowNames <- rownames(model$measures)
  ind <- which(rowNames == measureName)
  if (length(ind) == 0) {
    stop(paste0("measure not found (method = sur, measure = ", measureName, ")"))
  }
  r <- model$measures[ind[[1]], tarIndex]
  r
}

#' Summarizes Model Search Output
#'
#' This function estimates the reported models in the output of a \code{search.?} function and provides additional information.
#'
#' @param object Output from one of the \code{search.?} functions (see [search.sur], [search.varma], or [search.dc]).
#' @param y Data for dependent variables (Note that data is not saved in \code{object}).
#' @param x Data for exogenous variables (Note that data is not saved in \code{object}).
#' @param addModelBests If \code{TRUE} and search output contains best models, estimates them.
#' @param addModelAll If \code{TRUE} and search output contains all models, estimates them.
#' @param addItem1Bests If \code{TRUE} and search output contains \code{item1} part, estimates them.
#' @param printMsg If \code{TRUE}, prints details.
#' @param w Weight of observations in discrete choice estimation (Note that data is not saved in \code{object}).
#' @param newX New exogenous data in VARMA estimation (Note that data is not saved in \code{object}).
#' @param test If \code{TRUE}, verifies everything is working correctly. Please report errors.
#' @param ... Additional arguments
#'
#' @details
#' The output of \code{search.?} functions in this package only contains the information required
#' to re-estimate the models, not the actual estimation results. You can use [h.get.estim] to
#' get the estimation of specific indices. This function re-estimates everything in the \code{search.?}
#' output with a similar structure.
#'
#' @return A list with estimated models and other information, structured similarly to the given \code{search.?} output.
#'
#' @export
summary.ldtsearch <- function(object, y, x = NULL, addModelBests = TRUE,
                              addModelAll = FALSE, addItem1Bests = FALSE,
                              printMsg = FALSE, w = NULL, newX = NULL,
                              test = FALSE, ...) {
  test_perc <- 1e-8

  if (is.null(object)) {
    stop("argument is null.")
  }
  if (any(class(object) == "ldtsearch") == FALSE) {
    stop("Invalid class")
  }


  result <- list()
  result$method <- attr(object, "method")
  result$SearchPerc <- object$counts$searchedCount / object$counts$expectedCount * 100
  result$FailedPerc <- object$counts$failedCount / object$counts$searchedCount * 100
  result$MeasureNames <-
    names(object)[2:(length(object) - 1)] # the last element is the inputs

  i <- 0
  for (mea in result$MeasureNames) {
    i <- i + 1
    x_mea <- object[[mea]]
    targets <- names(x_mea)
    if (is.null(result$TargetNames)) {
      result$TargetNames <- c()
      for (targ in targets) {
        result$TargetNames <- append(result$TargetNames, x_mea[[targ]]$name)
      }
    }

    j <- 0
    for (targ in targets) {
      j <- j + 1
      x_tar <- x_mea[[targ]]
      x_names <- names(x_tar)


      x_mod <- x_tar$model
      if (is.null(x_mod) == FALSE) {
        if (addModelBests) {
          if (is.null(x_mod$bests) == FALSE) {
            k <- 0
            for (b in x_mod$bests) {
              k <- k + 1
              su_m <- h.get.estim(object, b$depIndices, b$exoIndices, y, x,
                printMsg = printMsg,
                params = b$parameters, newX = newX, w = w,
                distType = if (is.integer(b$dist) && b$dist == 0) "logit" else "probit"
              )
              if (test) {
                jt <- if (is.null(b$depIndices) || length(b$depIndices) == 0) {
                  1
                } else {
                  which(b$depIndices == j)[[1]]
                } # index of target in endogenous data
                testthat::expect_equal(
                  s.weight.from.measure(
                    getMeasureFrom(su_m, mea, jt, result$method), mea
                  ),
                  b$weight,
                  tolerance = test_perc
                )
              }
              result[[mea]][[targ]][[paste0("model")]][[paste0("bests")]][[paste0("best", k)]] <- su_m # nolint
            }
          }
        }

        if (addModelAll) {
          if (is.null(x_mod$all) == FALSE) {
            k <- 0
            for (b in x_mod$all) {
              k <- k + 1

              su_m <- h.get.estim(object, b$depIndices, b$exoIndices, y, x,
                printMsg = printMsg,
                params = b$parameters, newX = newX, w = w,
                distType = if (is.integer(b$dist) && b$dist == 0) "logit" else "probit"
              )
              if (test) {
                jt <- if (is.null(b$depIndices) || length(b$depIndices) == 0) {
                  1
                } else {
                  which(b$depIndices == j)[[1]]
                }
                testthat::expect_equal(
                  s.weight.from.measure(getMeasureFrom(su_m, mea, jt, result$method), mea),
                  b$weight,
                  tolerance = test_perc
                )
              }
              result[[mea]][[targ]][[paste0("model")]][[paste0("all")]][[paste0("model", k)]] <- su_m # nolint
            }
          }
        }
      }

      p <- if (is.null(x_mod)) 2 else 3
      if (length(x_names) >= p) {
        type1_name <- x_names[[p]]
        x_type1 <- x_tar[[type1_name]]
        if (is.null(x_type1) == FALSE) {
          if (addItem1Bests) {
            if (is.null(x_type1$bests) == FALSE) {
              items_names <- names(x_type1$bests)

              n <- 0
              for (item in x_type1$bests) {
                n <- n + 1
                k <- 0
                for (b in item) {
                  if (is.null(b)) {
                    next
                  }

                  k <- k + 1
                  if (k != 1) { # first one is name (note k-1 in the next line)

                    su_m <-
                      h.get.estim(object, b$depIndices, b$exoIndices, y, x,
                        printMsg = printMsg,
                        params = b$parameters, newX = newX, w = w,
                        distType = if (is.integer(b$dist) && b$dist == 0) "logit" else "probit"
                      )

                    if (test) {
                      jt <- if (is.null(b$depIndices) || length(b$depIndices) == 0) {
                        1
                      } else {
                        which(b$depIndices == j)[[1]]
                      }
                      testthat::expect_equal(
                        s.weight.from.measure(
                          getMeasureFrom(su_m, mea, jt, result$method), mea
                        ), b$weight,
                        tolerance = test_perc
                      )

                      if (result$method == "sur" || result$method == "dc") {
                        ind <- which(item$name == rownames(su_m$estimations$coefs))
                        testthat::expect_equal(su_m$estimations$coefs[ind, jt], b$mean,
                          tolerance = 1e-14
                        )
                        testthat::expect_equal(su_m$estimations$stds[ind, jt], sqrt(b$var),
                          tolerance = 1e-14
                        )
                      } else if (result$method == "varma") {
                        # item$name is Horizon1,...
                        ind <- as.integer(substr(item$name, 8, 8)) + su_m$prediction$startIndex - 1
                        testthat::expect_equal(as.numeric(su_m$prediction$means[jt, ind]),
                          as.numeric(b$mean),
                          tolerance = 1e-14
                        )
                        testthat::expect_equal(as.numeric(su_m$prediction$vars[jt, ind]),
                          as.numeric(b$var),
                          tolerance = 1e-14
                        )
                      }
                    }

                    result[[mea]][[targ]][[type1_name]][[paste0("bests")]][[items_names[[n]]]][[paste0("best", k - 1)]] <- su_m # nolint
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  class(result) <- c("summary.ldtsearch", "list")
  return(result)
}


#' Get Data from Model Search Output
#'
#' This function converts a specific part of the output from one of the \code{search.?} functions into a \code{data.frame}.
#'
#' @param x Output from one of the \code{search.?} functions (see [search.sur], [search.varma], or [search.dc]).
#' @param types One or more types of information to include in the data frame.
#' Can be \code{bestweights}, \code{allweights}, \code{inclusion}, \code{type1bests}, \code{cdf}, \code{extremebounds}, and/or \code{mixture}.
#' @param measures Indices or names of measures to use.
#' @param targets Indices or names of targets to use.
#' @param rows Indices or names of rows to use. If the requested object is a matrix
#' (or an array), determines the rows and cannot be \code{NULL}. For \code{type1bests}, this is the name of the variables.
#' @param columns Indices or names of columns to use. If the requested object is a matrix, determines the columns and cannot be \code{NULL}. For
#' \code{type1bests}, this is the name of the fields: \code{weight}, \code{mean}, and/or \code{var}.
#' @param itemIndices Indices of items such as \code{bests} to use.
#' @param colNameFun Function to determine column names. Argument is a list of names, i.e., one of the following items: \code{target},
#' \code{measure}, \code{row}, \code{column}, \code{item}. If \code{NULL}, uses \code{paste} function.
#' @param rowContent Character string determining type of information in rows of returned
#' data frame. Some items are not available for some \code{types}.
#' Use \code{row} for variables in rows of matrices such as \code{inclusion} or
#' \code{mixture}. Use \code{column} for columns of such matrices.
#' Use \code{item} for best models or models in the all field.
#' @param cdfIndex Integer for index of CDF if type is cdf
#' @param ... Additional arguments
#'
#' @details
#' There are five types of indices in this function: measures,
#' targets, bests, type1's items, equations.
#' Use \code{NULL} to use all available information or specify them.
#'
#'
#' @return A data frame containing requested data.
#'
#' @export
to.data.frame <- function(x, types = c(
                                      "bestweights", "allweights", "inclusion",
                                      "type1bests", "cdf", "extremebounds", "mixture"
                                    ),
                                    measures = NULL, targets = NULL,
                                    rows = NULL, columns = NULL, itemIndices = NULL,
                                    colNameFun = NULL,
                                    rowContent = c("measure", "target", "item", "row", "column"),
                                    cdfIndex = 0, ...) {
  if (is.null(x)) {
    stop("argument is null.")
  }
  if (is.null(colNameFun)){
    colNameFun <- function(ns) {
      colNameFun <- paste(ns[lengths(ns) > 0],  collapse = ".")
    }
  }

  getValueInMat <- function(mat, rowName, colName) {
    tryCatch(mat[[which(rownames(mat) == rowName), which(colnames(mat) == colName)]],
      error = function(e) NA
    )
  }

  c_fmi <- 2 # first measure index
  method <- attr(x, "method")

  type1name <- if (method == "sur") "coefs" else stop("Not implemented")

  if (is.vector(measures) == FALSE) measures <- c(measures)
  if (is.vector(targets) == FALSE) measures <- c(targets)
  if (is.vector(rows) == FALSE) measures <- c(rows)
  if (is.vector(columns) == FALSE) measures <- c(columns)
  if (is.vector(itemIndices) == FALSE) measures <- c(itemIndices)

  m_names <- names(x)
  t_names <- sapply(c(1:x$info$numTargets), function(j) x[[c_fmi]][[j]]$name)
  b_names <- paste0("best", itemIndices)

  # find indexes
  if (is.null(measures)) {
    measures <- c(c_fmi:(length(x) - 1))
  } else {
    if (is.numeric(measures)) {
      measures <- as.integer(measures) + c_fmi - 1
    } else {
      measures <- unlist(sapply(measures, function(m) which(m_names == tolower(m))))
    }
  }
  if (length(measures) == 0) {
    stop("no valid 'measure' is found.")
  }

  if (is.null(targets)) {
    targets <- c(1:x$info$numTargets)
  } else {
    if (is.numeric(targets)) {
      targets <- as.integer(targets)
    } else {
      targets <- unlist(sapply(targets, function(m) which(t_names == m)))
    }
  }

  if (length(targets) == 0) {
    stop("no valid 'target' is found.")
  }

  if (is.null(itemIndices)) {
    itemIndices <- c(1:x$info$searchItems$bestK)
  }

  m_names <- m_names[measures]
  t_names <- t_names[targets]

  types <- c(types)
  types <- sapply(types, function(ty) {
    tolower(match.arg(
      ty,
      c("bestweights", "allweights", "inclusion", "type1bests", "cdf", "extremebounds", "mixture")
    ))
  })
  if (is.vector(rowContent)) {
    rowContent <- rowContent[[1]]
  }
  rowContent <- match.arg(rowContent, c("measure", "target", "item", "row", "column"))

  if ((rowContent == "item" || rowContent == "column" || rowContent == "row") &&
    length(types) > 1) {
    warning("There are more than 1 element in 'types' and for this 'rowContent'
    the result (if no error occurs) might be misleading (row names are not generally the same).")
  }



  ws <- data.frame()

  for (type in types) {
    # for distinguishing type in multiple lists:
    if (length(types) > 1) {
      stype <- if (type == "bestweights") {
        "bw"
      } else if (type == "allweights") {
        "aw"
      } else if (type == "inclusion") {
        "inc"
      } else if (type == "type1bests") {
        "t1b"
      } else if (type == "cdf") {
        "cdf"
      } else if (type == "extremebounds") {
        "eb"
      } else if (type == "mixture") {
        "mix"
      } else {
        stop("Invalid type")
      }
    } else {
      stype <- NULL
    }

    if (type == "bestweights") {
      if (rowContent == "measure") {
        ws0 <- as.data.frame(lapply(itemIndices, function(i) {
          sapply(targets, function(j) {
            sapply(measures, function(m) x[[m]][[j]]$model$bests[[i]]$weight)
          })
        }))
        row.names(ws0) <- m_names
        colnames(ws0) <- c(sapply(
          paste0("best", itemIndices),
          function(b) sapply(t_names, function(t) colNameFun(list(type = stype, target = t, item = b)))
        ))
      } else if (rowContent == "target") {
        ws0 <- as.data.frame(lapply(itemIndices, function(i) {
          sapply(measures, function(m) {
            sapply(targets, function(j) x[[m]][[j]]$model$bests[[i]]$weight)
          })
        }))
        row.names(ws0) <- t_names
        colnames(ws0) <- c(sapply(
          paste0("best", itemIndices),
          function(b) {
            sapply(m_names, function(m) {
              colNameFun(list(
                type = stype,
                measure = m, item = b
              ))
            })
          }
        ))
      } else if (rowContent == "item") {
        ws0 <- as.data.frame(lapply(targets, function(j) {
          sapply(measures, function(m) {
            sapply(itemIndices, function(i) x[[m]][[j]]$model$bests[[i]]$weight)
          })
        }))
        row.names(ws0) <- paste0("best", itemIndices)
        colnames(ws0) <- c(sapply(t_names, function(t) {
          sapply(
            m_names,
            function(m) colNameFun(list(type = stype, measure = m, target = t))
          )
        }))
      } else {
        stop("invalid or not implemented 'rowContent'.")
      }

      ws <- tryCatch(as.data.frame(append(ws, ws0, length(ws)), row.names = row.names(ws0)),
        error = function(e) {
          warning(paste0("Could not append a type of infomation (type=", type, ")"))
          ws
        }
      )
    }

    if (type == "allweights") {
      allcount <- max(as.numeric(sapply(measures, function(m) {
        sapply(
          targets,
          function(t) length(x[[m]][[t]]$model$all)
        )
      })))
      if (allcount > 0) {
        if (rowContent == "measure") {
          ws0 <- as.data.frame(lapply(c(1:allcount), function(i) {
            sapply(targets, function(j) {
              sapply(measures, function(m) {
                tryCatch(x[[m]][[j]]$model$all[[i]]$weight,
                  error = function(e) NA
                )
              })
            })
          }))
          row.names(ws0) <- m_names
          colnames(ws0) <- c(sapply(paste0("model", c(1:allcount)), function(b) {
            sapply(
              t_names,
              function(t) colNameFun(list(type = stype, target = t, item = b))
            )
          }))
        } else if (rowContent == "target") {
          ws0 <- as.data.frame(lapply(c(1:allcount), function(i) {
            sapply(measures, function(m) {
              sapply(targets, function(j) {
                tryCatch(x[[m]][[j]]$model$all[[i]]$weight,
                  error = function(e) NA
                )
              })
            })
          }))
          row.names(ws0) <- t_names
          colnames(ws0) <- c(sapply(paste0("model", c(1:allcount)), function(b) {
            sapply(
              m_names,
              function(m) colNameFun(list(type = stype, measure = m, item = b))
            )
          }))
        } else if (rowContent == "item") {
          ws0 <- as.data.frame(lapply(targets, function(j) {
            sapply(measures, function(m) {
              sapply(c(1:allcount), function(i) {
                tryCatch(x[[m]][[j]]$model$all[[i]]$weight,
                  error = function(e) NA
                )
              })
            })
          }))
          row.names(ws0) <- paste0("model", c(1:allcount))
          colnames(ws0) <- c(sapply(t_names, function(t) {
            sapply(
              m_names,
              function(m) colNameFun(list(type = stype, measure = m, target = t))
            )
          }))
        } else {
          stop("invalid or not implemented 'rowContent'.")
        }

        ws <- tryCatch(as.data.frame(append(ws, ws0, length(ws)), row.names = row.names(ws0)),
          error = function(e) {
            warning(paste0("Could not append a type of infomation (type=", type, ")"))
            ws
          }
        )
      }
    }

    if (type == "type1bests") {
      if (is.null(rows)) { # get unique variable names
        rows <- unique(c(sapply(measures, function(m) {
          sapply(
            targets,
            function(t) sapply(x[[m]][[t]]$coefs$bests, function(c) c$name)
          )
        })))
      } else if (is.numeric(rows)) {
        rows <- tryCatch(
          sapply(measures, function(m) {
            sapply(
              targets,
              function(t) sapply(x[[m]][[t]]$coefs$bests, function(c) c$name)
            )
          })[rows],
          error = function(e) stop("Check the indices in the 'rows'.")
        )
      }

      if (is.null(columns)) { # get unique column.names
        columns <- c("weight", "mean", "var")
      } else if (is.numeric(rows)) {
        columns <- tryCatch(c("weight", "mean", "var")[columns],
          error = function(e) stop("Check the indices in the 'columns'.")
        )
      }


      if (rowContent == "measure") {
        ws0 <- as.data.frame(lapply(rows, function(r) {
          lapply(columns, function(c) {
            lapply(itemIndices, function(b) {
              lapply(targets, function(j) {
                sapply(measures, function(m) {
                  cind <- which(r == sapply(
                    x[[m]][[j]][[type1name]]$bests,
                    function(cb) cb$name
                  ))
                  if (length(cind) > 0) {
                    return(x[[m]][[j]][[type1name]]$bests[[cind[[1]]]][[b + 1]][[c]])
                  }
                  return(NA)
                })
              })
            })
          })
        }))
        row.names(ws0) <- m_names
        colnames(ws0) <- c(sapply(rows, function(r) {
          sapply(
            columns,
            function(c) {
              sapply(
                paste0("best", itemIndices),
                function(b) {
                  sapply(
                    t_names,
                    function(t) {
                      colNameFun(list(
                        type = stype, target = t,
                        column = c, row = r, item = b
                      ))
                    }
                  )
                }
              )
            }
          )
        }))
      } else if (rowContent == "target") {
        ws0 <- as.data.frame(lapply(rows, function(r) {
          lapply(columns, function(c) {
            lapply(itemIndices, function(b) {
              lapply(measures, function(m) {
                sapply(targets, function(j) {
                  cind <- which(r == sapply(
                    x[[m]][[j]][[type1name]]$bests,
                    function(cb) cb$name
                  ))
                  if (length(cind) > 0) {
                    return(x[[m]][[j]][[type1name]]$bests[[cind[[1]]]][[b + 1]][[c]])
                  }
                  return(NA)
                })
              })
            })
          })
        }))
        row.names(ws0) <- t_names
        colnames(ws0) <- c(sapply(rows, function(r) {
          sapply(
            columns,
            function(c) {
              sapply(
                paste0("best", itemIndices),
                function(b) {
                  sapply(
                    m_names,
                    function(m) {
                      colNameFun(list(
                        type = stype, measure = m,
                        column = c, row = r, item = b
                      ))
                    }
                  )
                }
              )
            }
          )
        }))
      } else if (rowContent == "item") {
        ws0 <- as.data.frame(lapply(rows, function(r) {
          lapply(columns, function(c) {
            lapply(measures, function(m) {
              lapply(targets, function(j) {
                sapply(itemIndices, function(b) {
                  cind <- which(r == sapply(
                    x[[m]][[j]][[type1name]]$bests,
                    function(cb) cb$name
                  ))
                  if (length(cind) > 0) {
                    return(x[[m]][[j]][[type1name]]$bests[[cind[[1]]]][[b + 1]][[c]])
                  }
                  return(NA)
                })
              })
            })
          })
        }))
        row.names(ws0) <- paste0("best", itemIndices)
        colnames(ws0) <- c(sapply(rows, function(r) {
          sapply(
            columns,
            function(c) {
              sapply(
                m_names,
                function(m) {
                  sapply(
                    t_names,
                    function(t) {
                      colNameFun(list(
                        type = stype, target = t,
                        column = c, row = r, measure = m
                      ))
                    }
                  )
                }
              )
            }
          )
        }))
      } else if (rowContent == "row") {
        ws0 <- as.data.frame(lapply(measures, function(m) {
          lapply(columns, function(c) {
            lapply(itemIndices, function(b) {
              lapply(targets, function(j) {
                sapply(rows, function(r) {
                  cind <- which(r == sapply(
                    x[[m]][[j]][[type1name]]$bests,
                    function(cb) cb$name
                  ))
                  if (length(cind) > 0) {
                    return(x[[m]][[j]][[type1name]]$bests[[cind[[1]]]][[b + 1]][[c]])
                  }
                  return(NA)
                })
              })
            })
          })
        }))
        row.names(ws0) <- rows
        colnames(ws0) <- c(sapply(m_names, function(m) {
          sapply(
            columns,
            function(c) {
              sapply(
                paste0("best", itemIndices),
                function(b) {
                  sapply(
                    t_names,
                    function(t) colNameFun(list(type = stype, target = t, column = c, measure = m, item = b))
                  )
                }
              )
            }
          )
        }))
      } else if (rowContent == "column") {
        ws0 <- as.data.frame(lapply(rows, function(r) {
          lapply(measures, function(m) {
            lapply(itemIndices, function(b) {
              lapply(targets, function(j) {
                sapply(columns, function(c) {
                  cind <- which(r == sapply(
                    x[[m]][[j]][[type1name]]$bests,
                    function(cb) cb$name
                  ))
                  if (length(cind) > 0) {
                    return(x[[m]][[j]][[type1name]]$bests[[cind[[1]]]][[b + 1]][[c]])
                  }
                  return(NA)
                })
              })
            })
          })
        }))
        row.names(ws0) <- columns
        colnames(ws0) <- c(sapply(rows, function(r) {
          sapply(
            m_names,
            function(m) {
              sapply(
                paste0("best", itemIndices),
                function(b) {
                  sapply(
                    t_names,
                    function(t) {
                      colNameFun(list(
                        type = stype, target = t,
                        measure = m, row = r, item = b
                      ))
                    }
                  )
                }
              )
            }
          )
        }))
      } else {
        stop("invalid or not implemented 'rowContent'.")
      }


      ws <- tryCatch(as.data.frame(append(ws, ws0, length(ws)), row.names = row.names(ws0)),
        error = function(e) {
          warning(paste0("Could not append a type of infomation (type=", type, ")"))
          ws
        }
      )
    }

    if (type == "inclusion" || type == "cdf" || type == "mixture" || type == "extremebounds") { # matrices
      cdfname <- paste0("cdf", cdfIndex)
      getMat <- function(m, t) {
        if (type == "inclusion") {
          x[[m]][[t]]$model$inclusion
        } else if (type == "cdf") {
          x[[m]][[t]][[type1name]]$cdfs[[cdfname]]
        } else if (type == "extremebounds") {
          x[[m]][[t]][[type1name]]$extremeBounds
        } else if (type == "mixture") {
          x[[m]][[t]][[type1name]]$mixture
        }
      }

      if (is.null(rows)) { # get unique row.names
        rows <- unique(c(sapply(
          measures,
          function(m) {
            sapply(
              targets,
              function(t) rownames(getMat(m, t))
            )
          }
        )))
      } else if (is.numeric(rows)) {
        rows <- tryCatch(
          sapply(
            measures,
            function(m) {
              sapply(
                targets,
                function(t) rownames(getMat(m, t))
              )
            }
          )[rows],
          error = function(e) stop("Check the indices in the 'rows'.")
        )
      }

      if (is.null(columns)) { # get unique column.names
        columns <- unique(c(sapply(
          measures,
          function(m) {
            sapply(
              targets,
              function(t) colnames(getMat(m, t))
            )
          }
        )))
      } else if (is.numeric(rows)) {
        columns <- tryCatch(
          sapply(measures, function(m) {
            sapply(
              targets,
              function(t) colnames(getMat(m, t))
            )
          })[columns],
          error = function(e) stop("Check the indices in the 'columns'.")
        )
      }


      if (rowContent == "measure") {
        ws0 <- as.data.frame(lapply(rows, function(r) {
          lapply(columns, function(c) {
            sapply(targets, function(j) {
              sapply(measures, function(m) getValueInMat(getMat(m, j), r, c))
            })
          })
        }))
        row.names(ws0) <- m_names
        colnames(ws0) <- c(sapply(
          rows,
          function(r) {
            sapply(
              columns,
              function(c) {
                sapply(
                  t_names,
                  function(t) colNameFun(list(type = stype, target = t, column = c, row = r))
                )
              }
            )
          }
        ))
      } else if (rowContent == "target") {
        ws0 <- as.data.frame(lapply(rows, function(r) {
          lapply(columns, function(c) {
            sapply(measures, function(m) {
              sapply(targets, function(j) getValueInMat(getMat(m, j), r, c))
            })
          })
        }))
        row.names(ws0) <- t_names
        colnames(ws0) <- c(sapply(
          rows,
          function(r) {
            sapply(
              columns,
              function(c) {
                sapply(
                  m_names,
                  function(m) colNameFun(list(type = stype, measure = m, column = c, row = r))
                )
              }
            )
          }
        ))
      } else if (rowContent == "row") {
        ws0 <- as.data.frame(lapply(measures, function(m) {
          lapply(columns, function(c) {
            sapply(targets, function(j) {
              sapply(rows, function(r) getValueInMat(getMat(m, j), r, c))
            })
          })
        }))
        row.names(ws0) <- rows
        colnames(ws0) <- c(sapply(
          m_names,
          function(m) {
            sapply(
              columns,
              function(c) {
                sapply(
                  t_names,
                  function(t) colNameFun(list(type = stype, target = t, column = c, measure = m))
                )
              }
            )
          }
        ))
      } else if (rowContent == "column") {
        ws0 <- as.data.frame(lapply(rows, function(r) {
          lapply(measures, function(m) {
            sapply(targets, function(j) {
              sapply(columns, function(c) getValueInMat(getMat(m, j), r, c))
            })
          })
        }))
        row.names(ws0) <- columns
        colnames(ws0) <- c(sapply(
          rows,
          function(r) {
            sapply(
              m_names,
              function(m) {
                sapply(
                  t_names,
                  function(t) colNameFun(list(type = stype, target = t, measure = m, row = r))
                )
              }
            )
          }
        ))
      } else {
        stop("invalid or not implemented 'rowContent'.")
      }

      ws <- tryCatch(as.data.frame(append(ws, ws0, length(ws)), row.names = row.names(ws0)),
        error = function(e) {
          warning(paste0("Could not append a type of infomation (type=", type, ")"))
          ws
        }
      )
    }
  }
  return(ws)
}




#' Plot Estimation Results
#'
#' This function plots estimated coefficients, which may be point estimates, bounds, intervals, or distributions.
#'
#' @param points A list where each element is a point estimate to be drawn
#' as a shape. Each element must also be a list with the following elements:
#' 1.\code{value}: estimated coefficient,
#' 2.\code{y}: vertical position (default=0),
#' 3.\code{shape}: shape of the point (default="circle"), and other usual attributes.
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
                       boundFun = function(b, type) {
                         if
                         (type == "xmin" || type == "ymin") {
                           0.9 * b
                         } else {
                           1.1 * b
                         }
                       },
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
        sdint <- if.not.null(g$sdmultiplier, 2) * sd
        x <- seq(g$mean - sdint, g$mean + sdint, length.out = numPoints)
        dists[[length(dists) + 1]] <- list(x = x, y = dnorm(x, g$mean, sd))
      } else if (type == "gld") {
        quantiles <- if.not.null(g$quantiles, seq(0.01, 0.99, length.out = numPoints))
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
      sapply(points, function(p) if.not.null(p$y, 0)),
      sapply(bounds, function(p) if.not.null(p$ymin, def_bound_ymin)),
      sapply(intervals, function(p) if.not.null(p$y, (def_bound_ymin + def_bound_ymax) / 2)),
      sapply(dists, function(p) min(p$y, na.rm = TRUE))
    )), na.rm = TRUE), "ymin")
  }

  if (is.na(ylim[[2]])) {
    ylim[[2]] <- boundFun(max(as.numeric(c(
      sapply(points, function(p) if.not.null(p$y, 0)),
      sapply(bounds, function(p) if.not.null(p$ymax, def_bound_ymax)),
      sapply(intervals, function(p) if.not.null(p$y, (def_bound_ymin + def_bound_ymax) / 2)),
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
      pch <- if.not.null(g$pch, 21)
      col <- if.not.null(g$col, "black")
      cex <- if.not.null(g$cex, 3)
      graphics::points(
        x = g$value, y = if.not.null(g$y, 0), type = "p",
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

      col <- if.not.null(g$col, rgb(0, 0, 0.0, alpha = 0.2))
      density <- if.not.null(g$density, NULL)
      border <- if.not.null(g$border, NA)
      graphics::rect(
        xleft = g$xmin, xright = g$xmax,
        ybottom = if.not.null(g$ymin, def_bound_ymin),
        ytop = if.not.null(g$ymax, def_bound_ymax), density = density,
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

      col <- if.not.null(g$col, "black")
      lty <- if.not.null(g$lty, 1)
      lwd <- if.not.null(g$lwy, 1)
      pch <- if.not.null(g$pch, 8)
      y <- if.not.null(g$y, 0)
      cex <- if.not.null(g$cex, c(1, 1, 1))

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
      col <- if.not.null(g$col, "black")
      lty <- if.not.null(g$lty, 1)
      lwd <- if.not.null(g$lwd, 1)

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



#' Create Table of Coefficients
#'
#' This function summarizes a list of estimated models (output of \code{estim.?} functions) and creates
#' a table of coefficients.
#'
#' @param list A named list where each element is output from a \code{estim.?} function, all belonging to a common analysis.
#' @param depInd Index of dependent variable in common analysis.
#' @param regInfo List of pairs of keys and names to determine information at bottom of table. Use "" (empty) for empty rows.
#' A list of available options is given in the details section.
#' Note that \code{num_eq} and \code{num_endo} (and \code{num_x} and \code{num_exo}) will differ with PCA analysis enabled.
#' @param hnameFun Function to change name of headers.
#' @param vnamesFun Function to change name of variables or codes in \code{regInfo}.
#' @param vnamesFun_sub List for replacing special character vectors in \code{vnamesFun}.
#' @param vnamesFun_max Maximum length for names in \code{vnamesFun}.
#' @param tableFun Function with arguments \code{coef}, \code{std}, \code{pvalue}, \code{minInColm}, \code{maxInCol} for formatting estimated coefficients.
#' Can also use one of the following character strings for predefined formatting: \code{sign}, \code{sign_star}, \code{coef}, \code{coef_star}, \code{coef_star_std}.
#' @param formatNumFun Function to format numbers if \code{tableFun} uses predefined formatting.
#' @param numCoefs Determines number of coefficients to insert in table. If \code{NA}, inserts all coefficients.
#' @param formatLatex If \code{TRUE}, default options are for 'latex', otherwise, 'html'.
#'
#' @details
#' The first part of the table is the header, followed by the coefficients. At the bottom, you can insert
#' the following items by specifying \code{regInfo}:
#' \itemize{
#' \item "" : empty line
#' \item \code{num_obs}: _No. Obs.; Number of observations.
#' \item \code{num_endo}: No. Eq. (orig.); Original number of equations or endogenous variables before being changed by PCA analysis.
#' \item \code{pca_y_exact}: PCA Count (y);
#' \item \code{pca_y_cutoff}: PCA Cutoff (y);
#' \item \code{pca_y_max}: PCA Max (y);
#' \item \code{num_eq}: No. Eq.; Number of equations after PCA analysis.
#' \item \code{num_exo}: No. Exo. (orig.);
#' \item \code{pca_x_exact}: PCA Count (x);
#' \item \code{pca_x_cutoff}: PCA Cutoff (x);
#' \item \code{pca_x_max}: PCA Max (x);
#' \item \code{num_x}: No. Exo.;
#' \item \code{num_x_all}: No. Exo. (all); Number of explanatory variables in all equations.
#' \item \code{num_rest}: No. Rest.; Number of restrictions in equation
#' \item \code{sigma2}: S.E. Reg.
#' \item \code{...}: An available measure name.
#' }
#'
#' @return A data frame with requested information.
#' @export
table.coefs <- function(list, depInd = 1,
                        regInfo = list(
                          c("", " "), c("num_obs", "No. Obs."), c("num_eq", "No. Eq."),
                          c("num_x", "No. Exo."), c("sigma2", "S.E. Reg."),
                          c("aic", "AIC"), c("sic", "SIC")
                        ),
                        hnameFun = function(x) x,
                        vnamesFun = function(x) x,
                        vnamesFun_sub = list(c("%", "\\\\%"), c("_", "\\\\_")), vnamesFun_max = 20,
                        tableFun = "coef_star", formatNumFun = function(colIndex, x) {
                          x
                        }, numCoefs = NA, formatLatex = TRUE) {
  get_stars <- function(pvalue) {
    if (is.nan(pvalue)) { # e.g., restricted to zero
      return(if (formatLatex) "\\textsuperscript{(r)}" else "<sup>(r)</sup>")
    }
    paste0(
      (if (formatLatex) "\\textsuperscript{" else "<sup>"),
      if (pvalue <= 0.01) {
        "***"
      } else if (pvalue <= 0.05) {
        "**"
      } else if (pvalue <= 0.1) {
        "*"
      } else {
        ""
      }, (if (formatLatex) "}" else "</sup>")
    )
  }

  if (tableFun == "sign") {
    tableFun <- function(j, coef, std, pvalue) {
      if (coef > 0) {
        "+"
      } else if (coef < 0) {
        "-"
      } else {
        "0"
      }
    }
  }
  if (tableFun == "sign_star") {
    tableFun <- function(j, coef, std, pvalue) {
      paste0(if (coef > 0) {
        "+"
      } else if (coef < 0) {
        "-"
      } else {
        "0"
      }, get_stars(pvalue))
    }
  } else if (tableFun == "coef") {
    tableFun <- function(j, coef, std, pvalue) {
      formatNumFun(j, coef)
    }
  } else if (tableFun == "coef_star") {
    tableFun <- function(j, coef, std, pvalue) {
      paste0(formatNumFun(j, coef), get_stars(pvalue))
    }
  } else if (tableFun == "coef_star_std") {
    tableFun <- function(j, coef, std, pvalue) {
      paste0(formatNumFun(j, coef), get_stars(pvalue), " (", formatNumFun(j, std), ")")
    }
  }

  list_names <- names(list)
  if (is.null(list_names)) {
    col_names <- paste0("m", rep(1:length(list)))
  } else {
    col_names <- lapply(list_names, function(n) hnameFun(n))
  }
  vnames <- unique(unlist(lapply(list, function(e) row.names(e$estimations$coefs))))
  vnames_0 <- sapply(vnames, function(n) {
    r <- vnamesFun(n)
    for (sg in vnamesFun_sub) {
      r <- gsub(sg[[1]], sg[[2]], r, fixed = TRUE)
    }
    if (nchar(r) > vnamesFun_max) {
      r <- paste0(substr(r, 1, vnamesFun_max - 3), "...")
    }
    r
  })
  regnam_0 <- sapply(regInfo, function(n) n[[2]])
  numCoefs <- as.integer(numCoefs)
  if (is.na(numCoefs)) {
    vvnames <- c(vnames_0, regnam_0)
  } else if (numCoefs > 0) {
    vvnames <- c(vnames_0[1:numCoefs], regnam_0)
  } else {
    vvnames <- regnam_0
  }

  r_table <- matrix("", (if (is.na(numCoefs)) {
    length(vnames) + length(regInfo)
  } else {
    length(regInfo) + numCoefs
  }), length(col_names))
  rownames(r_table) <- vvnames
  colnames(r_table) <- col_names

  j <- 0
  for (e in list) {
    j <- j + 1

    ns <- row.names(e$estimations$coefs)
    i <- 0
    # insert coefficients
    for (v in ns) {
      i <- i + 1
      if (is.na(numCoefs) == FALSE && numCoefs < i) {
        break
      }
      ind <- which(vnames == v)
      if (length(ind) == 1) {
        k <- ind[[1]]
        coef <- e$estimations$coefs[i, depInd]
        std <- e$estimations$stds[i, depInd]
        pv <- e$estimations$pValues[i, depInd]
        r <- tableFun(j, coef, std, pv)
        r_table[k, j] <- r
      }
    }


    # insert regression information
    simResNames <- names(e$simulation$results)
    i <- if (is.na(numCoefs)) length(vnames) else numCoefs
    for (rn in regInfo) {
      i <- i + 1
      v <- NULL
      r <- rn[[1]]

      if (r == "") {
        v <- "" # for an empty cell/line
      } else if (r == "num_obs") {
        v <- formatNumFun(j, as.integer(e$counts$obs))
      } else if (r == "num_eq") {
        v <- formatNumFun(j, as.integer(e$counts$eq))
      } else if (r == "num_endo") {
        v <- formatNumFun(j, as.integer(ncol(e$info$y)))
      } else if (r == "num_x") {
        v <- formatNumFun(j, as.integer(length(which(
          is.na(e$estimations$coefs[, depInd]) == FALSE
        ))))
      } else if (r == "num_exo") {
        v <- formatNumFun(j, as.integer(ncol(e$info$x)))
      } else if (r == "num_x_all") {
        v <- formatNumFun(j, as.integer(e$counts$exoAll))
      } else if (r == "num_rest") {
        v <- if (is.null(e$estimations$isRestricted)) {
          0L
        } else {
          as.integer(sum(e$estimations$isRestricted[, depInd]))
        }
      } else if (r == "pca_y_exact") {
        v <- formatNumFun(j, as.integer(e$info$pcaOptionsY$exactCount))
      } else if (r == "pca_x_exact") {
        v <- formatNumFun(j, as.integer(e$info$pcaOptionsX$exactCount))
      } else if (r == "pca_y_cutoff") {
        v <- formatNumFun(j, if (is.null(e$info$pcaOptionsY$exactCount) == F &&
                                 e$info$pcaOptionsY$exactCount == 0) {
          e$info$pcaOptionsY$cutoffRate
        } else {
          NA
        })
      } else if (r == "pca_x_cutoff") {
        v <- formatNumFun(j, if (is.null(e$info$pcaOptionsX$exactCount) == F &&
                                 e$info$pcaOptionsX$exactCount == 0) {
          e$info$pcaOptionsX$cutoffRate
        } else {
          NA
        })
      } else if (r == "pca_y_max") {
        v <- formatNumFun(j, as.integer(if (is.null(e$info$pcaOptionsY$exactCount) == F &&
                                            e$info$pcaOptionsY$exactCount == 0) {
          e$info$pcaOptionsY$max
        } else {
          NA
        }))
      } else if (r == "pca_x_max") {
        v <- formatNumFun(j, as.integer(if (is.null(e$info$pcaOptionsX$exactCount) == F &&
                                            e$info$pcaOptionsX$exactCount == 0) {
          e$info$pcaOptionsX$max
        } else {
          NA
        }))
      } else if (r == "sigma2") {
        v <- formatNumFun(j, e$estimations$sigma[depInd, depInd])
      } else { # must be a measure
        ind <- which(r == rownames(e$measures))
        if (length(ind) != 0) {
          v <- formatNumFun(j, e$measures[ind, depInd])

          if (r == "f") {
            ind0 <- which("fProb" == rownames(e$measures))
            fp <- e$measures[ind0, depInd]
            if (is.na(fp)) {
              warning("p-value of 'F' statistics is NA.")
            }
            v <- paste0(v, get_stars(fp))
          }
        }
      }

      if (is.null(v) || length(v) == 0) {
        warning(paste0("Invalid 'regInfo' member is requested. 'NA' is used. code=", r))
        v <- NA
      }
      r_table[i, j] <- v
    }
  }

  return(r_table)
}


