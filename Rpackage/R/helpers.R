
#' Combine More Than One \code{ldtsearch} Objects
#'
#' @param list a list with \code{ldtsearch} objects
#' @param type1Name the name of 'type1' in the object
#'
#' @return the combined \code{ldtsearch} object
combineSearch <- function(list, type1Name = "coefs") {
  if (length(list) == 1) {
    return(list[[1]])
  }

  cmb_best_ind <- function(item, oldYNames, oldXNames, newYNames, newXNames) {
    if (is.null(item)) {
      return(NULL)
    }
    # send the best one to this function along with the old and new names
    # we should just update the indexes
    item$depIndices <- sort(as.integer(sapply(
      item$depIndices,
      function(q) which(oldYNames == newYNames[q])
    )))
    item$exoIndices <- sort(as.integer(sapply(
      item$exoIndices,
      function(q) which(oldXNames == newXNames[q])
    )))

    return(item)
  }

  cmb_inclusion <- function(first, second) {
    # assuming first column is mean and second column counts
    firstNames <- rownames(first)
    secondNames <- rownames(second)
    first[, 1] <- first[, 1] * first[, 2] # convert first column to sum
    j <- 0
    for (n in secondNames) {
      j <- j + 1
      ind <- which(n == firstNames)
      if (length(ind) > 0) {
        ind <- ind[[1]]
        first[[ind, 1]] <- first[[ind, 1]] + second[[j, 1]] * second[[j, 2]]
        first[[ind, 2]] <- first[[ind, 2]] + second[[j, 2]]
      }
    }
    first[, 1] <- first[, 1] / first[, 2] # convert back to mean
    return(first)
  }

  cmb_cdf <- function(first, second) {
    # assuming columns are: mean, count, weights
    firstNames <- rownames(first)
    secondNames <- rownames(second)
    first[, 1] <- first[, 1] * first[, 3] # convert first column to weighted sum
    j <- 0
    for (n in secondNames) {
      j <- j + 1
      ind <- which(n == firstNames)
      if (length(ind) > 0) {
        ind <- ind[[1]]
        first[[ind, 1]] <- first[[ind, 1]] + second[[j, 1]] * second[[j, 3]]
        first[[ind, 2]] <- first[[ind, 2]] + second[[j, 2]]
        first[[ind, 3]] <- first[[ind, 3]] + second[[j, 3]]
      }
    }
    first[, 1] <- first[, 1] / first[, 3] # convert back to mean
    return(first)
  }

  cmb_extremebound <- function(first, second) {
    # assuming first is lower and second column is upper bounds
    firstNames <- rownames(first)
    secondNames <- rownames(second)
    j <- 0
    for (n in secondNames) {
      j <- j + 1
      ind <- which(n == firstNames)
      if (length(ind) > 0) {
        ind <- ind[[1]]
        first[[ind, 1]] <- min(first[[ind, 1]], second[[j, 1]])
        first[[ind, 2]] <- max(first[[ind, 2]], second[[j, 2]])
      }
    }
    return(first)
  }

  cmb_mixture <- function(first, second) {
    # assuming that columns are: mean, variance, skewness, kurtosis,count, sumWeights

    firstNames <- rownames(first)
    secondNames <- rownames(second)
    j <- 0
    for (n in secondNames) {
      j <- j + 1
      ind <- which(n == firstNames)
      if (length(ind) > 0) {
        ind <- ind[[1]]
        mix <- s.combine.by.moments4(
          list(
            mean = first[[ind, 1]], variance = first[[ind, 2]],
            skewness = first[[ind, 3]], kurtosis = first[[ind, 4]],
            count = first[[ind, 5]], sumWeights = first[[ind, 6]]
          ),
          list(
            mean = second[[j, 1]], variance = second[[j, 2]],
            skewness = second[[j, 3]], kurtosis = second[[j, 4]],
            count = second[[j, 5]], sumWeights = second[[j, 6]]
          )
        )
        first[ind, 1] <- mix$mean
        first[ind, 2] <- mix$variance
        first[ind, 3] <- mix$skewness
        first[ind, 4] <- mix$kurtosis
        first[ind, 5] <- mix$count
        first[ind, 6] <- mix$sumWeights
      }
    }

    return(first)
  }


  result <- list[[1]]
  for (i in c(2:length(list))) {
    newR <- list[[i]]

    result$info$endTime <- newR$info$endTime

    result$counts$expectedCount <- result$counts$expectedCount + newR$counts$expectedCount
    result$counts$searchedCount <- result$counts$searchedCount + newR$counts$searchedCount
    result$counts$failedCount <- result$counts$failedCount + newR$counts$failedCount
    if (is.null(result$counts$failedDetails))
      result$counts$failedDetails = list()

    for (f in newR$counts$failedDetails) {
      j <- 0
      added <- FALSE
      for (f0 in result$counts$failedDetails){
        j <- j + 1
        if (f$message == f0$message){
          added <- TRUE
          result$counts$failedDetails[[j]]$count = f0$count + f$count
          break
        }
      }

      if (added == FALSE){
        ind <- length(result$counts$failedDetails) + 1
      }
    }

    for (m in c(2:(length(result) - 1))) {
      mea0 <- result[[m]]
      mea1 <- newR[[m]]
      for (t in c(1:length(mea0))) {
        tar0 <- mea0[[t]]
        tar1 <- mea1[[t]]

        # BEST
        if (is.null(tar0$model$bests) == FALSE) {
          # update indexes
          news <- lapply(tar1$model$bests, function(a) {
            cmb_best_ind(
              a, result$info$yNames,
              result$info$xNames, newR$info$yNames, newR$info$xNames
            )
          })
          # combine, sort by weights and select a subset
          b_count <- length(tar0$model$bests)
          all <- as.list(append(tar0$model$bests, news))

          idx <- as.numeric(sapply(all, function(v) v$weight))
          all <- all[order(idx, decreasing = TRUE)][1:b_count]
          names(all) <- paste0("best", c(1:b_count))

          result[[m]][[t]]$model$bests <- all
        }

        # ALL
        if (is.null(tar0$model$all) == FALSE) {
          for (a in tar1$model$all) {
            result[[m]][[t]]$model$all[[length(result[[m]][[t]]$model$all) + 1]] <-
              cmb_best_ind(
                a, result$info$yNames,
                result$info$xNames, newR$info$yNames, newR$info$xNames
              )
          }
        }
        # INCLUSION
        if (is.null(tar0$model$inclusion) == FALSE) {
          result[[m]][[t]]$model$inclusion <- cmb_inclusion(
            tar0$model$inclusion,
            tar1$model$inclusion
          )
        }


        # EXTREME BOUND
        if (is.null(tar0[[type1Name]]$extremeBounds) == FALSE) {
          eb <- cmb_extremebound(tar0[[type1Name]]$extremeBounds, tar1[[type1Name]]$extremeBounds)
          result[[m]][[t]][[type1Name]]$extremeBounds <- eb
        }

        # CDFs
        if (is.null(tar0[[type1Name]]$cdfs) == FALSE) {
          j <- 0
          for (c0 in tar0[[type1Name]]$cdfs) {
            j <- j + 1
            result[[m]][[t]][[type1Name]]$cdfs[[j]] <- cmb_cdf(c0, tar1[[type1Name]]$cdfs[[j]])
          }
        }

        # MIXTURE
        if (is.null(tar0[[type1Name]]$mixture) == FALSE) {
          a <- cmb_mixture(tar0[[type1Name]]$mixture, tar1[[type1Name]]$mixture)
          result[[m]][[t]][[type1Name]]$mixture <- a
        }

        # TYPE1 BESTS
        if (is.null(tar0[[type1Name]]$bests) == FALSE) {
          # for each item in this list, there is a list of bests

          for (item1 in tar1[[type1Name]]$bests) {
            item1 <- item1[lengths(item1) != 0]
            # update its indexes

            news <- lapply(
              item1[2:length(item1)] # first one is the name
              , function(a) {
                cmb_best_ind(
                  a, result$info$yNames, result$info$xNames,
                  newR$info$yNames, newR$info$xNames
                )
              }
            )
            # we should whether append it or update an existing list
            found <- FALSE
            o <- 0
            for (item0 in tar0[[type1Name]]$bests) {
              item0 <- item0[lengths(item0) != 0]
              o <- o + 1
              if (item1$name == item0$name) {
                found <- TRUE
                # combine, sort by weights and select a subset
                b_count <- length(item0)
                all <- as.list(append(item0[2:length(item0)], news))
                idx <- as.numeric(sapply(all, function(v) v$weight))
                all <- all[order(idx, decreasing = TRUE)][1:b_count]
                names(all) <- paste0("best", c(1:b_count))

                all <- c(name = item1$name, all)

                result[[m]][[t]][[type1Name]]$bests[[o]] <- all
                break
              }
            }
            if (found == FALSE) {
              news <- c(name = item1$name, news)
              result[[m]][[t]][[type1Name]]$bests[[o]] <- news
            }
          }
        }
      }
    }
  }
  return(result)
}


#' Stepwise estimation
#'
#' @param method sur, dc or varma
#' @param data exogenous (for sur and dc) or endogenous (for varma)
#' @param sizes determines the steps
#' @param counts determines the size in each step
#' @param savePre if not \code{NULL}, it saves and tries to load the
#' progress of search step in a file (name=\code{paste0(savePre,i)} where
#' \code{i} is the index of the step).
#' @param ... Additional arguments
#'
#' @return the result
Search_s <- function(method, data, sizes = list(c(1, 2), c(3, 4), c(5), c(6:10)),
                     counts = c(NA, 40, 30, 20), savePre, ...) {
  dots <- list(...)
  printMsg = dots$searchOptions$printMsg
  if (is.null(printMsg))
    printMsg = FALSE

  # TODO: check search items here
  # dots <- list(...)
  # if (dots$searchItems)
  if (length(sizes) != length(counts)) {
    stop("Invalid number of elements in 'counts'.")
  }

  if (is.null(colnames(data))) {
    stop("'data' must have column names.") # don't set it here, due to the loading part
  }

  estims <- list()
  data_i <- NULL

  for (i in c(1:length(sizes))) {
    if (printMsg) {
      cat(paste0("Step ", i, ". estimation started...\n"))
    }

    if (is.null(savePre) == FALSE) {
      res <- suppressWarnings(tryCatch(readRDS(paste0(
        savePre, i,
        ".RData"
      )), error = function(e) NULL))
      if (is.null(res) == FALSE) {
        if (printMsg) {
          cat(paste0("  ... finished (read from file).\n"))
        }
        estims[[i]] <- res$estim
        data_i <- res$data
        next()
      }
    }

    size_i <- sizes[[i]]
    count_i <- counts[[i]]
    if (is.na(count_i)) {
      count_i <- ncol(data)
    }

    if (i == 1) {
      if (count_i < ncol(data)) {
        warning(paste0("First ", count_i, " variables are selected in the first step."))
        if (printMsg) {
          cat(paste0("Selected Variables: First", count_i, "Variables.\n"))
        }
      } else if (printMsg) {
        cat(paste0("Selected Variables: All Variables.\n"))
      }
      data_i <- data[, c(1:count_i)] # however, the first one should
      # be NA such that no arbitrary selection happens here
    } else { # use previous to select
      est <- estims[[i - 1]]
      x_i_names <- colnames(data_i)

      # find best models and sort them
      bests <- list()
      for (m in c(2:(length(est) - 1))) {
        mea <- est[[m]]
        for (t in c(1:length(m))) {
          tar <- mea[[t]]
          if (is.null(tar$model$bests) == FALSE && length(tar$model$bests) > 0) {
            for (b in c(1:length(tar$model$bests))) {
              bst <- tar$model$bests[[b]]
              names <- NULL
              if (method == "sur" || method == "dc")
                names <- x_i_names[bst$exoIndices]
              else
                names <- x_i_names[bst$depIndices]
              bests[[length(bests) + 1]] <- list(index = b, names = names)
            }
          }
        }
      }

      # sort by index
      idx <- as.numeric(sapply(bests, function(b) {
        b$index
      }))
      bests <- bests[order(idx, decreasing = FALSE)]

      vars <- c()
      for (bns in bests) {
        vars <- append(vars, bns$names) # we add a group and its length might exceed the count
        vars <- unique(vars)
        if (length(vars) >= count_i) {
          break
        }
      }

      # TODO: we can use or fill the gap with Inclusion weights. However,
      # first we should determine which measure or target to use ?!

      if (printMsg) {
        cat(paste0("Selected Variables: ", paste0(vars, sep = ";", collapse = ""), "\n"))
      }

      data_i <- data_i[, vars, drop = FALSE]
    }
    if (any(ncol(data_i) < size_i)) {
      warning("There is not enough variables in this step. Increase the value of
      'bestK' or if you have fix variables, adjust the sizes. Search is stoped.")
      break
    }

    if (printMsg)
      cat("\n=================\n")

    if (method == "sur") {
      estims[[i]] <- search.sur(x = data_i, xSizes = size_i, ...)
    } else if (method == "dc") {
      estims[[i]] <- search.dc(x = data_i, xSizes = size_i, ...)
    } else if (method == "varma") {
      estims[[i]] <- search.varma(y = data_i, ySizes = size_i, ...)
    } else {
      stop("invalid method")
    }

    if (printMsg)
      cat("\n=================\n")

    if (estims[[i]]$counts$searchedCount == estims[[i]]$counts$failedCount) {
      if (printMsg){
         print("......Failures.......")
         print(estims[[i]]$counts$failedDetails)
      }
      stop("all estimations failed")
    }

    if (is.null(savePre) == FALSE) {
      saveRDS(list(estim = estims[[i]], data = data_i), paste0(savePre, i, ".RData"))
    }

    if (printMsg) {
      cat(paste0("  ... finished.\n\n"))
    }
  }


  result <- combineSearch(estims, if (method == "varma") "predictions" else "coefs")



  if (method == "sur") {
    class(result) <- c("ldtsearchsur", "ldtsearch", "list")
    attr(result, "method") <- "sur"
  } else if (method == "dc") {
    class(result) <- c("ldtsearchdc", "ldtsearch", "list")
    attr(result, "method") <- "dc"
  } else if (method == "varma") {
    class(result) <- c("ldtsearchvarma", "ldtsearch", "list")
    attr(result, "method") <- "varma"
  }

  return(result)
}

is.empty <- function(arg) {
  if (is.null(arg) || is.na(arg)) {
    return(TRUE)
  }
  if (is.list(arg)) {
    if (length(arg) == 0) {
      return(TRUE)
    }
  } else if (is.character(arg)) {
    if (arg == "") {
      return(TRUE)
    }
    return(FALSE)
  } else {
    stop("invalid type")
  }
}

list.has.array <- function(List, arr) {
  for (a in List) {
    if (length(a) == length(arr)) {
      if (all(a == arr)) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

if.not.null <- function(arg, default = NULL) {
  if (is.null(arg)) {
    return(default)
  }
  return(arg)
}
