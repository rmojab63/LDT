


search.rfunc <- function(func,
                       length1,
                       isInnerExogenous,
                       data = get.data(),
                       combinations = get.combinations(),
                       metrics = get.search.metrics(),
                       modelChecks = get.search.modelchecks(),
                       items = get.search.items(),
                       options = get.search.options()){

  stopifnot(is.function(func))
  args <- formals(func)
  required_args <- c("columnIndices", "numEndo")

  for (arg in required_args) {
    if (!(arg %in% names(args)))
      stop(paste("The function 'func' does not have the required argument:", arg))
    #if (arg %in% required_args[3,5] && !is.null(args[[arg]]))
    #  print(paste("The argument", arg, "in function 'func' does not have a default value of 'NULL'."))
  }

  stopifnot(is.number(length1))
  length1 <- as.integer(length1)
  stopifnot(is.logical(isInnerExogenous) && length(isInnerExogenous) == 1)
  stopifnot(is.list(data))
  stopifnot(is.list(combinations))

  combinations <- get.indexation(combinations, data, isInnerExogenous) # it also check for inconsistencies, etc.

  if (options$parallel){
    options$parallel <- FALSE
    warning("Parallel computation is disabled. It is not implemented in this function.")
  }

  if (is.null(metrics))
    metrics = get.search.metrics()
  else
    stopifnot(is.list(metrics))
  if (is.null(modelChecks))
    modelChecks = get.search.modelchecks()
  else
    stopifnot(is.list(modelChecks))
  if (is.null(items))
    items = get.search.items()
  else
    stopifnot(is.list(items))
  if (is.null(options))
    options = get.search.options()
  else
    stopifnot(is.list(options))

  if (is.list(combinations$sizes)){ # use steps
    # steps will re-call this function with modified combinations in which sizes is no longer a list
    res <- search.steps("rfunc", isInnerExogenous = isInnerExogenous, data = data, combinations = combinations,
                        metrics = metrics, modelChecks = modelChecks, items = items, options = options,
                        func = func, length1 = length1)
    res
  }
  else {

    startTime <- Sys.time()
    res <- .SearchRFunc(data, combinations, metrics, modelChecks, items, options,
                      func, length1, isInnerExogenous)
    endTime <- Sys.time()

    res$info$data <- data
    res$info$combinations <- combinations
    res$info$metrics <- metrics
    res$info$options <- options
    res$info$modelChecks <- modelChecks
    res$info$items <- items
    res$info$startTime <- startTime
    res$info$endTime <- endTime

    class(res) <- c("ldt.search", "list")
    attr(res, "method") <- "rFunc"
    res
  }
}

