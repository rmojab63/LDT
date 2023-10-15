
indices <- c()
func <- function(columnIndices, numEndo,
                 data = NULL, items = NULL, metrics = NULL, modelChecks = NULL){
  d <- paste0(columnIndices, collapse = ",")
  indices[length(indices) + 1] <<- d
  return(list(error = "TEST",
              metrics = NULL,
              extra = NULL,
              type1means = NULL,
              type1vars = NULL))
}


test_that("Searcher general works with inner exogenous", {

  x <- matrix(runif(100, 1, 10),20,5)
  indices <<- c()
  res = search.rfunc(rFuncName = "func", length1 = 1, isInnerExogenous = TRUE,
                     data = get.data(x, endogenous = 2, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1,2),
                                                     innerGroups = list(c(1), c(2), c(1,2)),
                                                     numTargets = 2))
  expect_equal(9, length(indices))
  expect_equal(length(indices), length(unique(indices)))

  indices <<- c()
  res = search.rfunc(rFuncName = "func", length1 = 1, isInnerExogenous = TRUE,
                     data = get.data(x, endogenous = 2, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1,2),
                                                     innerGroups = list(c(1), c(2), c(1,2)),
                                                     numTargets = 1))
  expect_equal(6, length(indices))
  expect_equal(length(indices), length(unique(indices)))

  indices <<- c()
  res = search.rfunc(rFuncName = "func", length1 = 1, isInnerExogenous = TRUE,
                     data = get.data(x, endogenous = 3, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1,2),
                                                     innerGroups = list(c(1), c(2), c(1,2)),
                                                     numTargets = 2))
  expect_equal(15, length(indices))
  expect_equal(length(indices), length(unique(indices)))


})

test_that("Searcher general works with inner endogenous", {

  x <- matrix(runif(100, 1, 10),20,5)
  indices <<- c()
  res = search.rfunc(rFuncName = "func", length1 = 1, isInnerExogenous = FALSE,
                     data = get.data(x, endogenous = 2, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1,2),
                                                     innerGroups = list(c(1), c(2), c(1,2)),
                                                     numTargets = 2))
  expect_equal(18, length(indices))
  expect_equal(length(indices), length(unique(indices)))

  indices <<- c()
  res = search.rfunc(rFuncName = "func", length1 = 1, isInnerExogenous = FALSE,
                     data = get.data(x, endogenous = 2, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1,2),
                                                     innerGroups = list(c(1,2)),
                                                     numTargets = 1))
  expect_equal(6, length(indices))
  expect_equal(length(indices), length(unique(indices)))

  indices <<- c()
  res = search.rfunc(rFuncName = "func", length1 = 1, isInnerExogenous = FALSE,
                     data = get.data(x, endogenous = 3, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1,2),
                                                     innerGroups = list(c(1), c(2), c(1,2)),
                                                     numTargets = 2))
  expect_equal(9, length(indices))
  expect_equal(length(indices), length(unique(indices)))


})


test_that("Searcher general works with partitioning", {

  x <- matrix(runif(200, 1, 10),20,10)
  indices <<- c()
  res = search.rfunc(rFuncName = "func", length1 = 1, isInnerExogenous = TRUE,
                     data = get.data(x, endogenous = 5, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1:3),
                                                     partitions = list(c(1,2),c(3,4),c(5)),
                                                     innerGroups = list(c(1:3)),
                                                     numTargets = 3))
  expect_equal(14, length(indices))
  expect_equal(length(indices), length(unique(indices)))

})

test_that("Searcher general works with fixing partitions", {

  x <- matrix(runif(200, 1, 10), 20, 10)
  indices <<- c()
  res = search.rfunc(rFuncName = "func", length1 = 1, isInnerExogenous = TRUE,
                     data = get.data(x, endogenous = 5, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1:3),
                                                     numFixPartitions = 2,
                                                     innerGroups = list(c(1:3)),
                                                     numTargets = 3))
  expect_equal(4, length(indices))
  expect_equal(length(indices), length(unique(indices)))

})

test_that("Searcher general works with step-wise search", {
  skip_on_cran()

  x <- matrix(runif(400, 1, 10), 4, 100)
  indices <<- c()
  #res = search.rfunc(rFuncName = "func", length1 = 1, isInnerExogenous = TRUE,
  #                   data = get.data(x, endogenous = 50, addIntercept = FALSE),
  #                   combinations = get.combinations(sizes =  list(c(1,2),c(3),c(4)),
  #                                                   stepsNumVariables = c(NA, 30, 10),
  #                                                   innerGroups = list(c(1:3)),
  #                                                   numTargets = 3),
  #                   options = get.search.options(reportInterval = 0))

  #TODO

  #expect_equal(4, length(indices))
  #expect_equal(length(indices), length(unique(indices)))

})

