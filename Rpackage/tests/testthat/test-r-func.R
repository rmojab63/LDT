
# it seems that a global function does not work in most situations
indices <- c()
func <- function(columnIndices, numEndo){
  d <- paste0(columnIndices, collapse = ",")
  indices[length(indices) + 1] <<- d
  #print(d)
}


test_that("Searcher general works with inner exogenous", {

  x <- matrix(runif(100, 1, 10),20,5)
  indices <<- c()
  res = search.rfunc(func = func, length1 = 1, isInnerExogenous = TRUE,
                     data = get.data(x, endogenous = 2, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1,2),
                                                     innerGroups = list(c(1), c(2), c(1,2)),
                                                     numTargets = 2))
  expect_equal(9, length(indices))
  expect_equal(length(indices), length(unique(indices)))

  indices <<- c()
  res = search.rfunc(func = func, length1 = 1, isInnerExogenous = TRUE,
                     data = get.data(x, endogenous = 2, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1,2),
                                                     innerGroups = list(c(1), c(2), c(1,2)),
                                                     numTargets = 1))
  expect_equal(6, length(indices))
  expect_equal(length(indices), length(unique(indices)))

  indices <<- c()
  res = search.rfunc(func = func, length1 = 1, isInnerExogenous = TRUE,
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
  res = search.rfunc(func = func, length1 = 1, isInnerExogenous = FALSE,
                     data = get.data(x, endogenous = 2, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1,2),
                                                     innerGroups = list(c(1), c(2), c(1,2)),
                                                     numTargets = 2))
  expect_equal(18, length(indices))
  expect_equal(length(indices), length(unique(indices)))

  indices <<- c()
  res = search.rfunc(func = func, length1 = 1, isInnerExogenous = FALSE,
                     data = get.data(x, endogenous = 2, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1,2),
                                                     innerGroups = list(c(1,2)),
                                                     numTargets = 1))
  expect_equal(6, length(indices))
  expect_equal(length(indices), length(unique(indices)))

  indices <<- c()
  res = search.rfunc(func = func, length1 = 1, isInnerExogenous = FALSE,
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
  res = search.rfunc(func = func, length1 = 1, isInnerExogenous = TRUE,
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
  res = search.rfunc(func = func, length1 = 1, isInnerExogenous = TRUE,
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
  #using global function fails:
  func <- function(columnIndices, numEndo){
    d <- paste0(columnIndices, collapse = ",")
    #indices[length(indices) + 1] <<- d
    #print(d)
  }
  res = search.rfunc(func = func, length1 = 1, isInnerExogenous = TRUE,
                     data = get.data(x, endogenous = 50, addIntercept = FALSE),
                     combinations = get.combinations(sizes = c(1,2,3,4, 5),
                                                     stepsNumVariables = c(NA, 30, 20, 10),
                                                     innerGroups = list(c(1:3)),
                                                     numTargets = 3))
  #expect_equal(4, length(indices))
  #expect_equal(length(indices), length(unique(indices)))

})

