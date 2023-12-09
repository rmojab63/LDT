

test_that("Long-run growth works (discrete)", {
  g <- get.longrun.growth(c(1,2,1))
  expect_equal(0, g, tolerance = 1e-16)

  g <- get.longrun.growth(c(1,2))
  expect_equal(100, g, tolerance = 1e-16)

  y <- abs(rnorm(10,2,5))
  g <- get.longrun.growth(y)
  #dy <- (diff(y) / y[1:(length(y)-1)])
  x = y[[1]] # start from y_0 and reach y_T
  for (i in c(1:(length(y)-1)))
    x <- x * (1 + g/100)
  expect_equal(x, y[[length(y)]], tolerance = 1e-12)

})

test_that("Long-run growth works (continuous)", {

  y <- abs(rnorm(10,2,5))
  g <- get.longrun.growth(y, cont = TRUE)
  #dy <- log(y[2:(length(y))]) - log(y[1:(length(y)-1)])
  x = y[[1]] # start from y_0 and reach y_T
  for (i in c(1:(length(y)-1)))
    x <- x * exp(g/100)
  expect_equal(x, y[[length(y)]], tolerance = 1e-12)

})

test_that("Long-run growth works (discrete with NAs)", {

  y <- c(NA, 0, abs(rnorm(10,2,5)), 0, NA, NA)
  g <- get.longrun.growth(y, FALSE, FALSE, 2, 3, TRUE)
  x = y[[3]] # start from y_0 and reach y_T
  for (i in c(1:(length(y)-6)))
    x <- x * (1 + g/100)
  if (is.na(x))
    stop("x is NA")
  expect_equal(x, y[[length(y)-3]], tolerance = 1e-12)

})

test_that("Long-run growth works (%)", {

  y <- c(NA, 0, c(60,70,80,90), 0, NA, NA)
  g <- get.longrun.growth(y, TRUE, TRUE, 2, 3, TRUE)
  expect_equal(75, g, tolerance = 1e-12)

  g <- get.longrun.growth(y, FALSE, TRUE, 2, 3, TRUE)
  expect_equal(74.64202, g, tolerance = 1e-3) # TODO: check its validity

})



