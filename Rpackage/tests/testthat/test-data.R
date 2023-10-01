
test_that("get.data handles matrix without name correctly", {
  data <- matrix(1:24, ncol = 6)

  result <- get.data(data, endogenous = 1)
  expect_equal(colnames(result$data), c("Y1","X1","X2","X3","X4","X5"))
  expect_equal(result$numExo, ncol(data) - 1)
  expect_equal(result$obsCount, 4)

  result <- get.data(data, endogenous = 2)
  expect_equal(colnames(result$data), c("Y1", "Y2", "X1","X2","X3","X4"))
  expect_equal(result$numExo, ncol(data) - 2)

  result <- get.data(data, endogenous = 2, addIntercept = TRUE)
  expect_equal(colnames(result$data), c("Y1", "Y2", "(Intercept)", "X1","X2","X3","X4"))
  expect_equal(result$numExo, ncol(data) - 1) # + intercept

  result <- get.data(data, endogenous = 2, addIntercept = TRUE, weights = c(1:4))
  expect_equal(colnames(result$data), c("Y1", "Y2", "(weights)", "(Intercept)", "X1","X2","X3","X4"))

  result <- get.data(data, endogenous = 2, addIntercept = FALSE, weights = c(1:4))
  expect_equal(colnames(result$data), c("Y1", "Y2", "(weights)", "X1","X2","X3","X4"))

  result <- get.data(data, endogenous = 2, addIntercept = TRUE, weights = c(1:4), lambdas = NA)
  expect_equal(length(result$lambdas), 2)

  result <- get.data(data, endogenous = 2, addIntercept = TRUE, weights = c(1:4), lambdas = c(0.5, 0))
  expect_equal(length(result$lambdas), 2)

  newX <- matrix(51:62, ncol = 4)
  result <- get.data(data, endogenous = 2, newData = newX)
  expect_equal(c(result$data[,1]), c(1,2,3,4))
  expect_equal(c(result$data[,2]), c(5,6,7,8))
  expect_equal(c(result$data[,3]), c(9,10,11,12))
  result <- get.data.append.newX(result)
  expect_equal(c(result$data[,1]), c(1,2,3,4,NA,NA,NA))
  expect_equal(c(result$data[,2]), c(5,6,7,8,NA,NA,NA))
  expect_equal(c(result$data[,3]), c(9,10,11,12,51,52,53))

  result <- get.data(data, endogenous = 2, newData = newX, addIntercept = TRUE)
  expect_equal(c(result$data[,1]), c(1,2,3,4))
  expect_equal(c(result$data[,2]), c(5,6,7,8))
  expect_equal(c(result$data[,3]), rep(1,4))
  expect_equal(c(result$data[,4]), c(9,10,11,12))
  result <- get.data.append.newX(result)
  expect_equal(c(result$data[,1]), c(1,2,3,4,NA,NA,NA))
  expect_equal(c(result$data[,2]), c(5,6,7,8,NA,NA,NA))
  expect_equal(c(result$data[,3]), rep(1,7))
  expect_equal(c(result$data[,4]), c(9,10,11,12,51,52,53))

  result <- get.data(data, endogenous = 2, newData = newX, addIntercept = TRUE, weights = c(1:4))
  expect_equal(c(result$data[,1]), c(1,2,3,4))
  expect_equal(c(result$data[,2]), c(5,6,7,8))
  expect_equal(c(result$data[,3]), c(rep(1:4)))
  expect_equal(c(result$data[,4]), rep(1,4))
  expect_equal(c(result$data[,5]), c(9,10,11,12))
  expect_equal(result$newObsCount,3)
  result <- get.data.append.newX(result)
  expect_equal(c(result$data[,1]), c(1,2,3,4,NA,NA,NA))
  expect_equal(c(result$data[,2]), c(5,6,7,8,NA,NA,NA))
  expect_equal(c(result$data[,3]), c(rep(1:4), NA, NA, NA))
  expect_equal(c(result$data[,4]), rep(1,7))
  expect_equal(c(result$data[,5]), c(9,10,11,12,51,52,53))
  expect_equal(result$newObsCount,3)

  newX <- matrix(51:62, ncol = 3)
  result <- get.data(data, endogenous = 3, newData = newX, addIntercept = TRUE, weights = c(1:4))
  expect_equal(c(result$data[,1]), c(1,2,3,4))
  expect_equal(c(result$data[,8]), c(21,22,23,24))
  result <- get.data.append.newX(result)
  expect_equal(c(result$data[,1]), c(1,2,3,4,NA,NA,NA,NA))
  expect_equal(c(result$data[,8]), c(21,22,23,24,59,60,61,62))
  expect_equal(result$newObsCount,4)

})


test_that("get.data handles matrix with name correctly", {
  data <- matrix(1:24, ncol = 6,
                 dimnames = list(NULL,c("V1", "V2", "V3", "V4", "V5", "V6")))

  result <- get.data(data, endogenous = c("V6", "V1"))
  expect_equal(colnames(result$data), c("V6","V1","V2","V3","V4","V5"))
  expect_equal(result$numExo, ncol(data) - 2)
  expect_equal(result$obsCount, 4)
  expect_equal(result$data[,"V6"], data[,"V6"])
  expect_equal(result$data[,"V1"], data[,"V1"])
  expect_equal(result$data[,"V3"], data[,"V3"])

  result <- get.data(data, endogenous = c("V6", "V1"), addIntercept = TRUE)
  expect_equal(colnames(result$data), c("V6","V1", "(Intercept)", "V2","V3","V4","V5"))
  expect_equal(result$numExo, ncol(data) - 1) # + intercept

  result <- get.data(data, endogenous = c("V6", "V1"), addIntercept = TRUE, weights = c(1:4))
  expect_equal(colnames(result$data), c("V6","V1", "(weights)", "(Intercept)", "V2","V3","V4","V5"))

  newX <- matrix(51:62, ncol = 4, dimnames = list(NULL, c("V2","V3","V4","V5")))
  result <- get.data(data, endogenous = c("V6", "V1"), newData = newX)
  result <- get.data.append.newX(result)
  expect_equal(c(result$data[,1]), c(21,22,23,24,NA,NA,NA))
  expect_equal(c(result$data[,2]), c(1,2,3,4,NA,NA,NA))
  expect_equal(c(result$data[,3]), c(5,6,7,8,51,52,53))

  newX <- matrix(51:62, ncol = 4, dimnames = list(NULL, c("V2","V3","V0","V5")))
  result <- get.data(data, endogenous = c("V6", "V1", "V4"), newData = newX, addIntercept = TRUE, weights = c(1:4))
  result <- get.data.append.newX(result)
  expect_equal(c(result$data[,1]), c(21,22,23,24,NA,NA,NA))
  expect_equal(c(result$data[,2]), c(1,2,3,4,NA,NA,NA))
  expect_equal(c(result$data[,"V2"]), c(5,6,7,8,51,52,53))

})



test_that("get.data handles equations correctly", {
  data <- data.frame(matrix(1:24, ncol = 6))
  colnames(data) <- c("X1", "X2", "Y2", "X3", "Y1", "X4")

  equations <- list(
     Y1 ~ X2 + X1,
     Y2 ~ X4 + X3
  )

  result <- get.data(data, equations = equations)
  expect_equal(colnames(result$data), c("Y1","Y2","X2","X1","X4","X3"))
  expect_equal(result$numExo, ncol(data) - 2)
  expect_equal(result$obsCount, 4)
  expect_equal(result$data[,"Y1"], data[,"Y1"])
  expect_equal(result$data[,"Y2"], data[,"Y2"])
  expect_equal(result$data[,"X3"], data[,"X3"])


  result <- get.data(data, equations = equations, addIntercept = TRUE)
  expect_equal(colnames(result$data), c("Y1","Y2", "(Intercept)", "X2","X1","X4","X3"))
  expect_equal(result$numExo, ncol(data) - 1)
  expect_equal(result$obsCount, 4)
  expect_equal(result$data[,"Y1"], data[,"Y1"])
  expect_equal(result$data[,"Y2"], data[,"Y2"])
  expect_equal(result$data[,"X3"], data[,"X3"])


  newX <- matrix(51:62, ncol = 4, dimnames = list(NULL, c("X1","X2","X3","X4")))
  result <- get.data(data, equations = equations, newData = newX)
  result <- get.data.append.newX(result)
  expect_equal(c(result$data[,"Y1"]), c(data[,"Y1"],NA,NA,NA))
  expect_equal(c(result$data[,"X2"]), c(data[,"X2"],54, 55, 56))
  expect_equal(c(result$data[, "X1"]), c(data[,"X1"],51, 52, 53))

  newX <- matrix(51:65, ncol = 5, dimnames = list(NULL, c("X1","X2","X3","X4", "X0")))
  result <- get.data(data, equations = equations, newData = newX)
  result <- get.data.append.newX(result)
  expect_equal(c(result$data[,"Y1"]), c(data[,"Y1"],NA,NA,NA))
  expect_equal(c(result$data[,"X2"]), c(data[,"X2"],54, 55, 56))
  expect_equal(c(result$data[, "X1"]), c(data[,"X1"],51, 52, 53))

})
