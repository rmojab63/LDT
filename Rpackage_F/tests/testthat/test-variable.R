
test_that("Bind Variables works with lag/lead adjustments", {
  mat <- matrix(c(1:28),7,4)
  mat[1,1]=NaN
  mat[6,1]=NaN
  mat[7,1]=NaN
  mat[5,2]=NaN
  mat[6,2]=NaN
  mat[7,2]=NaN
  mat[2,4]=NaN
  mat[5,3]=NaN
  mat[5,4]=NaN
  mat[4,3]=NaN
  freq <- f.monthly(2022,12)

  v1 <- variable(mat[,1], freq, "V1")
  v2 <- variable(mat[,2], freq, "V2")
  v3 <- variable(mat[,3], freq, "V3")
  v4 <- variable(mat[,4], freq, "V4")

  res <- bind.variables(list(v1,v2,v3,v4),interpolate = FALSE, adjustLeadLags = FALSE,
                            numExo = 0, horizon = 0)
  expected <- matrix(c(NaN,2,3,4,5,NaN,NaN,
                       8,9,10,11,NaN,NaN,NaN,
                       15,16,17,NaN,NaN,20,21,
                       22,NaN,24,25,NaN,27,28),7,4)
  colnames(expected) <- c("V1","V2","V3","V4")
  rownames(expected) <- c('2022M12', '2023M1', '2023M2', '2023M3', '2023M4', '2023M5','2023M6')
  expect_equal(expected, res$data)
  expect_equal(c(2,1,1,1,5,4,7,7,0,0,1,1,0,0,0,0,0,-1,2,2), as.numeric(res$info))

  res <- bind.variables(list(v1,v2,v3,v4),interpolate = TRUE, adjustLeadLags = TRUE,
                       numExo = 0, horizon = 0)
  expected <- matrix(c(NaN,NaN,NaN,2,3,4,5,
                       NaN,NaN,NaN,8,9,10,11,
                       15,16,17,18,19,20,21,
                       22,23,24,25,26,27,28),7,4)
  colnames(expected) <- c("V1","V2(-1)","V3(+2)","V4(+2)")
  rownames(expected) <- c('2022M10', '2022M11', '2022M12', '2023M1', '2023M2', '2023M3', '2023M4')
  expect_equal(expected, res$data)
  expect_equal(c(4,4,1,1,7,7,7,7,0,0,0,0,0,0,2,2,0,0,0,0), as.numeric(res$info))


  res <- bind.variables(list(v1,v2,v3,v4),interpolate = FALSE, adjustLeadLags = TRUE,
                       numExo = 1, horizon = 2)
  expected <- matrix(c(NaN,NaN,NaN,2,3,4,5,NaN,NaN,
                       NaN,NaN,NaN,8,9,10,11,NaN,NaN,
                       15,16,17,NaN,NaN,20,21,NaN,NaN,
                       NaN,NaN,22,NaN,24,25,NaN,27,28),9,4)
  colnames(expected) <- c("V1","V2(-1)","V3(+2)","V4")
  rownames(expected) <- c('2022M10', '2022M11', '2022M12', '2023M1', '2023M2', '2023M3', '2023M4','2023M5','2023M6')
  expect_equal(expected, res$data)
  expect_equal(c(4,4,1,3,7,7,7,9,0,0,1,1,0,0,0,0,0,0,0,0), as.numeric(res$info))


})
