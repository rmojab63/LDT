

test_that("Frequency conversion (from DateList to Daily) works", {

  startFreq <- f.list.date(c("20220904","20220903","20220902","20220901"), "20220901")
  v <- Variable(c(4,3,2,1), "var", startFreq, list())

  w <- c.datelist.daily(v)
  expect_equal(c(1,2,3,4), w$data)
  expect_equal("20220901", as.character(w$startFrequency))

  # with NA

  startFreq <- f.list.date(c("20220904","20220901"), "20220901")
  v <- Variable(c(4,1), "var", startFreq, list())

  w <- c.datelist.daily(v)
  expect_equal(c(1,NaN,NaN,4), w$data)
  expect_equal("20220901", as.character(w$startFrequency))

})


test_that("Frequency conversion (from Daily to Multi-Day) works", {

  startFreq <- f.daily(2022,9,1)
  v <- Variable(c(1,2,3,4,5,6), "var", startFreq, list())
  w <- c.daily.multidaily(v,2,function(x)mean(x),FALSE)
  expect_equal(c(1.5,3.5,5.5), w$data)
  expect_equal("20220901", as.character(w$startFrequency))

  # not divisible
  v <- Variable(c(1,2,3,4,5,6,7), "var", startFreq, list())

  w <- c.daily.multidaily(v,2,function(x)mean(x),FALSE)
  expect_equal(c(1.5,3.5,5.5,7.0), w$data)
  expect_equal("20220901", as.character(w$startFrequency))

  # from end
  w <- c.daily.multidaily(v,2,function(x)mean(x),TRUE)
  expect_equal(c(1.0, 2.5, 4.5, 6.5), w$data)
  expect_equal("20220901", as.character(w$startFrequency))

  # with NA
  v <- Variable(c(1,NA,3,4,5,6,7), "var", startFreq, list())
  w <- c.daily.multidaily(v,2,function(x)mean(x),TRUE)
  expect_true(is.na(w$data[[2]]))

})

