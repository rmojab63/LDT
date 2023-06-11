

test_that("Frequency conversion (from DateList to Daily) works", {

  startFreq <- f.list.date(c("20220904","20220903","20220902","20220901"), "20220901")
  v <- variable(c(4,3,2,1), "V", startFreq, list())

  w <- convert.to.daily(v)
  expect_equal(c(1,2,3,4), w$data)
  expect_equal("20220901", as.character(w$startFrequency))

  # with NA

  startFreq <- f.list.date(c("20220904","20220901"), "20220901")
  v <- variable(c(4,1), "V", startFreq, list())

  w <- convert.to.daily(v)
  expect_equal(c(1,NaN,NaN,4), w$data)
  expect_equal("20220901", as.character(w$startFrequency))

})

test_that("Frequency conversion (from Daily to Multi-Day) works", {

  startFreq <- f.daily(2022,9,1)
  v <- variable(c(1,2,3,4,5,6), "V", startFreq, list())
  w <- convert.to.multidaily(v,2,function(x)mean(x),FALSE)
  expect_equal(c(1.5,3.5,5.5), w$data)
  expect_equal("20220901", as.character(w$startFrequency))

  # not divisible
  v <- variable(c(1,2,3,4,5,6,7), "V", startFreq, list())

  w <- convert.to.multidaily(v,2,function(x)mean(x),FALSE)
  expect_equal(c(1.5,3.5,5.5,7.0), w$data)
  expect_equal("20220901", as.character(w$startFrequency))

  # from end
  w <- convert.to.multidaily(v,2,function(x)mean(x),TRUE)
  expect_equal(c(1.0, 2.5, 4.5, 6.5), w$data)
  expect_equal("20220901", as.character(w$startFrequency))

  # with NA
  v <- variable(c(1,NA,3,4,5,6,7), "V", startFreq, list())
  w <- convert.to.multidaily(v,2,function(x)mean(x),TRUE)
  expect_true(is.na(w$data[[2]]))

  # with c functions
  w <- convert.to.multidaily(v,2,"mean",TRUE)
  expect_equal(c(1.0, 3.0, 4.5, 6.5), w$data)

})

test_that("Frequency conversion (from Daily to Weekly) works", {

  startFreq <- f.daily(2023,6,1) # It is Thursday
  v <- variable(c(1:14), "V", startFreq, list()) # two weeks

  w <- convert.to.weekly(v,"thu",function(x)mean(x)) # start from Thursday
  expect_equal(c(sum(c(1:7))/7, sum(c(8:14))/7), w$data)
  expect_equal("20230601", as.character(w$startFrequency))

  # more than two weeks
  v <- variable(c(1:17), "V", startFreq, list())
  w <- convert.to.weekly(v,"thu",function(x)mean(x))
  expect_equal(c(sum(c(1:7))/7, sum(c(8:14))/7, (15+16+17)/3), w$data)
  expect_equal("20230601", as.character(w$startFrequency))

  # start of week : Wednesday
  w <- convert.to.weekly(v,"wed",function(x)mean(x))
  expect_equal(c(sum(c(NA, 1:6))/7, sum(c(7:13))/7, sum(c(14:17))/4), w$data)
  expect_equal("20230531", as.character(w$startFrequency))

  # start of week : Monday
  w <- convert.to.weekly(v,"mon",function(x)mean(x))
  expect_equal(c(sum(c(NA,NA,NA, 1:4))/7, sum(c(5:11))/7, sum(c(12:17))/6), w$data)
  expect_equal("20230529", as.character(w$startFrequency))

})

test_that("Frequency conversion (from Daily to X-Times-A-Year) works", {

  startFreq <- f.daily(2023,1,1)
  v <- variable(c(1:(365*2)), "V", startFreq, list())

  # Monthly
  w <- convert.to.XxYear(v,12,function(x)mean(x))
  expect_equal(24, length(w$data))
  expect_equal("2023M1", as.character(w$startFrequency))
  expect_equal(sum(c(397:425))/29, w$data[[14]])
  expect_equal(sum(c(701:730))/30, w$data[[24]])

  # Yearly
  w <- convert.to.XxYear(v,1,function(x)mean(x))
  expect_equal(2, length(w$data))
  expect_equal("2023", as.character(w$startFrequency))
  expect_equal(sum(c(1:365))/365, w$data[[1]])
  expect_equal(sum(c(366:730))/365, w$data[[2]])


  # 24 times a year
  w <- convert.to.XxYear(v,24,function(x)mean(x))
  expect_equal(2*24, length(w$data))
  expect_equal("2023:1", as.character(w$startFrequency))
  expect_equal(sum(c(1:15))/15, w$data[[1]])
  expect_equal(sum(c(716:730))/15, w$data[[48]])

  # 5 times a year
  w <- convert.to.XxYear(v,5,function(x)mean(x))
  expect_equal(10, length(w$data))
  expect_equal("2023:1", as.character(w$startFrequency))
  expect_equal(sum(c(1:73))/73, w$data[[1]])
  expect_equal(sum(c(74:(2*73)))/73, w$data[[2]])
})
