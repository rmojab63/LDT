
YearBased <- list(
  list(A = f.cross.section(10), B = "10", C = "cs"),
  list(A = f.yearly(2000), B = "2000", C = "y"),
  list(A = f.quarterly(2000, 4), B = "2000Q4", C = "q"),
  list(A = f.monthly(2000, 1), B = "2000M1", C = "m"),
  list(A = f.multi.yearly(2000, 3), B = "2000", C = "z3"),
  list(A = f.x.times.a.year(2000, 3, 1), B = "2000:1", C = "y3"),
  list(A = f.x.times.z.years(2000, 4, 2, 1), B = "2000:1", C = "x4z2")
)

WeekBased <- list(
  list(A = f.weekly(2000, 10, 11), B = "20001011", C = "w"),
  list(A = f.daily(2000, 10, 11), B = "20001011", C = "d"),
  list(A = f.multi.weekly(2000, 10, 11, 4), B = "20001011", C = "w4"),
  list(A = f.multi.daily(2000, 10, 11, 5), B = "20001011", C = "d5"),
  list(A = f.daily.in.week(2000, 10, 11, "mon", "sun", TRUE), B = "20001011", C = "i:mon-sun")
)

DayBased_Weekly <- list(
  list(A = f.hourly(WeekBased[[1]]$A, 5), B = "20001011:5", C = "ho|w"),
  list(A = f.minutely(WeekBased[[1]]$A, 100), B = "20001011:100", C = "mi|w"),
  list(A = f.secondly(WeekBased[[1]]$A, 200), B = "20001011:200", C = "se|w"),
  list(A = f.x.times.a.day(WeekBased[[1]]$A, 12, 6), B = "20001011:6", C = "da12|w")
)

DayBased_Daily <- list(
  list(A = f.hourly(WeekBased[[2]]$A, 5), B = "20001011:5", C = "ho|d"),
  list(A = f.minutely(WeekBased[[2]]$A, 100), B = "20001011:100", C = "mi|d"),
  list(A = f.secondly(WeekBased[[2]]$A, 200), B = "20001011:200", C = "se|d"),
  list(A = f.x.times.a.day(WeekBased[[2]]$A, 12, 6), B = "20001011:6", C = "da12|d")
)

DayBased_MultiWeekly <- list(
  list(A = f.hourly(WeekBased[[3]]$A, 5), B = "20001011:5", C = "ho|w4"),
  list(A = f.minutely(WeekBased[[3]]$A, 100), B = "20001011:100", C = "mi|w4"),
  list(A = f.secondly(WeekBased[[3]]$A, 200), B = "20001011:200", C = "se|w4"),
  list(A = f.x.times.a.day(WeekBased[[3]]$A, 12, 6), B = "20001011:6", C = "da12|w4")
)

DayBased_MultiDaily <- list(
  list(A = f.hourly(WeekBased[[4]]$A, 5), B = "20001011:5", C = "ho|d5"),
  list(A = f.minutely(WeekBased[[4]]$A, 100), B = "20001011:100", C = "mi|d5"),
  list(A = f.secondly(WeekBased[[4]]$A, 200), B = "20001011:200", C = "se|d5"),
  list(A = f.x.times.a.day(WeekBased[[4]]$A, 12, 6), B = "20001011:6", C = "da12|d5")
)

DayBased_DailyInWeek <- list(
  list(A = f.hourly(WeekBased[[5]]$A, 5), B = "20001011:5", C = "ho|i:mon-sun"),
  list(A = f.minutely(WeekBased[[5]]$A, 100), B = "20001011:100", C = "mi|i:mon-sun"),
  list(A = f.secondly(WeekBased[[5]]$A, 200), B = "20001011:200", C = "se|i:mon-sun"),
  list(A = f.x.times.a.day(WeekBased[[5]]$A, 12, 6), B = "20001011:6", C = "da12|i:mon-sun")
)

List_String <- list(
  list(A = f.list.string(c("A", "B", "C", "D"), "C"), B = "C", C = "Ls:A;B;C;D"),
  list(A = f.list.string(c("A", "B", "C", "D"), "D"), B = "D", C = "Ls:A;B;C;D")
)

List_Date <- list(
  list(A = f.list.date(c("20001001", "20001101", "20001202", "20010104"), "20001001"), B = "20001001", C = "Ld:20001001;20001101;20001202;20010104")
)

test_that_tp <- function(L) {
  for (f in L) {
    s <- as.character(f$A)
    expect_equal(f$B, s)

    c <- get.class.id(f$A)
    expect_equal(f$C, c)

    g <- as.frequency(s, c)
    expect_equal(f$A, g)
  }
}

test_that("Year-Based ToString and Parse works", {
  test_that_tp(YearBased)
})

test_that("Week-Based ToString and Parse works", {
  test_that_tp(WeekBased)
})

test_that("DayBased_Weekly ToString and Parse works", {
  test_that_tp(DayBased_Weekly)
})

test_that("DayBased_Daily ToString and Parse works", {
  test_that_tp(DayBased_Daily)
})

test_that("DayBased_MultiWeekly ToString and Parse works", {
  test_that_tp(DayBased_MultiWeekly)
})

test_that("DayBased_MultiDaily ToString and Parse works", {
  test_that_tp(DayBased_MultiDaily)
})

test_that("DayBased_DailyInWeek ToString and Parse works", {
  test_that_tp(DayBased_DailyInWeek)
})

test_that("List_String ToString and Parse works", {
  test_that_tp(List_String)
})

test_that("List_Date ToString and Parse works", {
  test_that_tp(List_Date)
})
