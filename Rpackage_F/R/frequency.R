
#' Create a Cross-Section Frequency
#'
#' This frequency is typically used for indexed data. It is represented by an integer
#' that indicates the position of the observation.
#'
#' @param position An integer representing the position of the observation.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#"} (the number is the \code{position})
#' \item **Class Id** \code{"cs"}
#' }
#'
#' @return An object of class \code{ldtf} which is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{position}{Determines the \code{position}.}
#'
#' @export
#'
#' @examples
#'
#' cs0 <- f.cross.section(10) # this initializes a cross-section frequency
#'
#' cs0_value_str <-  as.character(cs0) # this will be '10'.
#' cs0_class_str <- get.class.id(cs0) # this will be 'cs'.
#'
#' cs_new <- as.frequency("20", "cs")
#' #      this is a cross-section frequency. It points to position 20.
#'
f.cross.section <- function(position) {
  position = as.integer(position)
  res <- .F_cross_section(position)
  res
}



#' Create an Annual Frequency
#'
#' Use this function to create a frequency for time-series data that occurs annually.
#'
#' @param year An integer representing the year of the observation.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#"} (the number is the \code{year})
#' \item **Class Id** \code{"y"}
#' }
#'
#' @return An object of class \code{ldtf} which is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{year}{Determines the \code{year}.}
#'
#' @export
#'
#' @examples
#'
#' y0 <- f.yearly(2020) # this initializes a 'yearly' frequency
#'
#' y0_value_str <-  as.character(y0) # this will be '2020'.
#' y0_class_str <- get.class.id(y0) # this will be 'y'.
#'
#' y_new <- as.frequency("2021", "y") # this is a yearly frequency. It points to year 2021.
#'
f.yearly <- function(year) {
  year = as.integer(year)
  res <- .F_yearly(year)
  res
}



#' Create a Quarterly Frequency
#'
#' Use this function to create a frequency for time-series data that occurs quarterly.
#'
#' @param year An integer representing the year of the observation.
#' @param quarter An integer representing the quarter of the observation (It should be between 1 and 4).
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#q#"} (first '#' is the \code{year}, second '#' is the
#'  \code{quarter}; e.g., 2010q3 or 2010q4. Note that 2000q0 or 2000q5 are invalid.
#' \item **Class Id** \code{"q"}
#' }
#'
#' @return An object of class \code{ldtf} which is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{year}{Determines the \code{year}.}
#' \item{quarter}{Determines the \code{quarter}.}
#'
#' @export
#'
#' @examples
#'
#' q0 <- f.quarterly(2020, 2)
#' #     this is a quarterly frequency that refers to the second quarter of the year 2021.
#'
#' q0_value_str <-  as.character(q0) # this will be '2020Q2'.
#' q0_class_str <- get.class.id(q0) # this will be 'q'.
#'
#' q_new <- as.frequency("2021q3", "q")
#' #      this is a quarterly frequency that refers to the third quarter of the year 2021.
#'
#' # Don't make the following mistakes:
#' \donttest{
#' q_invalid <- try(f.quarterly(2020, 0))
#' q_invalid <- try(f.quarterly(2020, 5))
#' q_invalid <- try(as.frequency("2021q0", "q"))
#' q_invalid <- try(as.frequency("2021q5", "q"))
#' q_invalid <- try(as.frequency("2021", "q"))
#' }
#'
f.quarterly <- function(year, quarter){
  year = as.integer(year)
  quarter = as.integer(quarter)
  if (quarter < 1 || quarter > 4)
    stop("Invalid 'quarter'. It should be between 1 and 4.")
  res <- .F_quarterly(year, quarter)
  res
}



#' Create a Monthly Frequency
#'
#' Use this function to create a frequency for time-series data that occurs monthly.
#'
#' @param year An integer representing the year of the observation.
#' @param month An integer representing the month of the observation (It should be between 1 to 12).
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#m#"} (first # is the \code{year}, second # is
#'  the \code{month} (1 to 12); e.g., 2010m8 or 2010m12. Note that 2000m0 or 2000m13 are invalid.
#' \item **Class Id** \code{"m"}
#' }
#'
#' @return An object of class \code{ldtf} which is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{year}{Determines the \code{year}.}
#' \item{month}{Determines the \code{month}.}
#'
#' @export
#'
#' @examples
#'
#' m0 <- f.monthly(2020, 2)
#' #     this is a monthly frequency that refers to the second month of the year 2020.
#'
#' m0_value_str <-  as.character(m0) # this will be '2020M2'.
#' m0_class_str <- get.class.id(m0) # this will be 'm'.
#'
#' m_new <- as.frequency("2021m3", "m")
#' #     this is a monthly frequency that refers to the third month of the year 2021.
#'
#' # Don't make the following mistakes:
#' \donttest{
#' m_invalid <- try(f.monthly(2020, 0))
#' m_invalid <- try(f.monthly(2020, 5))
#' m_invalid <- try(as.frequency("2021m0", "m"))
#' m_invalid <- try(as.frequency("2021m13", "m"))
#' m_invalid <- try(as.frequency("2021", "m"))
#' }
#'
f.monthly <- function(year, month){
  year = as.integer(year)
  month = as.integer(month)
  if (month < 1 || month > 12)
    stop("Invalid 'month'. It should be between 1 and 12.")
  res <- .F_monthly(year, month)
  res
}



#' Create a Multi-Year Frequency
#'
#' Use this function to create a frequency for time-series data that occurs every \code{z} years.
#'
#' @param year An integer representing the year of the observation.
#' @param z An integer representing the number of years. It should be larger than zero.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#"} (the number is the \code{year}, which means the string representation is the first year of the interval)
#' \item **Class Id** \code{"z#"} ('#' represents the value: \code{z}; e.g., z3 means every 3 years)
#' }
#'
#' @return An object of class \code{ldtf} which is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{year}{Determines the \code{year}.}
#' \item{z}{Determines the value: \code{z}.}
#' @export
#'
#' @examples
#'
#' my0 <- f.multi.yearly(2020, 2)
#' #      this is a multi-year frequency that refers to the year 2020.
#' #      The next observation is expected in 2022 (not the next year).
#'
#' my0_value_str <-  as.character(my0) # this will be '2020'.
#' my0_class_str <- get.class.id(my0) # this will be 'z2'.
#'
#' my_new <- as.frequency("2020", "z3")
#' #      this is a multi-year frequency that refers to the year 2020.
#' #      However, the next observation is expected in 2023.
#'
#' # Don't make the following mistakes:
#' \donttest{
#' my_invalid <- try(f.multi.yearly(2020, 0))
#' my_invalid <- try(f.multi.yearly(2020, -5))
#' my_invalid <- try(as.frequency("2021", "z"))
#' }
#'
f.multi.yearly <- function(year, z){
  year = as.integer(year)
  z = as.integer(z)
  if (z <= 0)
    stop("Invalid 'z'. It should be a positive number.")
  res <- .F_multi_yearly(year, z)
  res
}



#' Create an \code{X-Times-A-Year} Frequency
#'
#' Use this function to create a frequency for time-series data that occurs \code{x} times every year.
#'
#' @param year An integer representing the year of the observation.
#' @param x An integer representing the number of observations in each year. It should be a positive integer.
#' @param position An integer representing the position of the current observation. It should be a positive integer and cannot be larger than \code{x}.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#:#"} (first # is the \code{year} and
#' the second # is the \code{position}; e.g., 2010:8/12 or 2010:10/10. Note that 2000:0/2 or 2000:13/12 are invalid.
#' \item **Class Id** \code{"y#"} (the number is the value: \code{x})
#' }
#'
#' @return An object of class \code{ldtf} which is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{year}{Determines the \code{year}.}
#' \item{x}{Determines the value: \code{x}.}
#' \item{position}{Determines the \code{position}.}
#'
#' @export
#' @examples
#'
#' xty0 <- f.x.times.a.year(2020, 3, 1)
#' #      this frequency divides the year 2020 into 3 partitions
#' #      and refers to the first partition.
#'
#' xty_value_str <-  as.character(xty0) # this will be '2020:1'.
#' xty_class_str <- get.class.id(xty0) # this will be 'y3'.
#'
#' xty_new <- as.frequency("2021:24", "z24")
#' #     this frequency divides the year 2021 into 24 partitions
#' #     and refers to the last partition.
#'
#' # Don't make the following mistakes:
#' \donttest{
#' xty_invalid <- try(f.x.times.a.year(2020, 3, 0))
#' xty_invalid <- try(f.x.times.a.year(2020, 24, 25))
#' xty_invalid <- try(as.frequency("2021:13", "y12"))
#' xty_invalid <- try(as.frequency("2021:0", "y1"))
#' xty_invalid <- try(as.frequency("2021", "y1"))
#' }
#'
f.x.times.a.year <- function(year, x, position){
  year = as.integer(year)
  x = as.integer(x)
  position = as.integer(position)

  if (x <= 0)
    stop("Invalid 'x'. It should be a positive number.")
  if (position <= 0 || position > x)
    stop("Invalid 'position'. It should be a positive number and less than the value of 'x'.")

  res <- .F_x_times_a_year(year, x, position)
  res
}



#' Create an \code{X-Times-Z-Years} Frequency
#'
#' Use this function to create a frequency for time-series data that occurs \code{x} times every \code{z} years.
#'
#' @param year An integer representing the year of the observation.
#' @param x An integer representing the number of partitions in each \code{z} years. It should be a positive integer.
#' @param z An integer representing the number of years. It should be a positive integer.
#' @param position An integer representing the position of the current observation. It should be a positive integer and cannot be larger than \code{x}.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#:#"} (Similar to \code{X-Times-A-Year}. Note that the string representation refers to the first year of the interval.)
#' \item **Class Id** \code{"x#z#"} (first '#' is the value: \code{x},
#' second '#' is the value: \code{z}; e.g., x23z4 means 23 times every 4 years)
#' }
#'
#' @return An object of class \code{ldtf}, which is also a list with the following members:
#' \item{class}{The class of this frequency.}
#' \item{year}{The \code{year}.}
#' \item{z}{The value: \code{z}.}
#' \item{x}{The value: \code{x}.}
#' \item{position}{The \code{position}.}
#'
#' @export
#'
#' @examples
#'
#' xtzy0 <- f.x.times.z.years(2020, 3, 2, 3)
#' # This frequency divides the year 2020 into 3 partitions and
#' # refers to the last partition. The next observation
#' # belongs to 2022 (not the next year).
#'
#' xtzy_value_str <- as.character(xtzy0) # This will be '2020:3'.
#' xtzy_class_str <- get.class.id(xtzy0) # This will be 'x3z2'.
#'
#' xtzy_new <- as.frequency("2021:3", "x3z4")
#' # This frequency divides the year 2021 into 3 partitions
#' # and refers to the last partition. The next observation occurs after 4 years.
#'
#' # Don't make the following mistakes:
#' \donttest{
#' xtzy_invalid <- try(f.x.times.z.years(2020, 3, 5, 0))
#' xtzy_invalid <- try(f.x.times.z.years(2020, 3, 0, 1))
#' xtzy_invalid <- try(as.frequency("2021:25", "x24y2"))
#' }
#'
f.x.times.z.years <- function(year, x, z, position){
  year = as.integer(year)
  z = as.integer(z)
  x = as.integer(x)
  position = as.integer(position)

  if (z <= 0)
    stop("Invalid 'z'. It should be a positive number.")
  if (x <= 0)
    stop("Invalid 'x'. It should be a positive number.")
  if (position <= 0 || position > x)
    stop("Invalid 'position'. It should be a positive number and less than the value of 'x'.")

  res <- .F_x_times_z_years(year, x, z, position)
  res
}


# returns date as list with year, month, day
get_date <- function(date){
  res <- list()
  if (is.null(date))
    stop("Invalid 'date'. It is null.")

  if (is.list(date)){
    if (is.null(date$year) || is.null(date$month) || is.null(date$day))
      stop("Date must contain year, month and day elements.")
    res <- date
  }
  else if (is.numeric(date) || is.integer(date)){
    if (length(date) != 3)
      stop("Invalid 'date' array. Its length must be 3.")
    res <- list(year=as.integer(date[[1]]),
                month = as.integer(date[[2]]),
                day=as.integer(date[[3]]))
  }
  else {
    date <- as.Date(date, tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%Y%m%d"))
    res <- list(year = as.integer(format(date, "%Y")),
                month = as.integer(format(date, "%m")),
                day = as.integer(format(date, "%d")))
  }

  if (res$year < 1400)
    stop("Invalid 'year'. It should be larger than 1400.")
  if (res$month <= 0 || res$month > 12)
    stop("Invalid 'month'. It should be between 1 and 12.")
  if (res$day <= 0 || res$day > 31)
    stop("Invalid 'day'.")

  res
}

get_date_reformat <- function(date){
  date <- get_date(date)
  res <- sprintf("%04d%02d%02d", date$year, date$month, date$day)
  res
}


#' Create a Weekly Frequency
#'
#' Use this function to create a frequency for time-series data that occurs weekly. The first day of the week is used as the reference.
#'
#' @param date The date, which can be a list with \code{year}, \code{month}, and \code{day} elements.
#' It can also be an integer array with 3 elements for year, month, and day respectively,
#' or an object that can be used as an argument for the \code{base::as.Date} function.
#' This date determines the start of the week.
#'
#' @details
#' To use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** The first day of the week in \code{"YYYYMMDD"} format.
#' \item **Class Id** \code{"w"}
#' }
#'
#' @return An object of class \code{ldtf}, which is also a list with the following members:
#' \item{class}{The class of this frequency.}
#' \item{year}{The \code{year}.}
#' \item{month}{The \code{month}.}
#' \item{day}{The \code{day}.}
#'
#' @export
#'
#' @examples
#'
#' w0 <- f.weekly(c(2023, 1, 2)) # This is 2/1/2023, which is Monday.
#'#    The next observation belongs to 9/1/2023.
#'
#' w0_value_str <-  as.character(w0) # this will be '20230102'.
#' w0_class_str <- get.class.id(w0) # this will be 'w'.
#'
#' w_new <- as.frequency("20230109", "w") # This is 9/1/2023.
#'
#' # Don't use invalid or unsupported dates:
#' \donttest{
#' w_invalid <- try(as.frequency("1399109", "w")) # this is a too old date and unsupported
#' w_invalid <- try(as.frequency("20230132", "w")) # invalid day in month
#' w_invalid <- try(as.frequency("20231331", "w")) # invalid month
#' }
#'
f.weekly <- function(date){

  date <- get_date(date)


  res <- .F_weekly(date$year, date$month, date$day)
  res
}



#' Create a Multi-Week Frequency
#'
#' Use this function to create a frequency for time-series data that occurs every 'k' weeks. The first day of the first week is used as the reference.
#'
#' @param date The date, which can be a list with \code{year}, \code{month}, and \code{day} elements.
#' It can also be an integer array with 3 elements for year, month, and day respectively,
#' or an object that can be used as an argument for the \code{base::as.Date} function.
#' @param k The number of weeks.
#'
#' @details
#' To use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** The first day of the first week in \code{"YYYYMMDD"} format.
#' \item **Class Id** \code{"w#"} (the number is the value of \code{k}; e.g., w3 means every 3 weeks)
#' }
#'
#' @return An object of class \code{ldtf}, which is also a list with the following members:
#' \item{class}{The class of this frequency.}
#' \item{year}{The \code{year}.}
#' \item{month}{The \code{month}.}
#' \item{day}{The \code{day}.}
#' \item{k}{The value of \code{k}.}
#'
#' @export
#'
#' @examples
#'
#' mw0 <- f.multi.weekly(c(2023, 1, 2), 3)
#' # This is 2/1/2023, which is Monday. The next observation belongs to 23/1/2023.
#'
#' mw0_value_str <- as.character(mw0) # This will be '20230102'.
#' mw0_class_str <- get.class.id(mw0) # This will be 'w3'.
#'
#' mw_new <- as.frequency("20230109", "w4") # This is 9/1/2023.
#'
#' # Don't use invalid or unsupported dates:
#' \donttest{
#' mw_invalid <- try(as.frequency("1399109", "w4")) # this is a too old date and unsupported
#' mw_invalid <- try(as.frequency("20230132", "w5")) # invalid day in month
#' mw_invalid <- try(as.frequency("20231331", "w2")) # invalid month
#' mw_invalid <- try(as.frequency("20231012", "w0"))
#' }
#'
f.multi.weekly <- function(date, k){

  date <- get_date(date)

  k = as.integer(k)


  if (k <= 0)
    stop("Invalid 'k'. It must be a positive integer.")

  res <- .F_multi_weekly(date$year, date$month, date$day, k)
  res
}





#' Create a Daily Frequency
#'
#' Use this function to create a frequency for time-series data that occurs daily.
#'
#' @param date The date, which can be a list with \code{year}, \code{month}, and \code{day} elements.
#' It can also be an integer array with 3 elements for year, month, and day respectively,
#' or an object that can be used as an argument for the \code{base::as.Date} function.
#'
#' @details
#' To use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"YYYYMMDD"} (similar to \code{Weekly})
#' \item **Class Id** \code{"d"}
#' }
#'
#' @return An object of class \code{ldtf}, which is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{year}{Determines the \code{year}.}
#' \item{month}{Determines the \code{month}.}
#' \item{day}{Determines the \code{day}.}
#'
#' @export
#' @examples
#'
#' d0 <- f.daily(c(2023, 1, 2)) # This is 2/1/2023. Next observation belongs to 3/1/2023.
#'
#' d0_value_str <-  as.character(d0) # this will be '20230102'.
#' d0_class_str <- get.class.id(d0) # this will be 'd'.
#'
#' d_new <- as.frequency("20230109", "d") # This is 9/1/2023.
#'
#' # Don't use invalid or unsupported dates:
#' \donttest{
#' # d_invalid <- try(as.frequency("1399109", "d")) # this is a too old date and unsupported
#' # d_invalid <- try(as.frequency("20230132", "d")) # invalid day in month
#' d_invalid <- try(as.frequency("20231331", "d")) # invalid month
#' }
#'
f.daily <- function(date){

  date <- get_date(date)

  res <- .F_daily(date$year, date$month, date$day)
  res
}

#' Create a Multi-Day Frequency
#'
#' Use this function to create a frequency for time-series data that occurs every \code{k} days. The first day of the interval is used as the reference.
#'
#' @param date The date, which can be a list with \code{year}, \code{month}, and \code{day} elements.
#' It can also be an integer array with 3 elements for year, month, and day respectively,
#' or an object that can be used as an argument for the \code{base::as.Date} function.
#' @param k The number of days in the interval.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format**: The first day of the interval in \code{"YYYYMMDD"} format.
#' \item **Class Id**: \code{"d#"} (where # is the value of \code{k}; e.g., d3 means every 3 days)
#' }
#'
#' @return An object of class \code{ldtf}. It is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{year}{Determines the \code{year}.}
#' \item{month}{Determines the \code{month}.}
#' \item{day}{Determines the \code{day}.}
#' \item{k}{Determines the value: \code{k}.}
#'
#' @export
#' @examples
#'
#' md0 <- f.multi.daily(c(2023, 1, 2), 4) # This is 2/1/2023. Next observation belongs to 6/1/2023.
#'
#' md0_value_str <-  as.character(md0) # this will be '20230102'.
#' md0_class_str <- get.class.id(md0) # this will be 'd4'.
#'
#' md_new <- as.frequency("20230109", "d") # This is 9/1/2023.
#'
#' # Don't use invalid or unsupported dates:
#' \donttest{
#' md_invalid <- try(as.frequency("1399109", "d3")) # this is a too old date and unsupported
#' md_invalid <- try(as.frequency("20230132", "d4")) # invalid day in month
#' md_invalid <- try(as.frequency("20231331", "d5")) # invalid month
#' }
#'
f.multi.daily <- function(date, k){

  date <- get_date(date)

  k = as.integer(k)


  if (k <= 0)
    stop("Invalid 'k'. It must be a positive integer.")

  res <- .F_multi_daily(date$year, date$month, date$day, k)
  res
}



#' Create a \code{Daily-In-Week} Frequency
#'
#' Use this function to create a frequency for time-series data that occurs daily within a subset of a week. The first day of the interval is used as the reference.
#'
#' @param date The date, which can be a list with \code{year}, \code{month}, and \code{day} elements.
#' It can also be an integer array with 3 elements for year, month, and day respectively,
#' or an object that can be used as an argument for the \code{base::as.Date} function.
#' @param weekStart The first day of the week, which can be \code{sun}, \code{mon}, \code{tue}, \code{wed}, \code{thu}, \code{fri}, or \code{sat}.
#' @param weekEnd The last day of the week, which can be one of the values listed for \code{weekStart}. Together, they define the week.
#' @param forward If the current date is not in the week and this value is true, it moves forward to the first day of the week. If this value is false, it moves backward to the last day of the week.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format**: The first day of the interval in \code{"YYYYMMDD"} format.
#' \item **Class Id**: \code{"i:...-..."} (where the first '...' represents \code{weekStart} and the second '...' represents \code{weekEnd}; e.g., \code{i:mon-fri} means a week from Monday to Friday)
#' }
#'
#' @return An object of class \code{ldtf}. It is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{year}{Determines the \code{year}.}
#' \item{month}{Determines the \code{month}.}
#' \item{day}{Determines the \code{day}.}
#' \item{weekStart}{Determines the \code{weekStart}.}
#' \item{weekEnd}{Determines the \code{weekEnd}.}
#' @export
#'
#' @examples
#'
#' dw0 <- f.daily.in.week(c(2023, 5, 16), "mon", "fri") # This is 16/5/2023.
#' dw0_value_str <-  as.character(dw0) # this will be '20230516'.
#' dw0_class_str <- get.class.id(dw0) # this will be 'i:mon-fri'.
#'
#' # Let's use the same date with another week definition:
#' dw1 <- f.daily.in.week(c(2023, 5, 16), "wed", "sat")
#' #     This is NOT 16/5/2023. It is 17/5/2023.
#' #     Since it was outside the week, we moved it forward.
#' dw2 <- f.daily.in.week(c(2023, 5, 16), "wed", "sat", FALSE)
#' #     This is 13/5/2023. The original day was outside the
#' #     week, but we moved backward too the end of
#' #     the previous week (which is Saturday).
#'
#' dw_new <- as.frequency("20230519", "i:sat-wed")
#' #     This is 20/1/2023 (by default, it moves forward).
#'
#' # Don't use invalid or unsupported dates:
#' \donttest{
#' dw_invalid <- try(as.frequency("1399109", "d3")) # this is a too old date and unsupported
#' dw_invalid <- try(as.frequency("20230132", "d4")) # invalid day in month
#' dw_invalid <- try(as.frequency("20231331", "d5")) # invalid month
#'
#' # don't use invalid week definitions:
#' dw_invalid <- try(f.daily.in.week(c(2023, 5, 16), "Wednesday", "sat"))
#' }
#'
f.daily.in.week <- function(date, weekStart = "mon",
                          weekEnd = "fri", forward = TRUE){

  date <- get_date(date)

  weekStart = as.character(weekStart)
  weekEnd = as.character(weekEnd)
  forward = as.logical(forward)

  valids <- c("sun", "mon", "tue", "wed", "thu", "fri", "sat")
  if (any(weekStart == valids) == FALSE)
    stop(paste0("Invalid 'weekStart'. It should be one of ", paste(dQuote(valids), collapse = ", "),"."))
  if (any(weekEnd == valids) == FALSE)
    stop(paste0("Invalid 'weekEnd'. It should be one of ", paste(dQuote(valids), collapse = ", "),"."))


  res <- .F_daily_in_week(date$year, date$month, date$day, weekStart,
                        weekEnd, forward)
  res
}







#' Create a \code{List-String} Frequency
#'
#' This frequency is typically used for labeled data. It is generally a list, but it can also be used to label observations outside this list.
#'
#' @param items The items in the list.
#' @param value The current item.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format**: \code{"..."} (where '...' represents the \code{value})
#' \item **Class Id**: \code{Ls} or \code{Ls:...} (where '...' represents the semi-colon-separated \code{items})
#' }
#'
#' @return An object of class \code{ldtf}, which is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{items}{Determines the \code{items}.}
#' \item{value}{Determines the \code{value}.}
#'
#' @export
#' @examples
#'
#' L0 <- f.list.string(c("A","B","C","D"), "C")
#'
#' L0_value_str <-  as.character(L0) # this will be 'C'.
#' L0_class_str <- get.class.id(L0) # this will be 'Ls:A;B;C;D'.
#'
#' L_new <- as.frequency("A", "Ls:A;B;C;D")
#' L_new0 <- as.frequency("A", "Ls") # compared to the previous one, its items will be empty
#'
#' # Don't make the following mistakes:
#' \donttest{
#' L_invalid <- try(as.frequency("E", "Ls:A;B;C;D")) # 'E' is not a member of the list
#' L_invalid <- try(f.list.string(c("A","B","C","D"), "E"))
#' }
#'
f.list.string <- function(items, value){
  value = as.character(value)
  items = as.character(items)

  res <- .F_list_string(items, value)
  res
}


#' Create a \code{List-Date} Frequency
#'
#' Use this frequency for data with date labels. It is generally a list of dates, but it can also be used to label observations outside this list.
#'
#' @param items The items in the list in \code{YYYYMMDD} format.
#' @param value The current value in \code{YYYYMMDD} format. If null, the first value in \code{items} is used.
#' @param reformat If the elements of \code{items} are not in \code{YYYYMMDD} format, set this to be \code{TRUE}.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format**: \code{"YYYYMMDD"} (i.e., the \code{item})
#' \item **Class Id**: \code{Ld} or \code{Ld:...} (where '...' represents the semi-colon-separated \code{items})
#' }
#'
#' @return An object of class \code{ldtf}. It is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{items}{Determines the \code{items}.}
#' \item{value}{Determines the \code{value}.}
#'
#' @export
#' @examples
#'
#' Ld0 <- f.list.date(c("20231101","20220903","20200823","20230303"), "20200823")
#'
#' Ld0_value_str <-  as.character(Ld0) # this will be '20200823'.
#' Ld0_class_str <- get.class.id(Ld0)
#' #      this will be 'Ld:20231101;20220903;20200823;20230303'.
#'
#' Ld_new <- as.frequency("20231101", "Ld:20231101;20220903;20200823;20230303")
#' Ld_new0 <- as.frequency("20231101", "Ld")
#' #     compared to the previous one, its items will be empty
#'
#' # Don't make the following mistakes:
#' \donttest{
#' Ld_invalid <- try(as.frequency("20231102", "Ld:20231101;20220903;20200823;20230303"))
#'   # 'E' is not a member of the list
#' Ld_invalid <- try(f.list.date(c("20231101","20220903","20200823","20230303"), "20231102"))
#' }
#'
f.list.date <- function(items, value = NULL, reformat = TRUE){

  if (is.null(value) == FALSE)
     value = as.character(value)
  items = as.character(items)

  if (reformat){
    items <- as.character(sapply(items, function(d)get_date_reformat(d)))
    if (is.null(value) == FALSE)
       value <- get_date_reformat(value)
  }
  if (is.null(value))
    value <- items[[1]]

  res <- .F_list_date(items, value)
  res
}







#' Create an 'Hourly' Frequency
#'
#' Use this function to create a frequency for time-series data that occurs hourly in a day or a subset of a week.
#'
#' @param day A 'Day-based' object of class \code{ldtf}, such as \code{Daily} or \code{Daily-In-Week}.
#' @param hour The index of the hour in the day, which should be between 1 and 24.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format**: \code{"YYYYMMDD:#"} (where # represents the value of \code{hour})
#' \item **Class Id**: \code{ho|...} (where '...' represents the 'class id' of \code{day})
#' }
#'
#' @return An object of class \code{ldtf}. It is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{day}{Determines the \code{day}.}
#' \item{hour}{Determines the \code{hour}.}
#' 
#' @export
#' @examples
#'
#' ho0 <- f.hourly(f.daily(c(2023,5,16)),4)
#'
#' ho0_value_str <-  as.character(ho0) # this will be '20230516:4'.
#' ho0_class_str <- get.class.id(ho0)
#' #    this will be 'ho|d'. The second part (i.e., 'd')
#' #    shows that this frequency is defined in a 'Daily' frequency.
#'
#' ho_new <- as.frequency("20231101:3", "ho|i:wed-sat")
#'
#' # Don't make the following mistakes:
#' \donttest{
#' ho_invalid <- try(as.frequency("20231101:3", "ho|j:wed-sat"))
#' #  invalid format in day-based frequency
#' ho_invalid <- try(f.hourly(f.daily(c(2023,5,16)),25)) # invalid hour
#' }
#'
f.hourly <- function(day, hour){
  hour = as.integer(hour)

  if (hour < 1 || hour > 24)
    stop("Invalid 'hour'. It should be between 1 and 24.")

  res <- .F_hourly(day, hour)
  res
}

#' Create a \code{Minute-ly} Frequency
#'
#' Use this function to create a frequency for time-series data that occurs every minute in a day or a subset of a week.
#'
#' @param day A 'Day-based' object of class \code{ldtf}, such as \code{Daily} or \code{Daily-In-Week}.
#' @param minute The index of the minute in the day, which should be between 1 and 1440.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format**: \code{"YYYYMMDD:#"} (where # represents the value of \code{minute})
#' \item **Class Id**: \code{mi|...} (where '...' represents the 'class id' of \code{day})
#' }
#'
#' @return An object of class \code{ldtf}. It is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{day}{Determines the \code{day}.}
#' \item{minute}{Determines the \code{minute}.}
#' 
#' @export
#' @examples
#'
#' mi0 <- f.minutely(f.daily(c(2023,5,16)),1200)
#'
#' mi0_value_str <-  as.character(mi0) # this will be '20230516:1200'.
#' mi0_class_str <- get.class.id(mi0)
#' #     this will be 'mi|d'. The second part (i.e., 'd')
#' #     shows that this frequency is defined in a 'Daily' frequency.
#'
#' mi_new <- as.frequency("20231101:3", "mi|i:wed-sat")
#'
#' # Don't make the following mistakes:
#' \donttest{
#' mi_invalid <- try(as.frequency("20231101:3", "mi|j:wed-sat"))
#' #   invalid format in day-based frequency
#' mi_invalid <- try(f.minutely(f.daily(c(2023,5,16)),2000)) # invalid minute
#' }
#'
f.minutely <- function(day, minute){
  minute = as.integer(minute)

  if (minute < 1 || minute > 1440)
    stop("Invalid 'minute'. It should be between 1 and 1440.")

  res <- .F_minutely(day, minute)
  res
}



#' Create a \code{Second-ly} Frequency
#'
#' Use this function to create a frequency for time-series data that occurs every second in a day or a subset of a week.
#'
#' @param day A 'Day-based' object of class \code{ldtf}, such as \code{Daily} or \code{Daily-In-Week}.
#' @param second The index of the second in the day, which should be between 1 and 86400.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format**: \code{"YYYYMMDD:#"} (where # represents the value of \code{second})
#' \item **Class Id**: \code{se|...} (where '...' represents the 'class id' of \code{day})
#' }
#'
#' @return An object of class \code{ldtf}. It is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{day}{Determines the \code{day}.}
#' \item{second}{Determines the \code{second}.}
#' 
#' @export
#' @examples
#'
#' se0 <- f.secondly(f.daily(c(2023,5,16)),40032)
#'
#' se0_value_str <-  as.character(se0) # this will be '20230516:40032'.
#' se0_class_str <- get.class.id(se0)
#' #     this will be 'se|d'. The second part (i.e., 'd') shows
#' #     that this frequency is defined in a 'Daily' frequency.
#'
#' se_new <- as.frequency("20231101:3", "se|i:wed-sat")
#'
#' # Don't make the following mistakes:
#' \donttest{
#' mi_invalid <- try(as.frequency("20231101:3", "se|j:wed-sat"))
#' #  invalid format in day-based frequency
#' mi_invalid <- try(f.secondly(f.daily(c(2023,5,16)),100000)) # invalid second
#' }
#'
f.secondly <- function(day, second){

  second = as.integer(second)

  if (second < 1 || second > 86400)
    stop("Invalid 'minute'. It should be between 1 and 86400.")

  res <- .F_secondly(day, second)
  res
}



#' Create an \code{X-Times-A-Day} Frequency
#'
#' Use this function to create a frequency for time-series data that occurs \code{x} times in a day or a subset of a week.
#'
#' @param day A 'Day-based' object of class \code{ldtf}, such as \code{Daily} or \code{Daily-In-Week}.
#' @param x The number of observations in each day.
#' @param position The position of the current observation, which should be a positive integer and cannot be larger than \code{x}.
#'
#' @details
#' In order to use the \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format**: \code{"#"} (where '#' represents the value of \code{position})
#' \item **Class Id**: \code{"da#|..."} (where '#' represents the value of \code{x} and '...' represents the 'class id' of \code{day})
#' }
#'
#' @return An object of class \code{ldtf}. It is also a list with the following members:
#' \item{class}{Determines the class of this frequency.}
#' \item{day}{Determines the \code{day}.}
#' \item{second}{Determines the \code{second}.}
#' 
#' @export
#' @examples
#'
#' xd0 <- f.x.times.a.day(f.daily(c(2023,5,16)),13, 12)
#'
#' xd0_value_str <-  as.character(xd0) # this will be '20230516:12'.
#' xd0_class_str <- get.class.id(xd0)
#' #      this will be 'da13|d'. The second part (i.e., 'd')
#' #      shows that this frequency is defined in a 'Daily' frequency.
#'
#' xd_new <- as.frequency("20231101:3", "da3|i:wed-sat")
#'
#' # Don't make the following mistakes:
#' \donttest{
#' xd_invalid <- try(as.frequency("20231101:3", "da|i:wed-sat"))
#' #  invalid format in day-based frequency
#' xd_invalid <- try(f.x.times.a.day(f.daily(c(2023,5,16)),4,0)) # invalid position
#' }
#'
f.x.times.a.day <- function(day, x, position){

  x = as.integer(x)
  position = as.integer(position)

  if (x <= 0)
    stop("Invalid 'x'. It should be a positive number.")
  if (position <= 0 || position > x)
    stop("Invalid 'position'. It should be a positive number and less than the value of 'x'.")


  res <- .F_x_times_a_day(day, x, position)
  res
}

