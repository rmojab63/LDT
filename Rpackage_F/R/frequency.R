
#' Creates a Cross-Section Frequency
#'
#' This frequency is typically used for indexed data. It is represented by an integer
#' that indicates the position of the observation.
#'
#' @param position Position of the observation
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#"} (the number is the \code{position})
#' \item **Class Id** \code{"cs"}
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{position} \tab Determines the \code{position}.
#' }
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



#' Creates an Annual Frequency
#'
#' Use it to create a frequency for time-series data that occurs annually.
#'
#' @param year Year of the observation
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#"} (the number is the \code{year})
#' \item **Class Id** \code{"y"}
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{year} \tab Determines the \code{year}.
#' }
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



#' Creates a Quarterly Frequency
#'
#' Use it to create a frequency for time-series data that occurs quarterly.
#'
#' @param year Year of the observation.
#' @param quarter Quarter of the observation (It should be between 1 and 4).
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#q#"} (first '#' is the \code{year}, second '#' is the
#'  \code{quarter}; e.g., 2010q3 or 2010q4. Note that 2000q0 or 2000q5 are invalid.
#' \item **Class Id** \code{"q"}
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{year} \tab Determines the \code{year}.\cr
#' \code{quarter} \tab Determines the \code{quareter}.
#' }
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
#' #   q_invalid <- f.quarterly(2020, 0)
#' #   q_invalid <- f.quarterly(2020, 5)
#' #   q_invalid <- as.frequency("2021q0", "q")
#' #   q_invalid <- as.frequency("2021q5", "q")
#' #   q_invalid <- as.frequency("2021", "q")
#'
f.quarterly <- function(year, quarter){
  year = as.integer(year)
  quarter = as.integer(quarter)
  if (quarter < 1 || quarter > 4)
    stop("Invalid 'quarter'. It should be between 1 and 4.")
  res <- .F_quarterly(year, quarter)
  res
}



#' Creates a Monthly Frequency
#'
#' Use it to create a frequency for time-series data that occurs monthly.
#'
#' @param year Year of the observation.
#' @param month Month of the observation (It should be between 1 to 12).
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#m#"} (first # is the \code{year}, second # is
#'  the \code{month} (1 to 12); e.g., 2010m8 or 2010m12. Note that 2000m0 or 2000m13 are invalid.
#' \item **Class Id** \code{"m"}
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{year} \tab Determines the \code{year}.\cr
#' \code{month} \tab Determines the \code{month}.
#' }
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
#' #   m_invalid <- f.monthly(2020, 0)
#' #   m_invalid <- f.monthly(2020, 5)
#' #   m_invalid <- as.frequency("2021m0", "m")
#' #   m_invalid <- as.frequency("2021m13", "m")
#' #   m_invalid <- as.frequency("2021", "m")
#'
f.monthly <- function(year, month){
  year = as.integer(year)
  month = as.integer(month)
  if (month < 1 || month > 12)
    stop("Invalid 'month'. It should be between 1 and 12.")
  res <- .F_monthly(year, month)
  res
}



#' Creates a Multi-Year Frequency
#'
#' Use it to create a frequency for time-series data that occurs every \code{z} years.
#'
#' @param year Year of the observation.
#' @param z Number of years. It should be larger than zero.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#"} (the number is the \code{year}, which means the string representation is the first year of the interval)
#' \item **Class Id** \code{"z#"} ('#' represents the value: \code{z}; e.g., z3 means every 3 years)
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{year} \tab Determines the \code{year}.\cr
#' \code{z} \tab Determines the value: \code{z}.
#' }
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
#' #   my_invalid <- f.multi.yearly(2020, 0)
#' #   my_invalid <- f.multi.yearly(2020, -5)
#' #   my_invalid <- as.frequency("2021", "z")
#'
f.multi.yearly <- function(year, z){
  year = as.integer(year)
  z = as.integer(z)
  if (z <= 0)
    stop("Invalid 'z'. It should be a positive number.")
  res <- .F_multi_yearly(year, z)
  res
}



#' Creates an \code{X-Times-A-Year} Frequency
#'
#' Use it to create a frequency for time-series data that occurs \code{x} times every year.
#'
#' @param year Year of the observation.
#' @param x Number of the observation in each year. It should be a positive integer.
#' @param position Position of the current observation. It should be a positive integer. It cannot be larger than \code{x}.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#:#"} (first # is the \code{year} and
#' the second # is the \code{position}; e.g., 2010:8/12 or 2010:10/10. Note that 2000:0/2 or 2000:13/12 are invalid.
#' \item **Class Id** \code{"y#"} (the number is the value: \code{x})
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{year} \tab Determines the \code{year}.\cr
#' \code{x} \tab Determines the value: \code{x}.\cr
#' \code{position} \tab Determines the \code{position}.
#' }
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
#' #   xty_invalid <- f.x.times.a.year(2020, 3, 0)
#' #   xty_invalid <- f.x.times.a.year(2020, 24, 25)
#' #   xty_invalid <- as.frequency("2021:13", "y12")
#' #   xty_invalid <- as.frequency("2021:0", "y1")
#' #   xty_invalid <- as.frequency("2021", "y1")
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



#' Creates an \code{X-Times-Z-Years} Frequency
#'
#' Use it to create a frequency for time-series data that occurs \code{x} times every \code{z} years.
#'
#' @param year Year of the observation.
#' @param x Number of partitions in each \code{z} years. It should be a positive integer.
#' @param z Number of years. It should be a positive integer.
#' @param position Position of the current observation. It should be a positive integer. It cannot be larger than \code{x}.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#:#"} (Similar to \code{X-Times-A-Year}. Note that the string representation refers to the first year of the interval.)
#' \item **Class Id** \code{"x#z#"} (first '#' is the value: \code{x},
#' second '#' is the value: \code{z}; e.g., x23z4 means 23 times every 4 years)
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{year} \tab Determines the \code{year}.\cr
#' \code{z} \tab Determines the value: \code{z}.\cr
#' \code{x} \tab Determines the value: \code{x}.\cr
#' \code{position} \tab Determines the \code{position}.
#' }
#'
#' @export
#'
#' @examples
#'
#' xtzy0 <- f.x.times.z.years(2020, 3, 2, 3)
#' #      this frequency divides the year 2020 into 3 partitions and
#' #      refers to the last partition. The next observation
#' #      belongs to 2022 (not the next year).
#'
#' xtzy_value_str <-  as.character(xtzy0) # this will be '2020:3'.
#' xtzy_class_str <- get.class.id(xtzy0) # this will be 'x3z2'.
#'
#' xtzy_new <- as.frequency("2021:3", "x3z4")
#' #      this frequency divides the year 2021 into 3 partitions
#' #      and refers to the last partition. The next observation occurs after 4 years.
#'
#' # Don't make the following mistakes:
#' #   xtzy_invalid <- f.x.times.z.years(2020, 3, 5, 0)
#' #   xtzy_invalid <- f.x.times.z.years(2020, 3, 0, 1)
#' #   xtzy_invalid <- as.frequency("2021:25", "x24y2")
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





#' Creates a Weekly Frequency
#'
#' Use it to create a frequency for time-series data that occurs weekly. We use first day of the week as the reference.
#'
#' @param year Year of the observation. It should be a valid year as an integer.
#' @param month Month of the observation. It should be a valid month of year as an integer.
#' @param day Day of the observation which will determine the start of the week. It should be a valid day of month as an integer.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** First day of the week in \code{"YYYYMMDD"} format.
#' \item **Class Id** \code{"w"}
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{year} \tab Determines the \code{year}.\cr
#' \code{month} \tab Determines the \code{month}.\cr
#' \code{day} \tab Determines the \code{day}.
#' }
#'
#' @export
#'
#' @examples
#'
#' w0 <- f.weekly(2023, 1, 2) # This is 2/1/2023 which is Monday. Next observation belongs to 9/1/2023.
#'
#' w0_value_str <-  as.character(w0) # this will be '20230102'.
#' w0_class_str <- get.class.id(w0) # this will be 'w'.
#'
#' w_new <- as.frequency("20230109", "w") # This is 9/1/2023.
#'
#' # Don't use invalid or unsupported dates:
#' #   w_invalid <- as.frequency("1399109", "w") # this is a too old date and unsupported
#' #   w_invalid <- as.frequency("20230132", "w") # invalid day in month
#' #   w_invalid <- as.frequency("20231331", "w") # invalid month
#'
f.weekly <- function(year, month, day){
  year = as.integer(year)
  month = as.integer(month)
  day = as.integer(day)

  if (year < 1400)
    stop("Invalid 'year'. It should be larger than 1400.")
  if (month <= 0 || month > 12)
    stop("Invalid 'month'. It should be between 1 and 12.")
  if (day <= 0 || day > 31)
    stop("Invalid 'day'.")

  res <- .F_weekly(year, month, day)
  res
}



#' Creates a Multi-Week Frequency
#'
#' Use it to create a frequency for time-series data that occurs every 'k' weeks. We use first day of the first week as the reference.
#'
#' @param year Year of the observation. It should be a valid year as an integer.
#' @param month Month of the observation. It should be a valid month of year as an integer.
#' @param day Day of the observation which will determine the start of the week. It should be a valid day of month as an integer.
#' @param k Number of weeks.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** First day of the first week in \code{"YYYYMMDD"} format.
#' \item **Class Id** \code{"w#"} (the number is value: \code{k}; e.g., w3 means every 3 weeks)
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{year} \tab Determines the \code{year}.\cr
#' \code{month} \tab Determines the \code{month}.\cr
#' \code{day} \tab Determines the \code{day}.\cr
#' \code{k} \tab Determines the value: \code{k}.
#' }
#' @export
#'
#' @examples
#'
#' mw0 <- f.multi.weekly(2023, 1, 2,3)
#' #      This is 2/1/2023 which is Monday. Next observation belongs to 23/1/2023.
#'
#' mw0_value_str <-  as.character(mw0) # this will be '20230102'.
#' mw0_class_str <- get.class.id(mw0) # this will be 'w3'.
#'
#' mw_new <- as.frequency("20230109", "w4") # This is 9/1/2023.
#'
#' # Don't use invalid or unsupported dates:
#' #   mw_invalid <- as.frequency("1399109", "w4") # this is a too old date and unsupported
#' #   mw_invalid <- as.frequency("20230132", "w5") # invalid day in month
#' #   mw_invalid <- as.frequency("20231331", "w2") # invalid month
#' #   mw_invalid <- as.frequency("20231012", "w0")
#'
f.multi.weekly <- function(year, month, day, k){

  year = as.integer(year)
  month = as.integer(month)
  day = as.integer(day)
  k = as.integer(k)

  if (year < 1400)
    stop("Invalid 'year'. It should be larger than 1400.")
  if (month <= 0 || month > 12)
    stop("Invalid 'month'. It should be between 1 and 12.")
  if (day <= 0 || day > 31)
    stop("Invalid 'day'.")
  if (k <= 0)
    stop("Invalid 'k'. It must be a positive integer.")

  res <- .F_multi_weekly(year, month, day, k)
  res
}






#' Creates a Daily Frequency
#'
#' Use it to create a frequency for time-series data that occurs daily.
#'
#' @param year Year of the observation. It should be a valid year as an integer.
#' @param month Month of the observation. It should be a valid month of year as an integer.
#' @param day Day of the observation. It should be a valid day of month as an integer.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"YYYYMMDD"} (similar to \code{Weekly})
#' \item **Class Id** \code{"d"}
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{year} \tab Determines the \code{year}.\cr
#' \code{month} \tab Determines the \code{month}.\cr
#' \code{day} \tab Determines the \code{day}.
#' }
#'
#' @export
#' @examples
#'
#' d0 <- f.daily(2023, 1, 2) # This is 2/1/2023. Next observation belongs to 3/1/2023.
#'
#' d0_value_str <-  as.character(d0) # this will be '20230102'.
#' d0_class_str <- get.class.id(d0) # this will be 'd'.
#'
#' d_new <- as.frequency("20230109", "d") # This is 9/1/2023.
#'
#' # Don't use invalid or unsupported dates:
#' #   d_invalid <- as.frequency("1399109", "d") # this is a too old date and unsupported
#' #   d_invalid <- as.frequency("20230132", "d") # invalid day in month
#' #   d_invalid <- as.frequency("20231331", "d") # invalid month
#'
f.daily <- function(year, month, day){

  year = as.integer(year)
  month = as.integer(month)
  day = as.integer(day)

  if (year < 1400)
    stop("Invalid 'year'. It should be larger than 1400.")
  if (month <= 0 || month > 12)
    stop("Invalid 'month'. It should be between 1 and 12.")
  if (day <= 0 || day > 31)
    stop("Invalid 'day'.")

  res <- .F_daily(year, month, day)
  res
}

#' Creates an Multi-Day Frequency
#'
#' Use it to create a frequency for time-series data that occurs every \code{k} days. The first day of the interval is used as the reference.
#'
#' @param year Year of the observation. It should be a valid year as an integer.
#' @param month Month of the observation. It should be a valid month of year as an integer.
#' @param day Day of the observation. It should be a valid day of month as an integer.
#' @param k Number of days.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** First day of the interval in \code{"YYYYMMDD"} format.
#' \item **Class Id** \code{"d#"} (the number is value: \code{k}; e.g., d3 means every 3 days)
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{year} \tab Determines the \code{year}.\cr
#' \code{month} \tab Determines the \code{month}.\cr
#' \code{day} \tab Determines the \code{day}.\cr
#' \code{k} \tab Determines the value: \code{k}.
#' }
#'
#' @export
#' @examples
#'
#' md0 <- f.multi.daily(2023, 1, 2, 4) # This is 2/1/2023. Next observation belongs to 6/1/2023.
#'
#' md0_value_str <-  as.character(md0) # this will be '20230102'.
#' md0_class_str <- get.class.id(md0) # this will be 'd4'.
#'
#' md_new <- as.frequency("20230109", "d") # This is 9/1/2023.
#'
#' # Don't use invalid or unsupported dates:
#' #   md_invalid <- as.frequency("1399109", "d3") # this is a too old date and unsupported
#' #   md_invalid <- as.frequency("20230132", "d4") # invalid day in month
#' #   md_invalid <- as.frequency("20231331", "d5") # invalid month
#'
f.multi.daily <- function(year, month, day, k){

  year = as.integer(year)
  month = as.integer(month)
  day = as.integer(day)
  k = as.integer(k)

  if (year < 1400)
    stop("Invalid 'year'. It should be larger than 1400.")
  if (month <= 0 || month > 12)
    stop("Invalid 'month'. It should be between 1 and 12.")
  if (day <= 0 || day > 31)
    stop("Invalid 'day'.")
  if (k <= 0)
    stop("Invalid 'k'. It must be a positive integer.")

  res <- .F_multi_daily(year, month, day, k)
  res
}



#' Creates an \code{Daily-In-Week} Frequency
#'
#' Use it to create a frequency for time-series data that occurs daily within a subset of a week. The first day of the interval is used as the reference.
#'
#' @param year Year of the observation. It should be a valid year as an integer.
#' @param month Month of the observation. It should be a valid month of year as an integer.
#' @param day Day of the observation which will determine the start of the week. It should be a valid day of month as an integer.
#' @param weekStart First day of the week. It can be \code{sun}, \code{mon}, \code{tue}, \code{wed}, \code{thu}, \code{fri}, and \code{sat}.
#' @param weekEnd Last day of the week. See \code{weekStart} for possible values. Together, they define the week
#' @param forward If the current date is not in the week and this value is true, it moves forward to the first day of the week. If this value is false, it moves backward to the last day of the week.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"YYYYMMDD"} (First day of the interval in \code{"YYYYMMDD"} format.)
#' \item **Class Id** \code{"i:...-..."} (the first '...' is \code{weekStart} and the second '...' is \code{weekEnd}; e.g., \code{i:mon-fri} means a week from Monday to Friday)
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{year} \tab Determines the \code{year}.\cr
#' \code{month} \tab Determines the \code{month}.\cr
#' \code{day} \tab Determines the \code{day}.\cr
#' \code{weekStart} \tab Determines the \code{weekStart}.\cr
#' \code{weekEnd} \tab Determines the \code{weekEnd}.
#' }
#' @export
#'
#' @examples
#'
#' dw0 <- f.daily.in.week(2023, 5, 16, "mon", "fri") # This is 16/5/2023.
#' dw0_value_str <-  as.character(dw0) # this will be '20230516'.
#' dw0_class_str <- get.class.id(dw0) # this will be 'i:mon-fri'.
#'
#' # Let's use the same date with another week definition:
#' dw1 <- f.daily.in.week(2023, 5, 16, "wed", "sat")
#' #     This is NOT 16/5/2023. It is 17/5/2023.
#' #     Since it was outside the week, we moved it forward.
#' dw2 <- f.daily.in.week(2023, 5, 16, "wed", "sat", FALSE)
#' #     This is 13/5/2023. The original day was outside the
#' #     week, but we moved backward too the end of
#' #     the previous week (which is Saturday).
#'
#' dw_new <- as.frequency("20230519", "i:sat-wed")
#' #     This is 20/1/2023 (by default, it moves forward).
#'
#' # Don't use invalid or unsupported dates:
#' #   dw_invalid <- as.frequency("1399109", "d3") # this is a too old date and unsupported
#' #   dw_invalid <- as.frequency("20230132", "d4") # invalid day in month
#' #   dw_invalid <- as.frequency("20231331", "d5") # invalid month
#'
#' # don't use invalid week definitions:
#' #   f.daily.in.week(2023, 5, 16, "Wednesday", "sat")
#'
f.daily.in.week <- function(year, month, day, weekStart = "mon",
                          weekEnd = "fri", forward = TRUE){

  year = as.integer(year)
  month = as.integer(month)
  day = as.integer(day)

  weekStart = as.character(weekStart)
  weekEnd = as.character(weekEnd)
  forward = as.logical(forward)

  if (year < 1400)
    stop("Invalid 'year'. It should be larger than 1400.")
  if (month <= 0 || month > 12)
    stop("Invalid 'month'. It should be between 1 and 12.")
  if (day <= 0 || day > 31)
    stop("Invalid 'day'.")

  valids <- c("sun", "mon", "tue", "wed", "thu", "fri", "sat")
  if (any(weekStart == valids) == FALSE)
    stop(paste0("Invalid 'weekStart'. It should be one of ", paste(dQuote(valids), collapse = ", "),"."))
  if (any(weekEnd == valids) == FALSE)
    stop(paste0("Invalid 'weekEnd'. It should be one of ", paste(dQuote(valids), collapse = ", "),"."))


  res <- .F_daily_in_week(year, month, day, weekStart,
                        weekEnd, forward)
  res
}







#' Creates an \code{List-String} Frequency
#'
#' This frequency is typically used for labeled data. It is generally a list, but it can be used to label observations outside this list.
#'
#' @param items Items of the list.
#' @param value Current item.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"..."} (in which ... is the \code{value})
#' \item **Class Id** \code{Ls} or \code{Ls:...} (in which ... is the semi-colon separated \code{items})
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{items} \tab Determines the \code{items}.\cr
#' \code{value} \tab Determines the \code{value}.
#' }
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
#' #   L_invalid <- as.frequency("E", "Ls:A;B;C;D") # 'E' is not a member of the list
#' #   L_invalid <- f.list.string(c("A","B","C","D"), "E")
#'
f.list.string <- function(items, value){
  value = as.character(value)
  items = as.character(items)

  res <- .F_list_string(items, value)
  res
}


#' Creates an \code{List-Date} Frequency
#'
#' Use this frequency for data with date labels. It is generally a list of dates, but it can be used to label observations outside this list.
#'
#' @param items Items of the in \code{YYYYMMDD} format.
#' @param value Current value in \code{YYYYMMDD} format.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"YYYYMMDD"} (i.e., \code{item})
#' \item **Class Id** \code{Ld} or \code{Ld:...} (in which '...' is the semi-colon separated \code{items})
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{items} \tab Determines the \code{items}.\cr
#' \code{value} \tab Determines the \code{value}.
#' }
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
#' #   Ld_invalid <- as.frequency("20231102", "Ld:20231101;20220903;20200823;20230303")
#' #     'E' is not a member of the list
#' #   Ld_invalid <- f.list.date(c("20231101","20220903","20200823","20230303"), "20231102")
#'
f.list.date <- function(items, value){

  value = as.character(value)
  items = as.character(items)

  res <- .F_list_date(items, value)
  res
}







#' Creates an 'Hourly' Frequency
#'
#' Use it to create a frequency for time-series data that occurs hourly in day or a subset of a week.
#'
#' @param day A 'Day-based' object of class \code{ldtf}, such as \code{Daily} or \code{Daily-In-Week}.
#' @param hour Index of hour in the day. It should be between 1 and 24.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"YYYYMMDD:#"} (the number is the value: \code{hour})
#' \item **Class Id** \code{ho|...} (the '...' is the 'Class String' of \code{day})
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{day} \tab Determines the \code{day}.\cr
#' \code{hour} \tab Determines the \code{hour}.
#' }
#' @export
#' @examples
#'
#' ho0 <- f.hourly(f.daily(2023,5,16),4)
#'
#' ho0_value_str <-  as.character(ho0) # this will be '20230516:4'.
#' ho0_class_str <- get.class.id(ho0)
#' #    this will be 'ho|d'. The second part (i.e., 'd')
#' #    shows that this frequency is defined in a 'Daily' frequency.
#'
#' ho_new <- as.frequency("20231101:3", "ho|i:wed-sat")
#'
#' # Don't make the following mistakes:
#' #   ho_invalid <- as.frequency("20231101:3", "ho|j:wed-sat")
#' #     invalid format in day-based frequency
#' #   ho_invalid <- f.hourly(f.daily(2023,5,16),25) # invalid hour
#'
f.hourly <- function(day, hour){
  hour = as.integer(hour)

  if (hour < 1 || hour > 24)
    stop("Invalid 'hour'. It should be between 1 and 24.")

  res <- .F_hourly(day, hour)
  res
}

#' Creates an \code{Minute-ly} Frequency
#'
#' Use it to create a frequency for time-series data that occurs every minute in day or a subset of a week.
#'
#' @param day A 'Day-based' object of class \code{ldtf}, such as \code{Daily} or \code{Daily-In-Week}.
#' @param minute Index of hour in the day. It should be between 1 and 1440.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"YYYYMMDD:#"} (the number is the value: \code{minute})
#' \item **Class Id** \code{mi|...} (the '...' is the 'Class String' of \code{day})
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{day} \tab Determines the \code{day}.\cr
#' \code{minute} \tab Determines the \code{minute}.
#' }
#' @export
#' @examples
#'
#' mi0 <- f.minutely(f.daily(2023,5,16),1200)
#'
#' mi0_value_str <-  as.character(mi0) # this will be '20230516:1200'.
#' mi0_class_str <- get.class.id(mi0)
#' #     this will be 'mi|d'. The second part (i.e., 'd')
#' #     shows that this frequency is defined in a 'Daily' frequency.
#'
#' mi_new <- as.frequency("20231101:3", "mi|i:wed-sat")
#'
#' # Don't make the following mistakes:
#' #   mi_invalid <- as.frequency("20231101:3", "mi|j:wed-sat")
#' #      invalid format in day-based frequency
#' #   mi_invalid <- f.minutely(f.daily(2023,5,16),2000) # invalid minute
#'
f.minutely <- function(day, minute){
  minute = as.integer(minute)

  if (minute < 1 || minute > 1440)
    stop("Invalid 'minute'. It should be between 1 and 1440.")

  res <- .F_minutely(day, minute)
  res
}



#' Creates an \code{Second-ly} Frequency
#'
#' Use it to create a frequency for time-series data that occurs every second in day or a subset of a week.
#'
#' @param day A 'Day-based' object of class \code{ldtf}, such as \code{Daily} or \code{Daily-In-Week}.
#' @param second Index of hour in the day. It should be between 1 and 86400.
#'

#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"YYYYMMDD:#"} (the number is the value: \code{second})
#' \item **Class Id** \code{se|...} (the '...' is the 'Class String' of \code{day})
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{day} \tab Determines the \code{day}.\cr
#' \code{second} \tab Determines the \code{second}.
#' }
#' @export
#' @examples
#'
#' se0 <- f.secondly(f.daily(2023,5,16),40032)
#'
#' se0_value_str <-  as.character(se0) # this will be '20230516:40032'.
#' se0_class_str <- get.class.id(se0)
#' #     this will be 'se|d'. The second part (i.e., 'd') shows
#' #     that this frequency is defined in a 'Daily' frequency.
#'
#' se_new <- as.frequency("20231101:3", "se|i:wed-sat")
#'
#' # Don't make the following mistakes:
#' #   mi_invalid <- as.frequency("20231101:3", "se|j:wed-sat")
#' #      invalid format in day-based frequency
#' #   mi_invalid <- f.secondly(f.daily(2023,5,16),100000) # invalid second
#'
f.secondly <- function(day, second){

  second = as.integer(second)

  if (second < 1 || second > 86400)
    stop("Invalid 'minute'. It should be between 1 and 86400.")

  res <- .F_secondly(day, second)
  res
}



#' Creates an \code{X-Times-A-Day} Frequency
#'
#' Use it to create a frequency for time-series data that occurs \code{x} times in day or a subset of a week.
#'
#' @param day A 'Day-based' object of class \code{ldtf}, such as \code{Daily} or \code{Daily-In-Week}.
#' @param x Number of observations in each day.
#' @param position Position of the current observation. It should be a positive integer. It cannot be larger than \code{x}.
#'
#' @details
#' In order to use \code{\link{as.frequency}} function for this type of frequency,
#' you need the following information:
#' \itemize{
#' \item **Character Format** \code{"#"} (in which '#' is the value: \code{position})
#' \item **Class Id** \code{"da#|..."} (in which '#' is the value: \code{x} and '...' is the 'Class String' of \code{day}))
#' }
#'
#' @return An object of class 'ldtf'. It is also a list with the following members:
#' \tabular{ll}{
#' \code{class} \tab Determines the class of this frequency.\cr
#' \code{day} \tab Determines the \code{day}.\cr
#' \code{second} \tab Determines the \code{second}.
#' }
#' @export
#' @examples
#'
#' xd0 <- f.x.times.a.day(f.daily(2023,5,16),13, 12)
#'
#' xd0_value_str <-  as.character(xd0) # this will be '20230516:12'.
#' xd0_class_str <- get.class.id(xd0)
#' #      this will be 'da13|d'. The second part (i.e., 'd')
#' #      shows that this frequency is defined in a 'Daily' frequency.
#'
#' xd_new <- as.frequency("20231101:3", "da3|i:wed-sat")
#'
#' # Don't make the following mistakes:
#' #   xd_invalid <- as.frequency("20231101:3", "da|i:wed-sat")
#' #      invalid format in day-based frequency
#' #   xd_invalid <- f.x.times.a.day(f.daily(2023,5,16),4,0) # invalid position
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

