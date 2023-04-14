
#' Creates a Cross-Section Frequency
#'
#' This frequency is generally for indexed (or, non-time-series) data. It is an integer that represents the position of the observation.
#'
#' @param position Position of the observation
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"#"} (number is \code{position})
#' \item **Class String:** \code{"cs"}
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_CrossSection <- function(position) {
  res <- .F_CrossSection(position)
  res
}



#' Creates a \code{Yearly} Frequency
#'
#' Frequency for a series that happens every year
#'
#' @param year Year of the observation
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"#"} (number is \code{year})
#' \item **Class String:** \code{"y"}
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_Yearly <- function(year) {
  res <- .F_Yearly(year)
  res
}



#' Creates a \code{Quarterly} Frequency
#'
#' Frequency for a series that happens every quarter
#'
#' @param year Year of the observation
#' @param quarter Quarter of the observation (1 to 4)
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"#q#"} (first # is \code{year}, second # is
#'  \code{quarter}; e.g., 2010q3 or 2010q4.
#'  Note that 2000q0 or 2000q5 are invalid.
#' \item **Class String:** \code{"q"}
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_Quarterly <- function(year, quarter){
  res <- .F_Quarterly(year, quarter)
  res
}



#' Creates a \code{Monthly} Frequency
#'
#' Frequency for a series that happens every month
#'
#' @param year Year of the observation
#' @param month Month of the observation
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"#m#"} (first # is the \code{year}, second # is
#'  \code{month} (1 to 12); e.g., 2010m8 or 2010m12.
#'  Note that 2000m0 or 2000m13 are invalid.
#' \item **Class String:** \code{"m"}
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_Monthly <- function(year, month){
  res <- .F_Monthly(year, month)
  res
}



#' Creates a \code{Multi-Yearly} Frequency
#'
#' Frequency for a series that happens every \code{z} years
#'
#' @param year Year of the observation
#' @param z Number of years
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"#"} (similar to \code{Yearly})
#' \item **Class String:** \code{"z#"} (integer represents the
#' \code{z}; e.g., z3)
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_MultiYearly <- function(year, z){
  res <- .F_MultiYearly(year, z)
  res
}



#' Creates an \code{X-Times-A-Year} Frequency
#'
#' Frequency for a series that happens \code{x} times every year
#'
#' @param year Year of the observation
#' @param x Number of observation in each year
#' @param position Position of the current observation
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"#:#"} (first # is \code{year} and second # is
#'  \code{position};
#'  e.g., 2010:8/12 or 2010:10/10.
#'  Note that 2000:0/2 or 2000:13/12 are invalid.
#' \item **Class String:** \code{"y#"} (the number is \code{x})
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_XTimesAYear <- function(year, x, position){
  res <- .F_XTimesAYear(year, x, position)
  res
}



#' Creates an \code{X-Times-Z-Years} Frequency
#'
#' Frequency for a series that happens \code{x} times each \code{z} years
#'
#' @param year Year of the observation
#' @param x Number of partitons in each z years
#' @param z Number of years
#' @param position Position of the current observation
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"#:#"} (Similar to \code{X-Times-A-Year})
#' \item **Class String:** \code{"x#z#"} (first # is \code{x}, second # is \code{z};
#' e.g., x23z4 means 23 times every 4 years)
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_XTimesZYear <- function(year, x, z, position){
  res <- .F_XTimesZYear(year, x, z, position)
  res
}

#' Creates a \code{Weekly} Frequency
#'
#' Frequency for a series that happens every week
#'
#' @param year Year of the observation
#' @param month Month of the observation
#' @param day Day of the observation. It points to the first day of the week
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"YYYYMMDD"} (\code{YYYY} is the \code{year}, \code{MM} is \code{month} and \code{DD} is \code{day})
#' \item **Class String:** \code{"w"}
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_Weekly <- function(year, month, day){
  res <- .F_Weekly(year, month, day)
  res
}



#' Creates a \code{Multi-Weekly} Frequency
#'
#' Frequency for a series that happens every 'k' weeks
#'
#' @param year Year of the observation
#' @param month Month of the observation
#' @param day First day of the observation. It points to the first day of the week
#' @param k Number of weeks
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"YYYYMMDD"} (similar to \code{Weekly})
#' \item **Class String:** \code{"w#"} (the number is \code{k}; e.g., w3 means every 3 weeks)
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_MultiWeekly <- function(year, month, day, k){
  res <- .F_MultiWeekly(year, month, day, k)
  res
}



#' Creates a \code{Daily} Frequency
#'
#' Frequency for a series that happens every day
#'
#' @param year Year of the observation
#' @param month Month of the observation
#' @param day Day of the observation.
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"YYYYMMDD"} (similar to \code{Weekly})
#' \item **Class String:** \code{"d"}
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_Daily <- function(year, month, day){
  res <- .F_Daily(year, month, day)
  res
}

#' Creates an \code{Multi-Daily} Frequency
#'
#' Frequency for a series that happens every \code{k} days
#'
#' @param year Year of the observation
#' @param month Month of the observation
#' @param day First day of the observation
#' @param k Number of the days
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"YYYYMMDD"} (similar to \code{Weekly})
#' \item **Class String:** \code{"d#"} (the number is \code{k})
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_MultiDaily <- function(year, month, day, k){
  res <- .F_MultiDaily(year, month, day, k)
  res
}



#' Creates an \code{Daily-In-Week} Frequency
#'
#' Frequency for a series that happens every in the days of a week
#'
#' @param year Year of the observation
#' @param month Month of the observation
#' @param day First day of the observation
#' @param weekStart First day of the week. It can be \code{sun}, \code{mon},
#' \code{tue}, \code{wed}, \code{thu}, \code{fri}, and \code{sat}
#' @param weekEnd Last day of the week. See \code{weekStart}.
#' Together, they define the week
#' @param forward If current date in not in the week,
#' if true, it moves forward to the first day of the week.
#' Otherwise, it moves backward to the last day of the week.
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"YYYYMMDD"} (similar to \code{Weekly})
#' \item **Class String:** \code{"i:...-..."} (the first ... is \code{weekStart}
#' and the second ... is \code{weekEnd}; e.g., \code{i:mon-fri} means
#' a week that is from Monday to Friday)
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_DailyInWeek <- function(year, month, day, weekStart,
                          weekEnd, forward){
  res <- .F_DailyInWeek(year, month, day, weekStart,
                        weekEnd, forward)
  res
}



#' Creates an \code{List-String} Frequency
#'
#' Frequency for a series that is labeled by string
#'
#' @param items Items of the list
#' @param value Current item
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"..."} (in which ... is the \code{value})
#' \item **Class String:** \code{Ls} or \code{Ls:...} (in which ...
#' is the semi-colon separated \code{items})
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_ListString <- function(items, value){
  res <- .F_ListString(items, value)
  res
}



#' Creates an \code{List-Date} Frequency
#'
#' Frequency for a series that is labeled by dates
#'
#' @param items Items of the in string format: \code{YYYYMMDD}
#' @param value Current value in string format: \code{YYYYMMDD}
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"YYYYMMDD"} (i.e., \code{item})
#' \item **Class String:** \code{Ld} or \code{Ld:...} (in which ...
#' is the semi-colon separated \code{items})
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_ListDate <- function(items, value){
  res <- .F_ListDate(items, value)
  res
}

#' Creates an 'Hourly' Frequency
#'
#' Frequency for a series that happens every hour
#'
#' @param day A 'Day-based' frequency such as \code{Daily} or \code{Daily-In-Week}
#' @param hour Index of hour in the day (1 to 24)
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"YYYYMMDD:#"} (the number is \code{hour})
#' \item **Class String:** \code{ho|...} (the ... is the 'Class String' of \code{day})
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_Hourly <- function(day, hour){
  res <- .F_Hourly(day, hour)
  res
}

#' Creates an 'Minute-ly' Frequency
#'
#' Frequency for a series that happens every minute
#'
#' @param day A 'Day-based' frequency such as daily or daily-in-week
#' @param minute Index of Minute in the day (1 to 1440)
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"YYYYMMDD:#"} (the number is \code{minute})
#' \item **Class String:** \code{mi|...} (the ... is the 'Class String' of \code{day})
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_Minute_ly <- function(day, minute){
  res <- .F_Minute_ly(day, minute)
  res
}



#' Creates an 'Second-ly' Frequency
#'
#' Frequency for a series that happens every second
#'
#' @param day A 'Day-based' frequency such as daily or daily-in-week
#' @param second Index of second in the day (1 to 86400)
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"YYYYMMDD:#"} (the number is \code{second})
#' \item **Class String:** \code{se|...} (the ... is the 'Class String' of \code{day})
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_Second_ly <- function(day, second){
  res <- .F_Second_ly(day, second)
  res
}



#' Creates an 'X-Times-A-Day' Frequency
#'
#' Frequency for a series that happens x times in a day
#'
#' @param day A 'Day-based' frequency such as daily or daily-in-week
#' @param x Number of observations in a day
#' @param position Current position
#'
#' @details
#' \itemize{
#' \item **Value String:** \code{"#"} (the number is \code{hour})
#' \item **Class String:** \code{"da#|..."} (the number is \code{x}
#' and ... is the 'Class String' of \code{day}))
#' }
#'
#' @return An object of class 'ldtf'
#' @export
F_XTimesADay <- function(day, x, position){
  res <- .F_XTimesADay(day, x, position)
  res
}


#' Converts an \code{ldtf} Object to String
#'
#' The format is explained in \code{F_?} functions.
#'
#' @param value value of the frequency. It must be an \code{ldtf}
#' object returned from \code{F_?} functions.
#'
#' @return An object of class 'ldtf'
#' @export
ToString_F <- function(value){
  res <- .ToString_F(value)
  res
}


#' Converts an \code{ldtf} Object to String
#'
#' The format is explained in \code{F_?} functions.
#'
#' @param value value of the frequency. It must be an \code{ldtf}
#' object returned from \code{F_?} functions.
#'
#' @return An object of class 'ldtf'
#' @export
ToClassString_F <- function(value){
  res <- .ToClassString_F(value)
  res
}

#' Similar to \code{ToString_F} and Return Value and Class as String
#'
#' The format is explained in \code{F_?} functions.
#'
#' @param value value of the frequency. It must be an \code{ldtf}
#' object returned from \code{F_?} functions.
#'
#' @return An object of class 'ldtf'
#' @export
ToString_F0 <- function(value){
  res <- .ToString_F0(value)
  res
}

#' Converts back a String to \code{ldtf} Object
#'
#' The format is explained in \code{F_?} functions.
#'
#' @param str value of the frequency. It must be an \code{ldtf}
#' object returned from \code{F_?} functions.
#' @param classStr class of the frequency
#'
#' @return An object of class 'ldtf'
#' @export
Parse_F <- function(str, classStr){
  res <- .Parse_F(str, classStr)
  res
}

#' Generates a Sequence for a frequency
#'
#'
#' @param start first element of the sequence. It must be an \code{ldtf}
#' object returned from \code{F_?} functions.
#' @param length Length of the sequence
#'
#' @return A of strings
#' @export
Sequence_F <- function(start, length){
  res <- .Sequence_F(start, length)
  res
}
