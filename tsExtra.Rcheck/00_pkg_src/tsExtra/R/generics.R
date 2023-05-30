


#' Converts Frequency to Character
#'
#' The format is explained in \code{F_?} functions.
#'
#' @param frequency value of the frequency. It must be an \code{ldtf} object returned from \code{F_?} functions.
#'
#' @return String representation of the value of the frequency
#' @export
as.character.ldtf <- function(frequency){
  res <- .ToString_F(frequency)
  res
}

#' Gets Class 'Id' of a Frequency
#'
#' Use this function to get the 'id' of a frequency class.
#'
#' @param frequency The frequency. It must be an \code{ldtf} object returned from \code{F_?} functions.
#'
#' @details
#' You need this ‘id’ to convert the character back to the object. Some frequencies have a constant class id. For example, the class id for ‘monthly’ data is m. Some class ‘ids’ have parameters in them. Note that the format is explained in F_? functions.
#'
#'
#' @return A character that represents the class Id of this frequency.
#' @export
#' @examples
#'
#' freq <- F_XTimesADay(F_Daily(2023,5,16),13, 12)
#' freq_class_id <- get.class.id(freq) # this will be 'da13|d'.
#'
get.class.id <- function(frequency){
  res <- .ToClassString_F(frequency)
  res
}

#' Return Value and Class as String
#'
#' Similar to \code{as.character()} method. However, it returns the class id too.
#'
#' @param frequency value of the frequency. It must be an \code{ldtf}
#' object returned from \code{F_?} functions.
#'
#' @return A list with the following items:
#' \itemize{
#' \item **value:** String representation of the frequency. If you just want this, use \code{as.character()} function.
#' \item **day:** Class Id of this frequency. If you just want this, use \code{\link{get.class.id}} function.
#' \item **classType:** Type of the class.
#' }
#' @export
#'
#' @seealso \code{\link{get.class.id}}
#'
#' @examples
#'
#' freq <- F_XTimesADay(F_Daily(2023,5,16),13, 12)
#' freq_class_id <- get.character.info(freq) # this will be 'da13|d'.
#'
get.character.info <- function(frequency){
  res <- .ToString_F0(frequency)
  res
}

#' Converts Character to frequency
#'
#' Use this method to convert back your character to a frequency. You need the class id information.
#'
#' @param str value of the frequency as a valid character. You might have get it from \code{\link{as.character.ldtf}} function
#' @param classStr class id of the frequency as a valid character.
#'
#' @return An object of class 'ldtf'
#' @export
as.ldtf <-function(str, strclass){

  str= as.character(str)
  strclass = as.character(strclass)

  res <- .Parse_F(str, classStr)
  res
}


#' Generates a Sequence from a Range of Frequency
#'
#' Use it to generate a list of characters, where each element is a string representation of a frequency within the specified range.
#'
#' @param from First frequency of the range.
#' @param to Last frequency of the range.
#' @param by An integer that determines the increment of the sequence.
#'
#' @details
#' The two arguments \code{from} and \code{to} should be a valid frequencies (see \code{f.?} methods).
#' Also, they should be consistent. You cannot create a sequence in which one is e.g. monthly, and the other is yearly.
#'
#'
#' @return A list of characters that represents the sequence.
#' @export
#' @seealso \code{\link{seq0}}
#'
#' @examples
#' from <- F_Monthly(2020,1)
#' to <- F_Monthly(2021,12)
#' sequence1 <- seq(from, to, 1) # this will be '2020M1', '2020M2', ..., '2021M12'
#' sequence2 <- seq(from, to, 2) # this will be '2020M1', '2020M3', ..., '2021M11'
#' sequence3 <- seq(from, to, 3) # this will be '2020M1', '2020M4', ..., '2021M10'
#'
#' # backward:
#' sequence4 <- seq(to, from, -1) # this will be '2021M12', '2021M11', ..., '2020M1'
#'
seq <- function(from, to, by = 1){
  by = as.integer(by)
  if (by == 0)
    stop("Invalid 'by' value. It cannot be zero." )

  res <- .Sequence_F(from, to, by)
  res
}

#' Generates a Sequence from a Range of Frequency
#'
#' Use it to generate a list of characters, where each element is a string representation of a frequency within the specified range.
#'
#' @param start First frequency of the sequence .
#' @param length Length of the sequence
#' @param by increment of the sequence
#'
#' @return A list of characters that represents the sequence.
#' @export
#' @seealso \code{\link{seq}}
#'
#' @examples
#' start <- F_Monthly(2020,1)
#' sequence1 <- seq0(start, 24, 1) # this will be '2020M1', '2020M2', ..., '2021M12'
#' sequence2 <- seq0(start, 24, 2) # this will be '2020M1', '2020M3', ..., '2023M11'
#' sequence3 <- seq0(start, 24, 3) # this will be '2020M1', '2020M4', ..., '2025M10'
#'
#' # backward:
#' sequence4 <- seq0(start, 24, -1) # this will be '2020M1', '2019M12', ..., '2018M2'
#'
seq0 <- function(start, length, by = 1){

  length = as.integer(length)
  if (length < 0)
    stop("Invalid 'length' value. It cannot be negative." )

  by = as.integer(by)
  if (by == 0)
    stop("Invalid 'by' value. It cannot be zero." )

  res <- .Sequence_F0(start, length, by)
  res
}
