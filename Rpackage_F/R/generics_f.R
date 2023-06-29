
#' Print a Frequency
#'
#' @param x A frequency which is the output of \code{f.?} functions in this package.
#' @param ... Additional arguments
#'
#' @return \code{NULL}
#' @export
print.ldtf <- function(x, ...) {
  if (is.null(x)) {
    stop("argument is null.")
  }
  if (any(class(x) == "ldtf") == FALSE) {
    stop("Invalid class")
  }
  s <- .ToString_F0(x)

  cat("Frequency: ", s$value, " (", s$classType, ": ", s$class, ")", sep = "")
  return(NULL)
}

#' Convert Frequency to Character
#'
#' This function converts a frequency to its string representation. The format is explained in the \code{f.?} functions.
#'
#' @param x The value of the frequency, which must be an \code{ldtf} object returned from the \code{f.?} functions.
#' @param ... Additional arguments.
#'
#' @return A string representation of the value of the frequency.
#' @export
as.character.ldtf <- function(x, ...){
  res <- .ToString_F(x)
  res
}

#' Get the Class Id of a Frequency
#'
#' Use this function to get the 'id' of a frequency class.
#'
#' @param frequency The frequency, which must be an \code{ldtf} object returned from the \code{f.?} functions.
#'
#' @details
#' You need this 'id' to convert the character back to the object. Some frequencies have a constant class id, such as 'm' for 'monthly' data. Some class 'ids' have parameters in them. Note that the format is explained in the \code{f.?} functions.
#'
#' @return A character string that represents the class id of this frequency.
#' @export
#' @examples
#'
#' freq <- f.x.times.a.day(f.daily(c(2023,5,16)),13, 12)
#' freq_class_id <- get.class.id(freq) # this will be 'da13|d'.
#'
get.class.id <- function(frequency){
  res <- .ToClassString_F(frequency)
  res
}

#' Convert Frequency to Character and Class Id
#'
#' This function returns the output of the [as.character.ldtf] and [get.class.id] functions.
#'
#' @param frequency The value of the frequency, which must be an \code{ldtf} object returned from the \code{f.?} functions.
#'
#' @return A list with the following items:
#' \itemize{
#' \item **value**: The string representation of the frequency. If you only want this, use the \code{as.character()} function.
#' \item **day**: The class Id of this frequency. If you only want this, use the \code{\link{get.class.id}} function.
#' \item **classType**: The type of the class.
#' }
#' @export
#'
#' @seealso \code{\link{get.class.id}}
#'
#' @examples
#'
#' freq <- f.x.times.a.day(f.daily(c(2023,5,16)),13, 12)
#' freq_class_id <- get.class.id0(freq)
#'
#' freq1 <- f.monthly(2020,3)
#' freq1_class_id <- get.class.id0(freq1)
get.class.id0 <- function(frequency){
  res <- .ToString_F0(frequency)
  res
}

#' Convert Character String to Frequency
#'
#' Use this function to convert a character string back to a frequency. You need the class id information.
#'
#' @param str The value of the frequency as a valid character, which you might have obtained from the \code{\link{as.character.ldtf}} function.
#' @param classId The class id of the frequency. These are explained in \code{f.?} functions.
#'
#' @return A frequency, which is an object of class 'ldtf'. See the \code{f.?} functions.
#' @export
as.frequency <-function(str, classId){

  str= as.character(str)
  classId = as.character(classId)

  res <- .Parse_F(str, classId)
  res
}


#' Generate a Sequence from a Range of Frequencies
#'
#' Use this function to generate a list of character strings, where each element is a string representation of a frequency within the specified range.
#'
#' @param from The first frequency of the sequence.
#' @param to The last frequency of the sequence.
#' @param by An integer that determines the increment of the sequence.
#'
#' @details
#' The two arguments \code{from} and \code{to} should be valid frequencies (see the \code{f.?} functions).
#' They should also be consistent; you cannot create a sequence in which one is, for example, monthly and the other is yearly.
#'
#' @return A list of character strings that represents the sequence.
#' @export
#' @seealso \code{\link{get.seq0}}
#'
#' @examples
#' from <- f.monthly(2020,1)
#' to <- f.monthly(2021,12)
#' sequence1 <- get.seq(from, to, 1) # this will be '2020M1', '2020M2', ..., '2021M12'
#' sequence2 <- get.seq(from, to, 2) # this will be '2020M1', '2020M3', ..., '2021M11'
#' sequence3 <- get.seq(from, to, 3) # this will be '2020M1', '2020M4', ..., '2021M10'
#'
#' # backward:
#' sequence4 <- get.seq(to, from, -1) # this will be '2021M12', '2021M11', ..., '2020M1'
#'
get.seq <- function(from, to, by = 1){
  by = as.integer(by)
  if (by == 0)
    stop("Invalid 'by' value. It cannot be zero." )

  res <- .Sequence_F(from, to, by)
  res
}

#' Generate a Sequence from a Range of Frequencies
#'
#' Use this function to generate a list of character strings, where each element is a string representation of a frequency within the specified range.
#'
#' @param start The first frequency of the sequence.
#' @param length The length of the sequence.
#' @param by An integer that determines the increment of the sequence.
#'
#' @return A list of character strings that represents the sequence.
#' @export
#' @seealso \code{\link{get.seq}}
#'
#' @examples
#' start <- f.monthly(2020,1)
#' sequence1 <- get.seq0(start, 24, 1) # this will be '2020M1', '2020M2', ..., '2021M12'
#' sequence2 <- get.seq0(start, 24, 2) # this will be '2020M1', '2020M3', ..., '2023M11'
#' sequence3 <- get.seq0(start, 24, 3) # this will be '2020M1', '2020M4', ..., '2025M10'
#'
#' # backward:
#' sequence4 <- get.seq0(start, 24, -1) # this will be '2020M1', '2019M12', ..., '2018M2'
#'
#' # Lists are a little different:
#' start_l <- f.list.string(c("A","B","C","D"), "C")
#' sequence5 <- get.seq0(start_l, 5, 1) # this will be 'C', 'D', 'out_item:1', ..., 'out_item:3'
#'
get.seq0 <- function(start, length, by = 1){

  length = as.integer(length)
  if (length < 0)
    stop("Invalid 'length' value. It cannot be negative." )

  by = as.integer(by)
  if (by == 0)
    stop("Invalid 'by' value. It cannot be zero." )

  res <- .Sequence_F0(start, length, by)
  res
}


#' Get Next Frequency
#'
#' Use this function to get the next frequency.
#'
#' @param freq A frequency.
#' @param count Determines the number of steps. If negative,
#'   it returns the previous frequency.
#'
#' @return The next frequency after the given frequency.
#' @export
#'
#' @examples
#' f <- f.yearly(2000)
#' fn <- next.freq(f, 10) # this is 2010
next.freq <- function(freq, count){
  count = as.integer(count)
  res <- .F_Next(freq, count)
  res
}

#' Get Interval between two frequencies
#'
#' Use this function to get the number of intervals between two frequencies.
#'
#' @param freq1 The first frequency.
#' @param freq2 The second frequency.
#'
#' @return The number of intervals between the two frequencies (\code{freq1} - \code{freq2}).
#' @export
#'
#' @examples
#' f1 <- f.yearly(2000)
#' f2 <- f.yearly(2010)
#' count <- minus.freqs(f1, f2) # this is -10
#' count <- minus.freqs(f2, f1) # this is 10
minus.freqs <- function(freq1, freq2){
  res <- .F_Minus(freq1, freq2)
  res
}
