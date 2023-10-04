

is.number <- function(x){
  is.numeric(x) && length(x) == 1
}
is.positive.number <- function(x){
  is.number(x) && x > 0
}
is.negative.number <- function(x){
  is.number(x) && x < 0
}
is.zero.or.positive.number <- function(x){
  is.number(x) && x >= 0
}
is.zero.or.negative.number <- function(x){
  is.numeric(x) && length(x) == 1 && x <= 0
}


