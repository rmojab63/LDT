
#' Gets Distances Between Variables
#'
#' @param data (numeric matrix) Data with variables in the columns.
#' @param distance (string) Determines how distances are calculated. It can be \code{correlation}, \code{absCorrelation}, \code{euclidean}, \code{manhattan}, \code{maximum}.
#' @param correlation (string) If \code{distance} is correlation, it determines the type of the correlation. It can be \code{pearson}, \code{spearman}.
#' @param checkNan (bool) If false, \code{NAN}s are not omitted.
#'
#' @return A symmetric matrix (lower triangle as a vector).
#'
#' @export
GetDistance <- function(data,
                        distance = "correlation",
                        correlation = "pearson",
                        checkNan = TRUE) {
  res <- .GetDistance(data, distance, correlation, checkNan)
  res
}

#' Hierarchical Clustering
#'
#'
#' @param distances (numeric vector) Determines the distances. This must be the lower triangle of a (symmetric) distance matrix (without the diagonal).
#' @param numVariables (int) Determines the number of variables. This should hold: '2 * length(\code{distances}) = \code{numVariables}(\code{numVariables} - 1)'.
#' @param linkage (string) Determines how Distances are calculated in a left-right node merge. It can be \code{single}, \code{complete}, \code{uAverage}, \code{wAverage}, \code{ward}.
#'
#' @return A list:
#' \item{merge}{(integer matrix)}
#' \item{height}{(numeric vector)}
#' \item{order}{(integer vector)}
#'
#' @export
ClusterH <- function(distances, numVariables,
                     linkage = "single"){
  res <- .ClusterH(distances, numVariables, linkage)
  res
}


#' Groups Variables with Hierarchical Clustering
#'
#' @details The results might be different from R's 'cutree' function. I don't know how 'cutree' works, but here I iterate over the nodes and whenever a split occurs, I add a group until the required number of groups is reached.
#'
#' @param data (numeric matrix) Data with variables in the columns.
#' @param nGroups (int) Number of groups
#' @param threshold (double) A threshold for omitting variables. If distance between two variables in a group is less than this value, the second one will be omitted. Note that a change in the order of the columns might change the results.
#' @param distance (string)  Determines how distances are calculated. It can be \code{correlation}, \code{absCorrelation}, \code{euclidean}, \code{manhattan}, \code{maximum}.
#' @param linkage (string) Determines how Distances are calculated in a left-right node merge. It can be \code{single}, \code{complete}, \code{uAverage}, \code{wAverage}, \code{ward}.
#' @param correlation (string) If \code{distance} is correlation, it determines the type of the correlation. It can be \code{pearson}, \code{spearman}.
#'
#' @return A list:
#' \item{groups}{(of integer vectors) indexes of variables in each group.}
#' \item{removed}{(integer vector) indexes of removed variables.}
#'
#' @export
ClusterHGroup <- function(data, nGroups = 2, threshold = 0,
                          distance = "correlation",
                          linkage = "single",
                          correlation = "pearson")

{
  res <- .ClusterHGroup(data, nGroups, threshold,
                       distance, linkage, correlation)
  res
}
