#' @title Value Non-Matching (inverse of %in% )
#' @description returns a vector of the positions of (first) non-matches of its first argument in its second.
#' @param x vector or \code{NULL}: the values to be matched
#' @param y vector or \code{NULL}: the values to be matched against
#' @details inverse of %in%
#' @export
'%!in%' <- function(x,y)!('%in%'(x,y))
