#' @import stringr
#' @title Split vector
#' @description Splits a vector according to a pattern into fragments and returns the fragment at the fragPos position
#' @param vec vector to split
#' @param pattern character that indicates the pattern by which the vector will be divided
#' @param fragPos the position of the fragment to return
#' @export
split.vec <- function(vec, pattern, fragPos=1){
  lista <- stringr::str_split(vec, pattern)
  vectorFrags <- vector()

  for( i in 1:length(lista)){
    vectorFrags[i] <-  lista[[i]][fragPos]
  }
  vectorFrags
}
