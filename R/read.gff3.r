#' @title Read gff3 file
#' @description reads a gff3 file from TransDecoder output
#' @param file gff3 file
#' @details This function read an gff3 file from TransDecoder output and rename columns as:
#' \describe{
#' \item{seqid}{sequence ID}
#' \item{source}{algorithm or operating procedure that generated the feature}
#' \item{feature}{The type of the feature}
#' \item{strt}{The start of the feature}
#' \item{end}{The end of the feature}
#' \item{score}{The score of the feature}
#' \item{st}{strand (+ or -)}
#' \item{phase}{For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame.}
#' \item{att}{attributes}
#' }
#' @export
read.gff3 <- function(file){ read.table(file, col.names = c("seqid", "source","feature","strt","end", "score","st","phase","att"))
}
