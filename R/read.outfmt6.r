#' @title Read outfmt6 file
#' @description reads a annotations file in output format 6
#' @param file outfmt6 file
#' @details This function read a annotations file in output format 6 and rename columns as:
#' \describe{
#' \item{seqid}{query (e.g., unknown gene) sequence id}
#' \item{sseqid}{subject (e.g., reference genome) sequence id}
#' \item{pident}{percentage of identical matches}
#' \item{length}{alignment length (sequence overlap)}
#' \item{mismatch}{number of mismatches}
#' \item{gapopen}{ number of gap openings}
#' \item{qstart}{start of alignment in query}
#' \item{qend}{end of alignment in query}
#' \item{sstart}{start of alignment in subject}
#' \item{send}{end of alignment in subject}
#' \item{evalue}{expect value}
#' \item{bitscore}{bit score}
#' }
#' @export
read.outfmt6 <- function(file){ read.table(file, col.names = c("seqid","sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
}
