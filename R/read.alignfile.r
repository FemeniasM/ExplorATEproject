#' @title Read .align file from RepeatMasker
#' @description read .align file from RepeatMasker
#' @details reads the RepeatMasker .align file and returns a data.frame with the sequence ID, the repetitions class, and the sequence divergence respect to the reference sequence in Kimura 2-parameter distances. If \code{average = TRUE} return a vector with means by sequences (\code{by = "seqID"}) or by repetition class (\code{by = "classRep"}) as indicated.
#' @param file .align file
#' @param average logical vector indicating if the mean should be taken
#' @param by if \code{average = TRUE}, then it must be indicated if the average is taken by sequences, \code{by = "seqID"}, or by repetition class, \code{by = "classRep"}.
#' @export
read.alignfile <- function(file, average=F, by=c("seqID", "classRep")){
  a <- system(paste0("file=",file,"; ",readalin_sh), intern = T)
  ALN <- data.frame(seqID = as.character(split.vec(a[!grepl("Kimura",a)], "\\ ", 5)),
                    beg = as.character(split.vec(a[!grepl("Kimura",a)], "\\ ", 6)),
                    end = as.character(split.vec(a[!grepl("Kimura",a)], "\\ ", 7)),
                    classRep = as.character(split.vec(split.vec(a[grepl("#",a)], "#", 2),"\\ ",1)), #se deben aplicar los mismos filtros para los nombres de las repeticiones que en el archivo RM
                    K2P = as.numeric(split.vec(a[grepl("Kimura",a)], "= ", 2))
  )
  if(average==TRUE){
    K2P_mean <- tapply(ALN$K2P, ALN[,by], mean)
    K2P_mean
  }else{ALN}
}

