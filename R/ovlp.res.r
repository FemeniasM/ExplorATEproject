#' @import data.table
#' @import IRanges
#' @importFrom data.table :=
#' @title Overlappings resolution
#' @description This function resolves overlapping repeats assigned to the same transcript and returns a data frame of repeats with no overlaps. The user can define the criteria to solve the overlaps, either by higher score (HS), longer length (LE) or lower Kimura's distances (LD).
#' @param RepMask RepeatMasker output file. If rm.cotrans = F, them you must enter a RepeatMasker output file without co-transcripted repeats.
#' @param rm.cotrnas logical vector indicating whether co-transcribed repeats should be removed
#' @param align .align file
#' @param ignore.aln.pos The RepeatMasker alignments file may have discrepancies in the repeats positions with respect to the output file. If you selected over.res = "LD", then you can choose whether to take into account the positions of the alignment file or to take the average per repeats class (default).
#' @param over.res Indicates the method by which the repetition overlap will be resolved ("HS" by default).
#' HS: higher score, bases are assigned to the element with the highest score
#' LS: longer element, bases are assigned to the longest element
#' LD: lower divergence, bases are assigned to the element with the least divergence.
#' in all cases both elements have the same characteristics, the bases are assigned to the first element.
#' @param anot annotation file in outfmt6 format. It is necessary when the option rm.cotrans = T
#' @param gff3 gff3 file. It is necessary when the option rm.cotrans = T
#' @param stranded logical vector indicating if the library is strand specific
#' @param cleanTEsProt logical vector indicating whether the search for TEs-related proteins should be carried out (e.g.
#' transposases, integrases, env, reverse transcriptase, etc.). We recommend that users use a curated annotations file,
#' in which these genes have been excluded; therefore the default option is F. When T is selected, a search is performed
#' against a database obtained from UniProt, so we recommend that the annotations file have this format for the subject
#' sequence id (e.g. "CO1A2_MOUSE"/"sp|Q01149|CO1A2_MOUSE"/"tr|H9GLU4|H9GLU4_ANOCA")
#' @param featureSum Returns statistics related to the characteristics of the transcripts. Requires a gff3 file. If TRUE, returns a list of the
#' @param outdir Output directory
#' @param threads Number of cores to use in the processing. By default threads = 1
#' @param trpt.length A data.frame with two columns: the first column must contain the name of the transcripts, and the second the length corresponding to each transcript. The default is trpt.length=NULL, and the lengths for each transcript are taken from the RepeatMasker file.
#' @export
ovlp.res <- function(RepMask,anot,gff3,stranded=T, outdir, rm.cotrans=F, trpt.length=NULL, align, threads=1, ignore.aln.pos=T, over.res=c("HS","LS","LD"),
                       by="classRep", ...){
  #removiendo cotranscriptos
  if(rm.cotrans==T){
    message(paste("[",Sys.time(),"][INFO]","removing co-transcribed repeats ..."))

    RM.clean <- rm.cotransRep(RepMask=RepMask,anot=anot,gff3=gff3,stranded=stranded,outdir=outdir, ...)

  }else{ #leyendo repeatmasker
    namcol <- c("scoreSW", "PersubM","PerBasDel","PerBasIns","namSeq","SPMQuer","EPMQuer","NbasAfEQuer", "st",
                "namRep","classRep","NbasCompRep","SPMRepdb","EPMRepdb","nu","HSM")

    if(class(RepMask)=="data.frame"){

      if(length(RepMask)!=length(namcol)){stop("missing columns in RepeatMasker file")}
      err <- namcol[!(colnames(RepMask)==namcol)]
      if(length(err)!=0){ stop(paste("The RepeatMasker file must contain the following column names: ", err ))}

      RM.clean <- RepMask
    }else{
      RM.clean <- read.RepMask(RepMask)
    }


  }

message(paste("[",Sys.time(),"][INFO]","making references for overlaps resolution ..."))
#### Ahora la divergencia de la secuencia se toma desde el archivo RepMask, con el porcentaje de sustituciones
  # #leyendo el archivo de alineamientos
  #
  # if(over.res=="LD"){
  #   if(class(align)=="data.frame"){
  #     ALN <- align
  #     namcol.alin <- c("seqID", "beg", "end","classRep","K2P")
  #     if(length(ALN)!=length(namcol.alin)){stop("missing columns in .align file")}
  #     err.0 <- namcol[!(colnames(ALN)==namcol.alin)]
  #     if(length(err.0)!=0){ stop(paste("The .align file must contain the following column names: ", err ))}
  #
  #   }else{ALN <- read.alignfile(align)}
  #
  #   if(class(ALN)!="data.frame"){stop("check .aling file")}
  # }

  # Primero verificamos cuales son las repeticiones que tienen solapamientos y luego creamos una tabla con estas repeticiones.

data.table::setDTthreads(threads = threads)

if(by=="namRep"){

  RM.clean$byRep <- paste0(RM.clean$namRep,"::",RM.clean$classRep)

  }else if(by=="classRep"){

    RM.clean$byRep <- RM.clean$classRep

  } else{

    message(paste("[",Sys.time(),"][INFO]","error in 'by' argument: repeats must be defined by 'classRep' or 'namRep'"))
}


  RM1=RM2=data.table::data.table(RM.clean[,c( "namSeq", "SPMQuer", "EPMQuer", "byRep")])

  data.table::setkey(RM2, namSeq, SPMQuer, EPMQuer)

  RM_ovlp <- data.table::foverlaps(RM1, RM1, type="any", nomatch=NULL)


  #incorporando socres para cada repeat
  scores <- RM.clean[, c("PersubM","scoreSW","namSeq", "SPMQuer", "EPMQuer", "byRep")]
  scores$score.div <- 100-scores$PersubM
  scores$score.length <- scores$EPMQuer-scores$SPMQuer

  RM_ovlp_scores <- merge(RM_ovlp, scores, by=c("namSeq","SPMQuer", "EPMQuer","byRep"))

  #resolviendo solapamientos

message(paste("[",Sys.time(),"][INFO]","resolving overlaps"))

RM_ovlp_res <- RM_ovlp_scores[
  ,IDx := cumsum(data.table::fcoalesce(+(SPMQuer > (data.table::shift(cummax(EPMQuer), type = 'lag') + 1L)), 0L)), by = namSeq][
    , .(SPMQuer = min(SPMQuer), EPMQuer = max(EPMQuer),
        byRep = byRep[which.max(ifelse(over.res=="LS", score.length, ifelse(over.res=="LD",score.div,scoreSW)))],
        ovlp.res.score = max(ifelse(over.res=="LS", score.length, ifelse(over.res=="LD",score.div,scoreSW))),
    mean.score.div=mean(score.div),
    mean.scoreSW=mean(scoreSW),
    rep.length=(max(EPMQuer)-min(SPMQuer))),
    by = .(namSeq, IDx)][, IDx := NULL]
RM_ovlp_res <- data.table::data.table(RM_ovlp_res)

if(is.null(trpt.length)){
  message(paste("[",Sys.time(),"][INFO]","Taking sequence lengths from the RepeatMasker file"))
  RM.clean$trpt.length <- as.numeric(RM.clean$EPMQuer+as.numeric(gsub("\\(|\\)","",RM.clean$NbasAfEQuer)))
  trpt.length <- RM.clean[order(RM.clean$namSeq),c("namSeq","trpt.length")]
  trpt.length <- trpt.length[!(duplicated(trpt.length$namSeq)),]
  }else{
  message(paste("[",Sys.time(),"][INFO]","Taking sequence lengths from",trpt.length, "file"))
  colnames(trpt.length) <- c("namSeq","trpt.length")
}

#importando longitudes de secuencias
RM_ovlp_res <- merge(RM_ovlp_res,trpt.length, by="namSeq")

#incorporando las variables de anotacion
message(paste("[",Sys.time(),"][INFO]","Creating annotation table with non-overlapped repeats"))
RM_ovlp_res <- RM_ovlp_res[which(RM_ovlp_res$rep.length > 0)][,
            .(start=SPMQuer,end=EPMQuer,rep.length,
              total_repeat_length=sum(rep.length),
              total_transcript_length=trpt.length,
              ovlp.res.score,per_identity=mean.score.div,mean.scoreSW,
              perc.rep.trpt=(sum(rep.length)/trpt.length)
            ),by=c("namSeq","byRep")]

  message(paste("[",Sys.time(),"][INFO]","Overlaps resolution finished ..."))
  #escribiendo la tabla
  df <- as.data.frame(RM_ovlp_res)

  if(by=="namRep"){

  df$namRep <- split.vec(df$byRep,"::")
  df$classRep <- split.vec(df$byRep,"::",2)
  df$byRep <- NULL
  df <- df[,c(1:3,11:12,4:10)]
  }

  if(by=="classRep"){

    df$classRep <- df$byRep
    df$byRep <- NULL
      df <- df[,c(1:3,11,4:10)]
  }


  if(!(dir.exists(outdir))){

    dir.create(outdir)
  }

message(paste("[",Sys.time(),"][INFO]","Overlaps resolution finished: writing files in the output directory..."))
write.table(df, file=file.path(paste0(outdir,"/RM.withoutOvlp.out")), quote = FALSE, row.names = FALSE,sep = "\t")

  df
}
