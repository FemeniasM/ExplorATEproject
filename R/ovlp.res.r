#' @import stringr
#' @import GenomicRanges
#' @import IRanges
#' @import parallel
#' @import doParallel
#' @import foreach
#' @importFrom foreach %dopar%
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
#' @export
ovlp.res <- function(RepMask,anot,gff3,stranded=T, outdir, rm.cotrans=F, align, threads=1, ignore.aln.pos=T, over.res=c("HS","LS","LD"), ...){
#removiendo cotranscriptos
if(rm.cotrans==T){
message("removing co-transcribed repeats ...")

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
RM.clean$seqLength <- as.numeric(RM.clean$EPMQuer+as.numeric(gsub("\\(|\\)","",RM.clean$NbasAfEQuer)))

message("making references for overlaps resolution ...")
  #leyendo el archivo de alineamientos

  if(over.res=="LD"){
    if(class(align)=="data.frame"){
      ALN <- align
      namcol.alin <- c("seqID", "beg", "end","classRep","K2P")
      if(length(ALN)!=length(namcol.alin)){stop("missing columns in .align file")}
      err.0 <- namcol[!(colnames(ALN)==namcol.alin)]
      if(length(err.0)!=0){ stop(paste("The .align file must contain the following column names: ", err ))}

    }else{ALN <- read.alignfile(align)}

    if(class(ALN)!="data.frame"){stop("check .aling file")}
  }

# Primero verificamos cuales son las repeticiones que tienen solapamientos y luego creamos una tabla con estas repeticiones.

cl <- parallel::makeCluster(threads)
doParallel::registerDoParallel(cl)
suppressWarnings(
RM.ovlp <- foreach::foreach (n=unique(RM.clean$namSeq), .combine = rbind ) %dopar% {
    RM.ranges <- GenomicRanges::GRanges(n,
                                      IRanges::IRanges(
                                        as.numeric(RM.clean$SPMQuer[RM.clean$namSeq==n]),
                                        as.numeric(RM.clean$EPMQuer[RM.clean$namSeq==n]),
                                        names = as.character(RM.clean$classRep[RM.clean$namSeq==n])
                                      ),
                                      gsub("C", "-",RM.clean$st[RM.clean$namSeq==n]))

  x <- data.frame( IRanges::subsetByOverlaps(RM.ranges[0],RM.ranges[0]) )
  for(i in 1:length(RM.ranges)){

    df <- data.frame( IRanges::subsetByOverlaps(RM.ranges[i],RM.ranges[-i]) )
    if(nrow(df)>0){
      x <- rbind(x,df)

    }

  }
#ahora creamos una tabla RM a la que le agregamos las criterios de seleccion
    tbl <- RM.clean[1,][-1,]
  if(over.res=="LD"){tbl$alinrep <- vector()}
  tbl$width <- vector()
  for(l in 1:nrow(x)){
    j <- RM.clean[RM.clean$namSeq==x$seqnames[l] & RM.clean$SPMQuer==x$start[l] & RM.clean$EPMQuer==x$end[l], ]
    tbl[l,1:19] <- j
    if(over.res=="LD"){
      if(ignore.aln.pos==T){
        tbl$alinrep[l] <- mean(ALN$K2P[ALN$seqID==as.character(x$seqnames[l]) & ALN$classRep==as.character(j$classRep[l])], na.rm = T)
      }else{
        tbl$alinrep[l] <- ALN[ALN$seqID==as.character(x$seqnames[l]) & ALN$classRep==as.character(j$classRep) & ALN$beg==x$start[l] & ALN$end==x$end[l],]

      }
    }
    tbl$width[l] <- as.numeric(x$width[l])
  }
  if(!(is.na(tbl$namSeq))){
  return(tbl)
  }

}
)
message("resolving overlaps ...")
suppressWarnings(
RM.ovlp.res <- foreach::foreach (n=unique(RM.ovlp$namSeq), .combine = rbind ) %dopar% {
sbtbl <- RM.ovlp[RM.ovlp$namSeq==n,]
  for(i in 1:(nrow(sbtbl)-1)){
    if(is.na(sbtbl$SPMQuer[i]) | is.na(sbtbl$SPMQuer[i+1])){
      next}else{
        #si esta contenido dentro:
        if(sbtbl$SPMQuer[i+1]<=sbtbl$SPMQuer[i] & sbtbl$EPMQuer[i+1]>=sbtbl$EPMQuer[i] | sbtbl$SPMQuer[i]<=sbtbl$SPMQuer[i+1] & sbtbl$EPMQuer[i]>=sbtbl$EPMQuer[i+1]){
          if(sbtbl$SPMQuer[i]<=sbtbl$SPMQuer[i+1] & sbtbl$EPMQuer[i]>=sbtbl$EPMQuer[i+1]){
            if(over.res=="HS"){
              if(sbtbl$scoreSW[i]>=sbtbl$scoreSW[i+1]){
                sbtbl[i+1,] <- NA
              }else{
                sbtbl$SPMQuer[i+1] <- sbtbl$SPMQuer[i]
                sbtbl$EPMQuer[i+1] <- sbtbl$EPMQuer[i]
                sbtbl[i,] <- NA
              }}

            if(over.res=="LE"){
              if(sbtbl$width[i]>=sbtbl$width[i+1]){
                sbtbl[i+1,] <- NA
              }else{
                sbtbl$SPMQuer[i+1] <- sbtbl$SPMQuer[i]
                sbtbl$EPMQuer[i+1] <- sbtbl$EPMQuer[i]
                sbtbl[i,] <- NA
              }}
            if(over.res=="LD"){
              if(sbtbl$alinrep[i]<=sbtbl$alinrep[i+1]){
                sbtbl[i+1,] <- NA
              }else{
                sbtbl$SPMQuer[i+1] <- sbtbl$SPMQuer[i]
                sbtbl$EPMQuer[i+1] <- sbtbl$EPMQuer[i]
                sbtbl[i,] <- NA
              }}
          }else{

            if(sbtbl$SPMQuer[i+1] <= sbtbl$SPMQuer[i] & sbtbl$EPMQuer[i+1] >= sbtbl$EPMQuer[i]){
              if(over.res=="HS"){
                if(sbtbl$scoreSW[i]<=sbtbl$scoreSW[i+1]){
                  sbtbl[i,] <- NA
                }else{
                  sbtbl$SPMQuer[i] <- sbtbl$SPMQuer[i+1]
                  sbtbl$EPMQuer[i] <- sbtbl$EPMQuer[i+1]
                  sbtbl[i+1,] <- NA
                }}
              if(over.res=="LE"){
                if(sbtbl$width[i]<=sbtbl$width[i+1]){
                  sbtbl[i,] <- NA
                }else{
                  sbtbl$SPMQuer[i] <- sbtbl$SPMQuer[i+1]
                  sbtbl$EPMQuer[i] <- sbtbl$EPMQuer[i+1]
                  sbtbl[i+1,] <- NA
                }}
              if(over.res=="LD"){
                if(sbtbl$alinrep[i]>sbtbl$alinrep[i+1]){
                  sbtbl[i,] <- NA
                }else{
                  sbtbl$SPMQuer[i] <- sbtbl$SPMQuer[i+1]
                  sbtbl$EPMQuer[i] <- sbtbl$EPMQuer[i+1]
                  sbtbl[i+1,] <- NA
                }}
            }
          }
          #Si está solapado parcialmente:
        }else{
          if(sbtbl$EPMQuer[i]>=sbtbl$SPMQuer[i+1]){
            if(sbtbl$classRep[i]==sbtbl$classRep[i+1]){
              sbtbl$EPMQuer[i] <-sbtbl$EPMQuer[i+1]
              if(over.res=="HS"){
                sbtbl$scoreSW[i] <- max(c(sbtbl$scoreSW[i],sbtbl$scoreSW[i+1]))
              }
              if(over.res=="LE"){
                sbtbl$width[i] <- max(c(sbtbl$width[i],sbtbl$width[i+1]))
              }
              if(over.res=="LD"){
                sbtbl$alinrep[i] <- min(c(sbtbl$alinrep[i],sbtbl$alinrep[i+1]))
              }
              sbtbl[i+1,] <- NA
            }else{
              if(over.res=="HS"){
                if(sbtbl$scoreSW[i]>=sbtbl$scoreSW[i+1]){
                  sbtbl$SPMQuer[i+1] <- (sbtbl$EPMQuer[i])+1
                }else{
                  sbtbl$EPMQuer[i] <- (sbtbl$SPMQuer[i+1])-1
                }
              }

              if(over.res=="LE"){
                if(sbtbl$width[i]>=sbtbl$width[i+1]){
                  sbtbl$SPMQuer[i+1] <- (sbtbl$EPMQuer[i])+1
                }else{
                  sbtbl$EPMQuer[i] <- (sbtbl$SPMQuer[i+1])-1
                }
              }

              if(over.res=="LD"){
                if(sbtbl$alinrep[i]<=sbtbl$alinrep[i+1]){
                  sbtbl$SPMQuer[i+1] <- (sbtbl$EPMQuer[i])+1
                }else{
                  sbtbl$EPMQuer[i] <- (sbtbl$SPMQuer[i+1])-1
                }
              }
            }
          }
        }
      }
  }
  return(sbtbl[!(is.na(sbtbl$namSeq)),])
}
)
parallel::stopCluster(cl)
foreach::registerDoSEQ()

for(n in 1:nrow(RM.ovlp.res)){
  RM.ovlp.res$seqLength[n] <- unique(RM.clean$seqLength[RM.clean$namSeq==RM.ovlp.res$namSeq[n]])
}

#combinamos los resueltos con los no solapados:
RM.clean.no_ovlp <- RM.clean[RM.clean$nu%!in%RM.ovlp$nu,]#no solapados
RM.clean.all <- rbind(RM.clean.no_ovlp,RM.ovlp.res[,colnames(RM.clean.no_ovlp)])

RM.clean.all <- RM.clean.all[order(RM.clean.all$SPMQuer),]
RM.clean.all <- RM.clean.all[order(RM.clean.all$namSeq),]
message("overlaps resolution finished ...")
write.table(RM.clean.all, paste0(outdir, "/RM.withoutOvlp.out"), quote = FALSE, row.names = FALSE,sep = "\t")
RM.clean.all
}
