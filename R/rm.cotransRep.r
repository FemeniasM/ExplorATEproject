#' @import stringr
#' @import GenomicRanges
#' @import IRanges
#' @title Excludes co-transcripted repeats
#' @description excludes repeats that are co-transcripted with gene coding sequences and return a clean RepeatMasker file for downstream analysis. If featureSum is TRUE, writes a table with the repeats that overlap either the 3'UTR, 5'UTR, or CDS regions.
#' @param RepMask RepeatMasker output file
#' @param anot annotation file in outfmt6 format
#' @param gff3 gff3 file
#' @param stranded logical vector indicating if the library is strand specific
#' @param cleanTEsProt logical vector indicating whether the search for TEs-related proteins should be carried out (e.g. transposases, integrases, env, reverse transcriptase, etc.). We recommend that users work with a curated annotations file, in which these genes have been excluded; therefore the default option is F. When T is selected, a search is performed against a database obtained from UniProt, so we recommend that the annotations file have this format for the subject sequence id (e.g. "CO1A2_MOUSE" or "sp|Q01149|CO1A2_MOUSE" or "tr|H9GLU4|H9GLU4_ANOCA")
#' @param featureSum The function returns a summary of the protein coding transcripts that contain repeats. Three files are created in the output directory:
#' features.summary.csv file, with the transcripts and their characteristics,
#' features.summary.pdf file, with a barplot graph with the number of repetitions of each TE family found in the 5'-UTR, 3 '-UTR and CDS regions
#' and RepeatMasker RM.clean.out file.
#' @param outdir Output directory
#' @export
rm.cotransRep <- function(RepMask,anot,gff3,stranded=T, cleanTEsProt=F, featureSum=F,outdir){
  if(dir.exists(outdir)==F){
      dir.create(outdir)
    }

    if(class(RepMask)=="data.frame"){
    RM <- RepMask
    namcol <- c("scoreSW", "PersubM","PerBasDel","PerBasIns","namSeq","SPMQuer","EPMQuer","NbasAfEQuer", "st",
                "namRep","classRep","NbasCompRep","SPMRepdb","EPMRepdb","nu","HSM")

    if(length(RM)!=length(namcol)){stop(paste("[",Sys.time(),"][ERROR]","missing columns"))}
    err <- namcol[!(colnames(RM)==namcol)]
    if(length(err)!=0){ stop(paste("[",Sys.time(),"][ERROR]","The RepeatMasker file must contain the following column names: ", err ))}

  }else{RM <- read.RepMask(RepMask)}
  message("Please check if the names of the subclass, superfamilies and/or families are correct:\n")
  message(message=for(i in sort(unique(RM$classRep))){cat(i, "\n")})

  {continue <- readline(prompt = "please press 'y' if the names are correct:")
    if(continue=="y"){

      RM <- RM[order(RM$namSeq),]
      RM <- RM[!(RM$scoreSW==""),]
      RM$j <- paste0(RM$namSeq,RM$SPMQuer,RM$EPMQuer,RM$classRep)
      RM <- RM[!duplicated(RM$j),]
      RM$j <-NULL

      ANOT <- read.outfmt6(anot)
      ANOT <- ANOT[order(ANOT$evalue,decreasing = F),]
      ANOT <- ANOT[!duplicated(ANOT$seqid),]
      if (cleanTEsProt==T) {
        lista <- stringr::str_split(ANOT$sseqid,"\\|")
        for( i in 1:length(lista)){
          ANOT$rename[i] <-  ifelse(lista[[i]][1]=="sp",lista[[i]][3],ifelse(lista[[i]][1]=="tr",lista[[i]][3],lista[[i]][1]))
        }
        ANOT <- ANOT[ANOT$rename%!in%ExplorATE::idTEsProt,]
      }

      ANOT$originalid <- ANOT$seqid

      CDS.withRep <- ANOT[as.character(ANOT$originalid)%in%as.character(RM$namSeq), c("seqid","sseqid", "originalid")]
      message(paste("[",Sys.time(),"][INFO]",length(unique(CDS.withRep$seqid)), " coding genes transcripts contain RTEs repeats: ", round(length(unique(CDS.withRep$seqid))/(length(unique(ANOT$seqid)))*100, 2), "% coding genes transcripts"))
      write.csv(CDS.withRep, paste0(outdir, "/cotranscriptRep.csv"), quote = FALSE, row.names = FALSE)


      if (featureSum==T) {

        GFF <- read.gff3(gff3)
        if (stranded==T) {
          GFF <- GFF[GFF$st=="+",]
        }
        GFF.RTEs <- GFF[GFF$seqid%in%unique(CDS.withRep$originalid),]
        GFF.RTEs <- GFF.RTEs[GFF.RTEs$feature=="CDS" | GFF.RTEs$feature=="three_prime_UTR" | GFF.RTEs$feature=="five_prime_UTR",]
        RM.RTEs <- RM[RM$namSeq%in%GFF.RTEs$seqid,]

        GFF.RTEs.CDS <- GFF.RTEs[GFF.RTEs$feature=="CDS",]
        GFF.RTEs.3p <- GFF.RTEs[GFF.RTEs$feature=="three_prime_UTR",]
        GFF.RTEs.5p <- GFF.RTEs[GFF.RTEs$feature=="five_prime_UTR",]


        GFF.RTEs.CDS.gr <- GenomicRanges::GRanges(GFF.RTEs.CDS$seqid, IRanges::IRanges(as.numeric(GFF.RTEs.CDS$strt), as.numeric(GFF.RTEs.CDS$end),
                                                                                       names = as.character(GFF.RTEs.CDS$feature)), GFF.RTEs.CDS$st)

        GFF.RTEs.3p.gr <- GenomicRanges::GRanges(GFF.RTEs.3p$seqid, IRanges::IRanges(as.numeric(GFF.RTEs.3p$strt), as.numeric(GFF.RTEs.3p$end),
                                                                                     names = as.character(GFF.RTEs.3p$feature)), GFF.RTEs.3p$st)

        GFF.RTEs.5p.gr <- GenomicRanges::GRanges(GFF.RTEs.5p$seqid, IRanges::IRanges(as.numeric(GFF.RTEs.5p$strt), as.numeric(GFF.RTEs.5p$end),
                                                                                     names = as.character(GFF.RTEs.5p$feature)), GFF.RTEs.5p$st)

        RM.RTEs.gr <- GenomicRanges::GRanges(RM.RTEs$namSeq, IRanges::IRanges(as.numeric(RM.RTEs$SPMQuer), as.numeric(RM.RTEs$EPMQuer),
                                                                              names = as.character(RM.RTEs$classRep)), gsub("C", "-",RM.RTEs$st))
        suppressWarnings({
          cds.a <- IRanges::subsetByOverlaps(RM.RTEs.gr,GFF.RTEs.CDS.gr)
          cds <- cbind(data.frame(cds.a)[,c(1,5)],
                       IRanges::ranges(cds.a), rep("CDS", nrow(data.frame(cds.a))))
          names(cds) <- c("seqid", "st.Rep", "strtRep", "endRep", "widRep", "classRep", "Feature")

          utr3.a <- IRanges::subsetByOverlaps(RM.RTEs.gr,GFF.RTEs.3p.gr)
          utr3p <- cbind(data.frame(utr3.a)[,c(1,5)],
                         IRanges::ranges(utr3.a),rep("3'-UTR", nrow(data.frame(utr3.a))))
          names(utr3p) <- c("seqid", "st.Rep", "strtRep", "endRep", "widRep", "classRep", "Feature")

          utr5.a <- IRanges::subsetByOverlaps(RM.RTEs.gr,GFF.RTEs.5p.gr)
          utr5p <- cbind(data.frame(utr5.a)[,c(1,5)],
                         IRanges::ranges(utr5.a),rep("5'-UTR", nrow(data.frame(utr5.a))))
          names(utr5p) <- c("seqid", "st.Rep", "strtRep", "endRep", "widRep", "classRep", "Feature")
        })
        features.summary <- rbind(utr5p, cds, utr3p)
        write.csv(features.summary, paste0(outdir, "/features.summary.csv"), quote = FALSE, row.names = FALSE)
        tbl <- table(features.summary[,6:7])
        {   pdf(paste0(outdir, "/features.summary.pdf"),4,4)
          barplot(t(tbl),las=2,beside = T, legend.text = colnames(tbl))
          dev.off()
        }
        message(paste("[",Sys.time(),"][INFO]","Features analysis finished, the results are found in the output directory"))
      }

      RM.clean <- RM[RM$namSeq%!in%CDS.withRep$originalid,]
      RM.clean$class <- split.vec(RM.clean$classRep, "/",1)
      RM.clean$supFam <- split.vec(RM.clean$classRep, "/",2)
      RM.clean$Fam <- split.vec(RM.clean$classRep, "/",3)

      if(!(dir.exists(outdir))){

      dir.create(outdir)
      }


      write.table(RM.clean, paste0(outdir, "/RM.cotransRepClean.out"), quote = FALSE, row.names = FALSE,sep = "\t")

      RM.clean
    }}
}
