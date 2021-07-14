#' @import stringr
#' @import GenomicRanges
#' @import IRanges
#' @title Makes a reference file for Salmon
#' @description This function creates decoys and a transcriptome that will be used by Salmon. It also creates a reference file to import the estimates after the Salmon run.
#' The user can enter a RepMask file without deleting co-transcribed or overlapping repeats with the RepMask argument, or enter a RepMask file without co-transcribed but overlapping repeats with the RepMask.clean argument, or a file free of co-transcribed or overlapping repeats with the RepMask.ovlp.clean argument. When the file contains co-transcribed repeats, it must indicate rm.cotrans = T and when the file contains overlaps it must indicate overlapping = T.
#' @param overlapping Indicates whether the RepMask file contains overlapping repetitions (TRUE) or not (FALSE). When the RepMask file contains overlapping repetitions, the ovlp.res() function will be used to solve them and the resolution criteria must be indicated (higher score (HS), longer length (LE) or lower Kimura distances (LD))
#' @param rule A numerical vector respectively indicating the minimum percentage of identity, length (in base pairs) of the repeat to be analyzed and the percentage of the length of class/family repeat with respect to the length of the transcript. Example: c(80, 80, 80). Default is c(0,0,0)
#' @param trme transcriptome in fasta format
#' @param RepMask RepeatMasker output file. If rm.cotrans = F it is assumed that the file does not contain cotranscribed repeats. If overlapping = F it is assumed that the file does not contain overlapping.
#' @param rm.cotrnas logical vector indicating whether co-transcribed repeats should be removed
#' @param align .align file
#' @param over.res Indicates the method by which the repetition overlap will be resolved.
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
#' @param by The column by which the repeats will be classified
#' @param annot_by A character vector indicating whether the annotations should be made by "transcripts" or by "fragments". When annot_by = "transcripts", the proportion of each transposon class/family in each transcript is calculated and the transcript is annotated with the class/family with the highest coverage.
#' @param ignore.aln.pos The RepeatMasker alignments file may have discrepancies in the repeats positions with respect to the output file. If you selected over.res = "LD", then you can choose whether to take into account the positions of the alignment file or to take the average per repeats class (default).
#' @param threads Number of cores to use in the processing. By default threads = 1
#' @param bedtools bedtools binary path
#' @export
mk.reference <- function(RepMask,overlapping=F, by=c("namRep","classRep", "class", "supFam", "Fam"), rule=c(0,0,0), trme, threads=1, annot_by="transcripts", bedtools="bedtools", outdir, over.res="HS", ...){

if(overlapping==T){
    RM <- ovlp.res(RepMask=RepMask, threads=threads, outdir=outdir,over.res=over.res,...)
  }else{
    RM <- RepMask
  }
RM <- RM[RM$classRep%!in%c("Unknown", "rRNA", "Satellite", "Simple_repeat","Low_complexity","RNA","scRNA","snRNA","srpRNA", "tRNA","Other"),]
RM$width <- (RM$EPMQuer - RM$SPMQuer)

if(sum(rule)==0 && annot_by=="fragments"){
    message("making references to fragments annotations ...")
    message("rule is 0-0-0 ...")
    BED <- cbind(RM[,c(5,6,7)],RM[,by])
    }else{

message("applying rule ...")
refSeqs <- RM[RM$PersubM<(100-rule[1]) & RM$width>rule[2],]
refSeqs$namSeqL <- paste0(refSeqs$namSeq,".",refSeqs$seqLength)
suppressWarnings(
  a <- as.data.frame(
    unlist(
      apply(
        t(
          tapply(
            refSeqs$width,
            list(refSeqs[,by],
                 refSeqs$namSeqL),
            sum)
        ),
        1,
        function(x){
          return(
            na.omit(x)
          )
        }
      )
    )
  )
)
  sm_df <- data.frame(
    namSeq=split.vec(row.names(a),"\\.",1),
    repNam=split.vec(row.names(a),"\\.",3),
    seqLength=as.numeric(split.vec(row.names(a),"\\.",2)),
    Sum=a[,1]
  )

sm_df$perRepOvTran <- (sm_df$Sum/sm_df$seqLength*100)

    if(annot_by=="transcripts"){
      BED <- sm_df[sm_df$perRepOvTran>rule[3],]
      BED <- BED[order(BED$perRepOvTran, decreasing = T),]
      BED <- BED[!duplicated(BED$namSeq),]
      BED <- data.frame(namSeq=BED$namSeq,start=1,end=as.numeric(BED$seqLength),repNam=BED$repNam)
      message("making references to transcripts annotations ...")
    }else{
      refSeqs <- refSeqs[refSeqs$namSeq%in% sm_df$namSeq[sm_df$perRepOvTran>rule[3]],]
      BED <- data.frame(namSeq=refSeqs$namSeq,start=refSeqs$SPMQuer,end=refSeqs$EPMQuer,repNam=refSeqs[,by])
      message("making references to fragments annotations ...")
    }
}

write.table(BED,"decoys_head.bed", quote = F, row.names = F, col.names=F, sep="\t")

allIDfasta <- system(paste("grep '>'",trme,"|sed 's/>//' | awk '{print $1}'"),intern = T)
decoys <- allIDfasta[allIDfasta%!in%BED$refSeqs.namSeq]
write.table(decoys,"decoys.txt", col.names = F, row.names = F, quote = F)

system(paste0("awk", ' \'BEGIN{while((getline<"decoys.txt")>0)l[">"$1]=1}/^>/{f=l[$1]}f\' ',trme," > decoy.fa"))

message("making trmeSalmon.fasta file")
  system("cat decoys_head.bed | sort -k1,1 -k2,2n > RM_or.bed")
  system(paste(bedtools,"merge -i RM_or.bed -c 4 -o collapse > RM_or_merged.bed"))
  system(paste(bedtools,"getfasta -fi",trme,"-bed RM_or_merged.bed -fo Rep.fa"))
  system("cat Rep.fa decoy.fa > trmeSalmon.fasta")

system(paste("rm decoys_head.bed Rep.fa RM_or_merged.bed RM_or.bed decoy.fa decoys.txt", paste0(trme,".fai")))

Ref.salmon <- data.frame(paste0(BED[,1],":",BED[,2],"-",BED[,3]),BED[,4])
names(Ref.salmon) <- c("seqID","repID")

message("writing files in the output directory...")
system(paste("mv trmeSalmon.fasta",outdir))
write.table(Ref.salmon,paste0(outdir,"/references.csv"), col.names = F, row.names = F, quote = F, sep = ";")
write.table(decoys,paste0(outdir,"/decoys.txt"), col.names = F, row.names = F, quote = F)

message(paste("The reference.csv, decoys.txt and trmeSalmon.fasta files are in", outdir, "directory"))

  Ref.salmon
}

