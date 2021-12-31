#' @importFrom dplyr %>% group_by summarize
#' @title Makes a reference file for Salmon
#' @description This function creates decoys and a transcriptome that will be used by Salmon. It also creates a reference file to import the estimates after the Salmon run.
#' The user can enter a RepMask file without deleting co-transcribed or overlapping repeats with the RepMask argument, or enter a RepMask file without co-transcribed but overlapping repeats with the RepMask.clean argument, or a file free of co-transcribed or overlapping repeats with the RepMask.ovlp.clean argument. When the file contains co-transcribed repeats, it must indicate rm.cotrans = T and when the file contains overlaps it must indicate overlapping = T.
#' @param overlapping Indicates whether the RepMask file contains overlapping repetitions (TRUE) or not (FALSE). When the RepMask file contains overlapping repetitions, the ovlp.res() function will be used to solve them and the resolution criteria must be indicated (higher score (HS), longer length (LE) or lower Kimura distances (LD))
#' @param rule A numerical vector respectively indicating the minimum percentage of identity, the percentage of the length of class/family repeat with respect to the length of the transcript, and the length (in base pairs) of the repeat to be analyzed.
#' #The position of the numbers indicates respectively:
#Percent identity. Is taken as 100 minus the defined percentage of divergence. 80% percent identity support up to 20% divergence
#Ratio of the repeats to the transcript. (number of base pairs of the repeat / total number of bases of the transcript * 100)
#Minimum transcipt length
#' Example: c(80, 60, 100) indicates that those repeats with 80% identity or more in at least 60% of the transcript, and are at least 100 bp in length will be annotated as target TEs.  Default is c(80,80,80)
#' @param by The column by which the repeats will be classified: 'classRep' (default) or 'namRep'.
#' @param best.by Defines if only the best match of each transcript/sequence id should be returned (by default best.by = NULL which shows all matches for the sequence). The user can choose whether to be based on the longest repeat length ('total_repeat_length') or the highest percent identity ('per_divergence'). The mk.reference() function uses the best.by argument when references are annotated by transcripts (annot_by = 'transcripts')
#' A logical vector indicating its only the longest repeats for each transcript is reported. By default best = TRUE
#' @param custom.lengths A a numerical vector indicating the lengths for each transcript. By default lengths are calculated from the RepeatMasker file.
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
#' @param annot_by A character vector indicating whether the annotations should be made by "transcripts" or by "fragments". When annot_by = "transcripts", the proportion of each transposon class/family in each transcript is calculated and the transcript is annotated with the class/family with the highest coverage.
#' @param ignore.aln.pos The RepeatMasker alignments file may have discrepancies in the repeats positions with respect to the output file. If you selected over.res = "LD", then you can choose whether to take into account the positions of the alignment file or to take the average per repeats class (default).
#' @param threads Number of cores to use in the processing. By default threads = 1
#' @param bedtools bedtools binary path
#' @export
mk.reference <- function(RepMask,overlapping=F, by="classRep", trme, threads=1, annot_by="transcripts", rule=c(80,80,80), best.by="total_repeat_length", bedtools="bedtools", outdir, over.res="HS",custom.lengths=NULL, ...){

if(overlapping==T){
    RM <- ovlp.res(RepMask=RepMask, threads=threads, outdir=outdir,over.res=over.res,...)
  }else{
    RM <- RepMask
  }
RM <- RM[RM$classRep%!in%c("Unknown", "rRNA", "Satellite", "Simple_repeat","Low_complexity","RNA","scRNA","snRNA","srpRNA", "tRNA","Other"),]

if(annot_by=="transcripts"){
message("applying the Wicker-like rule ...")
WRout <- Wickerlike.rule(RM,by=by, rule = rule, best.by = best.by)
write.table(unique(WRout$SeqID),"targetTEs.txt", quote = F, row.names = F, col.names=F)

message("making references from transcripts annotations ...")
system(paste0("awk", ' \'BEGIN{while((getline<"targetTEs.txt")>0)l[">"$1]=1}/^>/{f=l[$1]}f\' ',trme," > targetTEs.fa"))

Ref.salmon <- data.frame(seqID=WRout$SeqID,repID=WRout[,4])
} else {

if(annot_by=="fragments"){
message("applying the Wicker-like rule ...")
WRout <- Wickerlike.rule(RM,by=by, rule = rule)
BED <- WRout[,1:4]

message("making references from fragments annotations ...")
write.table(BED,"WRout.bed", quote = F, row.names = F, col.names=F, sep="\t")

system("cat WRout.bed | sort -k1,1 -k2,2n > RM_sort.bed")
system(paste(bedtools,"merge -i RM_sort.bed -c 4 -o collapse > RM_sort_merged.bed"))
system(paste(bedtools,"getfasta -fi",trme,"-bed RM_sort_merged.bed -fo targetTEs.fa"))

Ref.salmon <- data.frame(seqID=paste0(BED[,1],":",BED[,2],"-",BED[,3]),repID=BED[,4])
} else {
  message("error in annot_by argument")
  }
}

message("taking decoy sequences")
allIDfasta <- system(paste("grep '>'",trme,"|sed 's/>//' | awk '{print $1}'"),intern = T)
decoys <- allIDfasta[allIDfasta%!in%unique(WRout$SeqID)]
write.table(decoys,"decoys.temp", col.names = F, row.names = F, quote = F)

system(paste0("awk", ' \'BEGIN{while((getline<"decoys.temp")>0)l[">"$1]=1}/^>/{f=l[$1]}f\' ',trme," > decoy.fa"))

message("making trmeSalmon.fasta file")
system("cat targetTEs.fa decoy.fa > trmeSalmon.fasta")

if(annot_by=="transcripts"){
system("rm targetTEs.fa decoy.fa decoys.temp")
}

if(annot_by=="fragments"){
system(paste("rm WRout.bed targetTEs.fa RM_sort_merged.bed RM_sort.bed decoy.fa decoys.temp", paste0(trme,".fai")))
}

message("writing files in the output directory...")
system(paste("mv trmeSalmon.fasta",outdir))
write.table(Ref.salmon,paste0(outdir,"/references.csv"), col.names = F, row.names = F, quote = F, sep = ";")
write.table(decoys,paste0(outdir,"/decoys.txt"), col.names = F, row.names = F, quote = F)

message(paste("The reference.csv, decoys.txt and trmeSalmon.fasta files are in", outdir, "directory"))

  Ref.salmon

}
