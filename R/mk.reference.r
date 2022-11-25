#' @import seqinr
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
#' @param trpt.length A data.frame with two columns: the first column must contain the name of the transcripts, and the second the length corresponding to each transcript. The default is trpt.length=NULL, and the lengths for each transcript are taken from the RepeatMasker file.
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
#' @export
mk.reference <- function(RepMask,overlapping=T, by="classRep", trme, threads=1, annot_by="transcripts", rule=c(80,80,80), best.by="total_repeat_length", outdir, over.res="HS",trpt.length=NULL, ...){

  if(overlapping==T){
    RM <- ovlp.res(RepMask=RepMask, threads=threads, outdir=outdir,over.res=over.res,trpt.length=trpt.length,by=by,...)

    }else{

      RM <- RepMask

      }

RM <- RM[RM$classRep%!in%c("Unknown", "rRNA", "Satellite", "Simple_repeat","Low_complexity","RNA","scRNA","snRNA","srpRNA", "tRNA","Other"),]

message(paste("[",Sys.time(),"][INFO]","applying the Wicker-like rule ..."))
WRout <- RM[RM$per_identity>rule[1] & (RM$perc.rep.trpt*100)>rule[2] & RM$total_transcript_length>rule[3],]

WRout_uniq_seqNam <- unique(sort(WRout$namSeq))

message(paste("[",Sys.time(),"][INFO]","reading fasta files ..."))
fasta <- seqinr::read.fasta(trme)

message(paste("[",Sys.time(),"][INFO]","defining decoys sequences ..."))
decoys_fasta <- fasta[names(fasta)%!in%WRout_uniq_seqNam]
decoys <- names(decoys_fasta)

  if(annot_by=="transcripts"){
    #anotando los transcriptos: como pueden haber mas de una repeticion no solapada en cada transcripto, se usa best.by para anotar solo una rep por secuencia
    #por defecto "total_repeat_length" selecciona la repeticion que suma mas pb en todo el transcripto

    message(paste("[",Sys.time(),"][INFO]","annotating by transcripts ..."))

    WRout <- WRout[order(as.numeric(WRout[,best.by]), decreasing = T),]
    WRout  <- WRout[!duplicated(WRout$namSeq),]

    target_fasta <- fasta[names(fasta)%in%WRout$namSeq]

    message(paste("[",Sys.time(),"][INFO]","Generating salmon references with by =", paste0("'",by,"'")))

    Ref.salmon <- WRout[,c("namSeq",by)]


  } else if(annot_by=="fragments"){

    if(by=="namRep"){
      WRout$SeqID <- paste0(WRout$namSeq,":",WRout$start,":",WRout$end,":",WRout$namRep,"/",WRout$classRep)
    }
    if(by=="classRep"){
      WRout$SeqID <- paste0(WRout$namSeq,":",WRout$start,":",WRout$end,":",WRout$classRep)
    }

      BED <- WRout[,c("namSeq","start","end","SeqID")]

      target_fasta <- list()

      for(i in 1:nrow(BED)){
        target_fasta[[i]] <- seqinr::as.SeqFastadna(fasta[[BED$namSeq[i]]][BED$start[i]:BED$end[i]],
                                                         name = BED$SeqID[i],
                                                         Annot = BED$SeqID[i])

      }
      message(paste("[",Sys.time(),"][INFO]","Generating salmon references with by =", paste0("'",by,"'")))

      Ref.salmon <- WRout[,c("SeqID",by)]

  }

all.fasta <- c(target_fasta,decoys_fasta)

  if(!(dir.exists(outdir))){

    dir.create(outdir)
  }

message(paste("[",Sys.time(),"][INFO]","writing references files in the output directory..."))

seqinr::write.fasta(sequences = seqinr::getSequence(all.fasta), names = seqinr::getName(all.fasta), file.out = paste0(outdir,"/trmeSalmon.fasta"))

write.table(Ref.salmon,paste0(outdir,"/references.csv"), col.names = F, row.names = F, quote = F, sep = ";")

write.table(decoys,paste0(outdir,"/decoys.txt"), col.names = F, row.names = F, quote = F)


message(paste("[",Sys.time(),"][INFO]","The reference.csv, decoys.txt and trmeSalmon.fasta files are in", outdir, "directory"))

Ref.salmon
}

