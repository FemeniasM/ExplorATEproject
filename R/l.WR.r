#' @importFrom dplyr %>% group_by summarize
#' @title Applies a rule-based sequence selection criterion similar to Wicker's rule on processed transcripts/sequences.
#' @description Applies a rule-based sequence selection criterion similar to Wicker's rule on processed transcripts / sequences. For each sequence identification, take the percent identity (100 minus the percent divergence), the ratio of the TE identified in the sequence respect to the total sequence length, and the minimum transcript length. The function returns a table with the sequence IDs that meet the criteria, the repeats name and its positions within the sequence, and other auxiliary variables (percentage of identity, number of total bases for each repetition and length of the transcript). The user can define whether to summarize all the matches found in each sequence or only the best match.
#' @param RepMask RepeatMasker output file.
#' @param rule A numerical vector respectively indicating the minimum percentage of identity, the percentage of the length of class/family repeat with respect to the length of the transcript, and the length (in base pairs) of the repeat to be analyzed.
#' #The position of the numbers indicates respectively:
#Percent identity. Is taken as 100 minus the defined percentage of divergence. 80% percent identity support up to 20% divergence
#Ratio of the repeats to the transcript. (number of base pairs of the repeat / total number of bases of the transcript * 100)
#Minimum transcipt length
#' Example: c(80, 60, 100) indicates that those repeats with 80% identity or more in at least 60% of the transcript, and are at least 100 bp in length will be annotated as target TEs.  Default is c(80,80,80)
#' @param by The column by which the repeats will be classified: 'classRep' (default) or 'namRep'.
#' @param best.by Defines if only the best match of each transcript/sequence id should be returned (by default best.by = NULL which shows all matches for the sequence). The user can choose whether to be based on the longest repeat length ('total_repeat_length') or the highest percent identity ('per_identity').
#' A logical vector indicating its only the longest repeats for each transcript is reported. By default best = TRUE
#' @param custom.lengths A a numerical vector indicating the lengths for each transcript. By default lengths are calculated from the RepeatMasker file.
#' @export
Wickerlike.rule <- function(RepMask, by="classRep", rule=c(80,80,80), best.by=NULL, custom.lengths=NULL){

if(is.data.frame(RepMask)){
  RM <- RepMask
} else {
  if(is.character(RepMask)){
  RM <- read.RepMask(RepMask)
  } else { message(paste("error loading file"),RepMask) }
}

if (is.null(custom.lengths)){
  RM$total_transcript_length <- as.numeric(RM$EPMQuer)+as.numeric(gsub("\\(|\\)","",RM$NbasAfEQuer))
} else {
    RM$total_transcript_length <- custom.lengths
  }

  colnames(RM)[colnames(RM)==by] <- "REP"
repBytr <- RM %>%
  dplyr::group_by(namSeq, REP) %>%
  dplyr::summarize(total_repeat_length=sum(EPMQuer- SPMQuer), SeqID=namSeq,end=EPMQuer, start=SPMQuer, REP, per_identity=(100-PersubM),total_transcript_length)

repBytr <- as.data.frame(repBytr)
colnames(repBytr)[colnames(repBytr)=="REP"] <- by
keep <- apply(repBytr, 1, function(x) length(unique(x[!is.na(x)])) != 1)
repBytr <- repBytr[keep, ]

rule.repBytr <- repBytr[repBytr$per_identity>rule[1] & (repBytr$total_repeat_length/repBytr$total_transcript_length*100)>rule[2] & repBytr$total_transcript_length>rule[3],]

rule.repBytr <- rule.repBytr[,c("SeqID","start","end",by,"per_identity","total_repeat_length", "total_transcript_length")]
out.file <- rule.repBytr

if (is.null(best.by)==FALSE) {
rule.repBytr.sort <- rule.repBytr[order(as.numeric(rule.repBytr[,best.by]), decreasing = T),]
out.file <- rule.repBytr.sort[!duplicated(rule.repBytr.sort$SeqID),]
}

out.file

}
