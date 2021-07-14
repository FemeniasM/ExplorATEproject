#' @title Makes a reference file for Salmon in model organisms
#' @description This function creates decoys and a transcriptome that will be used by Salmon when the reference genome is available.
#' @param bedtools bedtools binary path
#' @param GFFgen gtf file to reference genome
#' @param genome reference genome in fasta format
#' @param trme transcriptome in fasta format
#' @param outdir output directory
#' @param RMtr RepeatMasker output file from transcriptome
#' @param RMgen RepeatMasker output file from genome
#' @export
mk.reference_mo <- function(bedtools="bedtools", GFFgen,genome, trme, RMgen,
                            RMtr,overlapping=T, by=c("namRep","classRep", "class", "supFam", "Fam"),
                            threads=1, outdir, over.res=c("HS","LE","LD"), ...){
  message("Identifying repeats in the genome ...")

  system(paste("awk -v OFS='\t' '{if ($3!=exon) {print $1,$4,$5}}'",GFFgen,"| sort -k1,1 -k2,2n > genome.bed"))
  system(paste("awk -v OFS='\t' 'NR>3 {print $5,$6,$7}'",RMgen,"| sort -k1,1 -k2,2n > RMgen.bed"))

  message("Extracting decoys from genome sequences")
  system(paste(bedtools,"intersect -a genome.bed -b RMgen.bed -u > RepGen.bed"))
  system(paste(bedtools,"merge -i RepGen.bed > RepGen_merged.bed"))
  system(paste(bedtools,"getfasta -fi", genome, "-bed RepGen_merged.bed -fo genome_found.fa"))
  system(paste0("awk '{a=$0; getline;split(a, b, ",'":"',");  r[b[1]] = r[b[1]]",'""',"$0} END { for (k in r) { print ", "k\"\\n\"r[k]"," } }' genome_found.fa > decoy.fa"))

  if(overlapping==T){
    RM_tr <- read.RepMask(RMtr)
    RM_tr <- RM_tr[RM_tr$classRep%!in%c("Unknown", "rRNA", "Satellite", "Simple_repeat","Low_complexity","RNA","scRNA","snRNA","srpRNA", "tRNA","Other"),]
    RM.ovlp.res <- ovlp.res(RepMask= RM_tr,
                            cleanTEsProt = F,
                            featureSum = F,
                            outdir = outdir,
                            rm.cotrans = F,
                            over.res=over.res,
                            threads=threads)
  }else{
    RM_tr <- read.RepMask(RMtr)
    RM.ovlp.res <- RM_tr[RM_tr$classRep%!in%c("Unknown", "rRNA", "Satellite", "Simple_repeat","Low_complexity","RNA","scRNA","snRNA","srpRNA", "tRNA","Other"),]
  }


  BED <- cbind(RM.ovlp.res[,c(5,6,7)],RM.ovlp.res[,by])
  BED <- BED[order(as.factor(BED$namSeq),as.factor(BED$SPMQuer)),]

  write.table(BED,"RM.ovres.bed", quote = F, row.names = F, col.names=F, sep="\t")
  system("pwd")
  message("Extracting repeats from transcriptome sequences")
  Ref.salmon <- data.frame(paste0(BED[,1],":",BED[,2],"-",BED[,3]),BED[,4])
  names(Ref.salmon) <- c("seqID","repID")
  write.table(Ref.salmon,paste0(outdir,"/references.csv"), col.names = F, row.names = F, quote = F, sep = ";")

  message("Making trmeSalmon.fasta file")
  system("cat RM.ovres.bed | sort -k1,1 -k2,2n > RM_or.bed")
  system(paste(bedtools,"merge -i RM_or.bed > RM_or_merged.bed"))
  system(paste(bedtools,"getfasta -fi",trme,"-bed RM_or_merged.bed -fo Rep.fa"))
  system("cat Rep.fa decoy.fa > trmeSalmon.fasta")

  message("Making decoys.txt file")
  system("grep '>' decoy.fa | awk '{print substr($1,2); }' > decoys.txt")
  if(!(normalizePath(outdir)==getwd())){
  system(paste("mv trmeSalmon.fasta",outdir))
  system(paste("mv decoys.txt",outdir))
  }
  system("rm genome.bed RMgen.bed RepGen.bed RepGen_merged.bed genome_found.fa Rep.fa RM.ovres.bed RM_or_merged.bed RM_or.bed decoy.fa *.fai")

  Ref.salmon
}
