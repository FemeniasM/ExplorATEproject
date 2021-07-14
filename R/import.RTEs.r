#' @import tximport
#' @import edgeR
#' @import SummarizedExperiment
#' @import csaw
#' @import DESeq2
#' @title Loads retrotransposons counts and return a DGEList object for use with edgeR.
#' @description Imports count estimates and calculates an offset that corrects for changes to the average transcript length across samples (details in https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#edgeR). The function returns a DGEList for use in downstream analysis with edgeR
#' @param path.sal output salmon directory
#' @param ref.sal salmon reference build with mk.reference function
#' @param conditions A string vector indicating the experimental conditions of each library. Example: c("S1", "S1", "S1", "S2", "S2", "S2", "S3", "S3", "S3")
#' @param import_to A string object indicating whether to import for edgeR or DESeq2. If none is selected, transcript-level estimates are returned.
#' @export
import.RTEs <- function(path.sal, ref.sal, conditions, import_to = c("edgeR, DESeq2")){
  dir <- list.dirs(path.sal, recursive = F, full.names = F)
  files <- file.path(path.sal, dir,'quant.sf')
  names(files) <- dir

  if(all(file.exists(files))==F){stop("error loading files")}else{

    group.0 <- factor(conditions)
    des <- cbind(names(files),conditions)
    message("The following experimental design has been created:\n")
    message(message=for(i in 1:nrow(des)){print(paste(des[i,1],des[i,2]))})

    txi_rte <- tximport::tximport(files, type = "salmon",tx2gene = ref.sal)
    cts <- txi_rte$counts
    normMat <- txi_rte$length
    tpm <- txi_rte$abundance

    message("calculating offsets")

    if(import_to=="edgeR"){
    normMat <- normMat/exp(rowMeans(log(normMat)))
    normCts <- cts/normMat

    eff.lib <- edgeR::calcNormFactors(normCts) * colSums(normCts)

    normMat <- sweep(normMat, 2, eff.lib, "*")
    normMat <- log(normMat)

    y <- edgeR::DGEList(cts,group = group.0)
    y <- edgeR::scaleOffset(y, normMat)

    keep <- edgeR::filterByExpr(y)
    y <- y[keep, ]

    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = y$counts, offset = y$offset))
    se$totals <- y$samples$lib.size
    y$log.CPM <- csaw::calculateCPM(se, use.offsets = TRUE, log = TRUE)

    message("DGEList object is ready")
    y
    }else{
      if(import_to=="DESeq2"){
      sampleTable <- data.frame(condition = factor(conditions))
      rownames(sampleTable) <- colnames(txi_rte$counts)

      dds <- DESeqDataSetFromTximport(txi_rte, sampleTable, ~condition)
      dds
      }else{tximport::tximport(files, type = "salmon",txOut=T)}
    }
  }
}
