ExplorATE - *Explore Active Transposable Elements* -
====
## Overview

ExplorATE uses the Selective Alignment strategy (SA) to filter co-transcribed transposons with genes, based on alignment scores. ExplorATE first identifies target TEs and decoy TE sequences, and second performs the quantification of the target TEs using the SA algorithm in Salmon.
ExplorATE allows the TE analysis in multiple organisms (with or without reference genome). If a reference genome is provided, the user could (1) use a set of target TEs from the intergenic regions of the genome, or (2) use a *de novo* transcriptome (and its RepeatMasker annotations) to define the target TEs. Further, ExplorATE uses the reference genome and a genome-derived RepeatMasker file to extract TEs and define decoy sequences. Users can resolve overlapping repeats before extract target TEs.
If the reference genome is not available, ExplorATE uses the repeats identified by RepeatMasker from a *de novo* transcriptome. This transcriptome-derived RepeatMasker file is processed to resolve overlapping repeats and define target TEs and decoy sequences.
Target TEs can be defined at the fragment level or at the transcript level. At the fragment level, all non-co-transcribed TEs are used as target. At the transcript level, a rule similar to Wicker's is used to establish the identity of the transcript. The algorithm uses this criterion to assign target transcripts based on the percentage of identity for a class/family of TEs (calculated as 100 - the percentage of divergence), the percentage for each TE class/family in the transcript (calculated as the ratio between the TE length and the transcript length), and the minimum transcript length. For example, ‘80-80-80’ Wicker-like rule is a selection criterion where a transcript is annotated as "target" f it contains a TE class/family with percentage of identity >80%, this TE class/family represents >80% of the transcript length, and the transcript is at least 80bp in length.
ExplorATE creates files for Salmon execution with Selective Alignment algorithm. The counts estimated by Salmon can be imported into R with specific functions from the R package.
The functions from the ExplorATE package allow users to explore the files and run the pipeline step by step. In addition, ExplorATE includes a shell [script](https://github.com/FemeniasM/ExplorATE_shell_script) that helps users to generate input files or run the pipeline for either model or non-model organisms. For analysis of human or mouse data, users can use the [TESSA](https://github.com/FemeniasM/TESSA) repository that applies Salmon's Selective Alignment on pre-built indexes of references derived from RepeatMasker and GENCODE. The general workflow for model and non-model organisms is shown below.

Flowchart with model organisms
![alt text](./user_guide/flowchartMO.png)

Flowchart with non-model organisms
![alt text](./user_guide/flowchart_nonMO.png)


### Requirements

ExplorATE requires `bedtools` to be installed on your system. If `bedtools` is not installed on your OS, please refer to this [link](https://bedtools.readthedocs.io/en/latest/content/installation.html) and follow the instructions for your OS.

ExplorATE requires some previously installed R packages, select those packages that are not available in your environment:

```{r eval=FALSE}

install.packages(c("stringr","foreach","doParallel"))
install.packages(c("BiocManager","devtools")) 

BiocManager::install(c("readr","GenomicRanges", "IRanges","csaw", "edgeR","SummarizedExperiment","DESeq2", "tximport"))

```

### Install ExplorATE

Install ExplorATE from GitHub

```{r eval=FALSE}

devtools::install_github("FemeniasM/ExplorATEproject")

```

### User guide and vignette

Check the [vignette](https://femeniasm.github.io/ExplorATE_vignette/) and the [user_guide](https://femeniasm.github.io/ExplorATE_user_guide/)

