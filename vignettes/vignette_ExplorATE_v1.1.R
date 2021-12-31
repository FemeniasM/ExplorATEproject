#===================#
# ExporATE vignette #
#===================#

# This vignette implements the ExplorATE package on different data sets, both for model and non-model organisms. Further, this vignette helps users to obtain files necessary for the execution of ExplorATE, such as gene models and reference genomes for model organisms. User can find detailed information on each function in the users' guide.

### 1. Installation ####
# To install the R package ExplorATE, make sure you meet the requirements detailed on the program page and type the following command at the R console:

devtools::install_github("FemeniasM/ExplorATEproject")

# More information about the installation can be found in the user guide

### 2 Obtaining the data test ####

# ExplorATE includes a test data set for model and non-model organisms. Users can create a folder `ExplorATE_data_test` containing the data sets by typing in the linux terminal:

# git clone https://github.com/FemeniasM/ExplorATE_data_test

# Inside the `ExplorATE_data_test` folder three subdirectories are created: `inputs_hs`, `inputs_dm` and `inputs_lp` that contain the data for humans, *Drosophila melanogaster* and *Lioleamus parthenos* respectively. Each of the data sets are described throughout the vignette.
# We will analyze the data set for *Liolaemus parthenos*, the only parthenogenetic lizard of the entire Pleurodonta (Iguanidae) clade. The FASTQ files are available from Gene Expression Omnibus (accession no.`GSE173261`). However in this vignette, we will use a data set with FASTQ files randomly sampled from the original file in order to reduce library sizes and generate replicates for differential expression analysis. Similarly, we subset the RepeatMasker output file and reduced the transcriptome size to run the pipeline in a reasonable time for the vignette. First clone the repository to your local directory and unzip all files.

# The `inputs_lp` folder contains the output files for `BLAST`, `TransDecoder`, `RepeatMasker`, and the *de novo* transcriptome. You can use these files directly (continuing with step 3) or generate them as shown below. To generate the input files with the `nmo_in` mode of the ExplorATE shell script, users can modify the following code:

# bash path/to/ExplorATE_shell_script/ExplorATE nmo_in -p 12 -n <blastp binary path> \\
# -m <hmmerscan binary path> -r <RepeatMasker binary path> -d <TransDecoder directory path> \\
# -u <SwissProt database> -f <Pfam database> -i <user-defined TE library> -t trme_lp.fa \\
# -o inputs_lp


# The script above generates the input files to run the ExplorATE pipeline. For the next steps in the vignette it is assumed that all input files are in a folder named `inputs_lp`

# The script will generate the corresponding input files, except for the transcriptome which must be *de novo* assembled by the user or use the provided transcriptome (`Lp_trme.fasta`). For the execution of the next steps in the vignette we assume that all the input files are in the `inputs_lp` folder.

### 3 Running the ExplorATE pipeline for non-model organisms with functions from the R package ####
# To run the pipeline from the R ExplorATE package, assign the `inputs_lp` folder as the working directory and create the reference files for the salmon run.

setwd("path/to/folder/inputs_lp")
Lp.references <- ExplorATE::mk.reference(
                                      RepMask = "RM_lp.out",
                                      gff3 = "geneModel_lp.gff3",
                                      anot = "blastAnot_lp.outfmt6",
                                      cleanTEsProt = T,
                                      featureSum = T,
                                      outdir = "out_lp",
                                      rm.cotrans = T,
                                      overlapping = T,
                                      trme = "trme_lp.fa",
                                      stranded = F,
                                      by = "classRep",
                                      threads = 12,
                                      rule = c(80,80,80),
                                      over.res= "HS",
                                      annot_by = "transcripts"
                                      )
# Next, the message `"please press 'y' if the names are correct"` is printed on the console to verify that the names assigned to each repetition do not contain ambiguities (consult the users' guide. Press `'y'+ENTER` and ExplorATE continues with the exclusion of elements co-transcribed with genes and the resolution of overlaps by the highest score. Finally, the transcripts with repeats are annotated by class , and the files `decoys.txt` and `trmeSalmon.fasta` are written to the output directory. These files will be used in the next step for abundance estimation in Salmon. With the `featureSum = T` argument, a file `features.summary.csv` is created with the transcripts containing co-transcribed TEs.

# To run Salmon, users can run the following code

ExplorATE::run.salmon(index = "Lparthenos_index",
                      decoys = "out_lp/decoys.txt",
                      salmon_path = "~/programs/salmon-latest_linux_x86_64/bin/salmon",
                      pe_se = "pe",
                      kmer = 31,
                      trme = "out_lp/trmeSalmon.fasta",
                      lib_dir = "reads/",
                      threads = 12
)

# This function creates a folder `quant_out` in the working directory that contains the Salmon estimates and a reference file `references.csv` that are used to import the estimates as shown in the next section.

### 3.1 Importing estimates to R and estimating differential expression ####

# The last step is to import the estimates into the R environment for subsequent analyzes such as differential expression analysis. The `import.RTEs()` function allows to import the estimates and create ready-to-use objects  `DGEList` or `DESeqDataSet` for edgeR or DESeq respectively. Below is an example of how to import the estimates and perform the expression analysis in edgeR.
Lp.references <- read.csv("out_lp/references.csv", sep = ";", header = F)

y <- ExplorATE::import.RTEs(
                      path.sal = "quant_out",
                      ref.sal = Lp.references,
                      import_to = "edgeR",
                      conditions = rep(c("Brain", "Liver", "Ovary"), each = 3)
                      )
# Now "y" is a `DGEList` object that is ready to perform the dispersion estimation.Only estimates for the ERV3 and CR1 families are included in the reduced data set. Users can download the complete set and perform differential expression analysis as follows.

# Users can explore the dispersion of the data through an MDS plot from the limma package:

limma::plotMDS(y)
abline(h=0, v=0)

# The GLM approach require build a design matrix to compare tissues and then the dispersion is estimated with the `estimateDisp()` function:

group.0 <- factor(rep(c("Brain", "Liver","Ovary"), each = 3))
des <- model.matrix( ~ 0 + group.0, data = y$samples)
colnames(des) <- c("Brain", "Liver","Ovary")
y <- edgeR::estimateDisp(y,des)

# Finally the quasi-likelihood (QL) F-test is used to determine the differential expression:

fit_qlf <- edgeR::glmQLFit(y, des)
my.contrasts <- limma::makeContrasts(OvsB=Ovary-Brain,
                                     OvsL=Ovary-Liver,
                                     BvsL=Brain-Liver,
                                     levels=des)

qlf_OvsB <- edgeR::glmQLFTest(fit_qlf, contrast=my.contrasts[,"OvsB"])
qlf_OvsL <- edgeR::glmQLFTest(fit_qlf, contrast=my.contrasts[,"OvsL"])
qlf_BvsL <- edgeR::glmQLFTest(fit_qlf, contrast=my.contrasts[,"BvsL"])

#The TEs differentially expressed in each contrast with FDR <.05 can be explored with the following code:
topTags_OvsB <- edgeR::topTags(qlf_OvsB, n = nrow(qlf_OvsB$table),sort.by = "none")
topTags_OvsB[topTags_OvsB$table$FDR<.05,]

topTags_OvsL <- edgeR::topTags(qlf_OvsL, n = nrow(qlf_OvsL$table),sort.by = "none")
topTags_OvsL[topTags_OvsL$table$FDR<.05,]

topTags_BvsL <- edgeR::topTags(qlf_BvsL, n = nrow(qlf_BvsL$table),sort.by = "none")
topTags_BvsL[topTags_BvsL$table$FDR<.05,]
