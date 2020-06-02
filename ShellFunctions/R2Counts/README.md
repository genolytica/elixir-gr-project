# generate-signal-tracks

A shell script for obtaining rna-seq data counts.

## Description

This shell script ```runMetaseqR2Counts.R```, calls the metaseqr2() function
of the metaseqR2 package.
The script can be invoked from the system shell instead of the R shell, using ```Rscript```.

The script performs the following tasks:


## Prerequisite R/Bioconductor packages
The following packages are required to run the script:
- metaseqR2

## Basic Example
```
Rscript run_metaseqR2.R \
  --samplelist=my_targets.txt \
  --contrast=A_vs_B \
  --org=hg19 \
  --counttype=exon \
  --normalization=edger \
  --statistics=edger \
  --figformat=png \
  --xprtwhere=. \
  --rc=0.5 
```

## List of arguments

The following table presents the input arguments in detail:

|Parameter        |Description                                                                                                                                          |
|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|
|samplelist       |A small tab-delimited file with the experiment description. The first line of the external tab delimited file should contain column names (names are not important). The first column MUST contain UNIQUE sample names. The second column MUST contain the raw BAM/BED files WITH their full path. Alternatively, the path argument should be provided. If path is not provided and if the files in the second column of the targets file do not contain a path to a directory, the current directory is assumed to be the BAM/BED file container. The third column MUST contain the biological condition where each of the samples in the first column should belong to.|
|excludelist      |A list of samples to exclude, in the same (list) format as sampleList above.|
|path             |An optional path where all the BED/BAM files are placed, to be prepended to the BAM/BED file names in the targets file.|
|filetype         |The type of raw input files. It can be "auto" for auto-guessing, "bed" for BED files, "sam" for SAM files or "bam" for BAM files.|
|contrast         |A character vector of contrasts to be tested in the statistical testing step(s) of the pipeline. Each element of contrast should STRICTLY have the format "ConditionA_vs_ConditionB_vs_...". Special attention is needed as fold change calculations are based on this argument.|
|org              |For human genomes "hg18", "hg19" or "hg38", for mouse genomes "mm9", "mm10", for rat genomes "rn5" or "rn6", for drosophila genome "dm3" or "dm6", for zebrafish genome "danrer7", "danrer10" or "danrer11", for chimpanzee genome "pantro4", "pantro5", for pig genome "susscr3", "susscr11", for Arabidopsis thaliana genome "tair10" and for Equus caballus genome "equcab2".|
|refdb            |The reference annotation repository from which to retrieve annotation elements to use with metaseqr2. It can be one of "ensembl" (default), "ucsc" or "refseq" or a user based one (similar to the org argument).|
|version          |An integer denoting the version of the annotation to use from the local annotation database or fetch on the fly. For Ensembl, it corresponds to Ensembl releases, while for UCSC/RefSeq, it is the date of creation (locally). Defaults to "auto".|
|translevel       |Perform differential expression analysis at which transcriptional unit, can be one of "gene"(default),"transcript" for reporting differential expression at the transcript level or "exon" for exon level."|
|counttype        |The type of reads inside the counts file. It can be one of "gene", "exon" or "utr" for Quant-Seq (Lexogen) protocol. This is a very important and mandatory parameter as it defines the course of the workflow.|
|utrOpts_frac     |For Quant-Seq (Lexogen) protocol: the fraction (0-1) of the 3' UTR region to count reads in.|
|utrOpts_minlen   |For Quant-Seq (Lexogen) protocol: the minimum acceptable 3'UTR length irrespective of utrOpts_frac argument.|
|utrOpts_dnstrm   |For Quant-Seq (Lexogen) protocol: the number of base pairs to flank the end of the 3' UTR of transcripts when analyzing Quant-Seq data.|
|exonFilters      |
|geneFilters      |
|whenApplyFilter  |A character string determining when to apply the exon and/or gene filters, relative to normalization. It can be "prenorm" to apply the filters and exclude genes from further processing before normalization, or "postnorm" to apply the filters after normalization (default)."|
|normalization    |A normalization algorithm to be applied on the count data. It can be one of "edaseq" for EDASeq normalization, "deseq" for the normalization algorithm in the DESq package (default), "edger" for the normalization algorithms present in the edgeR package "noiseq" for the normalization algorithms present in the NOISeq package "nbpseq" for the normalization algorithms present in the NBPSeq package or "none" to not normalize the data (highly unrecommended). Algorithm specific arguments can be passed through the normArgs argument).|
|statistics       |
|qcPlots          |
|figFormat        |
|outList          |
|exportWhere      |
|exportWhat       |
|exportScale      |
|exportValues     |
|exportStats      |
|exportCountsTable|
|restrictCores    |
|report           |
|reportTop        |
|reportTemplate   |
|saveGeneModel    |
|verbose          |
|runLog           |
|reportDb         |
|localDb          |
|progressFun      |
|offlineReport    |
|.exportR2C       |


## Indicative runtimes