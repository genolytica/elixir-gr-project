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

|Parameter              |Description                                                                                                                                          |
|-----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|
|samplelist             |A small tab-delimited file with the experiment description. The first line of the external tab delimited file should contain column names (names are not important). The first column MUST contain UNIQUE sample names. The second column MUST contain the raw BAM/BED files WITH their full path. Alternatively, the path argument should be provided. If path is not provided and if the files in the second column of the targets file do not contain a path to a directory, the current directory is assumed to be the BAM/BED file container. The third column MUST contain the biological condition where each of the samples in the first column should belong to.|
|excludelist            |A list of samples to exclude, in the same (list) format as sampleList above.|
|path                   |An optional path where all the BED/BAM files are placed, to be prepended to the BAM/BED file names in the targets file.|
|filetype               |The type of raw input files. It can be "auto" for auto-guessing, "bed" for BED files, "sam" for SAM files or "bam" for BAM files.|
|contrast               |A character vector of contrasts to be tested in the statistical testing step(s) of the pipeline. Each element of contrast should STRICTLY have the format "ConditionA_vs_ConditionB_vs_...". Special attention is needed as fold change calculations are based on this argument.|
|libsizelist            ||
|org                    |For human genomes "hg18", "hg19" or "hg38", for mouse genomes "mm9", "mm10", for rat genomes "rn5" or "rn6", for drosophila genome "dm3" or "dm6", for zebrafish genome "danrer7", "danrer10" or "danrer11", for chimpanzee genome "pantro4", "pantro5", for pig genome "susscr3", "susscr11", for Arabidopsis thaliana genome "tair10" and for Equus caballus genome "equcab2".|
|refdb                  |The reference annotation repository from which to retrieve annotation elements to use with metaseqr2. It can be one of "ensembl" (default), "ucsc" or "refseq" or a user based one (similar to the org argument).|
|version                |!An integer denoting the version of the annotation to use from the local annotation database or fetch on the fly. For Ensembl, it corresponds to Ensembl releases, while for UCSC/RefSeq, it is the date of creation (locally). Defaults to "auto".|
|translevel             |Perform differential expression analysis at which transcriptional unit, can be one of "gene"(default),"transcript" for reporting differential expression at the transcript level or "exon" for exon level."|
|counttype              |The type of reads inside the counts file. It can be one of "gene", "exon" or "utr" for Quant-Seq (Lexogen) protocol. This is a very important and mandatory parameter as it defines the course of the workflow.|
|utrOpts_frac           |For Quant-Seq (Lexogen) protocol: the fraction (0-1) of the 3' UTR region to count reads in.|
|utrOpts_minlen         |For Quant-Seq (Lexogen) protocol: the minimum acceptable 3'UTR length irrespective of utrOpts_frac argument.|
|utrOpts_dnstrm         |For Quant-Seq (Lexogen) protocol: the number of base pairs to flank the end of the 3' UTR of transcripts when analyzing Quant-Seq data.|
|exonfltr               |Set as FALSE for no exon filter to be applied.|
|exonfltr_exonsprgene   |minActiveExons filter: Exons per gene. Defaults to 5.|
|exonfltr_minexons      |minActiveExons filter: read presence is required in at least exonfltr_minexons of exonfltr_exonsprgene. Defaults to 2.|
|exonfltr_frac          |minActiveExons filter: if read presence in at least exonfltr_minexons of exonfltr_exonsprgene not true, read presence is required in a exonfltr_frac fraction of the total exons. Defaults to 1/5.|
|genefltr               |Set as FALSE for no gene filters to be applied.|
|genefltr1_length       |Gene filter1: length filter where genes are accepted for further analysis if they are above genefltr_length. Defaults to 500.|
|genefltr2_avgperbp     |Gene filter2: a gene is accepted for further analysis if it has more average reads than the genefltr_avgquantile of the average count distribution per genefltr_avgperbp base pairs. Defaults to 100.|
|genefltr2_avgquantile  |Gene filter2: a gene is accepted for further analysis if it has more average reads than the genefltr_avgquantile of the average count distribution per genefltr_avgperbp base pairs. Defaults to 0.25.|
|genefltr3_expmedian    |Gene filter3: based on the overall expression of a gene. Genes below the median of the overall count distribution are not accepted for further analysis. Defaults to TRUE.|
|genefltr3_expmean      |Gene filter3: based on the overall expression of a gene. Genes below the mean of the overall count distribution are not accepted for further analysis. Defaults to FALSE.|
|genefltr3_expquantile  |Gene filter3: based on the overall expression of a gene. Genes below the specified quantile of the total counts distribution are not accepted for further analysis.|
|genefltr3_expknown     |Gene filter3: based on the overall expression of a gene. A set of known not-expressed genes in the system under investigation are used to estimate an expression cutoff. The value of this filter is a character vector of HUGO gene symbols (MUST be contained in the annotation, thus it's better to use annotation='download') whose counts are used to build a 'null' expression distribution. The 90th quantile of this distribution is then.|
|genefltr4_biotype      |!Gene filter4: genes with a certain biotype (MUST be contained in the annotation, thus it's better to use annotation='download') are excluded from the analysis. This filter is a named list of logical, where names are the biotypes in each genome and values are TRUE or FALSE. If the biotype should be excluded, the value should be TRUE else FALSE.|
|genefltr5_frac         |Gene filter5: a gene is further considered for statistical testing if genefltr5_frac (x100 for a percentage value) have more than genefltr5_mincount reads across all samples, (genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE). Defaults to 0.25.|
|genefltr5_mincount     |Gene filter5: a gene is further considered for statistical testing if genefltr5_frac (x100 for a percentage value) have more than genefltr5_mincount reads across all samples (genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE). Defaults to 10.|
|genefltr5_percon       |Gene filter5: a gene is further considered for statistical testing if genefltr5_frac (x100 for a percentage value) have more than genefltr5_mincount reads across all samples (genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE). Defaults to FALSE.|
|whenApplyFilter        |!A character string determining when to apply the exon and/or gene filters, relative to normalization. It can be "prenorm" to apply the filters and exclude genes from further processing before normalization, or "postnorm" to apply the filters after normalization (default)."|
|normalization          |!A normalization algorithm to be applied on the count data. It can be one of "edaseq" for EDASeq normalization, "deseq" for the normalization algorithm in the DESq package (default), "edger" for the normalization algorithms present in the edgeR package "noiseq" for the normalization algorithms present in the NOISeq package "nbpseq" for the normalization algorithms present in the NBPSeq package or "none" to not normalize the data (highly unrecommended). Algorithm specific arguments can be passed through the normArgs argument).|
|statistics             |!...|
|qcplots                ||
|figformat              |The format of the output diagnostic plots. It can be one or more of "png", "jpg", "tiff", "bmp", "pdf", "ps". The native format "x11" (for direct display) is not provided as an option as it may not render the proper display of some diagnostic plots in some devices.|
|outlist                |!A logical controlling whether to export a list with the results in the running environment. Defaults to FALSE.|
|xprtwhere              |An output directory for the project results (report, lists, diagnostic plots etc).|
|xprtwhat               ||
|xprtscale              ||
|xprtvalues             ||
|xprtstats              ||
|xprtcountstbl          ||
|rc                     ||
|report                 ||
|topreport              ||
|templatereport         ||
|genemodel              ||
|verbose                ||
|runlog                 ||
|reportdb               ||
|localdb                ||
|progressfun            ||
|offlinereport          ||
|exportr2c              ||

## Indicative runtimes