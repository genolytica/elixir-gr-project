#!/usr/local/bin/Rscript

# Wrapper for the metaseqR2 script.

suppressPackageStartupMessages(library("optparse"))
source("metaseqR2Counts.R")

option_list <- list(
	make_option(
		opt_str="--targets",
		action="store",
		default=NULL,
		help=paste0(
			"A tab-delimited file with the experimental description.\n",
			"The file should be text tab-delimited and structured as follows:\n",
			"the first line of the external tab delimited file should contain\n",
			"column names (names are not important). The first column MUST contain\n",
			"UNIQUE sample names. The second column MUST contain the raw BAM/BED files\n",
			"WITH their full path. Alternatively, the path argument should be provided.\n",
			"If path is not provided and if the files in the second column of the targets\n",
			"file do not contain a path to a directory, the current directory is assumed\n",
			"to be the BAM/BED file container. The third column MUST contain the biological\n",
			"condition where each of the samples in the first column should belong to."
		)	
	),
	make_option(
		opt_str="--excludelist",
		action="store",
		default=NULL,
		help="A list of samples to exclude, in the same format as the input targets.txt."
	),
#	make_option(
#		opt_str="--filetype",
#		action="store",
#		default="auto",
#		help=paste0(
#			"the type of raw input files. It can be'auto' for auto-guessing,\n",
#			"'bed' for BED files,'sam' for SAM files or 'bam' for BAM files."
#		)
#	),
	make_option(
		opt_str="--path",
		action="store",
		default=NULL,
		help=paste0(
			"An optional path where all the BED/BAM files are placed,\n",
			"to be prepended to the BAM/BED file names in the targets file."
		)
	),
	make_option(
		opt_str="--contrast",
		action="store",
		default=NULL,
		help=paste0(
			"Example: condition1_vs_condition2. Condition1 MUST be the control\n",
			"condition or any reference that condition2 is checked against."
		)
	),
	make_option(
		opt_str="--libsizelist",
		action="store",
		default=NULL,
		help=paste0(
			"An optional named list where names represent samples (MUST be the same\n",
			"as the samples in the input targets.txt in samplelist) and members are\n",
			"the library sizes (sequencing depth) for each sample."
		)
	),
#	make_option(
#		opt_str="--embedCols_id",
#		action="store",
#		default=4,
#		help=paste0(
#			"Regarding the delimited file (or data frame) provided with the counts argument,\n",
#			"with embedded annotation: Unique gene accessions column number.\n",
#			"Default to 4 which is the standard feature name column in a BED file."
#		)
#	),
#	make_option(
#		opt_str="--embedCols_gc",
#		action="store",
#		default=NA,
#		help=paste0(
#			"Regarding the delimited file (or data frame) provided with the counts argument,\n",
#			"with embedded annotation: GC-content column number.\n",
#			"If not provided, GC content normalization provided by EDASeq will not be available."
#		)
#	),
#	make_option(
#		opt_str="--embedCols_name",
#		action="store",
#		default=NA,
#		help=paste0(
#			"Regarding the delimited file (or data frame) provided with the counts argument,\n",
#			"with embedded annotation: HUGO gene symbols column number.\n",
#			"If not provided, it will not be available when reporting results. In addition,\n",
#			"the 'known' gene filter will not be available for application."
#		)
#	),
#	make_option(
#		opt_str="--embedCols_bt",
#		action="store",
#		default=NA,
#		help=paste0(
#			"Regarding the delimited file (or data frame) provided with the counts argument,\n",
#			"with embedded annotation: Gene biotype column number.\n",
#			"If not provided, the 'biodetection','countsbio','saturation','filtered' & 'biodist'\n",
#			"plots will not be available."
#		)
#	),
#	make_option(
#		opt_str="--annotation",
#		action="store",
#		default=NULL,
#		help=paste0(
#			"It can be one of i) NULL (default) to use the existing annotation database or fetch on the fly\n",
#			"ii) 'embedded' if the annotation elements are embedded in the read counts file (restrictions apply)\n",
#			"iii) a list with a path to a GTF file and certain required metadata."
#		)
#	),
	make_option(
		opt_str="--org",
		action="store",
		default="hg19",
		help=paste0("For human genomes 'hg18', 'hg19' or 'hg38'\n",
							 "for mouse genomes 'mm9', 'mm10'\n", 
							 "for rat genomes 'rn5' or 'rn6'\n",
							 "for drosophila genome 'dm3' or 'dm6'\n",
							 "for zebrafish genome 'danrer7', 'danrer10' or 'danrer11'\n", 
							 "for chimpanzee genome 'pantro4', 'pantro5'\n", 
							 "for pig genome 'susscr3', 'susscr11\n'", 
							 "for Arabidopsis thaliana genome 'tair10'\n" ,
							 "for Equus caballus genome 'equcab2'."
		)
	),
	make_option(
		opt_str="--refdb",
		action="store",
		default="ensembl",
		help=paste0(
			"The reference annotation repository from which to retrieve annotation elements\n",
			"to use with metaseqr2. It can be one of'ensembl'(default), 'ucsc' or 'refseq'\n",
			"or a user based one (similar to the org argument)."
		)
	),
	make_option(
		opt_str="--version",
		action="store",
		default="auto",
		help=paste0(
			"An integer denoting the version of the annotation to use from the local annotation\n",
			"database or fetch on the fly. For Ensembl, it corresponds to Ensembl releases, while\n",
			"for UCSC/RefSeq, it is the date of creation (locally)."
		)
	),
	make_option(
		opt_str="--translevel",
		action="store",
		default="gene",
		help=paste0(
			"Perform differential expression analysis at which transcriptional unit, can be one of 'gene'(default),\n",
			"'transcript' for reporting differential expression at the transcript level or 'exon' for exon level."
		)
	),
	make_option(
		opt_str=c("--counttype"),
		action="store",
		default="gene",  #MUST FIX
		help=paste0(
			"The type of reads inside the counts file. It can be one of 'gene', 'exon' or 'utr'\n",
			"for Quant-Seq (Lexogen) protocol. This is a very important and mandatory parameter\n",
			"as it defines the course of the workflow."
		)
	),
	make_option(
		opt_str=c("--utrOpts_frac"),
		action="store",
		default=1,
		help="The fraction (0-1) of the 3' UTR region to count reads in."
	),
	make_option(
		opt_str=c("--utrOpts_minlen"),
		action="store",
		default=300,
		help="The minimum acceptable 3'UTR length irrespective of utrOpts_frac argument."
	),
	make_option(
		opt_str=c("--utrOpts_dnstrm"),
		action="store",
		default=50,
		help=paste0(
			"The number of base pairs to flank the end of the 3' UTR of transcripts\n",
			"when analyzing Quant-Seq data."
		)
	),
#------------------------------------------------------------------------------------------	
#	make_option(
#		opt_str="--exonfltr",
#		action="store_true",
#		default="list",
#		help=paste0(
#			"The supported exon filter in the current version is minActiveExons which\n",
#			"implements a filter for demanding m out of n exons of a gene to have a certain\n",
#			"read presence with parameters exonfltr_exonsprgene, exonfltr_minexons and exonfltr_frac.\n",
#			"The filter is described as follows: if a gene has up to exonfltr_exonsprgene exons,\n",
#			"then read presence is required in at least exonfltr_minexons of them, else read presence\n",
#			"is required in a exonfltr_frac fraction of the total exons. With the default values, the filter\n",
#			"instructs that if a gene has up to 5 exons, read presence is required in at least 2, else\n",
#			"in at least 20 exons, in order to be accepted. Set exonfltr=NULL to not apply any exon filtering.\n",
#			"To apply your own filter parameters, change arguments exonfltr_exonsprgene, exonfltr_minexons, exonfltr_minexons."  
#		)
#	),
	make_option(
		opt_str="--exonfltr_exonsprgene",
		action="store",
		default=5,
		help="Exon filter: Exons per gene. Defaults to 5."
	),
	make_option(
		opt_str=c("--exonfltr_minexons"),
		action="store",
		default=2,
		help=paste0(
			"Exon filter: read presence is required in at least exonfltr_minexons of\n",
			"exonfltr_exonsprgene. Defaults to 2."
		)
	),
	make_option(
		opt_str=c("--exonfltr_frac"),
		action="store",
		default=1/5,
		help=paste0(
			"Exon filter: if read presence in at least exonfltr_minexons of exonfltr_exonsprgene\n",
			"not true, read presence is required in a exonfltr_frac fraction of the total exons.\n",
			"Defaults to 1/5."
		)
	),
#------------------------------------------------------------------------------------------	
	make_option(
		opt_str="--genefltr1_length",
		action="store",
		default=500,
		help=paste0(
			"Gene filter1: length filter where genes are accepted for further analysis if\n",
			"they are above genefltr_length. Defaults to 500."
		)
	),
	make_option(
		opt_str=c("--genefltr2_avgperbp"),
		action="store",
		default=100,
		help=paste0(
			"Gene filter2: a gene is accepted for further analysis if it has more average\n",
			"reads than the genefltr_avgquantile of the average count distribution per\n",
			"genefltr_avgperbp base pairs. Defaults to 100."
		)
	),
	make_option(
		opt_str=c("--genefltr2_avgquantile"),
		action="store",
		default=0.25,
		help=paste0(
			"Gene filter2: length filter where genes are accepted for further analysis if\n",
			"they are above genefltr_length. Defaults to 500."
		)
	),
	make_option(
		opt_str=c("--genefltr3_expmedian"),
		action="store",
		default=TRUE,
		help=paste0(
			"Gene filter3: based on the overall expression of a gene. Genes below the median\n",
			"of the overall count distribution are not accepted for further analysis. Defaults to TRUE."
		)
	),
	make_option(
		opt_str=c("--genefltr3_expmean"),
		action="store",
		default=FALSE,
		help=paste0(
			"Gene filter3: based on the overall expression of a gene. Genes below the mean\n",
			"of the overall count distribution are not accepted for further analysis. Defaults to FALSE."
		)
	),
	make_option(
		opt_str=c("--genefltr3_expquantile"),
		action="store",
		default=NA,
		help=paste0(
			"Gene filter3: based on the overall expression of a gene. Genes below the specified\n",
			"quantile of the total counts distribution are not accepted for further analysis."
		)
	),
	make_option(
		opt_str=c("--genefltr3_expknown"),
		action="store",
		default=NA,
		help=paste0(
			"Gene filter3: based on the overall expression of a gene. A set of known not-expressed\n",
			"genes in the system under investigation are used to estimate an expression cutoff.\n",
			"The value of this filter is a character vector of HUGO gene symbols (MUST be contained\n",
			"in the annotation, thus it's better to use annotation='download') whose counts are used\n",
			"to build a 'null' expression distribution. The 90th quantile of this distribution is then."
		)
	),
	make_option(
		opt_str=c("--genefltr4_biotype"),
		action="store",
		default="getDefaults('biotype.filter',org[1])",
		help=paste0(
			"Gene filter4: genes with a certain biotype (MUST be contained in the annotation,\n",
			"thus it's better to use annotation='download') are excluded from the analysis.\n",
			"This filter is a named list of logical, where names are the biotypes in each genome\n",
			"and values are TRUE or FALSE. If the biotype should be excluded, the value should be TRUE else FALSE."
		)
	),
	make_option(
		opt_str=c("--genefltr5_frac"),
		action="store",
		default=0.25,
		help=paste0(
			"Gene filter5: a gene is further considered for statistical testing if genefltr5_frac\n",
			"(x100 for a percentage value) have more than genefltr5_mincount reads across all samples\n",
			"(genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE)\n",
			"Defaults to 0.25."
		)
	),
	make_option(
		opt_str=c("--genefltr5_mincount"),
		action="store",
		default=10,
		help=paste0(
			"Gene filter5: a gene is further considered for statistical testing if genefltr5_frac\n",
			"(x100 for a percentage value) have more than genefltr5_mincount reads across all samples\n",
			"(genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE)\n",
			"Defaults to 10."
		)
	),
	make_option(
		opt_str=c("--genefltr5_percon"),
		action="store",
		default=FALSE,
		help=paste0(
			"Gene filter5: a gene is further considered for statistical testing if genefltr5_frac\n",
			"(x100 for a percentage value) have more than genefltr5_mincount reads across all samples\n",
			"(genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE)\n",
			"Defaults to FALSE."
		)
	),	
#------------------------------------------------------------------------------------------	
	make_option(
		opt_str=c("--whenapplyfilter"),
		action="store",
		default="postnorm",
		help=paste0(
			"A character string determining when to apply the exon and/or gene filters,\n",
			"relative to normalization. It can be 'prenorm' to apply apply the filters\n",
			"and exclude genes from further processing before normalization, or 'postnorm'\n",
			"to apply the filters after normalization (default)."
		)
	),
#	make_option(
#		opt_str=c("--normalization"),
#		action="store",
#		default="deseq",
#		help=paste0(
#			"The normalization algorithm to be applied on the count data.\n",
#			"It can be one of 'edaseq' for EDASeq normalization, 'deseq' for the normalization\n",
#			"algorithm in the DESq package (default), 'edger' for the normalization algorithms\n",
#			"present in the edgeR package 'noiseq' for the normalization algorithms present in\n",
#			"the NOISeq package 'nbpseq' for the normalization algorithms present in the NBPSeq\n",
#			"package or 'none' to not normalize the data (highly unrecommended). Algorithm specific\n",
#			"arguments can be passed through the normargs argument)."
#		)
#	),
#	make_option(
#		opt_str=c("--normargs"),
#		action="store",
#		default=NULL,
#		help=paste0(
#			"A named list whose names are the names of the normalization algorithm parameters and\n",
#			"its members parameter values. You should check the documentation of the packages\n",
#			"EDASeq, DESeq, edgeR, NOISeq and NBPSeq for the parameter names and parameter values."
#		)
#	),
#	make_option(
#		opt_str=c("--statistics"),
#		action="store",
#		default="deseq",
#		help=paste0(
#			"One or more statistical analyses to be performed by the metaseqr2 pipeline.\n",
#			"It can be one or more of 'deseq' (default) to conduct statistical test(s)\n",
#			"implemented in the DESeq package, 'edger' to conduct statistical test(s)\n",
#			"implemented in the edgeR package, 'limma' to conduct the RNA-Seq version of\n",
#			"statistical test(s) implemented in the limma package, 'noiseq' to conduct\n",
#			"statistical test(s) implemented in the NOISeq package, 'bayseq' to conduct\n",
#			"statistical test(s) implemented in the baySeq package, 'nbpseq' to conduct\n",
#			"statistical test(s) implemented in the NBPSeq package, 'deseq2' to conduct\n",
#			"statistical test(s) implemented in the DESeq2 package, 'dss' to conduct\n",
#			"statistical test(s) implemented in the DSS package and 'absseq' to conduct\n",
#			"statistical test(s) implemented in the ABSSeq package."
#		)
#	),
#	make_option(
#		opt_str=c("--statargs"),
#		action="store",
#		default=NULL,
#		help=paste0(
#			"A named list whose names are the names of the statistical algorithms used in the pipeline.\n",
#			"Each member is another named list whose names are the algorithm parameters and its members\n",
#			"are the parameter values. You should check the documentations of DESeq, edgeR, NOISeq, baySeq,\n",
#			"limma and NBPSeq for these parameters."
#		)
#	),	
#	make_option(
#		opt_str=c("--adjmethod"),
#		action="store",
#		default="BH",
#		help=paste0(
#			"The multiple testing p-value adjustment method. It can be one of p.adjust.methods or\n",
#			"'qvalue' from the qvalue Bioconductor package. Defaults to 'BH' for Benjamini-Hochberg correction."
#		)
#	),
#	make_option(
#		opt_str=c("--metap"),
#		action="store",
#		default="simes",
#		help=paste0(
#			"The meta-analysis method to combine p-values from multiple statistical tests.\n",
#			"It can be one of 'simes'(default), 'bonferroni', 'minp', 'maxp', 'weight, 'pandora',\n",
#			"'dperm_min', 'dperm_max', 'dperm_weight', 'fisher', 'fperm', 'whitlock' or 'none'."
#		)
#	),
#	make_option(
#		opt_str=c("--weight"),
#		action="store",
#		default="rep(1/length(statistics),length(statistics))",
#		help=paste0(
#			"A vector of weights with the same length as the statistics vector containing a weight\n",
#			"for each statistical test. It should sum to 1."
#		)
#	),
#	make_option(
#		opt_str=c("--nperm"),
#		action="store",
#		default=10000,
#		help=paste0(
#			"The number of permutations performed to derive the meta p-value when metap='fperm'\n",
#			"or metaP='dperm'. It defaults to 10000."
#		)
#	),
#	make_option(
#		opt_str=c("--pcut"),
#		action="store",
#		default=NA,
#		help=paste0(
#			"a p-value cutoff for exporting differentially genes, default is\n",
#			"to export all the non-filtered genes."
#		)
#	),	
#	make_option(
#		opt_str=c("--logoffset"),
#		action="store",
#		default=1,
#		help=paste0(
#			"An offset to be added to values during logarithmic transformations\n",
#			"in order to avoid Infinity (default is 1)."
#		)
#	),
#	make_option(
#		opt_str=c("--preset"),
#		action="store",
#		default=NULL,
#		help=paste0(
#			"An analysis strictness preset. preset can be one of 'all_basic',\n",
#			"'all_normal', 'all_full', 'medium_basic', 'medium_normal', 'medium_full',\n",
#			"'strict_basic', 'strict_normal' or 'strict_full', each of which control\n",
#			"the strictness of the analysis and the amount of data to be exported."
#		)
#	),
#	make_option(
#		opt_str=c("--qcplots"),
#		action="store",
#		default=NULL,
#		help=paste0(
#			"A set of diagnostic plots to show/create. It can be one or more of\n",
#			"'mds', 'biodetection', 'rnacomp', 'countsbio', 'saturation', 'readnoise',\n",
#			"'filtered', 'boxplot', 'gcbias', 'lengthbias', 'meandiff', 'meanvar', 'deheatmap',\n",
#			"'volcano', 'mastat', 'biodist', 'statvenn', 'foldvenn'.\n",
#			"Defaults to qcPlots=NULL -no diagnostic plots will be created."
#		)
#	),	
	make_option(
		opt_str=c("--figformat"),
		action="store",
#		default=0,
		help=paste0(
			"The format of the output diagnostic plots. It can be one or more of 'png', 'jpg',\n",
			"'tiff', 'bmp', 'pdf', 'ps'. The native format 'x11' (for direct display) is not provided\n",
			"as an option as it may not render the proper display of some diagnostic plots in some devices."
		)
	),
	make_option(
		opt_str=c("--outlist"),
		action="store",
		default=FALSE,
		help=paste0(
			"A logical controlling whether to export a list with the results\n",
			"in the running environment."
		)
	),
	make_option(
		opt_str=c("--xprtwhere"),
		action="store",
		default=NA,
		help="An output directory for the project results."
	),
	make_option(
		opt_str=c("--xprtwhat"),
		action="store",
#		default="counts",
		help=paste0(
			"The content of the final lists.\n",
			"It can be one or more of 'annotation',to bind the annotation elements for each gene,\n",
			"'p_value', to bind the p-values of each method,\n",
			"'adj_p_value', to bind the multiple testing adjusted p-values,\n",
			"'meta_p_value', to bind the combined p-value from the meta-analysis,\n",
			"'adj_meta_p_value', to bind the corrected combined p-value from the meta-analysis,\n",
			"'fold_change', to bind the fold changes of each requested contrast,\n",
			"'stats', to bind several statistics calclulated on raw and normalized counts(xprtstats argument),\n",
			"'counts', to bind the raw and normalized counts for each sample."
		)
	),
	make_option(
		opt_str=c("--xprtscale"),
		action="store",
#		default=,
		help=paste0(
			"Export values from one or more transformations applied to the data. It can be one or more of\n",
			"'natural', 'log2', 'log10', 'vst' (Variance Stabilizing Transormation, documentation in DESeq package)\n",
			"and 'rpgm' which is ratio of mapped reads per gene model (either the gene length or the sum of\n",
			"exon lengths, depending on countType argument). Note that this is not RPKM as reads are already\n",
			"normalized for library size using one of the supported normalization methods. Also, 'rpgm' might\n",
			"be misleading when normalization is other than 'deseq'."
		)
	),
	make_option(
		opt_str=c("--xprtvalues"),
		action="store",
#		default=normalized,
		help=paste0(
			"It can be one or more of 'raw' to export raw values (counts etc.)\n",
			"and 'normalized' to export normalized counts."
		)
	),
	make_option(
		opt_str=c("--xprtstats"),
		action="store",
#		default=0,
		help=paste0(
			"Calculate and export several statistics on raw and normalized counts, condition-wise.\n",
			"It can be one or more of 'mean', 'median', 'sd', 'mad', 'cv' for the Coefficient of Variation,\n",
			"'rcv' for a robust version of CV where the median and the MAD are used instead of the mean and\n",
			"the standard deviation."
		)
	),
	make_option(
		opt_str=c("--xprtcountstbl"),
		action="store",
		default=FALSE,
		help=paste0(
			"Exports the calculated read counts table when input is read from bam files and\n",
			"the normalized count table in all cases. Defaults to FALSE."
		)
	),	
	make_option(
		opt_str=c("--rc"),
		action="store",
		default=0.6,
		help=paste0(
			"In case of parallel execution of several subfunctions, the fraction of the\n",
			"available cores to use. In some cases if all available cores are used (rc=1)\n",
			"and the system does not have sufficient RAM, the pipeline running machine\n",
			"might significantly slow down."
		)
	),
	make_option(
		opt_str=c("--report"),
		action="store",
		default=TRUE,
		help=paste0(
			"A logical value controlling whether to produce a summary report or not.\n",
			"Defaults to TRUE."
		)
	),
	make_option(
		opt_str=c("--topreport"),
		action="store",
		default=0.1,
		help=paste0(
			"A fraction of top statistically significant genes to append to the HTML report.\n",
			"This helps in keeping the size of the report as small as possible, as appending\n",
			"the total gene list might create a huge HTML file. Users can always retrieve the\n",
			"whole gene lists from the report links. Defaults to 0.1 (top 10 genes).\n",
			"Set to NA or NULL to append all the statistically significant genes to the HTML report."
		)
	),
	make_option(
		opt_str=c("--templatereport"),
		action="store",
		default="default",
		help=paste0(
			"An HTML template to use for the report.\n",
			"Do not change this unless you know what you are doing."
		)
	),
	make_option(
		opt_str=c("--genemodel"),
		action="store",
		default=TRUE,
		help=paste0(
			"In case of exon analysis, a list with exon counts for each gene will be saved\n",
			"to the file xprtwhere/data/gene_model.RData. This file can be used as input\n",
			"to metaseqR for exon count based analysis, in order to avoid the time consuming\n",
			"step of assembling the counts for each gene from its exons"
		)
	),
	make_option(
		opt_str=c("--verbose"),
		action="store",
		default=TRUE,
		help="Print informative messages during execution? Defaults to TRUE."
	),
	make_option(
		opt_str=c("--runlog"),
		action="store",
		default=TRUE,
		help=paste0(
			"Write a log file of the metaseqr2 run using package log4r. Defaults to TRUE.\n",
			"The filename will be auto-generated."
		)
	),
	make_option(
		opt_str=c("--reportdb"),
		action="store",
		default="sqlite",
		help=paste0(
			"Database system to use for storing the report intereactive graphs.\n",
			"Can be 'sqlite' (default) or 'dexie'."
		)
	),
	make_option(
		opt_str=c("--localdb"),
		action="store",
		default="file.path(system.file(package='metaseqR2'),'annotation.sqlite')",
		help="The metaseqR2 annotation database location."
	),
	make_option(
		opt_str=c("--progressfun"),
		action="store",
		default=NULL,
		help=paste0(
			"A function which updates a Progress object from shiny. This function\n",
			"must accept a detail argument. See http://shiny.rstudio.com/articles/progress.html"
		)
	),
	make_option(
		opt_str=c("--offlinereport"),
		action="store",
		default=TRUE,
		help=paste0(
			"TRUE (default) to download and include the required JavaScript libraries\n",
			"to properly view the report offline. Ignored if report=FALSE"
		)
	),
#	make_option(
#		opt_str=c("--createtracks"),
#		action="store",
#		default=FALSE,
#		help=paste0(
#			"Option to create normalized bigWig files to display in a genome browser\n",
#			"(e.g. UCSC). Defaults to FALSE."
#		)
#	),
#	make_option(
#		opt_str=c("--overwritetracks"),
#		action="store",
#		default=FALSE,
#		help="Overwrite tracks if they already exist? Defaults to FALSE."
#	),
#	make_option(
#		opt_str=c("--trackspath"),
#		action="store",
#		default="file.path(exportWhere,'tracks')",
#		help=paste0(
#			"where to export the bigWig files,\n",
#			"defaults to file.path(exportWhere,'tracks')."
#		)
#	),
	make_option(
		opt_str=c("--exportr2c"),
		action="store",
		default=TRUE,
		help=paste0(
			"EXPORTS READ TO COUNTS\n",
			"Defaults to TRUE."
		)
	),	
)

opt <- parse_args(OptionParser(option_list=option_list))

# TODO: more checks
if (!(opt$org %in% c("hg18", "hg19", "hg38", "mm9","mm10", "rn5", "rn6", "dm3", "dm6", "danrer7","pantro4", "susscr3", "tair10", "equcab2" )))
	stop("The organism must be one of \"hg18\", \"hg19\", \"hg38\", \"mm9\",\"mm10\", \"rn5\", \"rn6\", \"dm3\", \"dm6\", \"danrer7\",\"pantro4\", \"susscr3\", \"tair10\", \"equcab2\"!")

#counts=opt$counts,
#createTracks=opt$createtracks,overwriteTracks=opt$overwritetracks,trackExportPath=opt$trackspath,
#sampleList=opt$targets,

metaseqR2(targets=opt$targets,excludeList=opt$excludelist,fileType=opt$filetype,
	path=opt$path,contrast=opt$contrast,libsizeList=opt$libsizelist,embedCols=list(idCol=opt$embedCols_id,
	gcCol=opt$embedCols_gc,nameCol=opt$embedCols_name,btCol=opt$embedCols_bt),annotation=opt$annotation,
	org=opt$org,refdb=opt$refdb,version=opt$version,transLevel=opt$translevel,countType=opt$counttype,
	utrOpts=list(frac=opt$utrOpts_frac,minLength=opt$utrOpts_minlen,downstream=opt$utrOpts_dnstrm),
#	minActiveExons=opt$exonfltr
	exonFilters=list(minActiveExons=list(exonsPerGene=opt$exonfltr_exonsprgene,minExons=opt$exonfltr_minexons,frac=opt$exonfltr_frac)),
	geneFilters=list(length=list(length=opt$genefltr1_length),avgReads=list(averagePerBp=opt$genefltr2_avgperbp,quantile=opt$genefltr2_avgquantile),
	expression=list(median=opt$genefltr3_expmedian,mean=opt$genefltr3_expmean,quantile=opt$genefltr3_expquantile,known=opt$genefltr3_expknown),
	biotype=opt$genefltr4_biotype,	presence=list(frac=opt$genefltr5_frac,minCount=opt$genefltr5_mincount,perCondition=opt$genefltr5_percon)),
	whenApplyFilter=opt$whenapplyfilter,normalization=opt$normalization,normArgs=opt$normargs,statistics=opt$statistics,
	statArgs=opt$statargs,adjustMethod=opt$adjmethod,metaP=opt$metap,weight=opt$weight,nperm=opt$nperm,pcut=opt$pcut,
	logOffset=opt$logoffset,preset=opt$preset,qcPlots=opt$qcplots,figFormat=opt$figformat,outList=opt$outlist,
	exportWhere=opt$xprtwhere,exportWhat=opt$xprtwhat,exportScale=opt$xprtscale,exportValues=opt$xprtvalues,
	exportStats=opt$xprtstats,exportCountsTable=opt$xprtcountstbl,restrictCores=opt$rc,report=opt$report,
	reportTop=opt$topreport,reportTemplate=opt$templatereport,saveGeneModel=opt$genemodel,verbose=opt$verbose,
	runLog=opt$runlog,reportDb=opt$reportdb,localDb=opt$localdb,progressFun=opt$progressfun,offlineReport=opt$offlinereport,
	.exportR2C=opt$exportr2c)
