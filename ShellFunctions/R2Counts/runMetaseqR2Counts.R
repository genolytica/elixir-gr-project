#! /usr/local/bin/Rscript

# Wrapper for the metaseqr2() function

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(metaseqR2))

option_list <- list(
	make_option(
		opt_str="--samples",
		action="store",
		help="Sample IDs - Enter as comma-separated list (no space)."
	),
	make_option(
		opt_str="--files",
		action="store",
		help="File names - Enter as comma-separated list (no space)."
	),
	make_option(
		opt_str="--conditions",
		action="store",
		help="Sample Conditions - Enter as comma-separated list (no space)."
	),
	make_option(
		opt_str="--paired",
		action="store",
		default=NULL,
		help=paste0(
			"'single' for single-end reads, 'paired' for paired-end reads or 'mixed' for BAMs\n",
			"that contain both single- and paired-end reads. If this column is not provided,\n",
			"single-end reads will be assumed. - Enter as comma-separated list (no space)."
		)
	),
	make_option(
		opt_str="--strandp",
		action="store",
		default=NULL,
		help=paste0(
			"'forward' for a forward (5'->3') strand library construction protocol,\n",
			"'reverse' for a reverse (3'->5') strand library construction protocol, or\n",
			"'no' for unstranded/unknown protocol. If this column is not provided,\n",
			"unstranded reads will be assumed - Enter as comma-separated list (no space)."
		)
	),
#	make_option(
#		opt_str="--counts",
#		action="store",
#		help=paste0(
#			"A text tab-delimited file containing gene, exon or 3'UTR counts in one of the following formats:\n",
#			"i) the first column contains unique gene or exon identifiers and the rest of the columns contain\n",
#			"the read counts for each sample.\n",
#			"ii) The first n columns should contain only **gene** annotation elements like chromosomal locations,\n",
#			"gene accessions, exon accessions, GC content etc. and the rest columns should contain gene read counts.\n",
#			"iii) counts can also be an .RData file with previous analysis elements\n",
#			"iv) counts can be a list representing the gene model."
#		)
#	),
#	make_option(
#		opt_str="--samplelist",
#		action="store",
#		help=paste0(
#			"A tab-delimited file with the experimental description.\n",
#			"The file should be text tab-delimited and structured as follows:\n",
#			"the first line of the external tab delimited file should contain\n",
#			"column names (names are not important). The first column MUST contain\n",
#			"UNIQUE sample names. The second column MUST contain the raw BAM/BED files\n",
#			"WITH their full path. Alternatively, the path argument should be provided.\n",
#			"If path is not provided and if the files in the second column of the targets\n",
#			"file do not contain a path to a directory, the current directory is assumed\n",
#			"to be the BAM/BED file container. The third column MUST contain the biological\n",
#			"condition where each of the samples in the first column should belong to."
#		)	
#	),
	make_option(
		opt_str="--excludelist",
		action="store",
		default=NULL,
		help="A list of samples to exclude, in the same format as the input targets.txt."
	),
	make_option(
		opt_str="--path",
		action="store",
#		default=NULL,
		help=paste0(
			"An optional path where all the BED/BAM files are placed,\n",
			"to be prepended to the BAM/BED file names in the targets file."
		)
	),
	make_option(
		opt_str="--filetype",
		action="store",
		default="auto",
		help=paste0(
			"The type of raw input files. It can be'auto' for auto-guessing,\n",
			"'bed' for BED files,'sam' for SAM files or 'bam' for BAM files."
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
	make_option(
		opt_str="--idcol",
		action="store",
		default=4,
		help=paste0(
			"Integer denoting the column number in the file (or data frame) provided with\n",
			"the counts argument, where the unique gene accessions are. Default to 4 which\n",
			"is the standard feature name column in a BED file."
		)
	),
	make_option(
		opt_str="--gccol",
		type="integer",
		action="store",
		default=NA,
		help=paste0(
			"Integer denoting the column number in the file (or data frame) provided with the counts\n",
			"argument, where each gene's GC content is given. If not provided, GC content normalization\n",
			"provided by EDASeq will not be available."
		)
	),
	make_option(
		opt_str="--namecol",
		type="integer",
		action="store",
		default=NA,
		help=paste0(
			"Integer denoting the column number in the file (or data frame) provided with the counts argument,\n",
			"where the HUGO gene symbols are given. If not provided, it will not be available when reporting results.\n",
			"In addition, the 'known' gene filter will not be available for application."
		)
	),	
	make_option(
		opt_str="--btcol",
		type="integer",
		action="store",
		default=NA,
		help=paste0(
			"Integer denoting the column number in the file (or data frame) provided with the counts argument,\n",
			"where the gene biotypes are given. If not provided, the 'biodetection', 'countsbio', 'saturation',\n",
			"'filtered' and 'biodist' plots will not be available."
		)
	),
	make_option(
		opt_str="--annotation",
		action="store",
		default=NULL,
		help=paste0(
			"It can be one of i) NULL (default) to use the existing annotation database or fetch on the fly\n",
			"ii) 'embedded' if the annotation elements are embedded in the read counts file (restrictions apply)\n",
			"iii) a list with a path to a GTF file and certain required metadata."
		)
	),
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
		opt_str="--counttype",
		action="store", 
		help=paste0(
			"The type of reads inside the counts file. It can be one of 'gene', 'exon' or 'utr'\n",
			"for Quant-Seq (Lexogen) protocol. This is a very important and mandatory parameter\n",
			"as it defines the course of the workflow."
		)
	),
	make_option(
		opt_str="--utrOpts_frac",
		action="store",
		default=1,
		help="The fraction (0-1) of the 3' UTR region to count reads in."
	),
	make_option(
		opt_str="--utrOpts_minlen",
		action="store",
		default=300,
		help="The minimum acceptable 3'UTR length irrespective of utrOpts_frac argument."
	),
	make_option(
		opt_str="--utrOpts_dnstrm",
		action="store",
		default=50,
		help=paste0(
			"The number of base pairs to flank the end of the 3' UTR of transcripts\n",
			"when analyzing Quant-Seq data."
		)
	),
	make_option(
		opt_str="--exonfltr",
		action="store_false",
		default=TRUE,
		help=paste0(
			"The supported exon filter in the current version is minActiveExons which\n",
			"implements a filter for demanding m out of n exons of a gene to have a certain\n",
			"read presence with parameters exonfltr_exonsprgene, exonfltr_minexons and exonfltr_frac.\n",
			"The filter is described as follows: if a gene has up to exonfltr_exonsprgene exons,\n",
			"then read presence is required in at least exonfltr_minexons of them, else read presence\n",
			"is required in a exonfltr_frac fraction of the total exons. With the default values, the filter\n",
			"instructs that if a gene has up to 5 exons, read presence is required in at least 2, else\n",
			"in at least 20 exons, in order to be accepted. Set exonfltr=FALSE to NOT apply\n",
			"the minActiveExons filter. To apply your own filter parameters, change arguments\n",
			"exonfltr_exonsprgene, exonfltr_minexons, exonfltr_minexons."  
		)
	),
	make_option(
		opt_str="--exonfltr_exonsprgene",
		action="store",
		default=5,
		help="minActiveExons filter: Exons per gene. Defaults to 5."
	),
	make_option(
		opt_str="--exonfltr_minexons",
		action="store",
		default=2,
		help=paste0(
			"minActiveExons filter: read presence is required in at least exonfltr_minexons of\n",
			"exonfltr_exonsprgene. Defaults to 2."
		)
	),
	make_option(
		opt_str="--exonfltr_frac",
		action="store",
		default=1/5,
		help=paste0(
			"minActiveExons filter: if read presence in at least exonfltr_minexons of exonfltr_exonsprgene\n",
			"not true, read presence is required in a exonfltr_frac fraction of the total exons.\n",
			"Defaults to 1/5."
		)
	),
	make_option(
		opt_str="--genefltr",
		action="store_false",
		default=TRUE,
		help=paste0("Set genefltr=FALSE to NOT apply any gene filtering. To apply your own\n",
			"filters parameters, change arguments:\n",
			"genefltr1_length (default 500) for length filter, where genes are accepted for further analysis\n",    
			"genefltr2_avgperbp (default 100) for average reads filter, where a gene is accepted for further analysis\n",
			"if it has more average reads than the genefltr_avgquantile (default 0.25) of the average count\n",
			"distribution per genefltr_avgperbp base pairs\n",
			"genefltr2_avgquantile (default 0.25)\n",
			"genefltr3_expmedian (default TRUE) for a filter based on the overall expression of a gene.\n",
			"Genes below the median of the overall count distribution are not accepted for further analysis.\n", 
			"genefltr3_expmean (default FALSE) for a filter based on the overall expression of a gene.\n",
			"Genes below the mean of the overall count distribution are not accepted for further analysis.\n",      
			"genefltr3_expquantile (default NA) for a filter based on the overall expression of a gene.\n",
			"Genes below the the specified quantile of the total counts distribution are not accepted\n",
			"for further analysis.\n",  
			"genefltr3_expknown (default NA) for a filter based on the overall expression of a gene.\n",
			"A set of known not-expressed genes in the system under investigation are used to estimate\n",
			"an expression cutoff.\n",   
			"genefltr4_biotype!!! genes with a certain biotype (MUST be contained in the annotation,\n",
			"thus it's better to use annotation='download') are excluded from the analysis.\n",    
			"genefltr5_frac (default 0.25) where a gene is further considered for statistical testing if genefltr5_frac\n",
			"(x100 for a percentage value) have more than genefltr5_mincount reads across all samples\n",
			"(genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE)\n",      
			"genefltr5_mincount (default 10) where a gene is further considered for statistical testing if genefltr5_frac\n",
			"(x100 for a percentage value) have more than genefltr5_mincount reads across all samples\n",
			"(genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE)\n",    
			"genefltr5_percon (default FALSE) where a gene is further considered for statistical testing if genefltr5_frac\n",
			"(x100 for a percentage value) have more than genefltr5_mincount reads across all samples\n",
			"(genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE)"
		)
	),
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
		opt_str="--genefltr2_avgperbp",
		action="store",
		default=100,
		help=paste0(
			"Gene filter2: a gene is accepted for further analysis if it has more average\n",
			"reads than the genefltr_avgquantile of the average count distribution per\n",
			"genefltr_avgperbp base pairs. Defaults to 100."
		)
	),
	make_option(
		opt_str="--genefltr2_avgquantile",
		action="store",
		default=0.25,
		help=paste0(
			"Gene filter2: a gene is accepted for further analysis if it has more average\n",
			"reads than the genefltr_avgquantile of the average count distribution per\n",
			"genefltr_avgperbp base pairs. Defaults to 0.25."
		)
	),
	make_option(
		opt_str="--genefltr3_expmedian",
		action="store_false",
		default=TRUE,
		help=paste0(
			"Gene filter3: based on the overall expression of a gene. Genes below the median\n",
			"of the overall count distribution are not accepted for further analysis. Defaults to TRUE."
		)
	),
	make_option(
		opt_str="--genefltr3_expmean",
		action="store_true",
		default=FALSE,
		help=paste0(
			"Gene filter3: based on the overall expression of a gene. Genes below the mean\n",
			"of the overall count distribution are not accepted for further analysis. Defaults to FALSE."
		)
	),
	make_option(
		opt_str="--genefltr3_expquantile",
		action="store",
		default=NA,
		help=paste0(
			"Gene filter3: based on the overall expression of a gene. Genes below the specified\n",
			"quantile of the total counts distribution are not accepted for further analysis."
		)
	),
	make_option(
		opt_str="--genefltr3_expknown",
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
		opt_str="--genefltr4_biotype",
		action="store",
#		default="getDefaults('biotypeFilter',org[1])",
		help=paste0(
			"Gene filter4: genes with a certain biotype (MUST be contained in the annotation,\n",
			"thus it's better to use annotation='download') are excluded from the analysis.\n",
			"This filter is a named list of logical, where names are the biotypes in each genome\n",
			"and values are TRUE or FALSE. If the biotype should be excluded, the value should be TRUE else FALSE."
		)
	),
	make_option(
		opt_str="--genefltr5_frac",
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
		opt_str="--genefltr5_mincount",
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
		opt_str="--genefltr5_percon",
		action="store_true",
		default=FALSE,
		help=paste0(
			"Gene filter5: a gene is further considered for statistical testing if genefltr5_frac\n",
			"(x100 for a percentage value) have more than genefltr5_mincount reads across all samples\n",
			"(genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE)\n",
			"Defaults to FALSE."
		)
	),		
	make_option(
		opt_str="--whenapplyfilter",
		action="store",
		default="postnorm",
		help=paste0(
			"A character string determining when to apply the exon and/or gene filters,\n",
			"relative to normalization. It can be 'prenorm' to apply the filters and exclude genes\n",
			"from further processing before normalization, or 'postnorm' to apply the filters\n",
			"after normalization (default)."
		)
	),
	make_option(
		opt_str="--normalization",
		action="store",
		default="deseq",
		help=paste0(
			"The normalization algorithm to be applied on the count data.\n",
			"It can be one of 'edaseq' for EDASeq normalization, 'deseq' for the normalization\n",
			"algorithm in the DESq package (default), 'edger' for the normalization algorithms\n",
			"present in the edgeR package 'noiseq' for the normalization algorithms present in\n",
			"the NOISeq package 'nbpseq' for the normalization algorithms present in the NBPSeq\n",
			"package or 'none' to not normalize the data (highly unrecommended). Algorithm specific\n",
			"arguments can be passed through the normargs argument)."
		)
	),
	make_option(
		opt_str="--normargs",
		action="store",
		default=NULL,
		help=paste0(
			"A named list whose names are the names of the normalization algorithm parameters and\n",
			"its members parameter values. Leave NULL for the defaults of normalization."
		)
	),
	make_option(
		opt_str="--statistics",
		action="store",
		default="deseq",
		help=paste0(
			"One or more statistical analyses to be performed by the metaseqr2 pipeline.\n",
			"It can be one or more of 'deseq' (default) to conduct statistical test(s)\n",
			"implemented in the DESeq package, 'edger' to conduct statistical test(s)\n",
			"implemented in the edgeR package, 'limma' to conduct the RNA-Seq version of\n",
			"statistical test(s) implemented in the limma package, 'noiseq' to conduct\n",
			"statistical test(s) implemented in the NOISeq package, 'bayseq' to conduct\n",
			"statistical test(s) implemented in the baySeq package, 'nbpseq' to conduct\n",
			"statistical test(s) implemented in the NBPSeq package, 'deseq2' to conduct\n",
			"statistical test(s) implemented in the DESeq2 package, 'dss' to conduct\n",
			"statistical test(s) implemented in the DSS package and 'absseq' to conduct\n",
			"statistical test(s) implemented in the ABSSeq package."
		)
	),
	make_option(
		opt_str="--statargs",
		action="store",
		default=NULL,
		help=paste0(
			"A named list whose names are the names of the statistical algorithms used in the pipeline.\n",
			"Each member is another named list whose names are the algorithm parameters and its members\n",
			"are the parameter values. Leave NULL for the defaults of statistics."
		)
	),
	make_option(
		opt_str="--adjmethod",
		action="store",
		default="BH",
		help=paste0(
			"The multiple testing p-value adjustment method. It can be one of p.adjust.methods or\n",
			"'qvalue' from the qvalue Bioconductor package. Defaults to 'BH' for Benjamini-Hochberg correction."
		)
	),
	make_option(
		opt_str="--metap",
		action="store",
		default="simes",
		help=paste0(
			"The meta-analysis method to combine p-values from multiple statistical tests.\n",
			"It can be one of 'simes'(default), 'bonferroni', 'minp', 'maxp', 'weight, 'pandora',\n",
			"'dperm_min', 'dperm_max', 'dperm_weight', 'fisher', 'fperm', 'whitlock' or 'none'."
		)
	),
	make_option(
		opt_str="--weight",
		action="store",
#		default="rep(1/length(statistics),length(statistics))",
		help=paste0(
			"A vector of weights with the same length as the statistics vector containing a weight\n",
			"for each statistical test. It should sum to 1."
		)
	),
	make_option(
		opt_str="--nperm",
		action="store",
		default=10000,
		help=paste0(
			"The number of permutations performed to derive the meta p-value when metap='fperm'\n",
			"or metaP='dperm'. It defaults to 10000."
		)
	),
	make_option(
		opt_str="--pcut",
		action="store",
		default=NA,
		help=paste0(
			"A p-value cutoff for exporting differentially genes, default is\n",
			"to export all the non-filtered genes."
		)
	),	
	make_option(
		opt_str="--logoffset",
		action="store",
		default=1,
		help=paste0(
			"An offset to be added to values during logarithmic transformations\n",
			"in order to avoid Infinity (default is 1)."
		)
	),
	make_option(
		opt_str="--poffset",
		action="store",
		default=NULL,
		help=paste0(
			"A value between 0 and 1 to multiply potential zero p-values with for the combination methods\n",
			"including weighting or NULL (default)."
		)
	),
	make_option(
		opt_str="--preset",
		action="store",
		default=NULL,
		help=paste0(
			"An analysis strictness preset. preset can be one of 'all_basic',\n",
			"'all_normal', 'all_full', 'medium_basic', 'medium_normal', 'medium_full',\n",
			"'strict_basic', 'strict_normal' or 'strict_full', each of which control\n",
			"the strictness of the analysis and the amount of data to be exported."
		)
	),
	make_option(
		opt_str="--qcplots",
		action="store",
		default="mds,biodetection,countsbio,saturation,readnoise,filtered,correl,pairwise,boxplot,gcbias,lengthbias,meandiff,meanvar,rnacomp,deheatmap,volcano,biodist,mastat,statvenn,foldvenn,deregulogram",
		help=paste0(
			"A set of diagnostic plots to show/create. It can be one or more of\n",
			"'mds', 'biodetection', 'rnacomp', 'countsbio', 'saturation', 'readnoise',\n",
			"'filtered', 'boxplot', 'gcbias', 'lengthbias', 'meandiff', 'meanvar', 'deheatmap',\n",
			"'volcano', 'mastat', 'biodist', 'statvenn', 'foldvenn'.\n",
			"Enter as comma-separated list OR set as NULL for no diagnostic plots to be created."
		)
	),
	make_option(
		opt_str="--figformat",
		action="store",
		default="png,jpg,tiff,bmp,pdf,ps",
		help=paste0(
			"The format of the output diagnostic plots.\n",
			"It can be one or more of 'png','jpg','tiff','bmp','pdf','ps'.\n",
			"Enter as comma-separated list to select 1 or more of the above."
		)
	),
	make_option(
		opt_str="--outlist",
		action="store_true",
		default=FALSE,
		help=paste0(
			"A logical controlling whether to export a list with the results\n",
			"in the running environment."
		)
	),
	make_option(
		opt_str="--xprtwhere",
		action="store",
		type="character",
		default=NA,
		help="An output directory for the project results."
	),
	make_option(
		opt_str="--xprtwhat",
		action="store",
		default="annotation,p_value,adj_p_value,meta_p_value,adj_meta_p_value,fold_change,stats,counts,flags",
		help=paste0(
			"The content of the final lists.\n",
			"It can be one or more of 'annotation',to bind the annotation elements for each gene,\n",
			"'p_value', to bind the p-values of each method,\n",
			"'adj_p_value', to bind the multiple testing adjusted p-values,\n",
			"'meta_p_value', to bind the combined p-value from the meta-analysis,\n",
			"'adj_meta_p_value', to bind the corrected combined p-value from the meta-analysis,\n",
			"'fold_change', to bind the fold changes of each requested contrast,\n",
			"'stats', to bind several statistics calclulated on raw and normalized counts(xprtstats argument),\n",
			"'counts', to bind the raw and normalized counts for each sample.\n",
			"Enter as comma-separated list to select 1 or more of the above."
		)
	),
	make_option(
		opt_str="--xprtscale",
		action="store",
		default="natural,log2,log10,vst,rpgm",
		help=paste0(
			"Export values from one or more transformations applied to the data. It can be one or more of\n",
			"'natural', 'log2', 'log10', 'vst' (Variance Stabilizing Transormation, documentation in DESeq package)\n",
			"and 'rpgm' which is ratio of mapped reads per gene model (either the gene length or the sum of\n",
			"exon lengths, depending on countType argument). Note that this is not RPKM as reads are already\n",
			"normalized for library size using one of the supported normalization methods. Also, 'rpgm' might\n",
			"be misleading when normalization is other than 'deseq'.\n",
			"Enter as comma-separated list to select 1 or more of the above."
		)
	),
	make_option(
		opt_str="--xprtvalues",
		action="store",
		default="raw,normalized",
		help=paste0(
			"It can be one or more of 'raw' to export raw values (counts etc.) and 'normalized' to export normalized counts.\n",
			"Enter as comma-separated list to select 1 or more of the above."
		)
	),
	make_option(
		opt_str="--xprtstats",
		action="store",
		default="mean,median,sd,mad,cv,rcv",
		help=paste0(
			"Calculate and export several statistics on raw and normalized counts, condition-wise.\n",
			"It can be one or more of 'mean', 'median', 'sd', 'mad', 'cv' for the Coefficient of Variation,\n",
			"'rcv' for a robust version of CV where the median and the MAD are used instead of the mean and\n",
			"the standard deviation. Enter as comma-separated list to select 1 or more of the above."
		)
	),
	make_option(
		opt_str="--xprtcountstbl",
		action="store_true",
		default=FALSE,
		help=paste0(
			"Exports the calculated read counts table when input is read from bam files and\n",
			"the normalized count table in all cases. Defaults to FALSE."
		)
	),
	make_option(
		opt_str="--rc",
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
		opt_str="--report",
		action="store_false",
		default=TRUE,
		help=paste0(
			"A logical value controlling whether to produce a summary report or not.\n",
			"Defaults to TRUE."
		)
	),
	make_option(
		opt_str="--topreport",
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
		opt_str="--templatereport",
		action="store",
#		default="default",
		help=paste0(
			"An HTML template to use for the report.\n",
			"Do not change this unless you know what you are doing."
		)
	),
	make_option(
		opt_str="--genemodel",
		action="store_false",
		default=TRUE,
		help=paste0(
			"In case of exon analysis, a list with exon counts for each gene will be saved\n",
			"to the file xprtwhere/data/gene_model.RData. This file can be used as input\n",
			"to metaseqR for exon count based analysis, in order to avoid the time consuming\n",
			"step of assembling the counts for each gene from its exons"
		)
	),
	make_option(
		opt_str="--verbose",
		action="store_false",
		default=TRUE,
		help="Print informative messages during execution? Defaults to TRUE."
	),
	make_option(
		opt_str="--runlog",
		action="store_false",
		default=TRUE,
		help=paste0(
			"Write a log file of the metaseqr2 run using package log4r. Defaults to TRUE.\n",
			"The filename will be auto-generated."
		)
	),
	make_option(
		opt_str="--reportdb",
		action="store",
		default="sqlite",
		help=paste0(
			"Database system to use for storing the report intereactive graphs.\n",
			"Can be 'sqlite' (default) or 'dexie'."
		)
	),
	make_option(
		opt_str="--localdb",
		action="store",
#		default="file.path(system.file(package='metaseqR2'),'annotation.sqlite')",
		help="The metaseqR2 annotation database location."
	),
	make_option(
		opt_str="--progressfun",
		action="store",
		default=NULL,
		help=paste0(
			"A function which updates a Progress object from shiny. This function\n",
			"must accept a detail argument. See http://shiny.rstudio.com/articles/progress.html"
		)
	),
	make_option(
		opt_str="--offlinereport",
		action="store_false",
		default=TRUE,
		help=paste0(
			"TRUE (default) to download and include the required JavaScript libraries\n",
			"to properly view the report offline. Ignored if report=FALSE"
		)
	),
	make_option(
		opt_str="--exportr2c",
		action="store_false",
		default=TRUE,
		help=paste0(
			"EXPORTS READ TO COUNTS\n",
			"Defaults to TRUE."
		)
	)
)

opt <- parse_args(OptionParser(option_list=option_list))

#Create targets file
samplenames.v <- unlist(strsplit(opt$samples, split=","))
filenames.v <- unlist(strsplit(opt$files, split=","))
conditions.v <- unlist(strsplit(opt$conditions, split=","))

if (!is.null(opt$paired)){
  paired.v <- unlist(strsplit(opt$paired, split=","))
}else{paired.v <- rep(NA,times=length(samplenames.v))}
if (!is.null(opt$strandp)){
  stranded.v <- unlist(strsplit(opt$strandp, split=","))
}else{stranded.v <- rep(NA,times=length(samplenames.v))}

targets <- data.frame(samplename=samplenames.v,
  filename=filenames.v,
  condition=conditions.v,
  paired=paired.v,
  stranded=stranded.v)
  
targets <- t(na.omit(t(targets)))

write.table(targets,file=file.path(opt$path,"targets.txt"),sep="\t",row.names=FALSE,quote=FALSE)
#write.table(targets,file="targetsT.txt",sep="\t",row.names=FALSE,quote=FALSE)

# Check if version is numeric or "auto"
if (opt$version!="auto"){
	versionNo <- as.numeric(opt$version)
}else{versionNo <- "auto"} 

# Set EXON & GENE Filters
exonFilters <- list()
geneFilters <- list()

if (isTRUE(opt$exonfltr)){
			exonFilters=list(
				minActiveExons=list(
					exonsPerGene=opt$exonfltr_exonsprgene,
					minExons=opt$exonfltr_minexons,
					frac=opt$exonfltr_frac))
}else{assign("exonFilters", NULL, envir = .GlobalEnv)}

if (isTRUE(opt$genefltr)){
			geneFilters=list(
				length=list(
					length=opt$genefltr1_length),
			avgReads=list(
				averagePerBp=opt$genefltr2_avgperbp,
				quantile=opt$genefltr2_avgquantile),
			expression=list(
				median=opt$genefltr3_expmedian,
				mean=opt$genefltr3_expmean,
				quantile=opt$genefltr3_expquantile,
				known=opt$genefltr3_expknown),
#			biotype=opt$genefltr4_biotype,
			biotype=getDefaults("biotypeFilter",opt$org),
			presence=list(
				frac=opt$genefltr5_frac,
				minCount=opt$genefltr5_mincount,
				perCondition=opt$genefltr5_percon))
}else{assign("geneFilters", NULL, envir = .GlobalEnv)}

# Statistics -> VECTOR
statistics.v <- unlist(strsplit(opt$statistics, split=","))

# Default WEiGHT calculation
opt$weight <- rep(1/length(statistics.v),length(statistics.v))

# QC-plots to be exported -> VECTOR
qcplots.v <- vector()
#if (is.null(opt$qcplots)){
if (opt$qcplots == "NULL"){
	qcplots.v <- NULL
}else{qcplots.v <- unlist(strsplit(opt$qcplots, split=","))}

# Figure Formats -> VECTOR
figform.v <- unlist(strsplit(opt$figformat, split=","))

# Outdir
#if  (!is.empty(opt$xprtwhere)){
#	path <- unlist(strsplit(opt$xprtwhere, split="/"))
#xprtwhere <- file.path(tempdir(),opt$xprtwhere)
#}

# Create "export" VECTORS...
xprtwhat.v <- unlist(strsplit(opt$xprtwhat, split=","))
xprtscale.v <- unlist(strsplit(opt$xprtscale, split=","))
xprtvalues.v <- unlist(strsplit(opt$xprtvalues, split=","))
xprtstats.v <- unlist(strsplit(opt$xprtstats, split=","))

opt$localdb <- file.path(system.file(package='metaseqR2'),'annotation.sqlite')
opt$templatereport <- "default"
			
# TODO: more checks
if (!(opt$org %in% c("hg18", "hg19", "hg38", "mm9","mm10", "rn5", "rn6", "dm3", "dm6", "danrer7","pantro4", "susscr3", "tair10", "equcab2" )))
	stop("The organism must be one of \"hg18\", \"hg19\", \"hg38\", \"mm9\",\"mm10\", \"rn5\", \"rn6\", \"dm3\", \"dm6\", \"danrer7\",\"pantro4\", \"susscr3\", \"tair10\", \"equcab2\"!")

metaseqr2(
#    counts=opt$counts,
#    sampleList=opt$samplelist,
	sampleList=file.path(opt$path,"targets.txt"),
    excludeList=opt$excludelist,
    path=opt$path,
    fileType=opt$filetype,
    contrast=opt$contrast,
    libsizeList=opt$libsizelist,
    embedCols=list(
        idCol=opt$idcol,
        gcCol=opt$gccol,
        nameCol=opt$namecol,
        btCol=opt$btcol
    ),
    annotation=opt$annotation,
    org=opt$org,
    refdb=opt$refdb,
    version=versionNo,
    transLevel=opt$translevel,
    countType=opt$counttype,
    utrOpts=list(
       frac=opt$utrOpts_frac,
       minLength=opt$utrOpts_minlen,
       downstream=opt$utrOpts_dnstrm
    ),
    exonFilters=exonFilters,
    geneFilters=geneFilters,
    whenApplyFilter=opt$whenapplyfilter,
    normalization=opt$normalization,
    normArgs=opt$normargs,
    statistics=statistics.v,
    statArgs=opt$statargs,
    adjustMethod=opt$adjmethod,
    metaP=opt$metap,
    weight=opt$weight,
    nperm=opt$nperm,
    pcut=opt$pcut,
    logOffset=opt$logoffset,
    pOffset=opt$poffset,
    preset=opt$preset,
    qcPlots=qcplots.v,
    figFormat=figform.v,
    outList=opt$outlist,
    exportWhere=opt$xprtwhere,
    exportWhat=xprtwhat.v,
    exportScale=xprtscale.v,
    exportValues=xprtvalues.v,
    exportStats=xprtstats.v,
    exportCountsTable=opt$xprtcountstbl,
    restrictCores=opt$rc,
    report=opt$report,
    reportTop=opt$topreport,
    reportTemplate=opt$templatereport,
    saveGeneModel=opt$genemodel,
    verbose=opt$verbose,
    runLog=opt$runlog,
    reportDb=opt$reportdb,
    localDb=opt$localdb,
    .progressFun=opt$progressfun,
    offlineReport=opt$offlinereport,
    .exportR2C=opt$exportr2c
)