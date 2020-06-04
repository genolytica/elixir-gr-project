#! /usr/local/bin/Rscript

# Wrapper for the metaseqr2() function

suppressPackageStartupMessages(library("optparse"))
library(metaseqR2)

option_list <- list(
	make_option(
		opt_str="--samplelist",
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
#	make_option(
# 		opt_str="--version",
#		action="store",
#		type="numeric",
#		default="auto",
#		help=paste0(
#			"An integer denoting the version of the annotation to use from the local annotation\n",
#			"database or fetch on the fly. For Ensembl, it corresponds to Ensembl releases, while\n",
#			"for UCSC/RefSeq, it is the date of creation (locally)."
#		)
#	),
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
		action="store",
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
		action="store",
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
		action="store",
		default=TRUE,
		help=paste0(
			"Gene filter3: based on the overall expression of a gene. Genes below the median\n",
			"of the overall count distribution are not accepted for further analysis. Defaults to TRUE."
		)
	),
	make_option(
		opt_str="--genefltr3_expmean",
		action="store",
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
#	make_option(
#		opt_str="--genefltr4_biotype",
#		action="store",
#		default="getDefaults('biotype.filter',org[1])",
#		help=paste0(
#			"Gene filter4: genes with a certain biotype (MUST be contained in the annotation,\n",
#			"thus it's better to use annotation='download') are excluded from the analysis.\n",
#			"This filter is a named list of logical, where names are the biotypes in each genome\n",
#			"and values are TRUE or FALSE. If the biotype should be excluded, the value should be TRUE else FALSE."
#		)
#	),
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
		action="store",
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
#	make_option(
#		opt_str="--qcplots",
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
		opt_str="--figformat",
		action="store",
		help=paste0(
			"The format of the output diagnostic plots. It can be one or more of 'png', 'jpg',\n",
			"'tiff', 'bmp', 'pdf', 'ps'. The native format 'x11' (for direct display) is not provided\n",
			"as an option as it may not render the proper display of some diagnostic plots in some devices."
		)
	),
	make_option(
		opt_str="--outlist",
		action="store",
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
		opt_str="--rc",
		action="store",
		default=0.6,
		help=paste0(
			"In case of parallel execution of several subfunctions, the fraction of the\n",
			"available cores to use. In some cases if all available cores are used (rc=1)\n",
			"and the system does not have sufficient RAM, the pipeline running machine\n",
			"might significantly slow down."
		)
	)
)

opt <- parse_args(OptionParser(option_list=option_list))

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
			#biotype=opt$genefltr4_biotype,
			presence=list(
				frac=opt$genefltr5_frac,
				minCount=opt$genefltr5_mincount,
				perCondition=opt$genefltr5_percon))
}else{assign("geneFilters", NULL, envir = .GlobalEnv)}
			
# TODO: more checks
if (!(opt$org %in% c("hg18", "hg19", "hg38", "mm9","mm10", "rn5", "rn6", "dm3", "dm6", "danrer7","pantro4", "susscr3", "tair10", "equcab2" )))
	stop("The organism must be one of \"hg18\", \"hg19\", \"hg38\", \"mm9\",\"mm10\", \"rn5\", \"rn6\", \"dm3\", \"dm6\", \"danrer7\",\"pantro4\", \"susscr3\", \"tair10\", \"equcab2\"!")


metaseqr2(
    sampleList=opt$samplelist,
    excludeList=opt$excludelist,
    fileType=opt$filetype,
    path=opt$path,
    contrast=opt$contrast,
    libsizeList=opt$libsizelist,
    org=opt$org,
    refdb=opt$refdb,
    #version=opt$version,
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
    statistics=opt$statistics,
    #qcPlots=opt$qcplots,
    figFormat=opt$figformat,
    outList=opt$outlist,
    exportWhere=opt$xprtwhere,
    exportWhat=c("annotation","p_value","adj_p_value","fold_change","stats","counts","flags"),
    exportScale="natural",
    exportValues=c("raw","normalized"),
    exportStats="mean",
    restrictCores=opt$rc,
    .exportR2C=TRUE
)