# This R script calculates ...
#
# The script performs the following tasks:
#
# Usage:
#
# Parameters:
#
# In all cases, the program will run in parallel if multiple cores exist in the 
# running machine and the R package "parallel" is present and rc is not NULL. 
# Otherwise, a single core will be used (e.g. in a Windows machine where the
# parallel package is not supported).

metaseqR2Counts <- function(
	targets,
#	sampleList,
    excludeList=NULL,
#    fileType=c("auto","sam","bam","bed"),
    path=NULL,
    contrast=NULL,
    libsizeList=NULL,
#    embedCols=list(
#        idCol=4,
#        gcCol=NA,
#        nameCol=NA,
#        btCol=NA
#    ),
#    annotation=NULL,
    org=c("hg18","hg19","hg38","mm9","mm10","rn5","rn6","dm3","dm6",
        "danrer7","pantro4","susscr3","tair10","equcab2"),
    refdb=c("ensembl","ucsc","refseq"),
    version="auto",
    transLevel=c("gene","transcript","exon"),
    countType=c("gene","exon","utr"),
    utrOpts=list(
        frac=1,
        minLength=300,
        downstream=50
    ),
    exonFilters=list(
        minActiveExons=list(
            exonsPerGene=5,
            minExons=2,
            frac=1/5
        )
    ),
    geneFilters=list(
        length=list(
            length=500
        ),
        avgReads=list(
            averagePerBp=100,
            quantile=0.25
        ),
        expression=list(
            median=TRUE,
            mean=FALSE,
            quantile=NA,
            known=NA,
            custom=NA
        ),
        biotype=getDefaults("biotypeFilter",org[1]),
        presence=list(
            frac=0.25,
            minCount=10,
            perCondition=FALSE
        )
    ),
    whenApplyFilter=c("postnorm","prenorm"),
#    normalization=c("deseq","deseq2","edaseq","edger","noiseq","nbpseq",
#        "absseq","dss","each","none"),
#    normArgs=NULL,
#    statistics=c("deseq","deseq2","edger","noiseq","bayseq","limma","nbpseq",
#        "absseq","dss"),
#    statArgs=NULL,
#    adjustMethod=sort(c(p.adjust.methods,"qvalue")),
#    metaP=if (length(statistics)>1) c("simes","bonferroni","fisher",
#        "dperm_min","dperm_max","dperm_weight","fperm","whitlock","minp","maxp",
#        "weight","pandora","none") else "none",
#    weight=rep(1/length(statistics),length(statistics)),
#    nperm=10000,
#    pcut=NA,
#    logOffset=1,
#   preset=NULL, # An analysis strictness preset
#    qcPlots=c(
#        "mds","biodetection","countsbio","saturation","readnoise","filtered",
#        "correl","pairwise","boxplot","gcbias","lengthbias","meandiff",
#        "meanvar","rnacomp","deheatmap","volcano","biodist","mastat","statvenn",
#        "foldvenn","deregulogram"
#    ),
    figFormat=c("png","jpg","tiff","bmp","pdf","ps"),
    outList=FALSE,
    exportWhere=NA, # An output directory for the project
    exportWhat=c("annotation","p_value","adj_p_value","meta_p_value",
        "adj_meta_p_value","fold_change","stats","counts","flags"),
    exportScale=c("natural","log2","log10","vst","rpgm"),
    exportValues=c("raw","normalized"),
    exportStats=c("mean","median","sd","mad","cv","rcv"),
    exportCountsTable=FALSE,
    restrictCores=0.6,
    report=TRUE,
    reportTop=0.1,
    reportTemplate="default",
    saveGeneModel=TRUE,
    verbose=TRUE,
    runLog=TRUE,
    reportDb=c("dexie","sqlite"),
    localDb=file.path(system.file(package="metaseqR2"),"annotation.sqlite"),
    progressFun=NULL,
    offlineReport=TRUE,
	.exportR2C=TRUE){
	
	if (!require(metaseqR2))
		stop("Bioconductor package metaseqR2 is required!")
	
	annotation <- NULL
	fileType <- "auto"
	normalization <- "deseq"
	statistics <- "edger"
	
#	assign("VERBOSE",verbose,envir=metaseqR2:::metaEnv)	
#	sampleList <- file.path(getwd(),targets)
	sampleList <- targets
#	sampleList<- file.path(getwd(),"targets.txt")	
#	setwd(dirname(getwd()))
#	assign("VERBOSE",TRUE,envir=metaseqR:::meta.env)
#	assign("VERBOSE",verbose,envir=metaseqR2:::metaEnv)	
#	sampleList <- file.path(getwd(),targets)
#	sampleList <- targets
#	message("Sample List1 ",typeof(targets))
#	message("Sample List1 ",file.exists(targets))
#	message("Sample List1 ",View(targets))
#	sampleList <- targets
#	sampleList <- file.path(getwd(),targets)
#	message("Type ",typeof(sampleList))
#	message("Exists ",file.exists(sampleList))
#	message("Missing ",missing(sampleList))
#	message("Is List ",is.list(sampleList))
#	message("View ",View(sampleList))
#	message("...",View(read.delim(sampleList)))
#	message("Sample List2 ",View(sampleList))
#	targets <- file.path(getwd(),"targets.txt")	
#	sampleList <- read.delim(targets)		

	
	metaseqR2::metaseqr2(sampleList,excludeList,fileType,path,contrast,libsizeList,
	embedCols,annotation,org,refdb,version,transLevel,countType,utrOpts,exonFilters,
	geneFilters,whenApplyFilter,normalization,normArgs,statistics,statArgs,adjustMethod,
	metaP,weight,nperm,pcut,logOffset,preset,qcPlots,figFormat,outList,exportWhere,exportWhat,
	exportScale,exportValues,exportStats,exportCountsTable,restrictCores,report,reportTop,
	reportTemplate,saveGeneModel,verbose,runLog,reportDb,localDb,progressFun,offlineReport,
	.exportR2C)
}