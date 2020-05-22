# This R script creates bigWig files to be used for exploring RNA signal in genome browsers. 
# When strands are separated, a UCSC genome browser trackhub is created to group tracks
# for the same sample. A link to the created data is returned.
#
# The script performs the following tasks:
# - Reads the targets file in order to read in the BAM files and make the Rles.
# - Checks if bigWig files already exist (--overwrite option).
# - Acquires sequence information from the first BAM file. Chromosome information
#	is obtained through the scanBamHeader function of the Rsamtools package.
# - Gets coverage using Rsamtools package and assigns sequence information for bigwig
#	while creating a list of GRanges objects (one for each sample).
# - Calculates normalization parameters and performs normalization.
# - For unstranded mode, the bigWig files for every sample and the track lines file
#	(tracks.txt) are exported in --exportpath option or working directory by default.
# - For stranded mode, a hub is created in which the bigWig files for positive
#	and negative strand for each sample are exported, including the generated genomes,
#	hub (entry point) and trackDb files as well as the hub link.
#
# Usage:
# tracks <- createSignalTracks (
#						targets,
#						org,
#						urlBase=NULL,
#						stranded=FALSE,
#						normTo=1e+9,
#						exportPath=".",
#						hubInfo=list(name="MyHub",shortLabel="My hub", longLabel="My hub",email="someone@example.com"),
#						overwrite=FALSE,
#						rc=NULL
#				)					 
#
# Parameters:
# targets	: a tab-delimited file which contains the sample names, the BAM/BED file names and
#			  the biological conditions/groups for each sample/file. The file should be text
#			  tab-delimited and structured as follows: the first line of the external tab-delimited	
#			  file should contain column names (names are not important). The first column MUST contain
#			  UNIQUE sample names. The second column MUST contain the raw BAM/BED files WITH their full
#			  path. If the files in the second column of the targets file do not contain a path to a directory,
#			  the current directory is assumed to be the BAM/BED file container. The third column MUST contain
#			  the biological condition where each of the samples in the first column should belong to. 
# org		: for human genomes "hg18", "hg19" or "hg38",
#			  for mouse genomes "mm9", "mm10",
#			  for rat genomes "rn5" or "rn6",
#			  for drosophila genome "dm3" or "dm6",
#			  for zebrafish genome "danrer7", "danrer10" or "danrer11",
#			  for chimpanzee genome "pantro4", "pantro5",
#			  for pig genome "susscr3", "susscr11",
#			  for Arabidopsis thaliana genome "tair10" and
#			  for Equus caballus genome "equcab2"
# urlBase	: a valid URL which is prepended to the created bigWig files produced.
#			  The base path of the bigDataUrl in UCSC Genome Browser track lines.
# stranded	: can be TRUE or FALSE depending on whether you wish to create stranded tracks by
#			  separating + and - strand reads. In the case of stranded tracks, a UCSC Genome Brower
#			  trackhub is created. Individual tracks can be retrieved from the trackhub.
# normTo	: a large integer, denoting the total sum of signal to be used as the normalization target.
#			  It defaults to 1e+9. This means that if for a particular sample the sum of signal is 1.5e+9
#			  (sum(sapply(coverage(x),sum)) == 1.5e+9) then this is linearly scaled to 1e+9.
# exportPath: path to export tracks.
# hubInfo	: information regarding the track hub created when stranded=TRUE.
#			  A list with the track hub description in case of stranded tracks.
#			  Please see the track hub specifications at the UCSC Genome Browser site.
# overwrite	: overwrite tracks if they exist? Defaults to FALSE.
# rc		: fraction of cores to use.
#
# In all cases, the program will run in parallel if multiple cores exist in the 
# running machine and the R package "parallel" is present and rc is not NULL. 
# Otherwise, a single core will be used (e.g. in a Windows machine where the
# parallel package is not supported).

createSignalTracks <- function(targets,org,urlBase=NULL,stranded=FALSE,
    normTo=1e+9,exportPath=".",hubInfo=list(name="MyHub",shortLabel="My hub",
    longLabel="My hub",email="someone@example.com"),overwrite=FALSE,rc=NULL){
    if (!requireNamespace("rtracklayer"))
        stopwrap("Bioconductor package rtracklayer is required to build ",
            "tracks!")
    
    if (!is.list(targets)) {
        if (file.exists(targets))
            targets <- readTargets(targets)
        else
            stopwrap("targets must be a valid existing text file!")
    }
	
	if (!require(parallel))
        stop("R-core package parallel is required!")
	if (!require(rtracklayer))
        stop("R-core package rtracklayer is required!")
    if (!require(Rsamtools))
        stop("Bioconductor package Rsamtools is required!")
    if (!require(GenomeInfoDb))
        stop("Bioconductor package GenomeInfoDb is required!")
	if (!require(GenomicAlignments))
        stop("Bioconductor package GenomicAlignments is required!")
	if (!require(GenomicRanges))
        stop("Bioconductor package GenomicAlignments is required!")
    if (!require(BiocManager))
        stop("The Bioconductor package BiocManager is required!")
    
    # Read in BAMs and make Rles
    message("Creating tracks...")
    if (stranded) {
        message("  stranded mode")
        return(.createStrandedSignalTracks(targets,org,normTo,urlBase,
            exportPath,hubInfo,overwrite,rc=rc))
    }
    else {
        message("  unstranded mode")
        return(.createUnstrandedSignalTracks(targets,org,normTo,urlBase,
            exportPath,overwrite,rc=rc))
    }
}

.createStrandedSignalTracks <- function(targets,org,normTo,urlBase,
    exportPath,hubInfo,overwrite,rc=NULL) {
    # First check if bigWig files already exist
    trackDbPath <- file.path(exportPath, getUcscOrganism(org))
    posCheck <- file.path(trackDbPath,
        paste0(unlist(targets$samples,use.names=FALSE),"_plus.bigWig"))
    negCheck <- file.path(trackDbPath,
        paste0(unlist(targets$samples,use.names=FALSE),"_minus.bigWig"))
    posEx <- file.exists(posCheck)
    negEx <- file.exists(negCheck)
    if (all(posEx) && all(negEx) && !overwrite) {
        message("Requested tracks have already been created! ",
            "Use overwrite=TRUE to re-create them.")
        return("")
    }
    
    # Positive and negative color options
    posBaseColours <- .getPosBaseColors()
    negBaseColours <- .getNegBaseColors()
    posCol <- rep(posBaseColours,length.out=length(targets$samples))
    posCol <- rep(posCol,lengths(targets$samples))
    negCol <- rep(negBaseColours,length.out=length(targets$samples))
    negCol <- rep(negCol,lengths(targets$samples))
    names(posCol) <- names(negCol) <- unlist(targets$samples,use.names=FALSE)
    
    # The BAM files
    bams <- unlist(targets$files,use.names=FALSE)
    
    # Get seqinfo form first BAM
    preSf <- .chromInfoFromBAM(bams[1])
    vchrs <- getValidChrs(org)
    preSf <- preSf[intersect(vchrs,rownames(preSf)),,drop=FALSE]
    sf <- .chromInfoToSeqInfoDf(preSf,o=org,asSeqinfo=TRUE)
    
    # Get coverage and assign seqinfo for bigwig
    message("Reading positive strand reads from BAM files to Rle...")
    pbg <- cmclapply(bams,function(b,v,s) {
        message("  reading ",b)
        reads <- trim(unlist(grglist(readGAlignments(file=b,
            param=ScanBamParam(scanBamFlag(isMinusStrand=FALSE))))))
        cov <- coverage(reads)
        seqs <- as.character(seqlevels(reads))
        vv <- intersect(v,seqs)
        cov <- cov[vv]
        gr <- as(slice(cov,1),"GRanges")
        seqinfo(gr) <- s
        return(gr)
    },vchrs,sf,rc=rc)
    names(pbg) <- names(posCol)
    
    message("Reading negative strand reads from BAM files to Rle...")
    nbg <- cmclapply(bams,function(b,v,s) {
        message("  reading ",b)
        reads <- trim(unlist(grglist(readGAlignments(file=b,
            param=ScanBamParam(scanBamFlag(isMinusStrand=TRUE))))))
        cov <- coverage(reads)
        seqs <- as.character(seqlevels(reads))
        vv <- intersect(v,seqs)
        cov <- cov[vv]
        gr <- as(slice(cov,1),"GRanges")
        # Inverse the coverage
        gr$score <- -gr$score
        seqinfo(gr) <- s
        return(gr)
    },vchrs,sf,rc=rc)
    names(nbg) <- names(negCol)
    
    # Calculate normalization factors
    rawPosSums <- sapply(pbg,function(x) sum(x$score))
    rawNegSums <- sapply(nbg,function(x) sum(x$score))
    
    rat <- rawPosSums/-rawNegSums
    posNormTo <- rat*0.5*normTo
    negNormTo <- normTo - posNormTo
    
    posNormFacs <- 0.5*posNormTo/rawPosSums
    negNormFacs <- -0.5*negNormTo/rawNegSums
    
    names(posNormFacs) <- names(pbg)
    names(negNormFacs) <- names(nbg)
    
    # Normalize (should be quick)
    message("Normalizing positive strand...")
    npbg <- cmclapply(names(pbg),function(n,B,N) {
        B[[n]]$score <- B[[n]]$score*N[n]
        return(B[[n]])
    },pbg,posNormFacs)
    names(npbg) <- names(pbg)
    
    message("Normalizing negative strand...")
    nnbg <- cmclapply(names(nbg),function(n,B,N) {
        B[[n]]$score <- B[[n]]$score*N[n]
        return(B[[n]])
    },nbg,negNormFacs)
    names(nnbg) <- names(nbg)
    
    # Start creating the hub
    # First the dir structure
    #trackDbPath <- file.path(exportPath,metaseqR2:::getUcscOrganism(org))
    if (!dir.exists(trackDbPath))
        dir.create(trackDbPath,recursive=TRUE)
    
    # Put the bigWig files there
    message("Exporting bigWig files...")
    posBwFiles <- file.path(trackDbPath,paste0(names(npbg),"_plus.bigWig"))
    names(posBwFiles) <- names(npbg)
    for (n in names(npbg)) {
        message("  exporting + strand for sample ",n," to ",posBwFiles[n])
        export.bw(npbg[[n]],posBwFiles[n])
    }
    negBwFiles <- file.path(trackDbPath,paste0(names(nnbg),"_minus.bigWig"))
    names(negBwFiles) <- names(nnbg)
    for (n in names(nnbg)) {
        message("  exporting - strand for sample ",n," to ",negBwFiles[n])
        export.bw(nnbg[[n]],negBwFiles[n])
    }
    
    # Write the genomes file
    gf <- .makeTrackhubGenomesFile(org)
    writeLines(paste(names(gf)," ",gf,sep=""),
        file.path(exportPath,"genomes.txt"))
    
    # Write the hub file
    hub <- .makeTrackhubEntrypoint(hubInfo)
    writeLines(paste(names(hub)," ",hub,sep=""),
        file.path(exportPath,"hub.txt"))
    
    # Write the trackDb file...
    pretdb <- .makeTrackhubMultiwig(names(posBwFiles),posBwFiles,negBwFiles,
        posCol,negCol,urlBase,org)
    # Clear previous file if exists, otherwise append will cause problems
    if (file.exists(file.path(trackDbPath,"trackDb.txt")))
        unlink(file.path(trackDbPath,"trackDb.txt"))
    for (p in pretdb) {
        tracks <- p$tracks
        p$tracks <- NULL
        meta <- paste(paste(names(p)," ",p,sep=""),collapse="\n")
        post <- paste(names(tracks$positive)," ",tracks$positive,sep="")
        post <- paste(paste0("    ",post),sep="",collapse="\n")
        negt <- paste(names(tracks$negative)," ",tracks$negative,sep="")
        negt <- paste(paste0("    ",negt),sep="",collapse="\n")
        aTrack <- c(meta,"\n\n",post,"\n\n",negt,"\n\n")
        cat(aTrack,file=file.path(trackDbPath,"trackDb.txt"),sep="",append=TRUE)
    }
    
    # Generate the hub link
    return(paste0(urlBase,"/hub.txt"))
}

.createUnstrandedSignalTracks <- function(targets,org,normTo,urlBase,
    exportPath,overwrite,rc=NULL) {
    # First check if bigWig files already exist
    bCheck <- file.path(exportPath,paste0(unlist(targets$samples,
        use.names=FALSE),".bigWig"))
    bEx <- file.exists(bCheck)
    if (all(bEx) && !overwrite) {
        message("Requested tracks have already been created! ",
            "Use overwrite=TRUE to re-create them.")
        return("")
    }
    
    # Create color scheme
    posBaseColours <- .getPosBaseColors()
    posCol <- rep(posBaseColours,length.out=length(targets$samples))
    posCol <- rep(posCol,lengths(targets$samples))
    names(posCol) <- unlist(targets$samples,use.names=FALSE)
    
    # The BAM files
    bams <- unlist(targets$files,use.names=FALSE)
    
    # Get seqinfo form first BAM
    preSf <- .chromInfoFromBAM(bams[1])
    vchrs <- getValidChrs(org)
    preSf <- preSf[intersect(vchrs,rownames(preSf)),,drop=FALSE]
    sf <- .chromInfoToSeqInfoDf(preSf,o=org,asSeqinfo=TRUE)
    
    # Get standed coverage and assign seqinfo for bigwig
    message("Reading BAM files to Rle...")
    bg <- cmclapply(bams,function(b,v,s) {
        message("  reading ",b)
        reads <- trim(unlist(grglist(readGAlignments(file=b))))
        cov <- coverage(reads)
        seqs <- as.character(seqlevels(reads))
        vv <- intersect(v,seqs)
        cov <- cov[vv]
        gr <- as(slice(cov,1),"GRanges")
        seqinfo(gr) <- s
        return(gr)
    },vchrs,sf,rc=rc)
    names(bg) <- names(posCol)
    
    # Calculate normalization factors
    rawSums <- sapply(bg,function(x) sum(x$score))
    normFacs <- normTo/rawSums
    names(normFacs) <- names(bg)
    
    # Normalize (should be quick)
    message("Normalizing...")
    nbg <- cmclapply(names(bg),function(n,B,N) {
        B[[n]]$score <- B[[n]]$score*N[n]
        return(B[[n]])
    },bg,normFacs)
    names(nbg) <- names(bg)
    
    message("Exporting bigWig files...")
    bwFiles <- file.path(exportPath,paste0(names(nbg),".bigWig"))
    names(bwFiles) <- names(nbg)
    for (n in names(nbg)) {
        message("  exporting for sample ",n," to ",bwFiles[n])
        export.bw(nbg[[n]],bwFiles[n])
    }
    
    message("Creating track lines...")
    opts <- .makeUnstrandedTrackOpts(posCol)
    trackLines <- character(length(bwFiles))
    for (i in 1:length(bwFiles)) {
        opts[[i]]$name <- sub("^([^.]*).*","\\1",basename(bwFiles[i]))
        opts[[i]]$name <- gsub(" ","_",opts[[i]]$name)
        opts[[i]]$bigDataUrl <- paste0(urlBase,"/",basename(bwFiles[i]))
        popts <- paste(names(opts[[i]]),"=",opts[[i]],sep="")
        trackLines[i] <- paste(popts,collapse=" ")
        trackLines[i] <- paste("track",trackLines[i])
    }
    
    tracksFile <- file.path(exportPath,"tracks.txt")
    writeLines(trackLines,tracksFile)
    tracksLink <- paste0(urlBase,"/tracks.txt")
    return(tracksLink)
}

.chromInfoFromBAM <- function(bam) {
    # Danger of including non-canonical chromosomes
    b <- scanBamHeader(bam)
    ci <- as.data.frame(b[[bam]]$targets)
    names(ci) <- "length"
    return(ci)
}

.chromInfoToSeqInfoDf <- function(ci,o="custom",circ=FALSE,asSeqinfo=FALSE) {
    if (asSeqinfo)
        return(Seqinfo(seqnames=rownames(ci),seqlengths=ci[,1],
            isCircular=rep(circ,nrow(ci)),genome=o))
    else
        return(data.frame(chromosome=rownames(ci),length=as.integer(ci[,1])))
}

.makeUnstrandedTrackOpts <- function(cols) {
    opts <- vector("list",length(cols))
    for (i in 1:length(cols)) {
        opts[[i]]$type <- "bigWig"
        opts[[i]]$color <- paste(t(col2rgb(cols[i])),collapse=",")
        opts[[i]]$visibility <- "full"
        opts[[i]]$maxHeightPixels <- "128:64:16"
    }
    return(opts)
}

.makeTrackhubMultiwig <- function(tnames,pfiles,nfiles,pcols,ncols,
    url,org) {
    opts <- vector("list",length(pcols))
    for (i in 1:length(pfiles)) {
        opts[[i]]$track <- paste0(tnames[i],"_stranded")
        opts[[i]]$type <- "bigWig"
        opts[[i]]$container <- "multiWig"
        opts[[i]]$aggregate <- "transparentOverlay"
        opts[[i]]$showSubtrackColorOnUi <- "on"
        opts[[i]]$shortLabel <- tnames[i]
        opts[[i]]$longLabel <- paste0(tnames[i]," signal")
        opts[[i]]$boxedCfg <- "on"
        opts[[i]]$autoScale <- "on"
        opts[[i]]$visibility <- "full"
        opts[[i]]$maxHeightPixels <- "128:64:16"
        opts[[i]]$tracks <- list(
            positive=list(
                track=paste0(tnames[i],"_plus"),
                parent=paste0(tnames[i],"_stranded"),
                type="bigWig",
                color=paste(t(col2rgb(pcols[i])),collapse=","),
                bigDataUrl=paste0(url,"/",getUcscOrganism(org),"/",
                    paste0(tnames[i],"_plus.bigWig"))
            ),
            negative=list(
                track=paste0(tnames[i],"_minus"),
                parent=paste0(tnames[i],"_stranded"),
                type="bigWig",
                color=paste(t(col2rgb(ncols[i])),collapse=","),
                bigDataUrl=paste0(url,"/",getUcscOrganism(org),"/",
                    paste0(tnames[i],"_minus.bigWig"))
            )
        )
    }
    return(opts)
}

.makeTrackhubGenomesFile <- function(org) {
    return(list(
        genome=getUcscOrganism(org),
        trackDb=paste0(org,"/trackDb.txt")
    ))
}

.makeTrackhubEntrypoint <- function(h) {
    return(list(
        hub=h$name,
        shortLabel=h$shortLabel,
        longLabel=h$longLabel,
        genomesFile="genomes.txt",
        email=h$email
    ))
}

.getPosBaseColors <- function() {
     return(c("#B40000","#00B400","#0000B4","#B45200","#9B59B6","#21BCBF",
        "#BC4800","#135C34","#838F00","#4900B5"))
}

.getNegBaseColors <- function() {
     return(c("#FF7575","#79FF79","#8484FF","#FFB77C","#EBB7FF","#63E5E7",
        "#FFB88B","#5BFFA5","#F4FF75","#C69FFF"))
}

getUcscOrganism <- function(org) {
    switch(org,
        hg18 = { return("hg18") },
        hg19 = { return("hg19") },
        hg38 = { return("hg38") },
        mm9 = { return("mm9") },
        mm10 = { return("mm10") },
        rn5 = { return("rn5") },
        rn6 = { return("rn6") },
        dm3 = { return("dm3") },
        dm6 = { return("dm6") },
        danrer7 = { return("danRer7") },
        danrer10 = { return("danRer10") },
        danrer11 = { return("danRer11") },
        pantro4 = { return("panTro4") },
        pantro5 = { return("panTro5") },
        susscr3 = { return("susScr3") },
        susscr11 = { return("susScr11") },
        equcab2 = { return("equCab2") },
        tair10 = { return("TAIR10") }
    )
}

getValidChrs <- function(org) {
    switch(org,
        hg18 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        hg19 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        hg38 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        mm9 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY"
            ))
        },
        mm10 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY"
            ))
        },
        rn5 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX"
            ))
        },
        rn6 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX"
            ))
        },
        dm3 = {
            return(c(
                "chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet",
                "chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet",
                "chrYHet"
            ))
        },
        dm6 = {
            return(c(
                "chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet",
                "chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet",
                "chrYHet"
            ))
        },
        danrer7 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        danrer10 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        danrer11 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        pantro4 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        pantro5 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        pantro6 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        susscr3 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr2","chr3","chr4","chr5","chr6","chr7",
                "chr8","chr9","chrX","chrY"
            ))
        },
        susscr11 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr2","chr3","chr4","chr5","chr6","chr7",
                "chr8","chr9","chrX","chrY"
            ))
        },
        equcab2 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr26","chr27","chr28","chr29","chr3","chr30",
                "chr31","chr4","chr5","chr6","chr7","chr8","chr9","chrX"#,"chrY"
            ))
        },
        tair10 = {
            return(c(
                "chr1","chr2","chr3","chr4","chr5"
            ))
        }
    )
}

cmclapply <- function(...,rc) {
    if (suppressWarnings(!requireNamespace("parallel")) 
        || .Platform$OS.type!="unix")
        m <- FALSE
    else {
        m <- TRUE
        ncores <- parallel::detectCores()
        if (ncores==1) 
            m <- FALSE
        else {
            if (!missing(rc) && !is.null(rc))
                ncores <- ceiling(rc*ncores)
            else 
                m <- FALSE
        }
    }
	if (m)
        return(mclapply(...,mc.cores=ncores,mc.set.seed=FALSE))
    else
        return(lapply(...))
}

readTargets <- function(input,path=NULL) {
    if (missing(input) || !file.exists(input))
        stopwrap("The targets file should be a valid existing text file!")
    tab <- read.delim(input,strip.white=TRUE)
    samples <- as.character(tab[,1])
    conditions <- unique(as.character(tab[,3]))
    rawfiles <- as.character(tab[,2])
    if (!is.null(path)) {
        tmp <- dirname(rawfiles) # Test if there is already a path
        if (any(tmp=="."))
            rawfiles <- file.path(path,basename(rawfiles))
    }
    # Check if they exist!!!
    for (f in rawfiles) {
        if (!file.exists(f))
            stopwrap("Raw reads input file ",f," does not exist! Please check!")
    }
    if (length(samples) != length(unique(samples)))
        stopwrap("Sample names must be unique for each sample!")
    if (length(rawfiles) != length(unique(rawfiles)))
        stopwrap("File names must be unique for each sample!")
    sampleList <- vector("list",length(conditions))
    names(sampleList) <- conditions
    for (n in conditions)
        sampleList[[n]] <- samples[which(as.character(tab[,3])==n)]
    fileList <- vector("list",length(conditions))
    names(fileList) <- conditions
    for (n in conditions) {
        fileList[[n]] <- rawfiles[which(as.character(tab[,3])==n)]
        names(fileList[[n]]) <- samples[which(as.character(tab[,3])==n)]
    }
    if (ncol(tab)>3) { # Has info about single- or paired-end reads / strand
        if (ncol(tab)==4) { # Stranded or paired
            whats <- tolower(as.character(tab[,4]))
            if (!all(whats %in% c("yes","no","forward","reverse",
                "single","paired")))
                stopwrap("Unknown options for paired-end reads and/or ",
                    "strandedness in targets file.")
            what <- whats[1]
            if (what %in% c("single","paired")) {
                hasPairedInfo <- TRUE
                hasStrandedInfo <- FALSE
            }
            else {
                if (what %in% c("yes","no")) {
                    .deprecatedWarning("readTargets")
                    tmp <- as.character(tab[,4])
                    tmp[tmp=="yes"] <- "forward"
                    tab[,4] <- tmp
                    hasPairedInfo <- FALSE
                    hasStrandedInfo <- TRUE
                }
                if (what %in% c("forward","reverse","no")) {
                    hasPairedInfo <- FALSE
                    hasStrandedInfo <- TRUE
                }
            }
        }
        if (ncol(tab)==5) { # Both
            whatsPaired <- tolower(as.character(tab[,4]))
            if (!all(whatsPaired %in% c("single","paired","mixed")))
                stopwrap("Unknown option for type of reads (single, paired, ",
                    "mixed) in targets file.")
            whatsStrand <- tolower(as.character(tab[,5]))
            if (!all(whatsStrand %in% c("yes","no","forward","reverse")))
                stopwrap("Unknown option for read strandedness in targets file")
            if (any(whatsStrand=="yes")) {
                .deprecatedWarning("readTargets")
                tmp <- as.character(tab[,5])
                tmp[tmp=="yes"] <- "forward"
                tab[,5] <- tmp
            }
            hasPairedInfo <- TRUE
            hasStrandedInfo <- TRUE
        }
        if (hasPairedInfo && !hasStrandedInfo) {
            pairedList <- vector("list",length(conditions))
            names(pairedList) <- conditions
            for (n in conditions) {
                pairedList[[n]] <- character(length(sampleList[[n]]))
                names(pairedList[[n]]) <- sampleList[[n]]
                for (nn in names(pairedList[[n]]))
                    pairedList[[n]][nn] <- as.character(tab[which(as.character(
                        tab[,1])==nn),4])
            }
        }
        else
            pairedList <- NULL
        if (hasStrandedInfo && !hasPairedInfo) {
            strandedList <- vector("list",length(conditions))
            names(strandedList) <- conditions
            for (n in conditions) {
                strandedList[[n]] <- character(length(sampleList[[n]]))
                names(strandedList[[n]]) <- sampleList[[n]]
                for (nn in names(strandedList[[n]]))
                    strandedList[[n]][nn] <- 
                        as.character(tab[which(as.character(tab[,1])==nn),4])
            }
        }
        else
            strandedList <- NULL
        if (hasStrandedInfo && hasPairedInfo) {
            strandedList <- vector("list",length(conditions))
            names(strandedList) <- conditions
            for (n in conditions) {
                strandedList[[n]] <- character(length(sampleList[[n]]))
                names(strandedList[[n]]) <- sampleList[[n]]
                for (nn in names(strandedList[[n]]))
                    strandedList[[n]][nn] <- as.character(tab[which(
                        as.character(tab[,1])==nn),5])
            }
            pairedList <- vector("list",length(conditions))
            names(pairedList) <- conditions
            for (n in conditions) {
                pairedList[[n]] <- character(length(sampleList[[n]]))
                names(pairedList[[n]]) <- sampleList[[n]]
                for (nn in names(pairedList[[n]]))
                    pairedList[[n]][nn] <- as.character(tab[which(as.character(
                        tab[,1])==nn),4])
            }
        }
    }
    else
        pairedList <- strandedList <- NULL
    # Guess file type based on only one of them
    tmp <- fileList[[1]][1]
    if (length(grep("\\.bam$",tmp,ignore.case=TRUE,perl=TRUE))>0)
        type <- "bam"
    else if (length(grep("\\.sam$",tmp,ignore.case=TRUE,perl=TRUE))>0)
        type <- "sam"
    else if (length(grep("\\.bed$",tmp,ignore.case=TRUE,perl=TRUE)>0))
        type <- "bed"
    else
        type <- NULL
    return(list(samples=sampleList,files=fileList,paired=pairedList,
        stranded=strandedList,type=type))
}

.deprecatedWarning <- function(func) {
    switch(func,
        readTargets = {
            warnwrap("\"yes\" and \"no\" for read strandedness have been ",
                "deprecated. Please use \"forward\", \"forward\" or \"no\". ",
                "Replacing \"yes\" with \"forward\"...")
        }
    )
}