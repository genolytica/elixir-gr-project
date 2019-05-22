library(metaseqR)

args = commandArgs(trailingOnly=TRUE)

work <- "/media/galadriel/hybridstat/elixir_project/"
gse <-args[1]
genome <- args[2]
cores <- as.numeric(args[3])

targets.file <- file.path(work,"datasets",gse,"metaseqR_out", "targets.txt")
targets <- read.targets(targets.file)

#SET CORES FRACTION TO USE
multic <- check.parallel(cores)
assign("VERBOSE",TRUE,envir=metaseqR:::meta.env)

annotation <- get.annotation(genome,"exon")
gene.data <- get.annotation(genome,"gene")
count.obj <- read2count(targets,annotation,file.type=targets$type,multic=multic)
exon.counts <- count.obj$counts
exon.data <- count.obj$mergedann
libsize <- count.obj$libsize
exon.counts <- cbind(exon.data[rownames(exon.counts),c("start","end","exon_id",
    "gene_id")],exon.counts[,unlist(targets$samples,use.names=FALSE)])
the.counts <- construct.gene.model(exon.counts,targets$samples,gene.data,
    multic=multic)
the.gene.counts <- the.exon.lengths <- 
    vector("list",length(unlist(targets$samples)))

names(the.gene.counts) <- names(the.exon.lengths) <- names(the.counts)
for (n in names(the.gene.counts)) {
    the.gene.counts[[n]] <- wapply(multic,the.counts[[n]],
        function(x) return(sum(x$count)))
    the.exon.lengths[[n]] <- wapply(multic,the.counts[[n]],
        function(x) return(sum(x$length)))
    the.gene.counts[[n]] <- do.call("c",the.gene.counts[[n]])
    the.exon.lengths[[n]] <- do.call("c",the.exon.lengths[[n]])
}
gene.counts <- do.call("cbind",the.gene.counts)

gene.length <- the.exon.lengths[[1]]
names(gene.length) <- rownames(gene.data)

norm.counts <- normalize.deseq(gene.counts,targets$samples,output="matrix")

b2c.out <- list(
    counts=gene.counts,
    norm=norm.counts,
    libsize=libsize,
    length=gene.length
)
save(b2c.out,file=file.path(work,"datasets",gse,"metaseqR_out",
    paste0(gse,".rda")))
