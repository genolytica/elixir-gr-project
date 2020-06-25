#! /usr/local/bin/Rscript

# Wrapper for the createSignalTracks() function

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
	make_option(
		opt_str="--path",
		action="store",
#		default=NULL,
		help="Directory where input files are located"
	),
#	make_option(
#		opt_str="--targets",
#		action="store",
#		help="A tab-delimited file with the experimental description"
#	),
	make_option(
		opt_str="--org",
		action="store",
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
		opt_str="--urlbase",
		action="store",
		default=NULL,
		help="A valid URL which is prepended to the created bigWig files"
	),
	make_option(
		opt_str="--stranded",
		action="store_true",
		default=FALSE,
		help=paste0(
			"Separate + and - strands and create separate bigWig files.\n",
			"In the case of stranded tracks, a UCSC Genome Brower trackhub is created.\n", 
			"Individual tracks can be retrieved from the trackhub."
		)
	),
	make_option(
		opt_str="--normto",
		action="store",
		default=1e+9,
		help=paste0(
			" Large integer, denoting the total sum of signal to be used as the normalization target.\n",
			"It defaults to 1e+9. This means that if for a particular sample the sum of signal is 1.5e+9\n",
			"(sum(sapply(coverage(x),sum)) == 1.5e+9) then this is linearly scaled to 1e+9."
		)
	),
	make_option(
		opt_str="--exportpath",
		action="store",
		default= ".",
		help=" Path to export tracks."
	),
	make_option(
		opt_str="--hubinfo_name",
		action="store",
		help=paste0(
			"Name of the track hub created in case of stranded tracks.\n",
			"Please see the track hub specifications at the UCSC Genome Browser site."
		)
	),
	make_option(
		opt_str="--hubinfo_sl",
		action="store",
		help=paste0(
			"Short label for the track hub created in case of stranded tracks.\n",
			"Please see the track hub specifications at the UCSC Genome Browser site."
		)
	),
	make_option(
		opt_str="--hubinfo_ll",
		action="store",
		help=paste0(
			"Long label for the track hub created in case of stranded tracks.\n",
			"Please see the track hub specifications at the UCSC Genome Browser site."
		)
	),
	make_option(
		opt_str="--hubinfo_email",
		action="store",
		help=paste0(
			"Email associated with the track hub created in case of stranded tracks.\n",
			"Please see the track hub specifications at the UCSC Genome Browser site."
		)
	),
	make_option(
		opt_str="--overwrite",
		action="store_true",
		default=FALSE,
		help="Overwrite tracks if they exist? Defaults to FALSE."
	),
	make_option(
		opt_str="--rc",
		action="store",
		default=NULL,
		help="Fraction of cores to use"
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

# TODO: more checks
if (!(opt$org %in% c("hg18", "hg19", "hg38", "mm9","mm10", "rn5", "rn6", "dm3", "dm6", "danrer7","pantro4", "susscr3", "tair10", "equcab2" )))
	stop("The organism must be one of \"hg18\", \"hg19\", \"hg38\", \"mm9\",\"mm10\", \"rn5\", \"rn6\", \"dm3\", \"dm6\", \"danrer7\",\"pantro4\", \"susscr3\", \"tair10\", \"equcab2\"!")
if (!is.numeric(opt$normto) || opt$normto<0)
	stop("Fragment length must be a positive large integer!")

	
createSignalTracks(targets=file.path(opt$path,"targets.txt"),org=opt$org,urlBase=opt$urlbase,stranded=opt$stranded,normTo=opt$normto,
	exportPath=opt$exportpath,hubInfo=list(name=opt$hubinfo_name,shortLabel=opt$hubinfo_sl,longLabel=opt$hubinfo_ll,
	email=opt$hubinfo_email),overwrite=opt$overwrite,rc=opt$rc)
