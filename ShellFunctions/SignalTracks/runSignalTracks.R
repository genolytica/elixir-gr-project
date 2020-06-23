#! /usr/local/bin/Rscript

# Wrapper for the createSignalTracks() function

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(metaseqR2))

option_list <- list(
	make_option(
		opt_str="--path",
		action="store",
		default=NULL,
		help="Directory where input files are located"
	),
	make_option(
		opt_str="--targets",
		action="store",
		help="A tab-delimited file with the experimental description"
	),
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
#print(opt$path)
#print(file.path(getwd(),opt$path))
#print(dirname(opt$targets))
#print(opt$path)
#print(list.files(opt$path))
#print(list.files(dirname(opt$targets)))
targets <- readTargets(opt$targets, path=opt$path)

# TODO: more checks
if (!(opt$org %in% c("hg18", "hg19", "hg38", "mm9","mm10", "rn5", "rn6", "dm3", "dm6", "danrer7","pantro4", "susscr3", "tair10", "equcab2" )))
	stop("The organism must be one of \"hg18\", \"hg19\", \"hg38\", \"mm9\",\"mm10\", \"rn5\", \"rn6\", \"dm3\", \"dm6\", \"danrer7\",\"pantro4\", \"susscr3\", \"tair10\", \"equcab2\"!")
if (!is.numeric(opt$normto) || opt$normto<0)
	stop("Fragment length must be a positive large integer!")
	
createSignalTracks(targets=targets,org=opt$org,urlBase=opt$urlbase,stranded=opt$stranded,normTo=opt$normto,
	exportPath=opt$exportpath,hubInfo=list(name=opt$hubinfo_name,shortLabel=opt$hubinfo_sl,longLabel=opt$hubinfo_ll,
	email=opt$hubinfo_email),overwrite=opt$overwrite,rc=opt$rc)
