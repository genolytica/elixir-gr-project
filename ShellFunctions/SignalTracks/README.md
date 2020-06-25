# generate-signal-tracks

A shell script for creating signal tracks.

## Description

This shell script ```runSignalTracks.R```, after creating an appropriate targets file,
calls the createSignalTracks() function of the metaseqR2 package and creates bigWig files
to be used for exploring RNA signal in genome browsers. When strands are separated, a UCSC genome browser
trackhub is created to group tracks for the same sample. A link to the created data is returned.
The script can be invoked from the system shell instead of the R shell, using ```Rscript```.

## Prerequisite R/Bioconductor packages
The following packages are required to run the script:
- optparse
- metaseqR2

## Basic Example
```
Rscript runSignalTracks.R \
	--samples=sample1,sample2,sample3,sample4 \
	--files=sample1.bam,sample2.bam,sample3.bam,sample4.bam \
	--conditions=control,control,treatment,treatment \
	--org=hg19 \
	--overwrite

# OR

Rscript runSignalTracks.R \
	--samples=sample1,sample2 \
	--files=sample1.bam,sample2.bam \
	--conditions=control,treatment \
	--strandp=forward,forward \
	--org=hg19 \
	--stranded \
	--hubinfo_name="MyHub" \
	--hubinfo_sl="MyHub" \
	--hubinfo_ll="MyHub" \
	--hubinfo_email="someone@example.com"
```

## List of arguments

The following table presents the input arguments in detail:

|Parameter   |Description                                                                                                                                          |
|------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|
|samples     |Sample IDs - Enter as comma-separated list (no space).|
|files       |File names - Enter as comma-separated list (no space).|
|condition   |Sample Conditions - Enter as comma-separated list (no space).|
|paired      |Should be "single" for single-end reads, "paired" for paired-end reads or "mixed" for BAMs that contain both single- and paired-end reads. If this column is not provided, single-end reads will be assumed. Enter as comma-separated list (no space).|
|strandp     |Should be "forward" for a forward (5'->3') strand library construction protocol, "reverse" for a reverse (3'->5') strand library construction protocol, or "no" for unstranded/unknown protocol. If this column is not provided, unstranded reads will be assumed. Enter as comma-separated list (no space).|
|path        |Directory where input files are located and where the targets.txt file will be created.|
|targets	 |A tab-delimited file which contains the sample names, the BAM/BED file names and the biological conditions/groups for each sample/file. The file should be text tab-delimited and structured as follows: the first line of the external tab-delimited	 file should contain column names (names are not important). The first column MUST contain UNIQUE sample names. The second column MUST contain the raw BAM/BED files WITH their full path. If the files in the second column of the targets file do not contain a path to a directory, the current directory is assumed to be the BAM/BED file container. The third column MUST contain the biological condition where each of the samples in the first column should belong to.|
|org		 |For human genomes "hg18", "hg19" or "hg38", for mouse genomes "mm9", "mm10", for rat genomes "rn5" or "rn6", for drosophila genome "dm3" or "dm6", for zebrafish genome "danrer7", "danrer10" or "danrer11", for chimpanzee genome "pantro4", "pantro5", for pig genome "susscr3", "susscr11", for Arabidopsis thaliana genome "tair10" and for Equus caballus genome "equcab2".|
|urlBase	 |A valid URL which is prepended to the created bigWig files produced. The base path of the bigDataUrl in UCSC Genome Browser track lines.|
|stranded	 |By default its false. Use the flag if you wish to create stranded tracks by separating + and - strand reads. In the case of stranded tracks, a UCSC Genome Brower trackhub is created. Individual tracks can be retrieved from the trackhub.|
|normTo	     |A large integer, denoting the total sum of signal to be used as the normalization target. It defaults to 1e+9. This means that if for a particular sample the sum of signal is 1.5e+9 (sum(sapply(coverage(x),sum)) == 1.5e+9) then this is linearly scaled to 1e+9.|
|exportPath  |Path to export tracks.|                                                                        
|hubInfo	 |Information regarding the track hub created when stranded=TRUE. A list with the track hub description in case of stranded tracks. Please see the track hub specifications at the UCSC Genome Browser site.|
|overwrite	 |Overwrite tracks if they exist? Defaults to false.|
|rc		     |Fraction of cores to use.|

## Indicative runtimes