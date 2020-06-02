# generate-signal-tracks

A shell script for creating signal tracks.

## Description

This shell script ```runSignalTracks.R```, calls the createSignalTracks() function
of the metaseqR2 package and creates bigWig files to be used for exploring RNA signal
in genome browsers. When strands are separated, a UCSC genome browser trackhub is
created to group tracks for the same sample. A link to the created data is returned.
The script can be invoked from the system shell instead of the R shell, using ```Rscript```.

The script performs the following tasks:
1. Reads the targets file in order to read in the BAM files and make the Rles.
2. Checks if bigWig files already exist (--overwrite option).
3. Acquires sequence information from the first BAM file. Chromosome information
   is obtained through the scanBamHeader function of the Rsamtools package.
4. Gets coverage using Rsamtools package and assigns sequence information for bigwig
   while creating a list of GRanges objects (one for each sample).
5. Calculates normalization parameters and performs normalization.
6. For unstranded mode, the bigWig files for every sample and the track lines file
   (tracks.txt) are exported in --exportpath option or working directory by default.
7. For stranded mode, a hub is created in which the bigWig files for positive
   and negative strand for each sample are exported, including the generated genomes,
   hub (entry point) and trackDb files as well as the hub link.

## Prerequisite R/Bioconductor packages
The following packages are required to run the script:
- metaseqR2

## Basic Example
```
Rscript runSignalTracks.R \
	--targets=targets.txt \
	--org=hg19 \
	--overwrite=TRUE

# OR

Rscript runSignalTracks.R \
	--targets=targets.txt \
	--org=hg19 \
	--stranded=TRUE \
	--hubinfo_name="MyHub" \
	--hubinfo_sl="MyHub" \
	--hubinfo_ll="MyHub" \
	--hubinfo_email="someone@example.com"
```

## List of arguments

The following table presents the input arguments in detail:

|Parameter   |Description                                                                                                                                          |
|------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|
|targets	 |A tab-delimited file which contains the sample names, the BAM/BED file names and the biological conditions/groups for each sample/file. The file should be text tab-delimited and structured as follows: the first line of the external tab-delimited	 file should contain column names (names are not important). The first column MUST contain UNIQUE sample names. The second column MUST contain the raw BAM/BED files WITH their full path. If the files in the second column of the targets file do not contain a path to a directory, the current directory is assumed to be the BAM/BED file container. The third column MUST contain the biological condition where each of the samples in the first column should belong to.|
|org		 |For human genomes "hg18", "hg19" or "hg38", for mouse genomes "mm9", "mm10", for rat genomes "rn5" or "rn6", for drosophila genome "dm3" or "dm6", for zebrafish genome "danrer7", "danrer10" or "danrer11", for chimpanzee genome "pantro4", "pantro5", for pig genome "susscr3", "susscr11", for Arabidopsis thaliana genome "tair10" and for Equus caballus genome "equcab2"|
|urlBase	 |A valid URL which is prepended to the created bigWig files produced. The base path of the bigDataUrl in UCSC Genome Browser track lines.|
|stranded	 |Can be TRUE or FALSE depending on whether you wish to create stranded tracks by separating + and - strand reads. In the case of stranded tracks, a UCSC Genome Brower trackhub is created. Individual tracks can be retrieved from the trackhub.|
|normTo	     |A large integer, denoting the total sum of signal to be used as the normalization target. It defaults to 1e+9. This means that if for a particular sample the sum of signal is 1.5e+9 (sum(sapply(coverage(x),sum)) == 1.5e+9) then this is linearly scaled to 1e+9.|
|exportPath  |Path to export tracks.|                                                                        
|hubInfo	 |Information regarding the track hub created when stranded=TRUE. A list with the track hub description in case of stranded tracks. Please see the track hub specifications at the UCSC Genome Browser site.|
|overwrite	 |Overwrite tracks if they exist? Defaults to FALSE.|
|rc		     |Fraction of cores to use.|

## Indicative runtimes