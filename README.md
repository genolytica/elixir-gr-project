# Dataset setup

## Retrieve data from SRA/ENA/EGA

1. From GEO/ENA/EGA accession number, we find the SRA links. For example, for
   the dataset GSE79183, the link between GEO accession numbers and SRA ids is:

|GEO id     |SRA id     |
|-----------|-----------|
|GSM2087387 |SRR3224187 |
|GSM2087388 |SRR3224188	|
|GSM2087389 |SRR3224189	|
|GSM2087390 |SRR3224190 |
|GSM2087392 |SRR3224191 |
|GSM2087393 |SRR3224192 |

2. In the working directory ```$WORK``` which can be found in our Slack channel 
   and in the datasets sub-directory, create a directory named after the 
   GEO/ENA/EGA accession and several sub-directories that will host the data
   to be downloaded and the various formats.

```
$ mkdir -p $WORK/GSE79183/fastq
$ mkdir -p $WORK/GSE79183/hisat2_out
$ mkdir -p $WORK/GSE79183/metaseqR_out
# ...
```

3. Use the ```fastq-dump``` tool from the SRA Toolkit (or EGA/ENA corresponding
   tools or download links) to get the FASTQ files. If the dataset description
   mentions a paired-end sequencing protocol, the ```--split-files``` option
   in ```fastq-dump``` must be used. The latest SRA toolking is located in 
   ```$WORK/tools/sratoolkit.2.9.6```. Sometimes, things may go a bit unexpected.
   This [link](https://www.biostars.org/p/222122/) may be of help. For example, 
   for the above dataset:

```
for ACC in SRR3224187 SRR3224188 SRR3224189 SRR3224190 SRR3224191 SRR3224192
do
    $WORK/tools/sratoolkit.2.9.6/bin/fastq-dump \
        --outdir  $WORK/GSE79183/fastq \
        --accession $ACC
done
```

## Alignment to the reference genome

The output of step 3 above should be a set of FASTQ files (single or pairs). 
Now, we must use these files to align them to a reference genome. For the 
particular dataset we are examining, this is the GRCh37 human genome version. 
We will use the HISAT2 splice-aware aligner and the indexes recommended in the
HISAST2 [website](https://ccb.jhu.edu/software/hisat2/manual.shtm). A first
reference and index for the above dataset is located at ```$WORK/reference```.
Apart from the reference genome, for RNA-Seq we must have a set of known splice
sites which can guide HISAT2 for the provision of better results. A template for
the HISAT2 pipeline (initially for paired-end reads) is implemented in the file
```$WORK/scripts/map_hisat2_advanced_paired-end.sh```. Changes may be required
for the script to run (e.g. for single-end reads). After all internal variables
are set,  ```sh $WORK_DIR/scripts/map_hisat2_advanced_paired-end.sh``` should be
enough to generate the BAM files required for further analysis. The output files
will be in ```$WORK/GSE79183/hisat2_out```.

## Quantification and differential expression analysis with metaseqR

After the alignment step, the metaseqR pipeline must be run to produce the data
that will be the basic of the functionality of the SeqCVIBE version that will
power the Elixir project deliverable. This is a straightforward proceudre, see
the ```$WORK/tools/metaseqR-template.R``` script.

## Generation of UCSC Genome Browser tracks

_TBD_

