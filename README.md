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
$ mkdir -p $WORK/datasets/GSE79183/fastq
$ mkdir -p $WORK/datasets/GSE79183/hisat2_out
$ mkdir -p $WORK/datasets/GSE79183/metaseqR_out
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
    echo "Downloading $ACC..."
    $WORK/tools/sratoolkit.2.9.6/bin/fastq-dump \
        --outdir  $WORK/datasets/GSE79183/fastq \
        --accession $ACC -v
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

## Generation of UCSC Genome Browser BigWig tracks

After finishing with the alignment to the reference genome, we must visualize
the RNA signal and also normalize it across samples. To this end, there is a
Perl script already developed ```normalize_bedgraph.pl```. See its ```--help```
argument for more details. Briefly, it accepts a list of BedGraph files and 
outputs a set of normalized BedGraph files. The normalization is being performed
by linearly up- or down-scaling the total BedGraph signal to ```1e+9``` (this is
the default value but it can be customized). The script also exports the factors
with which each signal is multiplied. We need these values for usage in SeqCVIBE
so they should be exported and retained.

Prior to using this script, the BAM files must be converted to BedGraph files.
This can be easily achieved:

```
mkdir -p $WORK/datasets/GSE79183/bigwig
for DIR in $WORK/datasets/GSE79183/hisat2_out/*
do
    SAMPLE=`basename $DIR`
    $WORK/tools/bedtools-2.28.0/bin/bedtools genomecov -split -bg \
        -ibam $DIR/$SAMPLE.bam | \
        grep -vP 'chrM|chrU|chrG|chrJ|rand|hap|loc|cox' | sort -k1,1 -k2g,2 > \
        $WORK/datasets/GSE79183/bigwig/$SAMPLE".bedGraph" &
done
```

It is also useful to get a "genome size" file (that is simply chromosomal 
lengths) which is useful later for the conversion to BigWig:

```
# Just one sample will provide what we need
SAMPLES=($(ls $WORK/datasets/GSE79183/hisat2_out/))
SAMPLE=${SAMPLES[0]}
samtools view -H $WORK/datasets/GSE79183/hisat2_out/$SAMPLE/$SAMPLE.bam | \
    grep '^@SQ' | sed s/@SQ\\tSN:// | \
    grep -vP 'chrM|chrU|chrG|chrJ|rand|hap|loc|cox' | sed s/\\tLN:/\\t/ > \
    $WORK/datasets/GSE79183/bigwig/genome.size
```

Then the normalization part:

```
NCORES=8
CURDIR=`pwd`
cd $WORK/datasets/GSE79183/bigwig
perl $WORK/scripts/normalize_bedgraph.pl \
    --input $WORK/datasets/GSE79183/bigwig/*.bedGraph \
    --exportfactors $WORK/datasets/GSE79183/bigwig/GSE79183_factors.txt \
    --ncores $NCORES
cd $CURDIR
```

After the execution finishes, the normalized BedGraph files 
(```*_norm.bedGraph```) must be converted to BigWig format.

```
GENOME=$WORK/datasets/GSE79183/bigwig/genome.size
for FILE in `ls $WORK/datasets/GSE79183/bigwig/*_norm.bedGraph`
do
    SAMPLE=`basename $FILE | sed s/_norm\.bedGraph//`
    $WORK/tools/kenttools-1.04.0/bedGraphToBigWig $FILE $GENOME \
        $WORK/datasets/GSE79183/bigwig/$SAMPLE".bigWig" &
done
```

Finally, remove all the BedGraph files as they are not needed any longer:

```
rm $WORK/datasets/GSE79183/bigwig/*.bedGraph
```

## Quantification with metaseqR

After the alignment step and signal track generation step for visualization, 
the metaseqR pipeline must be run to produce the data that will be the basis of 
the functionality of the SeqCVIBE version that will power the Elixir project 
deliverable. This is a straightforward proceudre, see the 
```$WORK/tools/metaseqR-template.R``` script.

Although metaseqR is primarily a differential expression analysis tool, we will
use is to generate gene expression quantifications at this point.

### Construction of the targets file

metaseqR is quite flexible regarding its inputs. One way that can be relatively
automated (to some extent) is to construct a targets file to be used with the
```read.targets``` function. The targets file contents can be automates
regarding the sample name and sample path. The condition name, whether it's a
single or paired-end library and whether a stranded protocol has been used.
Pairing can maybe be inferred from the ```LAYOUT``` field in SRA metadata using
some ```curl``` call or API (to be researched).

The basis of the targets file can be constructed as follows:

```
printf "%s\t%s\t%s\t%s\t%s\n" \
    "samplename" "filename" "condition" "paired" "stranded" > \
    $WORK/datasets/GSE79183/metaseqR_out/targets.txt
for DIR in $WORK/datasets/GSE79183/hisat2_out/*
do
    SAMPLE=`basename $DIR`
    LINK=`readlink -f $DIR/$SAMPLE".bam"`
    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$SAMPLE" "$LINK" "condition_name" "single_or_paired" "strand_mode" >> \
        $WORK/datasets/GSE79183/metaseqR_out/targets.txt
done
```

Then you must fill manually (at this first stage) the ```condition```, 
```paired``` and ```stranded``` columns according to metaseqR's instruction. For
this particular example, ```paired``` is ```paired```, ```stranded``` is 
```forward``` and condition is
```
WT
WT
WT
KO
KO
KO
```

### Reading and normalizing data with metaseqR

Because metaseqR is designed as a differential gene expression analysis
pipeline, the main ```metaseqr``` function is not useful. Instead, we are going
to use the functions ```read.targets```, ```read2count``` and ```norm.deseq```
to produce the data we need to construct an ```.rda``` (or ```.RData```, as long
as we are consistent) file to be used by SeqCVIBE. The next steps are supposed
to take place within R console. We suppose that ```work``` corresponds to the
```$WORK``` variable used so far.

1. Read the targets file

```
library(metaseqR)

# work <- "$WORK"

targets.file <- file.path(work,"datasets","GSE79183","metaseqR_out",
    "targets.txt")
targets <- read.targets(targets)
```

2. Do the counting

```
# Set up the awkward messaging and parallel environment, is better in the
# working version!
multic <- check.parallel(0.5)
assign("VERBOSE",TRUE,envir=metaseqR:::meta.env)

# Retrieve annotation, this example is mm10 and we are retrieving "exon"
# annotation for counting as per standard RNA-Seq protocols. Also, we are
# working with Ensembl annotations.
annotation <- get.annotation("mm10","exon")

# Retrieve also gene-based annotations used later
gene.data <- get.annotation("mm10","gene")

# The "annotation" data.frame can be saved to a text file and then for each
# mouse case, be read and re-used (until this is also launched in the new
# metaseqR version). Same for other organisms.

# Count
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

# Based on the sum of their exon lengths
gene.length <- the.exon.lengths[[1]]
names(gene.length) <- rownames(gene.data)

```

3. Normalize read counts with DESeq

```
norm.counts <- normalize.deseq(gene.counts,targets$samples,output="matrix")
```

4. Export the required data

```
b2c.out <- list(
    counts=gene.counts,
    norm=norm.counts,
    libsize=libsize,
    length=gene.length
)
save(b2c.out,file=file.path(work,"datasets","GSE79183","metaseqR_out",
    "GSE79183.rda")
```

## Directory structure to host BAM and signal files

_TBD_

### TODO

- CWL workflow of the described procedures to be used for new data analyses
- Rscript command line version to be used with CWL



