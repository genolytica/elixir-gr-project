#!/bin/bash

# USAGE: elixir_pipeline.sh <ACCESSION_ID> <human/mouse> <forward/reverse>
# design.txt: line 1 contains all Sample IDs of the dataset
#             line 2 contains each sample's condition
#             line 3 contains each sample's sequencing protocol
# Tab separated

set -e 

###########################################################################################################
# STEP 0
########
#Set env variables
source ~/conf.sh
export GSE="$1"
export SPECIES="$2"
export STRAND="$3"

readarray -t design <$WORK/datasets/$GSE/design.txt
SRR=($(echo "${design[0]}" | tr ' ' "\n"))
CONDITIONS=($(echo "${design[1]}" | tr ' ' "\n"))
TYPE=($(echo "${design[2]}" | tr ' ' "\n"))

#Prepare targets.txt
mkdir -p $WORK/datasets/$GSE/metaseqR_out
printf "%s\t%s\t%s\t%s\t%s\n" \
    "samplename" "filename" "condition" "paired" "stranded" > \
    $WORK/datasets/$GSE/metaseqR_out/targets.txt

for i in "${!SRR[@]}"
do
    printf "%s\t%s\t%s\t%s\t%s\n" \
    "${SRR[i]}" \
    "/media/galadriel/hybridstat/elixir_project/datasets/$GSE/hisat2_out/${SRR[i]}/${SRR[i]}.bam" \
    "${CONDITIONS[$i]}" \
    "${TYPE[$i]}" \
    "$STRAND" >>  \
    $WORK/datasets/$GSE/metaseqR_out/targets.txt
done

###########################################################################################################
# STEP 1
########
#Download SRR .fastq files
mkdir -p $WORK/datasets/$GSE/fastq

#convert SRP to SRRs: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP041753&go=go

for ACC in "${!SRR[@]}"
do 
    if [[ "${TYPE[$ACC]}" = 'single' ]]; then 
        echo "Downloading ${TYPE[$ACC]}-end ${SRR[ACC]}...";
        $WORK/tools/sratoolkit.2.9.6/bin/fastq-dump \
        --outdir  $WORK/datasets/$GSE/fastq \
        --gzip \
        --accession ${SRR[ACC]}
    elif [[ "${TYPE[$ACC]}" = 'paired' ]]; then 
        echo "Downloading ${TYPE[$ACC]}-end ${SRR[ACC]}..."
        $WORK/tools/sratoolkit.2.9.6/bin/fastq-dump \
        --outdir  $WORK/datasets/$GSE/fastq \
        --gzip \
        --accession ${SRR[ACC]} \
        --split-files
    fi
done

###########################################################################################################
# STEP 2a
########
# Quality control FastQC
mkdir -p $WORK/datasets/$GSE/fastqc

cd $WORK/datasets/$GSE/fastq
$WORK/tools/fastqc-0.11.8/fastqc \
    --outdir $WORK/datasets/$GSE/fastqc *.fastq.gz
while true; do
    read -p "Inspect FastQC files. Do you need to perform seqtk trimming?" yn
    case $yn in
        [Yy]* ) export fq="fastq_trim"; export fqc_dir="fastqc_trim"; echo "seqtk trimming enabled"; break;;
        [Nn]* ) export fq="fastq"; export fqc_dir="fastqc"; echo "seqtk trimming will be skipped"; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

# ###########################################################################################################
# # STEP 2b
# ########
# Trim bp ends (if prompted) and re-run FastQC
if [[ "$fq" = "fastq_trim" ]]; then
    mkdir -p $WORK/datasets/$GSE/fastq_trim
    for FILE in $WORK/datasets/$GSE/fastq/*.fastq.gz
    do 
        SAMPLE=`basename $FILE`
        TRIMMED=`basename $SAMPLE .gz`
        $WORK/tools/seqtk-1.3/seqtk trimfq $FILE > $WORK/datasets/$GSE/fastq_trim/$TRIMMED
        gzip $WORK/datasets/$GSE/fastq_trim/$TRIMMED
    done
    
    cd $WORK/datasets/$GSE/fastq_trim
    mkdir -p $WORK/datasets/$GSE/fastqc_trim
    
    $WORK/tools/fastqc-0.11.8/fastqc \
        --outdir $WORK/datasets/$GSE/fastqc_trim *.fastq.gz
fi


###########################################################################################################
# STEP 2c
########
# MultiQC report
mkdir -p $WORK/datasets/$GSE/multiqc

cd $WORK/datasets/$GSE/$fqc_dir

source /home/makis/cwlenv/bin/activate

multiqc . \
    -fd \
    -i $GSE \
    -b "Aggregate quality control for $GSE" \
    -o $WORK/datasets/$GSE/multiqc

deactivate

###########################################################################################################
# STEP 3
########
# Directory variables
export FASTQ_PATH=$WORK/datasets/$GSE/$fq
export HISAT2_OUTPUT=$WORK/datasets/$GSE/hisat2_out

# Command tool variables
export SAMTOOLS_COMMAND=$WORK/tools/samtools-1.9/samtools
export HISAT2_COMMAND=$WORK/tools/hisat2-2.1.0/hisat2
export BOWTIE2_COMMAND=$WORK/tools/bowtie2-2.3.5.1/bowtie2
export BEDTOOLS_COMMAND=$WORK/tools/bedtools-2.28.0/bin/bedtools

# Reference Genome Index
if [[ "$SPECIES" = 'human' ]]; then
    export HISAT2_INDEXES=$WORK/reference/hisat2/GRCh37/grch37_tran/genome_tran
    export BOWTIE2_INDEX=/media/raid/resources/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
elif [[ "$SPECIES" = 'mouse' ]]; then
    export HISAT2_INDEXES=$WORK/reference/hisat2/GRCm38/grcm38_tran/genome_tran
    export BOWTIE2_INDEX=/media/raid/resources/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome
else
    exit 1
fi

mkdir -p $WORK/datasets/$GSE/hisat2_out

for ACC in "${!SRR[@]}"
do
    bash $WORK/scripts/map_hisat2_advanced_"${TYPE[$ACC]}"-end.sh "${SRR[$ACC]}"
done

###########################################################################################################
# STEP 4a
########
#Bam to BedGraph
echo "Generating BedGraph files"
mkdir -p $WORK/datasets/$GSE/bigwig
for DIR in $WORK/datasets/$GSE/hisat2_out/*
do
    SAMPLE=`basename $DIR`
    $WORK/tools/bedtools-2.28.0/bin/bedtools genomecov -split -bg \
        -ibam $DIR/$SAMPLE.bam | \
        grep -vP 'chrM|chrU|chrG|chrJ|rand|hap|loc|cox' | sort -k1,1 -k2g,2 > \
        $WORK/datasets/$GSE/bigwig/$SAMPLE".bedGraph"
done

###########################################################################################################
# STEP 4b
########
#Calculate genome size
echo "Calcuating genome size"
SAMPLES=($(ls $WORK/datasets/$GSE/hisat2_out/))
SAMPLE=${SAMPLES[0]}
$WORK/tools/samtools-1.9/samtools view -H $WORK/datasets/$GSE/hisat2_out/$SAMPLE/$SAMPLE.bam | \
    grep '^@SQ' | sed s/@SQ\\tSN:// | \
    grep -vP 'chrM|chrU|chrG|chrJ|rand|hap|loc|cox' | sed s/\\tLN:/\\t/ > \
    $WORK/datasets/$GSE/bigwig/genome.size

###########################################################################################################
# STEP 4c
########
#Normalization
NCORES=8
echo "Normalizing"
cd $WORK/datasets/$GSE/bigwig
perl $WORK/scripts/normalize_bedgraph.pl \
    --input $WORK/datasets/$GSE/bigwig/*.bedGraph \
    --exportfactors $WORK/datasets/$GSE/bigwig/$GSE_factors.txt \
    --ncores $NCORES

###########################################################################################################
# STEP 4d
########
#Normalized BedGraph to BigWig
echo "Generating BigWig files"
GENOME="$WORK/datasets/$GSE/bigwig/genome.size"
for FILE in `ls $WORK/datasets/$GSE/bigwig/*_norm.bedGraph`
do
    SAMPLE=`basename $FILE | sed s/_norm\.bedGraph//`
    $WORK/tools/kenttools-1.04.0/bedGraphToBigWig $FILE $GENOME \
        $WORK/datasets/$GSE/bigwig/$SAMPLE".bigWig"
done


###########################################################################################################
# STEP 5b
########

if [[ "$SPECIES" = 'human' ]]; then
    genome="hg19"
elif [[ "$SPECIES" = 'mouse' ]]; then
    genome="mm10"
fi

echo "$genome" > genome

echo "Generating Analysis .Rda file"
Rscript $WORK/scripts/metaseqR_workflow.R $GSE $genome 0.12
