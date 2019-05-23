#!/bin/bash

###########################################################################################################
# STEP 0
########
#Set env variables
source ~/cong.sh
GSE="$1"
SPECIES="$2"
STRAND="$3"
PAIRED="$4"


IFS=$'\n' read -d '' -r -a design < $WORK/datasets/$GSE/design.txt
SRR=($(echo "${design[0]}" | tr ' ' "\n"))
CONDITIONS=($(echo "${design[1]}" | tr ' ' "\n"))

#SRR=('SRR3185384' 'SRR3185385' 'SRR3185386' 'SRR3185387' 'SRR3185388' 'SRR3185389')
#CONDITIONS=('WT' 'WT' 'WT' 'ENG_KO' 'ENG_KO' 'ENG_KO')


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
    "${CONDITIONS[$i]}"  \
    "$PAIRED" \
    "$STRAND" >> \
    $WORK/datasets/$GSE/metaseqR_out/targets.txt
done 
###########################################################################################################
# STEP 1
########
#Download SRR .fastq files
mkdir -p $WORK/datasets/$GSE/fastq

#convert SRP to SRRs: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP041753&go=go

for ACC in "${SRR[@]}"
do
    echo "Downloading $ACC..."
    $WORK/tools/sratoolkit.2.9.6/bin/fastq-dump \
        --outdir  $WORK/datasets/$GSE/fastq \
        --gzip \
        --accession $ACC -v
done

###########################################################################################################
# STEP 2a
########
#Quality control FastQC
mkdir -p $WORK/datasets/$GSE/fastqc

cd $WORK/datasets/$GSE/fastq
$WORK/tools/fastqc-0.11.8/fastqc \
    --outdir $WORK/datasets/$GSE/fastqc *.fastq.gz
while true; do
    read -p "Inspect FastQC files. Do you need to perform seqtk trimming?" yn
    case $yn in
        [Yy]* ) fq="fastq_trim"; echo "seqtk trimming enabled"; break;;
        [Nn]* ) fq="fastq"; echo "seqtk trimming will be skipped"; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

###########################################################################################################
# STEP 2b
########
#Trim bp ends (if prompted)
if [[ "$fq" = "fastq_trim" ]]; then
    mkdir -p $WORK/datasets/$GSE/fastq_trim
    for FILE in $WORK/datasets/GSE57397/fastq/*.fastq.gz
    do 
        SAMPLE=`basename $FILE`
        TRIMMED=`basename $SAMPLE .gz`
        $WORK/tools/seqtk-1.3/seqtk trimfq $FILE > $WORK/datasets/GSE57397/fastq_trim/$TRIMMED
        gzip $WORK/datasets/GSE57397/fastq_trim/$TRIMMED
    done
fi


###########################################################################################################
# STEP 3
########
#Alignment to reference
mkdir -p $WORK/datasets/$GSE/hisat2_out
bash $WORK/scripts/map_hisat2_advanced_$PAIRED-end.sh $GSE $SPECIES

#### ADD FUNCTIONALITY FOR PAIRED-END

###########################################################################################################
# STEP 4a
########
#Bam to BedGraph
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
cd $WORK/datasets/$GSE/bigwig
perl $WORK/scripts/normalize_bedgraph.pl \
    --input $WORK/datasets/$GSE/bigwig/*.bedGraph \
    --exportfactors $WORK/datasets/$GSE/bigwig/GSE78271_factors.txt \
    --ncores $NCORES

###########################################################################################################
# STEP 4d
########
#Normalized BedGraph to BigWig
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

Rscript $WORK/scripts/metaseqR_workflow.R $GSE $genome 0.12
