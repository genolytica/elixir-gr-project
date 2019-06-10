#!/bin/bash

set -e

# Number of cores to use
CORES=8

#Execute HISAT2
SAMPLE="$1"
echo "===== Mapping with hisat2 for $SAMPLE..."
mkdir -p $HISAT2_OUTPUT/$SAMPLE
$HISAT2_COMMAND -p $CORES \
    --add-chrname \
    --un $HISAT2_OUTPUT/$SAMPLE/unmapped.fastq \
    -x $HISAT2_INDEXES \
	-U $FASTQ_PATH/$SAMPLE".fastq.gz" \
	-S $HISAT2_OUTPUT/$SAMPLE/hisat2.sam
echo " "
echo "===== Trying now to map unmapped reads with bowtie2 for $SAMPLE..."
$BOWTIE2_COMMAND --local --very-sensitive-local \
    -p $CORES \
    -x $BOWTIE2_INDEX \
    -U $HISAT2_OUTPUT/$SAMPLE/unmapped.fastq | \
    $SAMTOOLS_COMMAND view -bhS \
    -o $HISAT2_OUTPUT/$SAMPLE/unmapped_remap_tmp.bad -
echo "===== Merging all reads for $SAMPLE..." 
$SAMTOOLS_COMMAND view -@ $CORES -HS $HISAT2_OUTPUT/$SAMPLE/hisat2.sam > \
    $HISAT2_OUTPUT/$SAMPLE/header.sam
# A reheader MUST be done as Bowtie2 returns different header than 
# hisat2 ??!!?? (@&#^*&^$
$SAMTOOLS_COMMAND reheader $HISAT2_OUTPUT/$SAMPLE/header.sam \
    $HISAT2_OUTPUT/$SAMPLE/unmapped_remap_tmp.bad > \
    $HISAT2_OUTPUT/$SAMPLE/unmapped_remap_tmp.uns
$SAMTOOLS_COMMAND sort -@ $CORES \
    $HISAT2_OUTPUT/$SAMPLE/unmapped_remap_tmp.uns \
    -o $HISAT2_OUTPUT/$SAMPLE/unmapped_remap.bam
# SAM file created by HISAT2 must be transformed into BAM before merging 
# the files
$SAMTOOLS_COMMAND view -@ $CORES -bhS $HISAT2_OUTPUT/$SAMPLE/hisat2.sam > \
    $HISAT2_OUTPUT/$SAMPLE/hisat2.bam
$SAMTOOLS_COMMAND merge -@ $CORES -f $HISAT2_OUTPUT/$SAMPLE/tmp.bam \
    $HISAT2_OUTPUT/$SAMPLE/hisat2.bam \
    $HISAT2_OUTPUT/$SAMPLE/unmapped_remap.bam
echo "===== Coordinate sorting all reads for $SAMPLE..."
$SAMTOOLS_COMMAND sort -@ $CORES $HISAT2_OUTPUT/$SAMPLE/tmp.bam \
    -o $HISAT2_OUTPUT/$SAMPLE/$SAMPLE.bam
echo "===== Indexing all merged reads for $SAMPLE..."
$SAMTOOLS_COMMAND index -@ $CORES $HISAT2_OUTPUT/$SAMPLE/$SAMPLE.bam
echo "===== Removing intermediate garbage for $SAMPLE..."
rm $HISAT2_OUTPUT/$SAMPLE/hisat2.bam $HISAT2_OUTPUT/$SAMPLE/tmp.bam \
    $HISAT2_OUTPUT/$SAMPLE/unmapped*.fastq \
    $HISAT2_OUTPUT/$SAMPLE/unmapped_remap.bam \
    $HISAT2_OUTPUT/$SAMPLE/unmapped_remap_tmp.bad \
    $HISAT2_OUTPUT/$SAMPLE/unmapped_remap_tmp.uns
echo " "
