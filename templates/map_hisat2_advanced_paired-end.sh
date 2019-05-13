#!/bin/bash

# Directory variables
WORK_DIR=MAIN_PROJECT_DIR
ACCESSION=GSEXXXXX
FASTQ_PATH=$WORK_DIR/datasets/$ACCESSION/fastq
HISAT2_OUTPUT=$WORK_DIR/datasets/$ACCESSION/hisat2_out

# Command tool variables
SAMTOOLS_COMMAND=$WORK_DIR/tools/samtools-1.9/samtools
HISAT2_COMMAND=$WORK_DIR/tools/hisat2-2.1.0/hisat2
BOWTIE2_COMMAND=$WORK_DIR/tools/bowtie2-2.3.5.1/bowtie2
BEDTOOLS_COMMAND=$WORK_DIR/tools/bedtools-2.28.0/bin/bedtools

# Reference Genome Index
HISAT2_INDEX=$WORK_DIR/reference/hisat2/GRCm38/grcm38_tran/genome_tran
#Bowtie Index
BOWTIE2_INDEX=$WORK_DIR/reference/bowtie2/mm10/genome

# Number of cores to use
CORES=16

#Execute HISAT2
for FILE in $FASTQ_PATH/*_1.fastq.gz
do	
	SAMPLE=`basename $FILE | sed s/\_1.fastq\.gz//`
	echo "===== Mapping with hisat2 for $SAMPLE..."
	mkdir -p $HISAT2_OUTPUT/$SAMPLE
	$HISAT2_COMMAND -p $CORES \
        --add-chrname \
        --un-conc $HISAT2_OUTPUT/$SAMPLE/unmapped.fastq \
        -x $HISAT2_INDEX \
        -1 $FASTQ_PATH/$SAMPLE"_1.fastq.gz" \
        -2 $FASTQ_PATH/$SAMPLE"_2.fastq.gz" | \
        $SAMTOOLS_COMMAND view -bS > $HISAT2_OUTPUT/$SAMPLE/hisat2.bam
	echo " "
	echo "===== Trying now to map unmapped reads with bowtie2 for $SAMPLE..."
	$BOWTIE2_COMMAND --local --very-sensitive-local --dovetail \
        -p $CORES \
        -x $BOWTIE2_INDEX \
        -1 $HISAT2_OUTPUT/$SAMPLE/unmapped.1.fastq \
        -2 $HISAT2_OUTPUT/$SAMPLE/unmapped.2.fastq | \
        $SAMTOOLS_COMMAND view -bhS \
        -o $HISAT2_OUTPUT/$SAMPLE/unmapped_remap_tmp.bad -
    echo "===== Merging all reads for $SAMPLE..." 
	$SAMTOOLS_COMMAND view -@ $CORES -H $HISAT2_OUTPUT/$SAMPLE/hisat2.bam > \
        $HISAT2_OUTPUT/$SAMPLE/header.sam
	# A reheader MUST be done as Bowtie2 returns different header than 
	# hisat2 ??!!?? (@&#^*&^$
	$SAMTOOLS_COMMAND reheader $HISAT2_OUTPUT/$SAMPLE/header.sam \
        $HISAT2_OUTPUT/$SAMPLE/unmapped_remap_tmp.bad > \
        $HISAT2_OUTPUT/$SAMPLE/unmapped_remap_tmp.uns
    $SAMTOOLS_COMMAND sort -@ $CORES \
        -o $HISAT2_OUTPUT/$SAMPLE/unmapped_remap.bam \
        $HISAT2_OUTPUT/$SAMPLE/unmapped_remap_tmp.uns
    $SAMTOOLS_COMMAND merge -@ $CORES -f $HISAT2_OUTPUT/$SAMPLE/tmp.bam \
        $HISAT2_OUTPUT/$SAMPLE/hisat2.bam \
        $HISAT2_OUTPUT/$SAMPLE/unmapped_remap.bam
	echo "===== Coordinate sorting all reads for $SAMPLE..."
    $SAMTOOLS_COMMAND sort -@ $CORES -o $HISAT2_OUTPUT/$SAMPLE/$SAMPLE".bam" \
        $HISAT2_OUTPUT/$SAMPLE/tmp.bam
        
    echo "===== Indexing all merged reads for $SAMPLE..."
    $SAMTOOLS_COMMAND index -@ $CORES $HISAT2_OUTPUT/$SAMPLE/$SAMPLE.bam
    echo "===== Removing intermediate garbage for $SAMPLE..."
    rm $HISAT2_OUTPUT/$SAMPLE/hisat2.bam $HISAT2_OUTPUT/$SAMPLE/tmp.bam \
        $HISAT2_OUTPUT/$SAMPLE/unmapped*.fastq \
        $HISAT2_OUTPUT/$SAMPLE/unmapped_remap.bam \
        $HISAT2_OUTPUT/$SAMPLE/unmapped_remap_tmp.bad \
        $HISAT2_OUTPUT/$SAMPLE/unmapped_remap_tmp.uns
    echo " "
done
