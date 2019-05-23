#!/bin/bash
set -e 

source ~/cong.sh

# Directory variables
FASTQ_PATH=$WORK/datasets/$GSE/$fq
HISAT2_OUTPUT=$WORK/datasets/$GSE/hisat2_out

# Command tool variables
SAMTOOLS_COMMAND=$WORK/tools/samtools-1.9/samtools
HISAT2_COMMAND=$WORK/tools/hisat2-2.1.0/hisat2
BOWTIE2_COMMAND=$WORK/tools/bowtie2-2.3.5.1/bowtie2
BEDTOOLS_COMMAND=$WORK/tools/bedtools-2.28.0/bin/bedtools

# Reference Genome Index
if [[ "$SPECIES" = 'human' ]]; then
	HISAT2_INDEXES=$WORK/reference/hisat2/GRCh37/grch37_tran/genome_tran
	BOWTIE2_INDEX=/media/raid/resources/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
elif [[ "$SPECIES" = 'mouse' ]]; then
	HISAT2_INDEXES=$WORK/reference/hisat2/GRCm38/grcm38_tran/genome_tran
	BOWTIE2_INDEX=/media/raid/resources/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome
else
	exit 1
fi

# Number of cores to use
CORES=10

#Execute HISAT2
for FILE in $FASTQ_PATH/*_1.fastq.gz
do	
	SAMPLE=`basename $FILE | sed s/\_1.fastq\.gz//`
	echo "===== Mapping with hisat2 for $SAMPLE..."
	mkdir -p $HISAT2_OUTPUT/$SAMPLE
	$HISAT2_COMMAND -p $CORES \
	    --add-chrname \
	    --un-conc $HISAT2_OUTPUT/$SAMPLE/unmapped.fastq \
	    -x $HISAT2_INDEXES \
		-1 $FASTQ_PATH/$SAMPLE"_1.fastq.gz" \
		-2 $FASTQ_PATH/$SAMPLE"_2.fastq.gz" \
		-S $HISAT2_OUTPUT/$SAMPLE/hisat2.sam
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
done
