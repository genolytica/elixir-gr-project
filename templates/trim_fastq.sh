#!/bin/bash
HOME_PATH=/media/raid/data/brfaa/charonis
FASTQ_PATH=$HOME_PATH/fastq
FASTQ_PATTERN=_sequence.txt
SEQTK_COMMAND=/opt/NGSTools/seqtk/seqtk
TRIM_OUT=$HOME_PATH/fastq_proc

# Number of cores to use
CORES=12

if [ -d $TRIM_OUT ]
then
	break;
else
	mkdir -p $TRIM_OUT
fi

# Trim with multiple cores
parallel -j $CORES --gnu 'SAMPLE=`basename {} | sed s/_sequence\.txt//`; /opt/NGSTools/seqtk/seqtk trimfq -b 1 {} > /media/raid/data/brfaa/charonis/fastq_proc/$SAMPLE.fastq' ::: $FASTQ_PATH/*$FASTQ_PATTERN

## Trim with single core
#for FILE in $FASTQ_PATH/*$FASTQ_PATTERN
#do
#	SAMPLE=`basename $FILE | sed s/\$FASTQ_PATTERN//`
#	echo "Processing $SAMPLE "
#	$SEQTK_COMMAND trimfq -b 1 $FILE > $TRIM_OUT/$SAMPLE.fastq
#done
