#!/bin/bash

THREADS=12
INDIR="./0_fastq"
OUTDIR="./1_fastp"

mkdir -p $OUTDIR

for FILE in $INDIR/*.fastq.gz; do

SAMPLE_NAME=$(basename $FILE .fastq.gz)


echo "Filtering $SAMPLE_NAME..."

fastp \
-i $FILE \
-o $OUTDIR/${SAMPLE_NAME}.trimmed.fastq.gz \
--thread $THREADS \
-h $OUTDIR/${SAMPLE_NAME}_fastp.html
done
