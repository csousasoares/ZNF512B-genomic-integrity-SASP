#!/bin/bash

INDEX="/home/cdsou/bioinformatics/hisat2_indexes/hg38_tran/genome_tran"
THREADS=12
INDIR="./1_fastp"
OUTDIR="./2_hisat2"

mkdir -p $OUTDIR

for FILE in $INDIR/*.trimmed.fastq.gz; do

SAMPLE_NAME=$(basename $FILE .trimmed.fastq.gz)

echo "Aligning $SAMPLE_NAME ..."

hisat2 \
-x $INDEX \
-U $FILE \
-p $THREADS \
--dta \
-S $OUTDIR/${SAMPLE_NAME}.sam \
2> $OUTDIR/${SAMPLE_NAME}_alignment_stats.txt

samtools view -bS $OUTDIR/${SAMPLE_NAME}.sam | \
samtools sort -o $OUTDIR/${SAMPLE_NAME}.sorted.bam
samtools index $OUTDIR/${SAMPLE_NAME}.sorted.bam
rm $OUTDIR/${SAMPLE_NAME}.sam

done
