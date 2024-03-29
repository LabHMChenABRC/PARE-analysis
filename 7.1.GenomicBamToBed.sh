#!/usr/bin/env bash
wd=$PWD
core=16
STAR_alignments_dir=$wd/PARE/alignments
OutDir_Bed=$STAR_alignments_dir/bed
mkdir -p $OutDir_Bed
date +"%b %d %T ..... Covert Bam file to bed format by bedtools"
for Genomic_Bam in $STAR_alignments_dir/*_Aligned.sortedByCoord.out.bam
do
	ID=$(basename ${Genomic_Bam%_Aligned.sortedByCoord.out.bam})
	date +"%b %d %T ..... Process $ID"
	bedtools bamtobed -i $Genomic_Bam | bgzip -@ $core >$OutDir_Bed/$ID.bed.gz
done
date +"%b %d %T ..... Finish"
