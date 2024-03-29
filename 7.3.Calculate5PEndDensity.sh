#!/usr/bin/env bash
wd=$PWD
BedDir=$wd/PARE/alignments/bed
gtfFile=$wd/Reference_sequence/TAIR10_representative_gene_models.gtf
Rscript $wd/GuitarPlotFast.R -b $BedDir -g $gtfFile -o $wd/PARE/metagenePlot
