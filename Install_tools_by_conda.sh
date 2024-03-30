#!/usr/bin/env bash
# Here, we use conda to create a environment for tools and softwares install in linux system
# To install conda on Linux, see conda website: https://conda.io/projects/conda/en/latest/user-guide/install/linux.html

# Create a environment "PARE_analysis" and add two channels
conda create -n PARE_analysis
conda activate PARE_analysis
conda config --add channels bioconda
conda config --add channels conda-forge

# Install tools:
conda install -c bioconda htslib # Basic bio-tool
conda install -c bioconda samtools # Basic bio-tool
conda install -c bioconda entrez-direct # Download sequence from NCBI genebank
conda install -c bioconda gffread # deal with gff/gtf files
conda install -c bioconda cutadapt # modfied/trim NGS reads
conda install -c bioconda bowtie # short reads and non-splice aware aligner
conda install -c bioconda star # Splice-aware aligner
conda install -c bioconda deeptools # make bigwig file from bam/cram file using bamCoverage
conda install -c bioconda bedtools # make bed file for Metagene plot using bamToBed

# install R and R packages for metagene plot
conda install -c conda-forge r-base
conda install -c bioconda bioconductor-guitar
conda install -c conda-forge r-data.table
conda install -c conda-forge r-r.utils
conda install -c conda-forge r-ggplot2
conda install -c conda-forge r-argparser

# Download GuitarPlotFast.R for metagene plot
wget -O $CONDA_PREFIX/bin/GuitarPlotFast.R https://raw.githubusercontent.com/LabHMChenABRC/PARE-analysis/main/Metagene-plot/GuitarPlotFast.R
chmod +x $CONDA_PREFIX/bin/GuitarPlotFast.R
