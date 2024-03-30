# PARE analysis
 Bioinformatics workflow of parallel analysis of RNA ends (PARE) sequencing in Linux OS
 
## Required tools and install:
Run bash script `Install_tools_by_conda.sh` to build a conda environment 

## Pipelines 
This pipeline includes pre-processing of raw data (adapter trimming and quality control), mapping of clean reads (genome and transcriptome alignments), count quantification (assignment and normalization), visualization of RNA decay profiles (plot of PARE counts along transcripts or specific sites), and global analysis/comparison (metagene plots). The parameters or steps specifically are required for PARE data analysis.
* Prefix the number of bash script (.sh) means analysis order.
* A command line tool `GuitarPlotFast.R` is used to make mategene plot. Please see README.md in `PARE-analysis/Metagene-plot`.
