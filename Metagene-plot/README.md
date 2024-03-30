
# GuitarPlotFast.R 
It is used to make a metagene plot based on PARE/Degradome dataset.

## Installation

Download GuitarPlotFast.R in your $HOME/bin direction which should be added to $PATH variable
```
wget -O ~/bin/GuitarPlotFast.R https://github.com/LabHMChenABRC/PARE-analysis/Metagene-plot/raw/main/GuitarPlotFast.R
```
make GuitarPlotFast.R executable 
```
chmod +x $CONDA_PREFIX/bin/GuitarPlotFast.R
```
## Usage
1. Make density files:
``` shell
GuitarPlotFast.R -b <Directory_of_bed.gz> -g <gtf> -o <Output_folder>

-b <Directory_to_bed> contains "bed.gz" files created by bedtools based on genomic bam files which are produced by STAR.
                      bedtools bamtobed -i <filename.bam/cram> | bgzip -@ <number of core  >directory_to_bed/filename.bed.gz
-g <gtf>              gtf file, using gffread can convert gff to gtf format:
                      gffread -T -o <.gtf> <.gff>
-o <Output_folder>    program will save .mrna.density files in this direction
```
This command will extract coding genes with 5'UTR, CDS, 3'UTR, upstream sequence, and downstream sequence lengths of at least 100bp each. And count the 5'Ends of unique alignments (MAPQ=255) on the forward strand of selected genes. The density distribution of the upstream sequence, 5'UTR, CDS, 3'UTR, and downstream sequence is scaled to 1:2:4:2:1 ratio.
STAR assigns unique alignments with MAPQ 255. 
It takes less than 15 min to calculate the density for 20M library

2. Create metagene plot:
``` shell
GuitarPlotFast.R -m -p <output_file_prefix> -d <metagenePlot_info_file> -o <output_folder>
-m                          switch to plot metagene figure
-p <output_file_prefix>     output file prefix
-d <metagenePlot_info_file> contains a header and give density file location and they group for plot
-o <output_folder>          metagene plot is saved in pdf here.
```
The below example of <metagenePlot_info_file> will used to create a plot with group1 WT, MutA and group2 WT, MutB 
| file                                 | plotgroup |
| ------------------------------------ | --------- |
| Path-of-density/WT.mrna.density      | 1,2       |
| Path-of-density/xrn4-6.mrna.density  | 1         |
| Path-of-density/fry1-6.mrna.density  | 2         |

